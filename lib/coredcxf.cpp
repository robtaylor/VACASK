#include <iomanip>
#include <cmath>
#include <filesystem>
#include "coredcxf.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

// Default parameters
DCXFParameters::DCXFParameters() {
    opParams.write = 0;
}

template<> int Introspection<DCXFParameters>::setup() {
    registerMember(out);
    registerMember(writeop);
    registerMember(write);
    registerNamedMember(opParams.nodeset, "nodeset");
    registerNamedMember(opParams.store, "store");
    
    return 0;
}
instantiateIntrospection(DCXFParameters);


DCXFCore::DCXFCore(
    OutputDescriptorResolver& parentResolver, DCXFParameters& params, OperatingPointCore& opCore, 
    std::unordered_map<Id,size_t>& sourceIndex, Circuit& circuit, 
    KluRealMatrix& jacobian, Vector<double>& incrementalSolution, 
    std::vector<Instance*>& sources, Vector<double>& tf, 
    Vector<double>& yin, Vector<double>& zin
) : AnalysisCore(parentResolver, circuit), params(params), outfile(nullptr), opCore_(opCore), sourceIndex(sourceIndex), 
    jacobian(jacobian), incrementalSolution(incrementalSolution), 
    sources(sources), tf(tf), yin(yin), zin(zin) {
    
    // Set analysis type for the initial operating point analysis
    auto& elsSystem = opCore_.solver().evalSetup();
    elsSystem.staticAnalysis = true;
    elsSystem.dcAnalysis = false;
    elsSystem.acAnalysis = true;
}

DCXFCore::~DCXFCore() {
    delete outfile;
}

bool DCXFCore::resolveOutputDescriptors(bool strict) {
    clearError();
    // Clear output sources
    outputSources.clear();
    // Clear source instance pointers, initialize to nullptrs
    sources.clear();
    sources.resize(sourceIndex.size(), nullptr);
    // Resolve output descriptors
    bool ok = true; 
    for (auto it = outputDescriptors.cbegin(); it != outputDescriptors.cend(); ++it) {
        Id name;
        size_t ndx;
        Instance *inst;
        switch (it->type) {
            case OutdTf:
            case OutdZin:    
            case OutdYin:
                // Get instance name and index
                name = it->idNdx.id;
                ndx = it->idNdx.ndx;
                // Find instance
                inst = circuit.findInstance(name);
                sources[ndx] = inst;
                if (strict) {
                    if (!inst) {
                        setError(DCXFError::NotFound);
                        errorInstance = name;
                        ok = false;
                        break;
                    }
                }
                // Instance found, but is not a source... this is always an error
                if (inst && !inst->model()->device()->isSource()) {
                    setError(DCXFError::NotSource);
                    errorInstance = name;
                    ok = false;
                    break;
                }
                // No instance and we reached this point, create constant source
                if (!inst) {
                    outputSources.emplace_back();
                }
                break;
        }
        if (!ok) {
            break;
        }
        switch (it->type) {
        case OutdTf:
            outputSources.emplace_back(&tf, ndx);
            break;
        case OutdZin:
            outputSources.emplace_back(&zin, ndx);
            break;
        case OutdYin:
            outputSources.emplace_back(&yin, ndx);
            break; 
        default:
            // Delegate to parent
            ok = parentResolver.resolveOutputDescriptor(*it, outputSources, strict);
        }
        if (!ok) {
            break;
        }
    }
    return ok;
}

bool DCXFCore::addDefaultOutputDescriptors() {
    // If output is suppressed, skip all this work
    if (!params.write) {
        return true;
    }
    if (savesCount==0) {
        return addAllTfZin(PTSave(Loc::bad, "default", Id(), Id()), sourceIndex);
    }
    return true;
}

bool DCXFCore::initializeOutputs(Id name, Status& s) {
    if (!params.write) {
        return true;
    }
    // Create output file if not created yet
    if (!outfile) {
        outfile = new OutputRawfile(
            name, outputDescriptors, outputSources,
            (circuit.simulatorOptions().core().rawfile==SimulatorOptions::rawfileBinary ? OutputRawfile::Flags::Binary : OutputRawfile::Flags::None) |
                OutputRawfile::Flags::Padded);
        outfile->setTitle(circuit.title());
        outfile->setPlotname("DC Transfer Function Analysis");
    }
    outfile->prologue();

    return true;
}

bool DCXFCore::finalizeOutputs(Status& s) {
    if (outfile) {
        outfile->epilogue();
        delete outfile;
        outfile = nullptr;
    }
    return true;
}

bool DCXFCore::deleteOutputs(Id name, Status& s) {
    if (!params.write) {
        return true;
    }

    // Cannot assume outfile is available
    auto fname = std::string(name)+".raw";
    if (std::filesystem::exists(fname)) {
        std::filesystem::remove(fname);
    }
    return true;
}
    
bool DCXFCore::rebuild(Status& s) {
    return true;
}

// System of equations is 
//   G(x) dx  = dJ
CoreCoroutine DCXFCore::coroutine(bool continuePrevious) {
    initProgress(1, 0);

    jacobian.setAccounting(circuit.tables().accounting());

    clearError();
    // Make sure structures are large enough
    incrementalSolution.resize(circuit.unknownCount()+1);
    tf.resize(sources.size());
    yin.resize(sources.size());
    zin.resize(sources.size());

    // Get output unknowns
    auto [ok, up, un] = getOutput(params.out);
    if (!ok) {
        co_yield CoreState::Aborted;
    }
    
    // Compute operating point
    auto opOk = opCore_.run(continuePrevious);
    if (!opOk) {
        setError(DCXFError::OperatingPointError);
        co_yield CoreState::Aborted;
    }

    auto& options = circuit.simulatorOptions().core();
    Int debug = options.smsig_debug;

    if (debug>0) {
        Simulator::dbg() << "Starting DC transfer function analysis.\n";
    }

    // Jacobian is already factored (done by op core)

    // Get RHS vector
    auto rhsVec = incrementalSolution.data();

    // Go through all listed sources, compute TF, Zin, and Yin
    auto nSrc = sources.size();
    bool error = false;
    for(decltype(nSrc) i=0; i<nSrc; i++) {
        // Prepare RHS
        zero(incrementalSolution);

        // Get instance
        auto inst = sources[i];
        // No instance, continue to next
        if (!inst) {
            continue;
        }

        if (debug>0) {
            Simulator::dbg() << "Computing response to '" << std::string(inst->name()) << "'.\n";
        }

        // Get excitation equations and response unknowns
        auto [e1, e2] = inst->sourceExcitation(circuit);
        auto [r1, r2] = inst->sourceResponse(circuit);
        
        // Set RHS, load negated residual to get the true value of delta after solve()
        rhsVec[e1] += -inst->scaledUnityExcitation();
        rhsVec[e2] -= -inst->scaledUnityExcitation();

        if (debug>=100) {
            Simulator::dbg() << "Linear system for instance " << inst->name() << "\n";
            jacobian.dump(Simulator::dbg(), dataWithoutBucket(incrementalSolution)); 
            Simulator::dbg() << "\n";
        }

        // Solve
        if (!jacobian.solve(dataWithoutBucket(incrementalSolution))) {
            setError(DCXFError::MatrixError);
            error = true;
            break;
        }
        
        // Set bucket to 0
        rhsVec[0] = 0.0;

        if (options.solutioncheck && !jacobian.isFinite(dataWithoutBucket(incrementalSolution), true, true)) {
            setError(DCXFError::SolutionError);
            if (options.smsig_debug) {
                Simulator::dbg() << "A solution entry is not finite. Solver failed.\n";
            }
            error = true;
            break;
        }
        
        // Compute TF
        tf[i] = rhsVec[up] - rhsVec[un];

        // Compute Yin and Zin
        if (inst->model()->device()->isVoltageSource()) {
            // Voltage source excitation
            yin[i] = (rhsVec[r1] - rhsVec[r2])*inst->responseScalingFactor()/inst->scaledUnityExcitation();
            if (yin[i]!=0.0) {
                zin[i] = 1.0/yin[i];
            } else {
                // Infinity
                zin[i] = 1e20;
            }
        } else {
            // Current source excitation
            zin[i] = (rhsVec[r1] - rhsVec[r2])*inst->responseScalingFactor()/inst->scaledUnityExcitation();
            if (zin[i]!=0.0) {
                yin[i] = 1.0/zin[i];
            } else {
                // Infinity
                yin[i] = 1e20;
            }
        }
    }

    if (debug>0) {
        if (error) {
            Simulator::dbg() << "DC transfer function analysis " << (error ? "exited prematurely" : "completed") << ".\n";
        }
    }
    
    // Dump solution
    if (params.write && outfile) {
        outfile->addPoint();
    }
    
    setProgress(1);

    co_yield CoreState::Finished;
}

bool DCXFCore::run(bool continuePrevious) {
    auto c = coroutine(continuePrevious);
    bool ok = true;
    while (!c.done()) {
        if (c.resume()==CoreState::Aborted) {
            ok = false;
            break;
        };
    }
    return ok;
}

bool DCXFCore::formatError(Status& s) const {
    auto nr = UnknownNameResolver(circuit);
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4);
    
    // First, handle AnalysisCore errors
    if (lastError!=Error::OK) {
        AnalysisCore::formatError(s);
        return false;
    }
    
    // Then handle DCXFCore errors
    switch (lastDcTfError) {
        case DCXFError::NotFound:
            s.set(Status::Analysis, std::string("Source '")+std::string(errorInstance)+"' not found.");
            break;
        case DCXFError::NotSource:
            s.set(Status::Analysis, std::string("Instance '")+std::string(errorInstance)+"' is not a source.");
            break;
        case DCXFError::EvalAndLoad:
            s.set(Status::Analysis, "Jacobian evaluation failed.");
            break;
        case DCXFError::MatrixError:
            jacobian.formatError(s, &nr);
            break;
        case DCXFError::SolutionError:
            jacobian.formatError(s, &nr);
            s.extend("Solution component is not finite.");
            break;
        case DCXFError::OperatingPointError:
            opCore_.formatError(s);
            break;
        default:
            return true;
    }
    s.extend("Leaving DC incremental analysis.");
    return false;
}

void DCXFCore::dump(std::ostream& os) const {
    AnalysisCore::dump(os);
    os << "  Results\n";
    auto nSrc = sources.size();
    for(decltype(nSrc) i=0; i<nSrc; i++) {
        auto inst = sources[i];
        if (!inst) {
            continue;
        }
        os << "    tf(" << inst->name() << ") " << tf[i] << "\n";
        os << "    zin(" << inst->name() << ") " << zin[i] << "\n";
        os << "    yin(" << inst->name() << ") " << yin[i] << "\n";
    }
}


}
