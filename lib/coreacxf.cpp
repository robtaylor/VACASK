#include <iomanip>
#include <cmath>
#include <filesystem>
#include "coreacxf.h"
#include "simulator.h"
#include "answeep.h"
#include "context.h"
#include "common.h"
#include <numbers>

namespace NAMESPACE {

// Default parameters
ACXFParameters::ACXFParameters() {
    opParams.writeOutput = 0;
}

template<> int Introspection<ACXFParameters>::setup() {
    registerMember(out);
    registerMember(from);
    registerMember(to);
    registerMember(step);
    registerMember(mode);
    registerMember(points);
    registerMember(values);
    registerMember(dumpop);
    registerNamedMember(opParams.nodeset, "nodeset");
    registerNamedMember(opParams.store, "store");

    return 0;
}
instantiateIntrospection(ACXFParameters);


ACXFCore::ACXFCore(
    OutputDescriptorResolver& parentResolver, ACXFParameters& params, OperatingPointCore& opCore, std::unordered_map<Id,size_t>& sourceIndex, 
    Circuit& circuit, 
    KluRealMatrix& dcJacobian, VectorRepository<double>& dcSolution, VectorRepository<double>& dcStates, 
    KluComplexMatrix& acMatrix, Vector<Complex>& acSolution, 
    std::vector<Instance*>& sources, Vector<Complex>& tf, Vector<Complex>& yin, Vector<Complex>& zin
) : AnalysisCore(parentResolver, circuit), params(params), outfile(nullptr), opCore_(opCore), sourceIndex(sourceIndex), 
    dcSolution(dcSolution), dcStates(dcStates), dcJacobian(dcJacobian), 
    acMatrix(acMatrix), acSolution(acSolution), sources(sources), tf(tf), yin(yin), zin(zin) {
    
    // Set analysis type for the initial operating point analysis
    auto& elsSystem = opCore_.solver().evalSetup();
    elsSystem.staticAnalysis = true;
    elsSystem.dcAnalysis = false;
    elsSystem.acAnalysis = true;
}

ACXFCore::~ACXFCore() {
    delete outfile;
}

bool ACXFCore::resolveOutputDescriptors(bool strict) {
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
                        setError(ACXFError::NotFound);
                        errorInstance = name;
                        ok = false;
                        break;
                    }
                }
                // Instance found, but is not a source... this is always an error
                if (inst && !inst->model()->device()->isSource()) {
                    setError(ACXFError::NotSource);
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
        case OutdFrequency:
            outputSources.emplace_back(&frequency);
            break;
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

bool ACXFCore::addCoreOutputDescriptors() {
    clearError();
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    if (!addOutputDescriptor(OutputDescriptor(OutdFrequency, "frequency"))) {
        lastError = Error::Descriptor;
        errorId = "frequency";
        return false;
    }
    return true;
}

bool ACXFCore::addDefaultOutputDescriptors() {
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    if (savesCount==0) {
        return addAllTfZin(PTSave(Loc::bad, "default", Id(), Id()), sourceIndex);
    }
    return true;
}

bool ACXFCore::initializeOutputs(Id name, Status& s) {
    // Create output file if not created yet
    if (!outfile) {
        outfile = new OutputRawfile(
            name, outputDescriptors, outputSources,
            (circuit.simulatorOptions().core().rawfile==SimulatorOptions::rawfileBinary ? OutputRawfile::Flags::Binary : OutputRawfile::Flags::None) |
                OutputRawfile::Flags::Padded | OutputRawfile::Flags::Complex);
        outfile->setTitle(circuit.title());
        outfile->setPlotname("AC Small Signal Transfer Function Analysis");
    }
    outfile->prologue();

    return true;
}

bool ACXFCore::finalizeOutputs(Status& s) {
    outfile->epilogue();
    delete outfile;
    outfile = nullptr;
    return true;
}

bool ACXFCore::deleteOutputs(Id name, Status& s) {
    if (!params.writeOutput) {
        return true;
    }

    // Cannot assume outfile is available
    auto fname = std::string(name)+".raw";
    if (std::filesystem::exists(fname)) {
        std::filesystem::remove(fname);
    }
    return true;
}
    
bool ACXFCore::rebuild(Status& s) {
    // AC analysis matrix
    if (!acMatrix.rebuild(circuit.sparsityMap(), circuit.unknownCount())) {
        acMatrix.formatError(s);
        return false;
    }
    
    // Resistive Jacobian entries remain bound to OP Jacobian, 
    // reactive parts will be bound to imaginary entries of acMatrix
    if (!circuit.bind(nullptr, Component::Real, std::nullopt, &acMatrix, Component::Imaginary, std::nullopt, s)) {
        return false;
    }
    
    return true;
}

// System of equations is 
//   (G(x) + i C(x)) dx = dJ
CoreCoroutine ACXFCore::coroutine(bool continuePrevious) {
    acMatrix.setAccounting(circuit.tables().accounting());
    
    clearError();

    auto n = circuit.unknownCount(); 
    // Make sure structures are large enough
    acSolution.resize(n+1);
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
        setError(ACXFError::OperatingPointError);
        co_yield CoreState::Aborted;
    }

    auto& options = circuit.simulatorOptions().core();
    Int debug = options.smsig_debug;

    if (debug>0) {
        Simulator::dbg() << "Starting AC small-signal transfer function analysis.\n";
    }
    
    // Evaluate resistive and reactive Jacobian, bypass is not allowed
    EvalSetup esReactive { 
        // Inputs, can be set here (we do not rotate)
        .solution = &dcSolution, 
        .states = &dcStates, 
        
        // Evaluation type reported to the model
        .acAnalysis = true, 
        
        // Evaluation 
        .enableLimiting = false, 
        .evaluateResistiveJacobian = true, 
        .evaluateReactiveJacobian = true, 
    };

    LoadSetup lsReactive { 
        // Outputs
        .loadResistiveJacobian = false, 
        .loadReactiveJacobian = true, 
    };

    // Copy OP Jacobian to real part of acMatrix, zero out imaginary part. 
    // Because the real part is taken from OP Jacobian it includes
    // the shunt resistors. 
    auto nnz = dcJacobian.nnz();
    auto Jr = dcJacobian.data();
    auto M = acMatrix.data();
    for(decltype(nnz) i=0; i<nnz; i++) {
        M[i] = Jr[i];
    }
    
    // Evaluate Jacobians 
    // Actually we only need to evaluate the reactive Jacobian 
    // because the resistive part was evaluated by OP analysis
    // We do both here in case OpenVAF has bugs with this corner case :)
    if (!circuit.evalAndLoad(&esReactive, nullptr, nullptr)) {
        // Load error
        setError(ACXFError::EvalAndLoad);
        if (debug>0) {
            Simulator::dbg() << "Error in AC Jacobian evaluation.\n";
        }
        co_yield CoreState::Aborted;
    }

    // Handle Abort, Finish, Stop
    if (esReactive.requests.abort) {
        if (debug>0) {
            Simulator::dbg() << "Abort requested during AC Jacobian evaluation. Exiting.\n";
        }
        co_yield CoreState::Aborted;
    }
    if (esReactive.requests.finish) {
        if (debug>0) {
            Simulator::dbg() << "Finish requested during AC Jacobian evaluation. Exiting.\n";
        }
        co_yield CoreState::Finished;
    }
    if (esReactive.requests.stop) {
        if (debug>0) {
            Simulator::dbg() << "Stop requested during AC Jacobian evaluation. Exiting.\n";
        }
        co_yield CoreState::Stopped;
    }

    // Create sweeper, put it in unique ptr to free it when method returns
    ScalarSweep sweeper;
    if (!sweeper.setup(params, errorStatus)) {
        setError(ACXFError::Sweeper);
        co_yield CoreState::Aborted;
    }
    if (progressReporter) {
        progressReporter->setValueFormat(ProgressReporter::ValueFormat::Scientific, 6);
        progressReporter->setValueDecoration("", "Hz");    
    }
    initProgress(sweeper.count(), 0);
    
    // Frequency sweep
    sweeper.reset();
    bool finished = false;
    frequency = -1.0;
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4);
    bool error = false;
    do {
        // Compute should always succeed
        Value v;
        if (!sweeper.compute(v, errorStatus)) {
            setError(ACXFError::SweepCompute);
            error = true;
            break;
        }

        // The value, however, must be convertible to real
        if (!v.convertInPlace(Value::Type::Real, errorStatus)) {
            setError(ACXFError::BadFrequency);
            error = true;
            break;
        }
        frequency = v.val<Real>();
        double omega = 2*std::numbers::pi*frequency;

        if (debug>0) {
            ss.str(""); ss << frequency;
            Simulator::dbg() << "frequency=" << ss.str() << "\n";
        }

        // Load AC matrix, we must update the imaginary part only
        acMatrix.zero(Component::Imaginary);
        lsReactive.reactiveJacobianFactor = omega;
        if (!circuit.evalAndLoad(nullptr, &lsReactive, nullptr)) {
            // Load error
            setError(ACXFError::EvalAndLoad);
            if (debug>0) {
                Simulator::dbg() << "Error in AC Jacobian load.\n";
            }
            error = true;
            break;
        }
        
        // Check if matrix entries are finite, no need to check RHS 
        // since we loaded it without any computation (i.e. we only used mag and phase)
        if (options.matrixcheck && !acMatrix.isFinite(true, true)) {
            auto nr = UnknownNameResolver(circuit);
            setError(ACXFError::MatrixError);
            if (debug>2) {
                Simulator::dbg() << "A matrix entry is not finite.\n";
            }
            error = true;
            break;
        }

        // Factor
        bool forceFullFactorization = false;        
        if (acMatrix.isFactored()) {
            // Refactor (if possible)
            if (!acMatrix.refactor()) {
                // Failed, try again by fully factoring
                forceFullFactorization = true;
            } 
        }
        if (forceFullFactorization || !acMatrix.isFactored()) {
            // Full factorization
            if (!acMatrix.factor()) {
                // Failed, give up
                setError(ACXFError::MatrixError);
                if (debug>0) {
                    Simulator::dbg() << "LU factorization failed.\n";
                }
                error = true;
                break;
            }
        }
        // Check if matrix is singular
        double rcond;
        if (!acMatrix.rcond(rcond)) {
            setError(ACXFError::MatrixError);
            if (debug>0) {
                Simulator::dbg() << "Condition number estimation failed.\n";
            }
            error = true;
            break;
        }
        if (rcond<1e-15) {
            setError(ACXFError::MatrixError);
            if (debug>0) {
                Simulator::dbg() << "Matrix is close to singular.\n";
            }
            error = true;
            break;
        }

        // Go through all sources listed in sources
        auto nSrc = sources.size();
        for(decltype(nSrc) i=0; i<nSrc; i++) {
            // Prepare RHS
            zero(acSolution); 

            // Get instance
            auto inst = sources[i];
            // No instance, continue to next
            if (!inst) {
                continue;
            }

            if (debug>1) {
                Simulator::dbg() << "Computing response to '" << std::string(inst->name()) << "'.\n";
            }

            // Get excitation equations and response unknowns
            auto [e1, e2] = inst->sourceExcitation(circuit);
            auto [r1, r2] = inst->sourceResponse(circuit);
            
            // Set RHS, load negated residual to get the true value of delta after solve()
            acSolution[e1] += -inst->scaledUnityExcitation();
            acSolution[e2] -= -inst->scaledUnityExcitation();

            // Set RHS bucket
            acSolution[0] = 0.0;

            if (debug>=100) {
                Simulator::dbg() << "Linear system for instance " << inst->name() << "\n";
                acMatrix.dump(Simulator::dbg(), dataWithoutBucket(acSolution)); 
                Simulator::dbg() << "\n";
            }

            // Solve, set bucket to 0.0
            if (!acMatrix.solve(dataWithoutBucket(acSolution))) {
                setError(ACXFError::MatrixError);
                if (debug>2) {
                    Simulator::dbg() << "Failed to solve factored system.\n";
                }
                error = true;
                break;
            }
            acSolution[0] = 0.0;

            if (options.solutioncheck && !acMatrix.isFinite(dataWithoutBucket(acSolution), true, true)) {
                setError(ACXFError::SolutionError);
                if (options.smsig_debug) {
                    Simulator::dbg() << "A solution entry is not finite. Solver failed.\n";
                }
                error = true;
                break;
            }

            // Collect results
            tf[i] = acSolution[up] - acSolution[un];

            // Compute Yin and Zin
            if (inst->model()->device()->isVoltageSource()) {
                // Voltage source excitation
                yin[i] = (acSolution[r1] - acSolution[r2])*inst->responseScalingFactor()/inst->scaledUnityExcitation();
                if (yin[i]!=0.0) {
                    zin[i] = 1.0/yin[i];
                } else {
                    // Infinity
                    zin[i] = 1e20;
                }
            } else {
                // Current source excitation
                zin[i] = (acSolution[r1] - acSolution[r2])*inst->responseScalingFactor()/inst->scaledUnityExcitation();
                if (zin[i]!=0.0) {
                    yin[i] = 1.0/zin[i];
                } else {
                    // Infinity
                    yin[i] = 1e20;
                }
            }
        }
        
        if (error) {
            break;
        }

        // Dump solution
        if (params.writeOutput && outfile) {
            outfile->addPoint();
        }

        finished = sweeper.advance();

        setProgress(sweeper.at(), frequency);
    } while (!finished && !error);

    if (debug>0) {
        Simulator::dbg() << "AC transfer function frequency sweep " << (finished ? "completed" : "exited prematurely") << ".\n";
    }

    if (!finished) {
        errorFreq = frequency;
    }

    // No need to bind resistive Jacobian enatries. 
    // OP analysis will still work fine, even in sweep. 
    // We only changed the bindings of the reactive Jacobian entries. 
    
    if (finished) {
        co_yield CoreState::Finished;
    } else {
        co_yield CoreState::Aborted;
    }
}

bool ACXFCore::run(bool continuePrevious) {
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

bool ACXFCore::formatError(Status& s) const {
    auto nr = UnknownNameResolver(circuit);
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4);
    
    // First, handle AnalysisCore errors
    if (lastError!=Error::OK) {
        AnalysisCore::formatError(s);
        return false;
    }
    
    // Then handle ACXFCore errors
    switch (lastAcTfError) {
        case ACXFError::NotFound:
            s.set(Status::Analysis, std::string("Source '")+std::string(errorInstance)+"' not found.");
            break;
        case ACXFError::NotSource:
            s.set(Status::Analysis, std::string("Instance '")+std::string(errorInstance)+"' is not a source.");
            break;
        case ACXFError::Sweeper:
        case ACXFError::SweepCompute:
            s.set(errorStatus);
            break;
        case ACXFError::EvalAndLoad:
            s.set(Status::Analysis, "Jacobian evaluation failed.");
            break;
        case ACXFError::MatrixError:
            acMatrix.formatError(s, &nr);
            break;
        case ACXFError::SolutionError:
            acMatrix.formatError(s, &nr);
            s.extend("Solution component is not finite.");
            break;
        case ACXFError::OperatingPointError:
            opCore_.formatError(s);
            break;
        case ACXFError::SingularMatrix:
            s.set(Status::Analysis, "Matrix is close to singular.");
            break;
        case ACXFError::BadFrequency:
            s.set(Status::Analysis, "Frequency value cannot be converted to real.");
            break;
        default:
            return true;
    }
    if (errorFreq>=0) {
        ss.str(""); ss << errorFreq;
        s.extend(std::string("Leaving frequency sweep at frequency=")+ss.str()+".");
    } else {
        s.extend("Leaving frequency sweep.");
    }
    return false;
}

void ACXFCore::dump(std::ostream& os) const {
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
