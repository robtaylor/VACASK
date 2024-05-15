#include <iomanip>
#include <cmath>
#include <filesystem>
#include "coreactf.h"
#include "simulator.h"
#include "answeep.h"
#include "context.h"
#include "common.h"
#include <numbers>

namespace NAMESPACE {

// Default parameters
AcTfParameters::AcTfParameters() {
    opParams.writeOutput = 0;
}

template<> int Introspection<AcTfParameters>::setup() {
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
instantiateIntrospection(AcTfParameters);


AcTfCore::AcTfCore(
    Analysis& analysis, AcTfParameters& params, OperatingPointCore& opCore, std::unordered_map<Id,size_t>& sourceIndex, 
    Circuit& circuit, 
    KluRealMatrix& dcJacobian, VectorRepository<double>& dcSolution, VectorRepository<double>& dcStates, 
    KluComplexMatrix& acMatrix, Vector<Complex>& acSolution, 
    std::vector<Instance*>& sources, Vector<Complex>& tf, Vector<Complex>& yin, Vector<Complex>& zin
) : AnalysisCore(analysis, circuit), params(params), outfile(nullptr), opCore_(opCore), sourceIndex(sourceIndex), 
    dcSolution(dcSolution), dcStates(dcStates), dcJacobian(dcJacobian), 
    acMatrix(acMatrix), acSolution(acSolution), sources(sources), tf(tf), yin(yin), zin(zin) {
    
    // Set analysis type for the initial operating point analysis
    auto& elsSystem = opCore_.solver().evalSetupSystem();
    auto& elsResidual = opCore_.solver().evalSetupResidual();

    elsSystem.staticAnalysis = true;
    elsSystem.dcAnalysis = false;
    elsSystem.acAnalysis = true;

    elsResidual.staticAnalysis = true;
    elsResidual.dcAnalysis = false;
    elsResidual.acAnalysis = true;
}

AcTfCore::~AcTfCore() {
    delete outfile;
}

// Implement this in every derived class so that calls to 
// resolveOutputDescriptor() will be inlined. 
bool AcTfCore::resolveOutputDescriptors(bool strict, Status &s) {
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
                        s.set(Status::NotFound, std::string("Source '")+std::string(name)+"' not found.");
                        ok = false;
                        break;
                    }
                }
                // Instance found, but is not a source... this is always an error
                if (inst && !inst->model()->device()->isSource()) {
                    s.set(Status::NotFound, std::string("Instance '")+std::string(name)+"' is not a source.");
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
            outputSources.emplace_back(&(circuit.simulatorInternals().frequency));
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
            ok = analysis.resolveOutputDescriptor(*it, outputSources, strict, s);
        }
        if (!ok) {
            break;
        }
    }
    return ok;
}

bool AcTfCore::addCoreOutputDescriptors(Status& s) {
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    return addOutputDescriptor(OutputDescriptor(OutdFrequency, "frequency"));
}

bool AcTfCore::addDefaultOutputDescriptors(Status& s) {
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    if (savesCount==0) {
        return addAllTfZin(PTSave(Loc::bad, "default", Id(), Id()), false, sourceIndex, s);
    }
    return true;
}

bool AcTfCore::initializeOutputs(Id name, Status& s) {
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

bool AcTfCore::finalizeOutputs(Status &s) {
    outfile->epilogue();
    delete outfile;
    outfile = nullptr;
    return true;
}

bool AcTfCore::deleteOutputs(Id name, Status &s) {
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
    
bool AcTfCore::rebuild(Status& s) {
    // AC analysis matrix
    if (!acMatrix.rebuild(circuit.sparsityMap(), circuit.unknownCount())) {
        acMatrix.formatError(s);
        return false;
    }
    
    // Resistive Jacobian entries remain bound to OP Jacobian, 
    // reactive parts will be bound to imaginary entries of acMatrix
    if (!circuit.bind(nullptr, nullptr, Component::RealPart, nullptr, &acMatrix, Component::ImagPart, s)) {
        return false;
    }
    
    return true;
}

// System of equations is 
//   (G(x) + i C(x)) dx = dJ
bool AcTfCore::run(bool continuePrevious, Status& s) {
    auto n = circuit.unknownCount(); 
    // Make sure structures are large enough
    acSolution.resize(n+1);
    tf.resize(sources.size());
    yin.resize(sources.size());
    zin.resize(sources.size());

    // Get output unknowns
    auto [ok, up, un] = getOutput(params.out, s);
    if (!ok) {
        return false;
    }
    
    // Compute operating point
    auto opOk = opCore_.run(continuePrevious, s);
    if (!opOk) {
        return false;
    }

    auto& options = circuit.simulatorOptions().core();
    Int debug = options.smsig_debug;

    if (debug>0) {
        Simulator::dbg() << "Starting AC small-signal transfer function analysis.\n";
    }
    
    // Evaluate resistive and reactive Jacobian
    EvalAndLoadSetup elsEval { 
        // Inputs, can be set here (we do not rotate)
        .solution = &dcSolution, 
        .states = &dcStates, 
        
        // Evaluation 
        .enableLimiting = false, 
        .evaluateResistiveJacobian = true, 
        .evaluateReactiveJacobian = true, 
        
        // Evaluation type reported to the model
        .acAnalysis = true, 
        
        // Outputs - none
    };

    EvalAndLoadSetup elsLoad { 
        // Needs no inputs because we are not evaluating

        // Skip evaluation
        .skipEvaluation = true, 

        // Evaluation type reported to the model
        .acAnalysis = true, 
        
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
    if (!circuit.evalAndLoad(elsEval, nullptr, s)) {
        // Load error
        s.extend("Exiting AC small-signal transfer function analysis.");
        if (debug>0) {
            Simulator::dbg() << "Error in AC Jacobian evaluation.\n";
        }
        return false;
    }

    // Handle Abort, Finish, Stop
    circuit.updateEvalFlags(elsEval);
    if (circuit.checkFlags(Circuit::Flags::Abort)) {
        if (debug>0) {
            Simulator::dbg() << "Abort requested during AC Jacobian evaluation. Exiting.\n";
        }
        return false;
    }
    if (circuit.checkFlags(Circuit::Flags::Finish)) {
        if (debug>0) {
            Simulator::dbg() << "Finish requested during AC Jacobian evaluation. Exiting.\n";
        }
        return true;
    }
    if (circuit.checkFlags(Circuit::Flags::Stop)) {
        if (debug>0) {
            Simulator::dbg() << "Stop requested during AC Jacobian evaluation. Exiting.\n";
        }
        return true;
    }

    // Create sweeper, put it in unique ptr to free it when method returns
    auto sweeper = ScalarSweep::create(params, s);
    if (!sweeper) {
        return false;
    }
    std::unique_ptr<ScalarSweep> sweeperPtr(sweeper);

    // Frequency sweep
    sweeper->reset();
    bool finished = false;
    double freq = -1.0;
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4);
    bool error = false;
    do {
        // Compute should always succeed
        Value v;
        sweeper->compute(v, s);

        // The value, however, must be convertible to real
        if (!v.convertInPlace(Value::Type::Real, s)) {
            s.extend("Frequency value is not real.");
            if (debug>0) {
                Simulator::dbg() << "Frequency value is not real.\n";
            }
            error = true;
            break;
        }
        freq = v.val<Real>();
        double omega = 2*std::numbers::pi*freq;
        circuit.simulatorInternals().frequency = freq;

        if (debug>0) {
            ss.str(""); ss << freq;
            Simulator::dbg() << "frequency=" << ss.str() << "\n";
        }

        // Load AC matrix, we must update the imaginary part only
        acMatrix.zero(Component::ImagPart);
        elsLoad.reactiveJacobianFactor = omega;
        if (!circuit.evalAndLoad(elsLoad, nullptr, s)) {
            // Load error
            s.extend("Exiting AC small-signal transfer function analysis.");
            if (debug>0) {
                Simulator::dbg() << "Error in AC Jacobian load.\n";
            }
            error = true;
            break;
        }
        if (circuit.checkFlags(Circuit::Flags::Abort)) {
            if (debug>0) {
                Simulator::dbg() << "Abort requested during AC transfer function frequency sweep.\n";
            }
            return false;
        }
        
        // Check if matrix entries are finite, no need to check RHS 
        // since we loaded it without any computation (i.e. we only used mag and phase)
        if (!acMatrix.isFinite(options.infcheck, options.nancheck)) {
            auto nr = UnknownNameResolver(circuit);
            acMatrix.formatError(s, &nr);
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
                auto nr = UnknownNameResolver(circuit);
                acMatrix.formatError(s, &nr);
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
            acMatrix.formatError(s);
            if (debug>0) {
                Simulator::dbg() << "Condition number estimation failed.\n";
            }
            error = true;
            break;
        }
        if (rcond<1e-15) {
            if (debug>0) {
                Simulator::dbg() << "Matrix is close to singular.\n";
            }
            s.set(Status::LinearSolver, "Matrix is close to singular.");
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
                acMatrix.formatError(s);
                if (debug>2) {
                    Simulator::dbg() << "Failed to solve factored system.\n";
                }
                error = true;
                break;
            }
            acSolution[0] = 0.0;

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

        // Handle Finish and Stop
        if (circuit.checkFlags(Circuit::Flags::Finish)) {
            if (debug>0) {
                Simulator::dbg() << "Finish requested during AC transfer function frequency sweep.\n";
            }
            break;
        }
        if (circuit.checkFlags(Circuit::Flags::Stop)) {
            if (debug>0) {
                Simulator::dbg() << "Stop requested during AC transfer function frequency sweep.\n";
            }
            break;
        }
        
        finished = sweeper->advance();
    } while (!finished && !error);

    if (debug>0) {
        Simulator::dbg() << "AC transfer function frequency sweep " << (finished ? "completed" : "exited prematurely") << ".\n";
    }

    if (!finished) {
        if (freq>=0) {
            ss.str(""); ss << freq;
            s.set(Status::NotConverged, std::string("Leaving AC transfer function frequency sweep at frequency=")+ss.str()+".");
        } else {
            s.set(Status::NotConverged, "Leaving AC transfer function frequency sweep.");
        }
    }

    // No need to bind resistive Jacobian enatries. 
    // OP analysis will still work fine, even in sweep. 
    // We only changed the bindings of the reactive Jacobian entries. 
    
    return finished;
}


void AcTfCore::dump(std::ostream& os) const {
    AnalysisCore::dump(os);
    os << "  Results" << std::endl;
    auto nSrc = sources.size();
    for(decltype(nSrc) i=0; i<nSrc; i++) {
        auto inst = sources[i];
        if (!inst) {
            continue;
        }
        os << "    tf(" << inst->name() << ") " << tf[i] << std::endl;
        os << "    zin(" << inst->name() << ") " << zin[i] << std::endl;
        os << "    yin(" << inst->name() << ") " << yin[i] << std::endl;
    }
}

}
