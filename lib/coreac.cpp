#include <iomanip>
#include <cmath>
#include <filesystem>
#include "coreac.h"
#include "simulator.h"
#include "answeep.h"
#include "context.h"
#include "common.h"
#include <numbers>

namespace NAMESPACE {

// Default parameters
AcParameters::AcParameters() {
    opParams.writeOutput = 0;
}

template<> int Introspection<AcParameters>::setup() {
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
instantiateIntrospection(AcParameters);


AcCore::AcCore(
    Analysis& analysis, AcParameters& params, OperatingPointCore& opCore, Circuit& circuit, 
    KluRealMatrix& dcJacobian, VectorRepository<double>& dcSolution, VectorRepository<double>& dcStates, 
    KluComplexMatrix& acMatrix, Vector<Complex>& acSolution
) : AnalysisCore(analysis, circuit), params(params), outfile(nullptr), opCore_(opCore), 
    dcSolution(dcSolution), dcStates(dcStates), dcJacobian(dcJacobian), 
    acMatrix(acMatrix), acSolution(acSolution) {
    
    // Set analysis type for the initial operating point analysis
    auto& elsSystem = opCore_.solver().evalSetupSystem();
    elsSystem.staticAnalysis = true;
    elsSystem.dcAnalysis = false;
    elsSystem.acAnalysis = true;
}

AcCore::~AcCore() {
    delete outfile;
}

// Implement this in every derived class so that calls to 
// resolveOutputDescriptor() will be inlined. 
bool AcCore::resolveOutputDescriptors(bool strict) {
    // Clear output sources
    outputSources.clear();
    // Resolve output descriptors
    bool ok = true; 
    for (auto it = outputDescriptors.cbegin(); it != outputDescriptors.cend(); ++it) {
        Node *node;
        Instance *inst;
        switch (it->type) {
        case OutdSolComponent:
            ok = addComplexVarOutputSource(strict, it->id, acSolution);
            break;
        case OutdFrequency:
            outputSources.emplace_back(&(circuit.simulatorInternals().frequency));
            break;
        default:
            // Delegate to parent
            ok = analysis.resolveOutputDescriptor(*it, outputSources, strict);
            break;
        }
        if (!ok) {
            break;
        }
    }
    return ok;
}

bool AcCore::addCoreOutputDescriptors() {
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

bool AcCore::addDefaultOutputDescriptors() {
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    if (savesCount==0) {
        return addAllUnknowns(PTSave(Loc::bad, "default", Id(), Id()));
    }
    return true;
}

bool AcCore::initializeOutputs(Id name) {
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    // Create output file if not created yet
    if (!outfile) {
        outfile = new OutputRawfile(
            name, outputDescriptors, outputSources,
            (circuit.simulatorOptions().core().rawfile==SimulatorOptions::rawfileBinary ? OutputRawfile::Flags::Binary : OutputRawfile::Flags::None) |
                OutputRawfile::Flags::Padded | OutputRawfile::Flags::Complex);
        outfile->setTitle(circuit.title());
        outfile->setPlotname("AC Small Signal Analysis");
    }
    outfile->prologue();

    return true;
}

bool AcCore::finalizeOutputs() {
    if (outfile) {
        outfile->epilogue();
        delete outfile;
        outfile = nullptr;
    }
    return true;
}

bool AcCore::deleteOutputs(Id name) {
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
    
bool AcCore::rebuild(Status& s) {
    clearError();
    // AC analysis matrix
    if (!acMatrix.rebuild(circuit.sparsityMap(), circuit.unknownCount())) {
        setError(AcError::MatrixError);
        return false;
    }
    
    // Resistive Jacobian entries remain bound to OP Jacobian, 
    // reactive parts will be bound to imaginary entries of acMatrix
    if (!circuit.bind(nullptr, Component::Real, &acMatrix, Component::Imaginary, s)) {
        return false;
    }
    
    return true;
}

// System of equations is 
//   (G(x) + i C(x)) dx = dJ
CoreCoroutine AcCore::coroutine(bool continuePrevious) {
    acMatrix.setAccounting(circuit.tables().accounting());
    
    clearError();
    
    auto n = circuit.unknownCount(); 
    // Make sure structures are large enough
    acSolution.resize(n+1);
    
    // Compute operating point
    auto opOk = opCore_.run(continuePrevious);
    if (!opOk) {
        setError(AcError::OpError);
        co_yield CoreState::Aborted;
    }

    auto& options = circuit.simulatorOptions().core();
    Int debug = options.smsig_debug;

    if (debug>0) {
        Simulator::dbg() << "Starting AC small-signal analysis.\n";
    }
    
    // Evaluate resistive and reactive Jacobian, bypass is not allowed
    EvalSetup esReactive { 
        // Inputs
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
        .acResidual = acSolution.data()
    };

    // Copy OP Jacobian to real part of acMatrix, zero out imaginary part
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
        setError(AcError::EvalAndLoad);
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

    // Create sweeper
    ScalarSweep sweeper;
    if (!sweeper.setup(params, errorStatus)) {
        setError(AcError::Sweeper);
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
    double freq = -1.0;
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4);
    bool error = false;
    do {
        // Compute should always succeed
        Value v;
        if (!sweeper.compute(v, errorStatus)) {
            setError(AcError::SweepCompute);
            error = true;
            break;
        }

        // The value, however, must be convertible to real
        if (!v.convertInPlace(Value::Type::Real, errorStatus)) {
            setError(AcError::BadFrequency);
            if (debug>0) {
                Simulator::dbg() << "Frequency value cannot be converted to real.\n";
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

        // Zero out imaginary part, and RHS. 
        // Because the real part is taken from OP Jacobian it includes
        // the shunt resistors. 
        // Load imaginary part and AC residual. 
        acMatrix.zero(Component::Imaginary);
        zero(acSolution);
        lsReactive.reactiveJacobianFactor = omega;
        if (!circuit.evalAndLoad(nullptr, &lsReactive, nullptr)) {
            // Load error
            setError(AcError::EvalAndLoad);
            if (debug>0) {
                Simulator::dbg() << "Error in AC Jacobian load.\n";
            }
            error = true;
            break;
        }
        acSolution[0] = 0.0;

        // Change sign of residual because it is on the RHS 
        // and we need the small signal response with the correct sign
        for(decltype(n) i=0; i<=n; i++) {
            acSolution[i] = -acSolution[i];
        }
        
        if (debug>=100) {
            Simulator::dbg() << "Linear system at frequency " << freq << "\n";
            acMatrix.dump(Simulator::dbg(), dataWithoutBucket(acSolution)); 
            Simulator::dbg() << "\n";
        }

        // Check if matrix entries are finite, no need to check RHS 
        // since we loaded it without any computation (i.e. we only used mag and phase)
        if (options.matrixcheck && !acMatrix.isFinite(true, true)) {
            setError(AcError::MatrixError);
            if (debug>0) {
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
                setError(AcError::MatrixError);
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
            setError(AcError::MatrixError);
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
            setError(AcError::SingularMatrix);
            error = true;
            break;
        }

        // Solve, set bucket to 0.0
        if (!acMatrix.solve(dataWithoutBucket(acSolution))) {
            setError(AcError::MatrixError);
            if (debug>2) {
                Simulator::dbg() << "Failed to solve factored system.\n";
            }
            error = true;
            break;
        }
        acSolution[0] = 0.0;

        if (options.solutioncheck && !acMatrix.isFinite(dataWithoutBucket(acSolution), true, true)) {
            setError(AcError::SolutionError);
            if (options.smsig_debug) {
                Simulator::dbg() << "A solution entry is not finite. Solver failed.\n";
            }
            error = true;
            break;
        }
        
        // Dump solution point
        if (params.writeOutput && outfile) {
            outfile->addPoint();
        }
        
        finished = sweeper.advance();
        
        setProgress(sweeper.at(), freq);
    } while (!finished && !error);
    
    if (debug>0) {
        Simulator::dbg() << "AC frequency sweep " << (finished ? "completed" : "exited prematurely") << ".\n";
    }

    if (!finished) {
        errorFreq = freq;
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

bool AcCore::run(bool continuePrevious) {
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

bool AcCore::formatError(Status& s) const {
    auto nr = UnknownNameResolver(circuit);
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4);
    
    // First, handle AnalysisCore errors
    if (lastError!=Error::OK) {
        AnalysisCore::formatError(s);
        return false;
    }
    
    // Then handle AcCore errors
    switch (lastAcError) {
        case AcError::Sweeper:
        case AcError::SweepCompute:
            s.set(errorStatus);
            break;
        case AcError::EvalAndLoad:
            s.set(Status::Analysis, "Jacobian evaluation failed.");
            break;
        case AcError::MatrixError:
            acMatrix.formatError(s, &nr);
            break;
        case AcError::SolutionError:
            acMatrix.formatError(s, &nr);
            s.extend("Solution component is not finite.");
            break;
        case AcError::OpError:
            opCore_.formatError(s);
            break;
        case AcError::SingularMatrix:
            s.set(Status::Analysis, "Matrix is close to singular.");
            break;
        case AcError::BadFrequency:
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

void AcCore::dump(std::ostream& os) const {
    AnalysisCore::dump(os);
    os << "  Results\n";
    circuit.dumpSolution(os, acSolution.data(), "    ");
}

}
