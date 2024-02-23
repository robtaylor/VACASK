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
    auto& elsResidual = opCore_.solver().evalSetupResidual();

    elsSystem.staticAnalysis = true;
    elsSystem.dcAnalysis = false;
    elsSystem.acAnalysis = true;

    elsResidual.staticAnalysis = true;
    elsResidual.dcAnalysis = false;
    elsResidual.acAnalysis = true;
}

AcCore::~AcCore() {
    delete outfile;
}

// Implement this in every derived class so that calls to 
// resolveOutputDescriptor() will be inlined. 
bool AcCore::resolveOutputDescriptors(bool strict, Status &s) {
    // Clear output sources
    outputSources.clear();
    // Resolve output descriptors
    bool ok = true; 
    for (auto it = outputDescriptors.cbegin(); it != outputDescriptors.cend(); ++it) {
        Node *node;
        Instance *inst;
        switch (it->type) {
        case OutdSolComponent:
            ok = addComplexVarOutputSource(strict, it->id, acSolution, s);
            break;
        case OutdFrequency:
            outputSources.emplace_back(&(circuit.simulatorInternals().frequency));
            break;
        default:
            // Delegate to parent
            ok = analysis.resolveOutputDescriptor(*it, outputSources, strict, s);
            break;
        }
        if (!ok) {
            break;
        }
    }
    return ok;
}

bool AcCore::addCoreOutputDescriptors(Status& s) {
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    return addOutputDescriptor(OutputDescriptor(OutdFrequency, "frequency"));
}

bool AcCore::addDefaultOutputDescriptors(Status& s) {
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    if (savesCount==0) {
        return addAllUnknowns(PTSave(Loc::bad, "default", Id(), Id()), false, s);
    }
    return true;
}

bool AcCore::initializeOutputs(Id name, Status& s) {
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

bool AcCore::finalizeOutputs(Status &s) {
    if (outfile) {
        outfile->epilogue();
        delete outfile;
        outfile = nullptr;
    }
    return true;
}

bool AcCore::deleteOutputs(Id name, Status &s) {
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
    // AC analysis matrix
    if (!acMatrix.rebuild(circuit.sparsityMap(), circuit.unknownCount(), s)) {
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
bool AcCore::run(bool continuePrevious, Status& s) {
    auto n = circuit.unknownCount(); 
    // Make sure structures are large enough
    acSolution.resize(n+1);
    
    // Compute operating point
    auto opOk = opCore_.run(continuePrevious, s);
    if (!opOk) {
        return false;
    }

    auto& options = circuit.simulatorOptions().core();
    Int debug = options.smsig_debug;

    if (debug>0) {
        Simulator::dbg() << "Starting AC small-signal analysis.\n";
    }
    
    // Evaluate resistive and reactive Jacobian
    EvalAndLoadSetup elsEval { 
        // Inputs
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
    if (!circuit.evalAndLoad(elsEval, nullptr, s)) {
        // Load error
        s.extend("Exiting AC small-signal analysis.");
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

        // Zero out imaginary part, and RHS. 
        // Because the real part is taken from OP Jacobian it includes
        // the shunt resistors. 
        // Load imaginary part and AC residual. 
        acMatrix.zero(Component::ImagPart);
        zero(acSolution);
        elsLoad.reactiveJacobianFactor = omega;
        if (!circuit.evalAndLoad(elsLoad, nullptr, s)) {
            // Load error
            s.extend("Exiting AC small-signal analysis.");
            if (debug>0) {
                Simulator::dbg() << "Error in AC Jacobian load.\n";
            }
            error = true;
            break;
        }
        if (circuit.checkFlags(Circuit::Flags::Abort)) {
            if (debug>0) {
                Simulator::dbg() << "Abort requested during AC frequency sweep.\n";
            }
            return false;
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
        if (!acMatrix.isFinite(options.infcheck, options.nancheck, s)) {
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
            if (!acMatrix.factor(s)) {
                // Failed, give up
                if (debug>0) {
                    Simulator::dbg() << "LU factorization failed.\n";
                }
                error = true;
                break;
            }
        }
        // Check if matrix is singular
        double rcond;
        if (!acMatrix.rcond(rcond, s)) {
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
            s.set(Status::MatrixLU, "Matrix is close to singular.");
            error = true;
            break;
        }

        // Solve, set bucket to 0.0
        if (!acMatrix.solve(dataWithoutBucket(acSolution), s)) {
            if (debug>2) {
                Simulator::dbg() << "Failed to solve factored system.\n";
            }
            error = true;
            break;
        }
        acSolution[0] = 0.0;
        
        // Dump solution point
        if (params.writeOutput && outfile) {
            outfile->addPoint();
        }

        // Handle Finish and Stop
        if (circuit.checkFlags(Circuit::Flags::Finish)) {
            if (debug>0) {
                Simulator::dbg() << "Finish requested during AC frequency sweep.\n";
            }
            break;
        }
        if (circuit.checkFlags(Circuit::Flags::Stop)) {
            if (debug>0) {
                Simulator::dbg() << "Stop requested during AC frequency sweep.\n";
            }
            break;
        }
        
        finished = sweeper->advance();
    } while (!finished && !error);
    
    if (debug>0) {
        Simulator::dbg() << "AC frequency sweep " << (finished ? "completed" : "exited prematurely") << ".\n";
    }

    if (!finished) {
        if (freq>=0) {
            ss.str(""); ss << freq;
            s.set(Status::NotConverged, std::string("Leaving AC frequency sweep at frequency=")+ss.str()+".");
        } else {
            s.set(Status::NotConverged, "Leaving AC frequency sweep.");
        }
    }

    // No need to bind resistive Jacobian enatries. 
    // OP analysis will still work fine, even in sweep. 
    // We only changed the bindings of the reactive Jacobian entries. 
    
    return finished;
}


void AcCore::dump(std::ostream& os) const {
    AnalysisCore::dump(os);
    os << "  Results" << std::endl;
    circuit.dumpSolution(os, acSolution.data(), "    ");
}

}
