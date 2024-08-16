#include "coreopnr.h"
#include "simulator.h"
#include "common.h"
#include <iomanip>

namespace NAMESPACE {

// Equations and limiting
//
// SPICE formulation
//   G    .. Jacobian
//   x    .. unknowns
//   f(x) .. residual
// 
//   G_i x_{i+1} = - f(x_i) + G_i x_i
// 
// Delta formulation
//   G_i (x_{i+1} - x_i) = - f(x_i)
// 
// SPICE formulation with limiting
//   xl   .. limited unknowns
//   Gl   .. Jacobian at limited unknowns (G(xl))
//   
//   Gl_i x_{i+1} = -f(xl_i) + Gl_i xl_i
//
//   Gl_i x_{i+1} = -f(xl_i) + Gl_i xl_i + Gl_i x_i - Gl_i x_i
//
//   Gl_i (x_{i+1} - x_i) = -f(xl_i) + Gl_i (xl_i - x_i) 
// 
// In the spirit of nrsolver.cpp we rewrite this as
//
//   Gl_i (x_i - x_{i+1}) = f(xl_i) - Gl_i (xl_i - x_i) 
//                          -------   -----------------
//                          residual  RHS linearized residual contribution (storeb by OSDI models)
// 
// load_residual_resist()  ... loads residual (computed at the point of limiting)
// load_limit_rhs_resist() ... loads negated RHS linearized residual contribution
// 
// We first load residuals by calling load_residual_resist(), 
// then subtract RHS limited residual by calling load_limit_rhs_resist(). 
//
// The above equation can also be interpreted as
//   Gl_i (x_i - x_{i+1}) = f(xl_i) + Gl_i (x_i - xl_i)
//                          -----------------------------
//                          residual linearized at xl_i and computed at x_i

// Slots 0 (current) and -1 (future) are used for NR solver
// Slots 1, 2, ... correspond to past values (at t_{k}, t_{k-1}, ...)
// Therefore historyOffset needs to be set to 1
    
OpNRSolver::OpNRSolver(
    Circuit& circuit, KluRealMatrix& jac, 
    VectorRepository<double>& states, VectorRepository<double>& solution, 
    NRSettings& settings, Int forcesSize
) : NRSolver(circuit, jac, states, solution, settings) {
    resizeForces(forcesSize);

    // For constructing the linearized system in NR loop
    esSystem = EvalSetup {
        // Inputs
        .solution = &solution, 
        .states = &states, 

        // Signal this is static and DC analysis
        .staticAnalysis = true, 
        .dcAnalysis = true, 

        // Evaluation 
        .enableLimiting = true, 
        .evaluateResistiveJacobian = true, 
        .evaluateResistiveResidual = true, 
        .evaluateLinearizedResistiveRhsResidual = true, 
        .evaluateOpvars = true, 
    };

    lsSystem = LoadSetup {
        .states = &states, 
        .loadResistiveJacobian = true, 
    };

    csSystem = ConvSetup {
        .solution = &solution, 
        .states = &states
    };
}

void OpNRSolver::requestHighPrecision(bool f) {
    circuit.simulatorInternals().highPrecision = f;
}

bool OpNRSolver::rebuild() {
    // Call parent's rebuild
    if (!NRSolver::rebuild()) {
        // Assume parent has set the error flag
        return false;
    }

    // Allocate space in vetors
    auto n = circuit.unknownCount();
    maxResidualContribution_.resize(n+1);
    dummyStates.resize(circuit.statesCount());

    return true;
}

bool OpNRSolver::initialize(bool continuePrevious) {
    // This method is called once on entering run()
    // This is the right place to set up vectors

    // If bypass is enabled, prepare space for previous device states
    // Need to do this here because the user might sweep nr_bypass, but
    // the minimum requirement for calling rebuild() is that mapping 
    // changes. But nr_bypass does not affect mapping. 
    if (circuit.simulatorOptions().core().nr_bypass) {
        deviceStates.resize(circuit.deviceStatesCount());
    }
    
    // Set vectors for building linear system
    bool computeMaxResidualContribution = settings.residualCheck;

    lsSystem.resistiveResidual = delta.data();
    lsSystem.linearizedResistiveRhsResidual = delta.data();
    lsSystem.maxResistiveResidualContribution = computeMaxResidualContribution ? maxResidualContribution_.data() : nullptr;

    // Set up tolerance reference value for solution
    auto& options = circuit.simulatorOptions().core();
    if (options.relrefsol==SimulatorOptions::relrefPointLocal) {
        settings.globalSolRef = false;
        settings.historicSolRef = false;
    } else if (options.relrefsol==SimulatorOptions::relrefLocal) {
        settings.globalSolRef = false;
        settings.historicSolRef = true;
    } else if (options.relrefsol==SimulatorOptions::relrefPointGlobal) {
        settings.globalSolRef = true;
        settings.historicSolRef = false;
    } else if (options.relrefsol==SimulatorOptions::relrefGlobal) {
        settings.globalSolRef = true;
        settings.historicSolRef = true;
    } else if (options.relrefsol==SimulatorOptions::relrefRelref) {
        if (options.relref == SimulatorOptions::relrefAlllocal) {
            settings.globalSolRef = false;
            settings.historicSolRef = true;
        } else if (options.relref == SimulatorOptions::relrefSigglobal) {
            settings.globalSolRef = true;
            settings.historicSolRef = true;
        } else if (options.relref == SimulatorOptions::relrefAllglobal) {
            settings.globalSolRef = true;
            settings.historicSolRef = true;
        } else {
            lastError = Error::BadSolReference;
            return false;
        }
    } else {
        lastError = Error::BadSolReference;
        return false;
    }

    // Set up tolerance reference value for residual
    if (options.relrefres==SimulatorOptions::relrefPointLocal) {
        settings.globalResRef = false;
        settings.historicResRef = false;
    } else if (options.relrefres==SimulatorOptions::relrefLocal) {
        settings.globalResRef = false;
        settings.historicResRef = true;
    } else if (options.relrefres==SimulatorOptions::relrefPointGlobal) {
        settings.globalResRef = true;
        settings.historicResRef = false;
    } else if (options.relrefres==SimulatorOptions::relrefGlobal) {
        settings.globalResRef = true;
        settings.historicResRef = true;
    } else if (options.relrefres==SimulatorOptions::relrefRelref) {
        if (options.relref == SimulatorOptions::relrefAlllocal) {
            settings.globalResRef = false;
            settings.historicResRef = true;
        } else if (options.relref == SimulatorOptions::relrefSigglobal) {
            settings.globalResRef = false;
            settings.historicResRef = true;
        } else if (options.relref == SimulatorOptions::relrefAllglobal) {
            settings.globalResRef = true;
            settings.historicResRef = true;
        } else {
            lastError = Error::BadResReference;
            return false;
        }
    } else {
        lastError = Error::BadResReference;
        return false;
    }
    
    csSystem.inputDelta = delta.data();

    return true;
}

bool OpNRSolver::postSolve(bool continuePrevious) {
    // Check convergence if nr_bypass is enabled
    if (circuit.simulatorOptions().core().nr_bypass && !skipConvergenceCheck) {
        if (!circuit.converged(csSystem)) {
            lastError = Error::ConvergenceCheck;
            errorIteration = iteration;
            if (settings.debug>2) {
                Simulator::dbg() << "Instance convergence check error.\n";
            }
            return false;
        }
    }

    if (circuit.simulatorOptions().core().nr_bypass) {
        auto& acct = circuit.tables().accounting();
        acct.acctNew.bpiicount += esSystem.bypassableInstances;
        acct.acctNew.bpiiconv += esSystem.bypassableInstances-csSystem.nonConvergedInstances;
        acct.acctNew.bpiibypass += esSystem.bypassedInstances;
        acct.acctNew.bpiibpfailed += esSystem.failedBypassInstances;
        // Simulator::dbg() << "iter " << iteration << ", bypassable: " << esSystem.bypassableInstances 
        //     << ", unconverged " << csSystem.nonConvergedInstances 
        //     << ", bypassed " << esSystem.bypassedInstances 
        //     << ", bypass failed " << esSystem.failedBypassInstances << "\n";
    }
    
    return true;
}

void OpNRSolver::loadShunts(double gshunt, bool loadJacobian) {
    // Now load gshunt if it is greater than 0.0
    // Gshunt current (and its residual contribution) is
    //   gshunt * x
    double* xprev = solution.data();
    auto nUnknowns = circuit.unknownCount();
    for(decltype(nUnknowns) i=1; i<=nUnknowns; i++) {
        auto ptr = diagPtrs[i];
        if (!isFlow[i]) {
            // Jacobian
            if (loadJacobian) {
                *ptr += gshunt;
            }
            // Residual
            delta[i] += gshunt*xprev[i];
        }
    }
}

bool OpNRSolver::evalAndLoadWrapper(EvalSetup& evalSetup, LoadSetup& loadSetup) {
    lastError = Error::OK;
    if (!circuit.evalAndLoad(&evalSetup, &loadSetup, nullptr)) {
        // Load error
        lastError = Error::EvalAndLoad;
        if (settings.debug>2) {
            Simulator::dbg() << "Evaluation error.\n";
        }
        return false;
    }
    
    // Update circuit's flags (Abort, Finish, Stop)
    circuit.updateEvalFlags(evalSetup);

    // Handle abort right now, finish and stop are handled outside NR loop
    if (circuit.checkFlags(Circuit::Flags::Abort)) {
        if (settings.debug>2) {
            Simulator::dbg() << "Abort requested during evaluation.\n";
        }
        return false;
    }

    return true;
}

void OpNRSolver::setNodesetAndIcFlags(bool continuePrevious) {
    auto nsiter = circuit.simulatorOptions().core().op_nsiter;

    // Set nodesetEnabled flag in esSystem
    // Forces slot 0 (continuation nodesets) and 1 (user nodesets) 
    // are used in ordinary OP analysis
    // Nodesets are enabled in iterations 1..op_nsiter 
    // if continuePrevious is false. 
    // They are also enabled if slot 1 (user nodesets) is active. 
    esSystem.nodesetEnabled = (iteration<=nsiter) && (continuePrevious==false);
    
    // Set icEnabled flag in esSystem
    // Slot 2 holds permanent forces for computing initial conditions
    // when OP analysis is invoked from tran core. 
    // Whenever this slot is active transient forces are enables 
    // and we are applying initial conditions. 
    esSystem.icEnabled = forcesEnabled.size()>2 && forcesEnabled[2];
}

std::tuple<bool, bool> OpNRSolver::buildSystem(bool continuePrevious) {
    lastError = Error::OK;

    auto n = circuit.unknownCount();

    // Remove forces originating from nodesets after nsiter iterations
    auto nsiter = circuit.simulatorOptions().core().op_nsiter;
    // Do this only at nsiter+1 (first iteration has index 1)
    if (iteration==nsiter+1) {
        // Continuation nodesets
        enableForces(0, false);
        // User-specified nodesets
        enableForces(1, false);
    }

    // Set nodeset and IC flags
    setNodesetAndIcFlags(continuePrevious); 
    
    // Raw arrays
    double* xprev = solution.data();
    
    // Inject value for limiting debugging purposes
    // if (iteration==1) {
    //     xprev[1]=1.0;
    // }

    // Simulator::dbg() << "Before loading:\n";
    // Simulator::dbg() << "  Residual:\n";
    // circuit.dumpSolution(std::cout, delta.data(), "    ");
    // Simulator::dbg() << "\n";
    
    // Init limits if not in continue mode and iteration is 1
    esSystem.initializeLimiting = !continuePrevious && (iteration==1);
    // Write value to simulatorInternals
    circuit.simulatorInternals().initalizeLimiting = esSystem.initializeLimiting; 

    // Bypass enabled, take it into account at evaluation time
    if (circuit.simulatorOptions().core().nr_bypass) {
        // In continue mode bypass is allowed in first iteration
        // Otherwise we allow it starting with the second one 
        if (continuePrevious || iteration>1) {
            esSystem.allowBypass = true;
        } else {
            esSystem.allowBypass = false;
        }
        // For bypass check
        esSystem.deviceStates = deviceStates.data();
        // For convergence check
        csSystem.deviceStates = deviceStates.data();
    }

    // Evaluate and load
    auto evalSt = evalAndLoadWrapper(esSystem, lsSystem);
    // If bypass was forced for one iteration
    if (circuit.simulatorInternals().forceBypass) {
        // Turn forced bypass off after the system is built
        // (we are allowed to do it for one iteration only)
        // Skip device convergence checks for one iteration
        circuit.simulatorInternals().forceBypass = false;
        skipConvergenceCheck = true;
    } else {
        // This makes sure that the device convergence check is 
        // skipped only if bypass was forced. 
        skipConvergenceCheck = false;
    }
    if (!evalSt) {
        lastError = Error::EvalAndLoad;
        errorIteration = iteration;
        return std::make_tuple(false, esSystem.limitingApplied);
    }
    delta[0] = 0.0;

    // Simulator::dbg() << "After loading:\n";
    // Simulator::dbg() << "  Residual:\n";
    // circuit.dumpSolution(std::cout, delta.data(), "    ");
    // Simulator::dbg() << "\n";
    
    // Now load gshunt if it is greater than 0.0
    auto gshunt = circuit.simulatorInternals().gshunt;
    if (gshunt>0) {
        loadShunts(gshunt);
    }

    // Prevent convergence if limiting was applied
    return std::make_tuple(true, esSystem.limitingApplied); 
}

}
