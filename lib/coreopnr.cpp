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
) : circuit(circuit), states(states), 
    NRSolver(circuit.tables().accounting(), jac, solution, settings) {
    resizeForces(forcesSize);

    // For constructing the linearized system in NR loop
    evalSetup_ = EvalSetup {
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

    loadSetup_ = LoadSetup {
        .states = &states, 
        .loadResistiveJacobian = true, 
    };

    convSetup_ = ConvSetup {
        .solution = &solution, 
        .states = &states
    };
}

bool OpNRSolver::rebuild() {
    // Call parent's rebuild
    if (!NRSolver::rebuild()) {
        // Assume parent has set the error flag
        return false;
    }

    // Allocate space in vetors
    auto n = circuit.unknownCount();
    dummyStates.resize(circuit.statesCount());
    historicMaxSolution_.resize(n+1);
    historicMaxResidualContribution_.resize(n+1);
    maxResidualContribution_.resize(n+1);
    globalMaxResidualContribution_.resize(n+1);
    globalMaxSolution_.resize(n+1);
    resetMaxima();

    // Build flow node flags
    isFlow.resize(n+1);
    for(decltype(n) i=1; i<=n; i++) {
        auto* node = circuit.reprNode(i);
        isFlow[i] = node->maskedFlags(Node::Flags::NodeTypeMask)==Node::Flags::PotentialNode;
    }
    
    return true;
}

bool OpNRSolver::initialize(bool continuePrevious) {
    // This method is called once on entering run()
    // This is the right place to set up vectors

    // Clear flags
    clearFlags();

    // If not in continue mode set current states to 0
    if (!continuePrevious) {
        // Zero states
        states.zero();
    }

    // If bypass is enabled, prepare space for previous device states
    // Need to do this here because the user might sweep nr_bypass, but
    // the minimum requirement for calling rebuild() is that mapping 
    // changes. But nr_bypass does not affect mapping. 
    if (circuit.simulatorOptions().core().nr_bypass) {
        deviceStates.resize(circuit.deviceStatesCount());
    }
    
    // Set vectors for building linear system
    bool computeMaxResidualContribution = settings.residualCheck;

    loadSetup_.resistiveResidual = delta.data();
    loadSetup_.linearizedResistiveRhsResidual = delta.data();
    loadSetup_.maxResistiveResidualContribution = computeMaxResidualContribution ? maxResidualContribution_.data() : nullptr;

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
    
    convSetup_.inputDelta = delta.data();

    return true;
}

bool OpNRSolver::preIteration(bool continuePrevious) {
    // Clear maximal residual contribution
    zero(maxResidualContribution_);
    // Clear future states
    states.zeroFuture();
    // Pass iteration number to Verilog-A models
    circuit.simulatorInternals().iteration = iteration;
    return true;    
}

bool OpNRSolver::postSolve(bool continuePrevious) {
    // Check convergence if nr_bypass is enabled and convergence check is not to be skipped. 
    // The test and state storing is skipped if evaluation bypass was forced. 
    if (circuit.simulatorOptions().core().nr_bypass && !skipConvergenceCheck) {
        // When high precision is requested we only store instance state 
        // and assume instance is not converged. 
        convSetup_.storeStateOnly = highPrecision;

        if (!circuit.converged(convSetup_)) {
            lastError = Error::ConvergenceCheck;
            errorIteration = iteration;
            if (settings.debug>2) {
                Simulator::dbg() << "Instance convergence check error.\n";
            }
            return false;
        }
    }

    auto& acct = circuit.tables().accounting();
    acct.acctNew.bpinst += evalSetup_.bypassableInstances;
    acct.acctNew.bpopport += evalSetup_.bypassOpportunuties;
    acct.acctNew.bpbypassed += evalSetup_.bypassedInstances;
    acct.acctNew.bpiiconvcheck += convSetup_.instancesConvergenceChecks;
    acct.acctNew.bpiiconverged += convSetup_.convergedInstances;
    
    return true;
}

bool OpNRSolver::postConvergenceCheck(bool continuePrevious) {
    // If algorithm converged we are going to exit next and states must be 
    // rotated because the new state belongs to the current solution
    states.rotate();

    // Print debug information on convergence
    if (settings.debug) {
        std::stringstream ss;
        ss << std::scientific << std::setprecision(2);
        Simulator::dbg() << "Iteration " << std::to_string(iteration) << (preventedConvergence ? ", convergence not allowed" : "");
        if (!preventedConvergence) {
            Simulator::dbg() << (iterationConverged ? ", converged" : "");
            if (settings.residualCheck) {
                ss.str(""); ss << maxResidual;
                Simulator::dbg() << ", worst residual=" << ss.str() << " @ " << (maxResidualNode ? maxResidualNode->name() : "(unknown)");
            }
            if (iteration>1) {
                ss.str(""); ss << maxDelta;
                Simulator::dbg() << ", worst delta=" << ss.str() << " @ " << (maxDeltaNode ? maxDeltaNode->name() : "(unknown)");
            }
        }
        Simulator::dbg() << "\n";
    }
    return true;
}

bool OpNRSolver::postIteration(bool continuePrevious) {
    // Update pointMaxSolution_
    pointMaxSolution_ = 0;
    auto xnew = solution.futureData();
    auto n = solution.length()-1;
    for(decltype(n) i=1; i<=n; i++) {
        double c = std::fabs(xnew[i]);
        if (c>pointMaxSolution_) {
            pointMaxSolution_ = c;
        }
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
    evalSetup.requestHighPrecision = highPrecision;
    if (!circuit.evalAndLoad(&evalSetup, &loadSetup, nullptr)) {
        // Load error
        lastError = Error::EvalAndLoad;
        if (settings.debug>2) {
            Simulator::dbg() << "Evaluation error.\n";
        }
        return false;
    }

    // Store Abort, Finish, and Stop flag
    if (evalSetup.requests.abort) {
        setFlags(Flags::Abort);
    }
    if (evalSetup.requests.finish) {
        setFlags(Flags::Finish);
    }
    if (evalSetup.requests.stop) {
        setFlags(Flags::Stop);
    }
    
    // Handle abort right now, finish and stop are handled outside NR loop
    if (checkFlags(Flags::Abort)) {
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
    evalSetup_.nodesetEnabled = (iteration<=nsiter) && (continuePrevious==false);
    
    // Set icEnabled flag in esSystem
    // Slot 2 holds permanent forces for computing initial conditions
    // when OP analysis is invoked from tran core. 
    // Whenever this slot is active transient forces are enables 
    // and we are applying initial conditions. 
    evalSetup_.icEnabled = forcesEnabled.size()>2 && forcesEnabled[2];
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
    evalSetup_.initializeLimiting = !continuePrevious && (iteration==1);
    // Write value to simulatorInternals
    circuit.simulatorInternals().initalizeLimiting = evalSetup_.initializeLimiting; 

    // Bypass enabled, take it into account at evaluation time
    if (circuit.simulatorOptions().core().nr_bypass) {
        // In continue mode bypass is allowed in first iteration
        // Otherwise we allow it starting with the second one 
        if (continuePrevious || iteration>1) {
            evalSetup_.allowBypass = true;
        } else {
            evalSetup_.allowBypass = false;
        }
        // For bypass check
        evalSetup_.deviceStates = deviceStates.data();
        // For convergence check
        convSetup_.deviceStates = deviceStates.data();
    }

    
    // Force instance evaluation bypass if requested
    evalSetup_.forceBypass = circuit.simulatorInternals().requestForcedBypass;
    
    // Evaluate and load
    auto evalSt = evalAndLoadWrapper(evalSetup_, loadSetup_);
    // If bypass forcing was requested clear that request. 
    // It is allowed for one iteration only. 
    if (circuit.simulatorInternals().requestForcedBypass) {
        circuit.simulatorInternals().requestForcedBypass = false;
        // Skip device convergence checks for one iteration
        skipConvergenceCheck = true;
    } else {
        // This makes sure that the device convergence check is 
        // skipped only if bypass was forced. 
        skipConvergenceCheck = false;
    }
    if (!evalSt) {
        lastError = Error::EvalAndLoad;
        errorIteration = iteration;
        return std::make_tuple(false, evalSetup_.limitingApplied);
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
    return std::make_tuple(true, evalSetup_.limitingApplied); 
}

std::tuple<bool, bool> OpNRSolver::checkResidual() {
    // Compute norms only in debug mode
    bool computeNorms = settings.debug;

    // In residual we have the residual at previous solution
    // We are going to check that residual
    
    // Number of unknowns (vector length includes a bucket at index 0)
    auto n = circuit.unknownCount();

    // Results
    double maxResidual = 0.0;
    double maxNormResidual = 0.0;
    double l2normResidual2 = 0.0;
    Node* maxResidualNode = nullptr;
    
    // Assume residual is OK
    bool residualOk = true;
    
    // Get point maximum
    pointMaxResidualContribution_ = 0;
    for(decltype(n) i=1; i<=n; i++) {
        double c = std::fabs(maxResidualContribution_[i]);
        if (c>pointMaxResidualContribution_) {
            pointMaxResidualContribution_ = c;
        }
    }
    
    // Go through all variables (except ground)
    for(decltype(n) i=1; i<=n; i++) {
        // Representative node, associated flow nature index
        auto rn = circuit.reprNode(i);
        bool isPotential = ((rn->flags() & Node::Flags::PotentialNode) == Node::Flags::PotentialNode); 
        size_t ndx = isPotential ? 1 : 0;

        // Compute tolerance reference
        double tolref = std::fabs(maxResidualContribution_[i]);

        // Account for global and historic references
        if (settings.historicResRef) {
            if (settings.globalResRef) {
                tolref = std::max(tolref, globalMaxResidualContribution_[ndx]);
            } else {
                tolref = std::max(tolref, historicMaxResidualContribution_[i]);
            }
        } else if (settings.globalResRef) {
            tolref = std::max(tolref, pointMaxResidualContribution_);
        }

        // Residual tolerance (Designer's Guide to Spice and Spectre, chapter 2.2.2)
        auto tol = circuit.residualTolerance(rn, tolref);

        // Residual component
        double rescomp = fabs(delta[i]);
    
        // Normalized residual component
        double normResidual = rescomp/tol;

        if (computeNorms) {
            l2normResidual2 += normResidual*normResidual;
            // Update largest normalized component
            if (i==1 || normResidual>maxNormResidual) {
                maxResidual = rescomp;
                maxNormResidual = normResidual;
                maxResidualNode = rn;
            }
        }

        // See if residual component exceeds tolerance
        if (rescomp>tol) {
            residualOk = false;
            // Can exit if not computing norms
            if (!computeNorms) {
                break;
            }
        }
    }
    
    return std::make_tuple(true, residualOk); 
}

std::tuple<bool, bool> OpNRSolver::checkDelta() {
    // Compute norms only in debug mode
    bool computeNorms = settings.debug;

    // In delta we have the solution change
    // Check it for convergence
    
    // Number of unknowns (vector length includes a bucket at index 0)
    auto n = circuit.unknownCount();

    // Raw arrays
    double* xprev = solution.data();
    
    maxDelta = 0.0;
    maxNormDelta = 0.0;
    Node* maxDeltaNode = nullptr;
    
    // Check convergence (see if delta is small enough), 
    // but only if this is iteration 2 or later
    // In iteration 1 assume we did not converge
    
    // Assume we converged
    bool deltaOk = true;
    
    double* xdelta = delta.data();

    // Get point maximum
    pointMaxSolution_ = 0;
    for(decltype(n) i=1; i<=n; i++) {
        double c = std::fabs(xprev[i]);
        if (c>pointMaxSolution_) {
            pointMaxSolution_ = c;
        }
    }

    // Use 1-based index (with bucket) because same indexing is used for variables
    for(decltype(n) i=1; i<=n; i++) {
        // Representative node, associated potential nature index
        auto rn = circuit.reprNode(i);
        bool isPotential = ((rn->flags() & Node::Flags::PotentialNode) == Node::Flags::PotentialNode); 
        size_t ndx = isPotential ? 0 : 1;

        // Compute tolerance reference
        double tolref = std::fabs(xprev[i]);
        
        // Cannot account for new solution because damping has not been performed yet

        // Account for global and historic references
        if (settings.historicSolRef) {
            if (settings.globalSolRef) {
                tolref = std::max(tolref, globalMaxSolution_[ndx]);
            } else {
                tolref = std::max(tolref, historicMaxSolution_[i]);
            }
        } else if (settings.globalSolRef) {
            tolref = std::max(tolref, pointMaxSolution_);
        }
        
        // Compute tolerance
        double tol = circuit.solutionTolerance(rn, tolref);

        // Absolute solution change 
        double deltaAbs = fabs(delta[i]);
        
        if (computeNorms) {
            double normDelta = deltaAbs/tol;
            if (i==1 || normDelta>maxNormDelta) {
                maxDelta = deltaAbs;
                maxNormDelta = normDelta;
                maxDeltaNode = rn;
            }
        }

        // Check tolerance
        if (deltaAbs>tol) {
            // Did not converge
            deltaOk = false;
            
            // Can exit if not computing norms
            if (!computeNorms) {
                break;
            }
        }
    }
    
    return std::make_tuple(true, deltaOk);
}

void OpNRSolver::resetMaxima() {
    zero(historicMaxSolution_);
    zero(historicMaxResidualContribution_);
    zero(globalMaxSolution_);
    zero(globalMaxResidualContribution_);
    pointMaxSolution_ = 0;
    pointMaxResidualContribution_ = 0;
}  

void OpNRSolver::initializeMaxima(OpNRSolver& other) {
    historicMaxSolution_ = other.historicMaxSolution();
    globalMaxSolution_ = other.globalMaxSolution();
    historicMaxResidualContribution_ = other.historicMaxResidualContribution();
    globalMaxResidualContribution_ = other.globalMaxResidualContribution();
}

void OpNRSolver::updateMaxima() {
    auto n = circuit.unknownCount();
    auto* x = solution.data();
    auto* mrc = maxResidualContribution_.data();
    for(decltype(n) i=1; i<=n; i++) {
        bool isPotential = ((circuit.reprNode(i)->flags() & Node::Flags::PotentialNode) == Node::Flags::PotentialNode); 
        
        double c;
        size_t ndx;
        
        // Voltage nodes -> potential nature is voltage (0) 
        // Flow nodes -> potential nature is current (1) 
        ndx = isPotential ? 0 : 1;
        c = std::fabs(x[i]);
        if (c>historicMaxSolution_[i]) {
            historicMaxSolution_[i] = c;
        }
        if (c>globalMaxSolution_[ndx]) {
            globalMaxSolution_[ndx] = c;
        }
        
        // Voltage nodes -> flow nature is current (1) 
        // Flow nodes -> flow nature is voltage (0) 
        ndx = isPotential ? 1 : 0;
        c = std::fabs(mrc[i]);
        if (c>historicMaxResidualContribution_[i]) {
            historicMaxResidualContribution_[i] = c;
        }
        if (c>globalMaxResidualContribution_[ndx]) {
            globalMaxResidualContribution_[ndx] = c;
        }
    }
}

void OpNRSolver::dumpSolution(std::ostream& os, double* solution, const char* prefix) {
    circuit.dumpSolution(os, solution, prefix);
}

}
