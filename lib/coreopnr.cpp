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
    elsSystem = EvalAndLoadSetup { 
        // Inputs
        .solution = &solution, 
        .states = &states, 

        // Evaluation 
        .enableLimiting = true, 
        .evaluateResistiveJacobian = true, 
        .evaluateResistiveResidual = true, 
        .evaluateLinearizedResistiveRhsResidual = true, 
        .evaluateOpvars = true, 

        // Signal this is static and DC analysis
        .staticAnalysis = true, 
        .dcAnalysis = true, 
        
        // Outputs (vectors need to be set in the loop)
        .loadResistiveJacobian = true, 
        // .linearizedResidualLoadOnlyIfLimited = true, 
    };

    // For computing the residual in damped NR algorithm
    elsResidual = EvalAndLoadSetup { 
        // Inputs, set solution and states when residual is being computed
        .solution = &solution, 
        // Use future solution slot as basis for computing residual in damped NR
        .oldSolutionSlot = -1,
        .states = &states, 
        // New computed states will go to dummyStates to keep actual new states intact
        .dummyStates = &dummyStates,  

        // Evaluation
        .enableLimiting = true, 
        .initializeLimiting = false, 
        .evaluateResistiveResidual = true, 
        .evaluateLinearizedResistiveRhsResidual = true, 

        // Signal this is static and DC analysis
        .staticAnalysis = true, 
        .dcAnalysis = true, 
        
        // Outputs (vectors need to be set in the loop)
        // .linearizedResidualLoadOnlyIfLimited = true, 
    };
}

bool OpNRSolver::rebuild(Status& s) {
    // Call parent's rebuild
    if (!NRSolver::rebuild(s)) {
        return false;
    }

    // Allocate space in vetors
    auto n = circuit.unknownCount();
    maxResidualContribution.resize(n+1);
    dummyStates.resize(circuit.statesCount());
    
    return true;
}

bool OpNRSolver::initialize(bool continuePrevious, Status& s) {
    // This method is called once on entering run()
    // This is the right place to set up vectors
    
    // Call parent's initialize()
    if (!NRSolver::initialize(continuePrevious, s)) {
        return false;
    }

    // Set vectors for building linear system
    bool computeMaxResidualContribution = settings.residualCheck || settings.dampingSteps>0;
    elsSystem.resistiveResidual = delta.data();
    elsSystem.linearizedResistiveRhsResidual = delta.data();
    elsSystem.maxResistiveResidualContribution = computeMaxResidualContribution ? maxResidualContribution.data() : nullptr;
    
    // Set vectors for computing residual (for damping)
    elsResidual.resistiveResidual = delta.data(); 
    elsResidual.linearizedResistiveRhsResidual = delta.data();

    return true;
}

void OpNRSolver::loadShunts(bool loadJacobian) {
    // Now load gshunt if it is greater than 0.0
    // Gshunt current (and its residual contribution) is
    //   gshunt * x
    auto gshunt = circuit.simulatorInternals().gshunt;
    if (gshunt>0.0) {
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
}

bool OpNRSolver::evalAndLoadWrapper(EvalAndLoadSetup& els, Status& s) {
    if (!circuit.evalAndLoad(els, nullptr, s)) {
        // Load error
        if (settings.debug>2) {
            Simulator::dbg() << "Evaluation error.\n";
        }
        return false;
    }
    
    // Update circuit's flags (Abort, Finish, Stop)
    circuit.updateEvalFlags(els);

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

    // Set nodesetEnabled flag in elsSystem
    // Forces slot 0 (continuation nodesets) and 1 (user nodesets) 
    // are used in ordinary OP analysis
    // Nodesets are enabled in iterations 1..op_nsiter 
    // if continuePrevious is false. 
    // They are also enabled if slot 1 (user nodesets) is active. 
    elsSystem.nodesetEnabled = (iteration<=nsiter) && (continuePrevious==false);

    // Set icEnabled flag in elsSystem
    // Slot 2 holds permanent forces for computing initial conditions
    // when OP analysis is invoked from tran core. 
    // Whenever this slot is active transient forces are enables 
    // and we are applying initial conditions. 
    elsSystem.icEnabled = forcesEnabled.size()>2 && forcesEnabled[2];
}

std::tuple<bool, bool> OpNRSolver::buildSystem(bool continuePrevious, Status& s) {
    auto n = circuit.unknownCount();

    // Remove forces originating from nodesets after nsiter iterations
    auto nsiter = circuit.simulatorOptions().core().op_nsiter;
    if (iteration>nsiter) {
        // Continuation nodesets
        enableForces(0, false);
        // User-specified nodesets
        enableForces(1, false);
    }

    // Set nodeset and IC flags
    setNodesetAndIcFlags(continuePrevious); 
    
    // Clear linearized residual, and maximal residual contribution
    zero(maxResidualContribution);
    
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
    elsSystem.initializeLimiting = !continuePrevious && (iteration==1);
    // Write value to simulatorInternals
    circuit.simulatorInternals().initalizeLimiting = elsSystem.initializeLimiting; 

    // Evaluate and load
    if (!evalAndLoadWrapper(elsSystem, s)) {
        return std::make_tuple(false, elsSystem.limitingApplied);
    }
    delta[0] = 0.0;
    
    // Simulator::dbg() << "After loading:\n";
    // Simulator::dbg() << "  Residual:\n";
    // circuit.dumpSolution(std::cout, delta.data(), "    ");
    // Simulator::dbg() << "\n";
    
    // Now load gshunt if it is greater than 0.0
    loadShunts();

    // Prevent convergence if limiting was applied
    return std::make_tuple(true, elsSystem.limitingApplied); 
}

std::tuple<bool, double, double, double, Node*> OpNRSolver::checkResidual(bool* residualOk, bool computeNorms, Status& s) {
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
    if (residualOk) {
        *residualOk = true;
    }
    
    // Go through all variables (except ground)
    for(decltype(n) i=1; i<=n; i++) {
        // Get representative node for i-th variable
        auto rn = circuit.reprNode(i);
        // Residual tolerance (Designer's Guide to Spice and Spectre, chapter 2.2.2)
        auto tol = circuit.residualTolerance(rn, maxResidualContribution[i]);
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
            if (residualOk) {
                *residualOk = false;
                // Can exit if not computing norms
                if (!computeNorms) {
                    break;
                }
            }
        }
    }
    
    return std::make_tuple(true, maxResidual, maxNormResidual, l2normResidual2, maxResidualNode); 
}

std::tuple<bool, double, double, Node*> OpNRSolver::checkDelta(bool* deltaOk, bool computeNorms, Status& s) {
    // In delta we have the solution change
    // Check it for convergence
    
    // Number of unknowns (vector length includes a bucket at index 0)
    auto n = circuit.unknownCount();

    // Raw arrays
    double* xprev = solution.data();
    
    double maxDelta = 0.0;
    double maxNormDelta = 0.0;
    Node* maxDeltaNode = nullptr;
    
    // Check convergence (see if delta is small enough), 
    // but only if this is iteration 2 or later
    // In iteration 1 assume we did not converge
    
    // Assume we converged
    if (deltaOk) {
        *deltaOk = true;
    }

    // Use 1-based index (with bucket) because same indexing is used for variables
    for(decltype(n) i=1; i<=n; i++) {
        auto rn = circuit.reprNode(i);
        double tol = circuit.solutionTolerance(rn, xprev[i]);
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
            if (deltaOk) {
                *deltaOk = false;
            }
            // Can exit if not computing norms
            if (!computeNorms) {
                break;
            }
        }
    
    }
    
    return std::make_tuple(true, maxDelta, maxNormDelta, maxDeltaNode);
}

std::tuple<bool, bool> OpNRSolver::computeResidual(bool continuePrevious, Status& s) {
    // Number of unknowns (vector length includes a bucket at index 0)
    auto n = circuit.unknownCount();

    // Zero dummy state and residual
    zero(dummyStates);
    zero(delta);
    
    // Set nodeset and IC flags
    setNodesetAndIcFlags(continuePrevious); 
    
    // Set nodesetEnabled flag in elsSystem
    // Slots 0 and 1 are nodesets used in ordinary OP analysis
    elsSystem.nodesetEnabled = forcesEnabled[0] || forcesEnabled[1];
    // Set icEnabled flag in elsSystem
    // Slot 2 holds permanent forces for computing initial conditions
    // when OP analysis is invoked from tran core. 
    elsSystem.icEnabled = forcesEnabled.size()>2 && forcesEnabled[2];
    
    // This time do not initialize limiting, divert new state to dummy vector
    if (!evalAndLoadWrapper(elsResidual, s)) {
        return std::make_tuple(false, elsResidual.limitingApplied);
    }
    delta[0] = 0.0;
    
    // Now load gshunt if it is greater than 0.0
    // Do not load Jacobian
    loadShunts(false);

    return std::make_tuple(true, elsResidual.limitingApplied);
}

}
