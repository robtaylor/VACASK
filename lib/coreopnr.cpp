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
    auto gshunt = circuit.simulatorInternals().gshunt;
    if (gshunt>0) {
        loadShunts(gshunt);
    }

    // Prevent convergence if limiting was applied
    return std::make_tuple(true, elsSystem.limitingApplied); 
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
    auto gshunt = circuit.simulatorInternals().gshunt;
    if (gshunt>0) {
        loadShunts(gshunt, false);
    }

    return std::make_tuple(true, elsResidual.limitingApplied);
}

}
