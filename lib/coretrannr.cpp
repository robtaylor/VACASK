#include "coretrannr.h"
#include "simulator.h"
#include "common.h"
#include <iomanip>

namespace NAMESPACE {

TranNRSolver::TranNRSolver(
    Circuit& circuit, KluRealMatrix& jac, 
    VectorRepository<double>& states, VectorRepository<double>& solution, 
    NRSettings& settings, IntegratorCoeffs& integCoeffs
) : OpNRSolver(circuit, jac, states, solution, settings, 2), integCoeffs(&integCoeffs) {
    // TranNRSolver has 2 force slots
    // 0 .. continuation nodesets for sweep and homotopy
    //      cannot contain branch forces
    // 1 .. forces explicitly specified via nodeset analysis parameter
    //      can contain branch forces
    // Slots containing branch forces affect the circuit topology. 
    // They need to be set before rebuild() is called. 

    // Set analysis type
    evalSetup_.staticAnalysis = false;
    evalSetup_.dcAnalysis = false;
    evalSetup_.tranAnalysis = true;
    
    // For constructing the linearized system in NR loop
    // Add reactive Jacobian and residual evaluation
    evalSetup_.evaluateReactiveJacobian = true;
    evalSetup_.evaluateReactiveResidual = true;
    evalSetup_.evaluateLinearizedReactiveRhsResidual = true;

    // Make sure reactive residual is stored in the state vector at evaluation
    evalSetup_.storeReactiveState = true;
    
    // Need this for evaluation of residual derivative
    evalSetup_.integCoeffs = &integCoeffs;

    // Breakpoints, timestep limiting
    evalSetup_.computeNextBreakpoint = true;
    evalSetup_.computeBoundStep = true;

    // Also check reactive residual and Jacobian convergence
    convSetup_.checkReactiveConvergece = true;
    
    
    // Set up Jacobian loading
    loadSetup_.loadResistiveJacobian = false;
    loadSetup_.loadReactiveJacobian = false;
    loadSetup_.loadTransientJacobian = true;
    loadSetup_.integCoeffs = &integCoeffs;
}

bool TranNRSolver::initialize(bool continuePrevious) {
    // This method is called once on entering run()
    // This is the right place to set vectors
    
    // Call parent's initialize()
    if (!OpNRSolver::initialize(continuePrevious)) {
        // Assume parent has set lastError
        return false;
    }

    // Set output vector for building linear system (reactive residual derivative)
    loadSetup_.reactiveResidualDerivative = delta.data();

    // Compute reactive residual derivative contribution
    // Update it in the max resistive residual contribution vector 
    loadSetup_.maxReactiveResidualDerivativeContribution = loadSetup_.maxResistiveResidualContribution; 
    
    return true;
}  

}
