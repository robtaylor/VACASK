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
    // 0 .. continuation nodesets for sweep
    //      cannot contain branch forces
    // 1 .. forces explicitly specified via nodeset analysis parameter
    //      can contain branch forces
    // Slots containing branch forces affect the circuti topology. 
    // They need to be set before rebuild() is called. 

    // For constructing the linearized system in NR loop
    // Add reactive evaluation
    elsSystem.evaluateReactiveJacobian = true;
    elsSystem.evaluateReactiveResidual = true;
    elsSystem.evaluateLinearizedReactiveRhsResidual = true;

    // Set analysis type
    elsSystem.staticAnalysis = false;
    elsSystem.dcAnalysis = false;
    elsSystem.tranAnalysis = true;
    
    // Add reactive Jacobian loading
    elsSystem.integCoeffs = &integCoeffs;
    elsSystem.loadResistiveJacobian = false;
    elsSystem.loadReactiveJacobian = false;
    elsSystem.loadTransientJacobian = true;
    
    // Residual loading
    // elsSystem.linearizedResidualLoadOnlyIfLimited = true;
    
    // State storage
    elsSystem.storeTerminalReactiveResidualState = true;
    elsSystem.storeTerminalReactiveResidualDerivativeState = true;
    
    // Breakpoints, timestep limiting
    elsSystem.computeNextBreakpoint = true;
    elsSystem.computeBoundStep = true;
    
    // For computing the residual (damping)
    // Add reactive evaluation
    elsResidual.evaluateReactiveResidual = true;
    elsResidual.evaluateLinearizedReactiveRhsResidual = true;

    // Set up analysis type
    elsResidual.staticAnalysis = false;
    elsResidual.dcAnalysis = false;
    elsResidual.tranAnalysis = true;
    
    // Integrator
    elsSystem.integCoeffs = &integCoeffs;
    
    // Residual loading
    // elsResidual.linearizedResidualLoadOnlyIfLimited = true;
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
    elsSystem.reactiveResidualDerivative = delta.data();

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
    } else {
        lastError = Error::BadResReference;
        return false;
    }
    
    // Set output vector for computing residual (reactive residual derivative)
    elsResidual.reactiveResidualDerivative = delta.data();
    
    // Compute reactive residual derivative contribution to max residual contribution
    // if we are computing max resistive residual contribution to max residual contribution 
    // Update it in the same vector as max resistive residual contribution. 
    elsSystem.maxReactiveResidualDerivativeContribution = elsSystem.maxResistiveResidualContribution; 
    
    return true;
}  

}
