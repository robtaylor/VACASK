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

    // Set analysis type
    esSystem.staticAnalysis = false;
    esSystem.dcAnalysis = false;
    esSystem.tranAnalysis = true;
    
    // For constructing the linearized system in NR loop
    // Add reactive Jacobian and residual evaluation
    esSystem.evaluateReactiveJacobian = true;
    esSystem.evaluateReactiveResidual = true;
    esSystem.evaluateLinearizedReactiveRhsResidual = true;

    // Make sure reactive residual is stored in the state vector at evaluation
    esSystem.storeReactiveState = true;
    
    // Need this for evaluation of residual derivative
    esSystem.integCoeffs = &integCoeffs;

    // Breakpoints, timestep limiting
    esSystem.computeNextBreakpoint = true;
    esSystem.computeBoundStep = true;

    // Also check reactive residual and Jacobian convergence
    csSystem.checkReactiveConvergece = true;
    
    
    // Set up Jacobian loading
    lsSystem.loadResistiveJacobian = false;
    lsSystem.loadReactiveJacobian = false;
    lsSystem.loadTransientJacobian = true;
    lsSystem.integCoeffs = &integCoeffs;
    
    
    // For computing the residual (damping)
    // Set analysis type
    esResidual.staticAnalysis = false;
    esResidual.dcAnalysis = false;
    esResidual.tranAnalysis = true;

    // Add reactive residual evaluation
    esResidual.evaluateReactiveResidual = true;
    esResidual.evaluateLinearizedReactiveRhsResidual = true;

    // Make sure reactive residual is stored in the state vector at evaluation
    esResidual.storeReactiveState = true;
    
    // Need this for evaluation of residual derivative
    esResidual.integCoeffs = &integCoeffs;
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
    lsSystem.reactiveResidualDerivative = delta.data();

    // Compute reactive residual derivative contribution
    // Update it in the max resistive residual contribution vector 
    lsSystem.maxReactiveResidualDerivativeContribution = lsSystem.maxResistiveResidualContribution; 
    
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
    
    // Set output vector for reactive residual derivative
    lsResidual.reactiveResidualDerivative = delta.data();
    
    return true;
}  

}
