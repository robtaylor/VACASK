#ifndef __ELSETUP_DEFINED
#define __ELSETUP_DEFINED

#include <cmath>
#include <limits>
#include "ansupport.h"
#include "coretrancoef.h"
#include "common.h"


namespace NAMESPACE {

class Circuit;

class Model;
class Instance;

typedef struct DeviceRequests {
    // Return information on what happened during evaluation
    // Verilog-A abort/finish/stop
    bool abort {};
    bool finish {};
    bool stop {};

    // Methods
    void clear() {
        abort = false;
        finish = false;
        stop = false;
    };
} Requests;

typedef struct EvalSetup {
    // {} for default initialization
    // State and solution repository
    VectorRepository<double>* solution {};
    VectorRepository<double>::DepthIndexDelta oldSolutionSlot {0};
    VectorRepository<double>* states {}; 
    // For diverting new state output to a bucket (when not nullptr)
    Vector<double>* dummyStates {};

    // Data for instance bypass check (previous values)
    double* deviceStates {};
    
    // What mode are we running in - information for evaluator
    bool staticAnalysis {};
    bool dcAnalysis {};
    bool acAnalysis {};
    bool tranAnalysis {};
    bool noiseAnalysis {};
    bool nodesetEnabled {};
    bool icEnabled {};

    // Limiting control
    bool enableLimiting {}; 
    bool initializeLimiting {};

    // Core evaluations
    bool evaluateResistiveJacobian {};
    bool evaluateReactiveJacobian {};
    bool evaluateResistiveResidual {};
    bool evaluateReactiveResidual {};
    bool evaluateLinearizedResistiveRhsResidual {};
    bool evaluateLinearizedReactiveRhsResidual {};
    bool evaluateNoise {};
    bool evaluateOpvars {};

    // Allow bypassing core evaluation
    bool allowBypass {}; 
    
    // Store reactive residual in states or dummyStates
    bool storeReactiveState {};

    // .. what to evaluate beside core 
    bool computeBoundStep {};
    bool computeNextBreakpoint {};
    bool computeMaxFreq {};

    // Numerical differentiation of residual contributions after core evaluations
    // Results are written to states or dummyStates only if not nullptr
    IntegratorCoeffs* integCoeffs {};
    
    // Return information on what happened during evaluation
    // Verilog-A abort/finish/stop
    struct DeviceRequests requests;
    
    // Limiting applied (i.e. $discontinuity(-1))
    bool limitingApplied {};
    
    // Discontinuity signalled
    // Negative when no discontinuity, set first by an instance that calls $discontinuity with 
    // a nonnegative argument, updated by subsequent instances that call $discontinuity 
    // with a lower nonnegative argument. 
    Int discontinuity;
    
    // For setting the upper bound on the timestep
    // Infinite initially, set first by an instance that calls $bound_step with 
    // an argument greater than 0, updated by subsequent instances that call 
    // $bound_step with a lower argument that is greater than 0. 
    double boundStep {};

    // Next breakpoint
    // Infinite initially, set first by an instance if the set value is greater than current
    // time, updated by subsequent instances that set it to a value greater than current time. 
    double nextBreakPoint;
    
    // For setting maximal source frequency
    // Zero initially, increased by instances that generate a signal. 
    double maxFreq {};

    // Counter of instaces that are not converged, is reset by initialize()
    size_t bypassableInstances;
    size_t bypassedInstances;

    // 
    // Internals
    // 

    // Fast access pointers - do not set manually
    double* oldSolution; // with bucket
    double* oldStates; // states (current data)
    double* newStates; // can be either from states (future data) or dummyStates (current data)
    
    // Methods
    void clearFlags() {
        requests.clear();
        discontinuity = -1;
        limitingApplied = false;
    };

    bool initialize() {
        DBGCHECK(states && states->size()<2, "States history must have at least two slots.");
        DBGCHECK(solution && solution->size()<2, "Solution history must have at least two slots.");
        if (solution) {
            oldSolution = solution->data(oldSolutionSlot);
        }
        if (states) {
            oldStates = states->data();
        }
        if (dummyStates) {
            // Dummy states are given when we want to avoid tainting future states
            newStates = dummyStates->data();
        } else if (states) {
            newStates = states->futureData();
        }
        
        if (integCoeffs) {
            DBGCHECK(states->size()<integCoeffs->a().size()+1, "Integration method requires a state history with at least "+std::to_string(integCoeffs->a().size()+1)+" slots.");
            DBGCHECK(states->size()<integCoeffs->b().size()+1, "Integration method requires a state history with at least "+std::to_string(integCoeffs->b().size()+1)+" slots.");
        }

        nextBreakPoint = -1.0;
        boundStep = -1.0;
        maxFreq = 0.0;

        bypassableInstances = 0;
        bypassedInstances = 0;

        return true;
    };

    void clearBounds() { 
        boundStep=std::numeric_limits<double>::infinity(); 
        nextBreakPoint=std::numeric_limits<double>::infinity(); 
        discontinuity=-1; 
        maxFreq=0.0; 
    };
    void setBoundStep(double bound) { if (bound<boundStep) boundStep=bound; };
    void setDiscontinuity(Int i) { if (i<0) return; if (discontinuity<0 || i<discontinuity) discontinuity=i; };
    bool setBreakPoint(double t, SimulatorInternals& internals) {
        if (std::abs(t-internals.time) <= timeRelativeTolerance*internals.time) {
            // Breakpoint now or close to now, it is too late to take it into account. 
            // It should have been set earlier. 
            // Signal discontinuity
        } else if (t<internals.time) {
            // Breakpoint in past, ignore
        } else {
            // Set next breakpoint
            if (t<nextBreakPoint) {
                nextBreakPoint = t;
            }
        }
        return true;
    };
    void setMaxFreq(double freq) { if (freq>maxFreq) maxFreq=freq; }; 
} EvalSetup;


typedef struct LoadSetup {
    // {} for default initialization
    // States - need them whenever
    // - maxReactiveResidualContribution is not nullptr
    // - maxReactiveResidualDerivativeContribution is not nullptr
    // - reactiveResidualDerivative is not nullptr
    // From states we retrieve reactive residual and its derivative wrt time. 
    VectorRepository<double>* states {}; 
    
    // What part of Jacobian to bound locations
    
    // Add resistive Jacobian to bound locations
    bool loadResistiveJacobian {};
    
    // Add reactive Jacobian to bound locations
    // Multiplies Jacobian entries with reactiveJacobianFactor before adding them. 
    bool loadReactiveJacobian {};
    double reactiveJacobianFactor {};

    // Add to bound locations
    // - resistive Jacobian 
    // - reactive Jacobian scaled by integCoeffs->leadingCoeff()
    bool loadTransientJacobian {}; 
    IntegratorCoeffs* integCoeffs {};
    
    // Where to load resistive residual, skip if nullptr
    double* resistiveResidual {}; // with bucket

    // Where to load reactive residual, skip if nullptr
    double* reactiveResidual {}; // with bucket
    
    // Where to load linearized resistive residual, skip loading if nullptr
    double* linearizedResistiveRhsResidual {}; // with bucket

    // Where to load linearized reactive residual, skip loading if nullptr
    double* linearizedReactiveRhsResidual {}; // with bucket

    // Where to load reactive residual derivative, skip if nullptr
    // Assumes reactive residual derivative was computed at evaluation time 
    // and stored in the states vector (i.e. integCoeffs was not nullptr). 
    double* reactiveResidualDerivative {}; // with bucket

    // Maximal resistive residual contribution per node, skip if nullptr
    double* maxResistiveResidualContribution {}; // with bucket

    // Maximal reactive residual contribution per node, skip if nullptr
    // Assumes reactive residual was stored at evaluation time in the states vector
    // (i.e. storeReactiveState was set to true)
    double* maxReactiveResidualContribution {}; // with bucket

    // Maximal reactive residual derivative contribution per node, skip if nullptr
    // Assumes reactive residual derivative was computed at evaluation time 
    // and stored in the states vector (i.e. integCoeffs was not nullptr). 
    double* maxReactiveResidualDerivativeContribution {}; // with bucket

    // Where to load DC small-signal residual, skip if nullptr
    double* dcIncrementResidual {}; // with bucket

    // Where to load AC small-signal residual, skip if nullptr
    Complex* acResidual {}; // with bucket
    
    // 
    // Internals
    // 

    // Fast access pointers - do not set manually
    double* oldStates; // states (current data)
    double* newStates; // states (future data)
    
    // Methods
    bool initialize() {
        DBGCHECK(states && states->size()<2, "States history must have at least two slots.");
        if (states) {
            oldStates = states->data();
            newStates = states->futureData();
        } else {
            oldStates = newStates = nullptr;
        }
        
        return true;
    };
} LoadSetup;


typedef struct ConvSetup {
    // {} for default initialization
    // State and solution repository
    VectorRepository<double>* solution {};
    VectorRepository<double>::DepthIndexDelta oldSolutionSlot {0};
    VectorRepository<double>* states {}; 

    // Data for instance convergence check (previous values)
    double* deviceStates {};

    // Delta vector for inputs convergence test
    double* inputDelta {};

    // Check reactive residual and Jacobian for convergence
    bool checkReactiveConvergece {};

    // Counter of instaces that are not converged, is reset by initialize()
    size_t instancesConvergenceChecks;
    size_t convergedInstances;

    // 
    // Internals
    // 

    // Fast access pointers - do not set manually
    double* oldSolution; // with bucket
    double* oldStates; // states (current data)
    double* newStates; // can be either from states (future data) or dummyStates (current data)
    
    // Methods
    bool initialize() {
        instancesConvergenceChecks = 0;
        convergedInstances = 0;

        oldSolution = solution->data(oldSolutionSlot);
        
        oldStates = states->data();
        newStates = states->futureData();
        
        return true;
    };
} ConvSetup;

}

#endif
