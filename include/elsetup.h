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


typedef struct EvalSetup {
    // {} for default initialization
    // State and solution repository
    VectorRepository<double>* solution {};
    VectorRepository<double>::DepthIndexDelta oldSolutionSlot {0};
    VectorRepository<double>* states {}; 
    // For diverting new state output to a bucket (when not nullptr)
    Vector<double>* dummyStates {};
    
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

    // Master switch for skipping core evaluations
    bool skipCoreEvaluation {}; 

    // .. what to evaluate beside core 
    bool computeBoundStep {};
    bool computeNextBreakpoint {};
    bool computeMaxFreq {};

    // Store reactive residual in states or dummyStates
    bool storeReactiveState {};

    // Numerical differentiation of residual contributions after core evaluations
    // Results are written to states or dummyStates only if not nullptr
    IntegratorCoeffs* integCoeffs {};
    
    // Return information on what happened during evaluation
    // Verilog-A abort/finish/stop
    bool abortRequested {};
    bool finishRequested {};
    bool stopRequested {};

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

    // 
    // Internals
    // 

    // Fast access pointers - do not set manually
    double* oldSolution; // with bucket
    double* oldStates; // states (current data)
    double* newStates; // can be either from states (future data) or dummyStates (current data)
    
    // Methods
    void clearFlags() {
        discontinuity = -1;
        limitingApplied = false;
        abortRequested = false;
        finishRequested = false;
        stopRequested = false;
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
    
    // What part of Jacobian to load to bound locations
    bool loadResistiveJacobian {};
    bool loadReactiveJacobian {};
    double reactiveJacobianFactor {};
    bool loadTransientJacobian {}; // uses integCoeffs
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
    double* reactiveResidualDerivative {}; // with bucket

    // Maximal resistive residual contribution per node
    // Updated only if not nullptr
    double* maxResistiveResidualContribution {}; // with bucket

    // Maximal reactive residual contribution per node
    // Updated only if not nullptr
    double* maxReactiveResidualContribution {}; // with bucket

    // Maximal reactive residual derivative contribution per node
    // Updated only if not nullptr and integCoeffs is not nullptr
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

}

#endif
