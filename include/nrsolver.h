#ifndef __NRSOLVER_DEFINED
#define __NRSOLVER_DEFINED

#include "ansupport.h"
#include "options.h"
#include "klumatrix.h"
#include "ansolution.h"
#include "status.h"
#include "common.h"

#include "circuit.h"

namespace NAMESPACE {

typedef struct NRSettings {
    // Input
    int debug {0};
    Int itlim {100};
    Int itlimCont {50};
    Int convIter {1};
    bool residualCheck {true};
    Real dampingFactor {0.8};
    bool matrixCheck {};
    bool rhsCheck {};
    bool solutionCheck {};
    bool historicSolRef {false};
    bool globalSolRef {false};
    bool historicResRef {false};
    bool globalResRef {false};
    Real forceFactor {1e5}; 
} NRSettings;

enum class NRSolverFlags : uint8_t { 
    Abort = 1,  // Exit analysis immediately, even in the middle of computing a point
    Finish = 2, // Wait until current point is computed to the end, then exit simulation
                 // i.e. for multipoint analyses (sweep, frequency sweep, time sweep) 
                 // wait until current point is computed, then exit
                 // Do not exit sweep. 
    Stop = 4,   // Stop analysis to possibly continue it later
                 // Exit sweep. 
};
DEFINE_FLAG_OPERATORS(NRSolverFlags);

class NRSolver : public FlagBase<NRSolverFlags> {
public:
    enum class Error {
        OK, 
        ForcesIndex, 
        EvalAndLoad, 
        ConvergenceCheck, 
        LinearSolver, 
        SolutionError, 
        Convergence, 
        BadSolReference, 
        BadResReference, 
    };

    NRSolver(
        Circuit& circuit, KluRealMatrix& jac, 
        VectorRepository<double>& states, VectorRepository<double>& solution, 
        NRSettings& settings
    );

    // Clear error
    void clearError() { lastError = Error::OK; }; 

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore, NameResolver* resolver=nullptr) const; 

    // Return value: ok, prevent convergence
    virtual std::tuple<bool, bool> buildSystem(bool continuePrevious) = 0;

    // Return values: ok, magnitude of residual component with maximal relative magnitude, 
    //                maximal relative magnitude of component, squared L2 norm of relative magnitude vector, 
    //                corresponding node identifier
    // Relative magnitude is computed wrt. tolerance. 
    virtual std::tuple<bool, double, double, double, Id> checkResidual(bool* residualOk, bool computeNorms) = 0;
    
    // Return values: ok, magnitude of delta component with maximal relative magnitude, 
    //                maximal relative magnitude of component, corresponding node identifier
    // Relative magnitude is computed wrt. tolerance. 
    virtual std::tuple<bool, double, double, Id> checkDelta(bool* deltaOk, bool computeNorms) = 0;

    // Rebuild internal structures that depend on topology
    virtual bool rebuild();

    // Initialize run (upsize internal structures)
    // Called once at the beginning of NRSolver::run() 
    virtual bool initialize(bool continuePrevious) = 0;

    // Pre-iteration tasks, called at the beginning of iteration
    virtual bool preIteration(bool continuePrevious) { return true; };

    // Post-solve tasks
    virtual bool postSolve(bool continuePrevious) { return true; };

    // Post-iteration tasks, called at the end of iteration
    virtual bool postIteration(bool continuePrevious) { return true; };

    // Resize forces repository
    void resizeForces(Int n);

    // Get forces from a forces slot
    Forces& forces(Int ndx);

    // Enable/disable forces slot
    bool enableForces(Int ndx, bool enable) {
        lastError = Error::OK;

        if (ndx<0 || ndx>=forcesList.size()) {
            lastError = Error::ForcesIndex;
            return false;
        }
        forcesEnabled[ndx] = enable; 
        return true;
    };

    // Do we have any forces
    bool haveForces() {
        auto nf = forcesList.size();
        for(decltype(nf) iForce=0; iForce<nf; iForce++) {
            if (forcesEnabled[iForce]) {
                return true;
            }
        }
        return false;
    };

    // Clear error, run solver
    // Return values: ok, number of iterations
    bool run(bool continuePrevious);

    // Return number of iterations spent in NR loop
    Int iterations() const { return iteration; }; 

protected:
    // Load forces
    bool loadForces(bool loadJacobian=true); 
    
    // High precision requested
    bool highPrecision;

    // Passed from outside
    Circuit& circuit;
    KluRealMatrix& jac;
    VectorRepository<double>& states;
    VectorRepository<double>& solution;
    NRSettings& settings;

    // Solution natures and residual natures are currently limited to 
    //   0 .. voltage
    //   1 .. current
    
    Vector<double*> diagPtrs;
    std::vector<std::vector<std::tuple<double*, double*>>> extraDiags;
    Vector<bool> isFlow;

    std::vector<Forces> forcesList;
    std::vector<bool> forcesEnabled;
    
    Vector<double> rowNorm;

    // Stores residual (before solving), and negative delta (after solving)
    Vector<double> delta;
    
    Int iteration;

    Error lastError;
    Int errorIteration;
};

}

#endif
