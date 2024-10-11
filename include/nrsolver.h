#ifndef __NRSOLVER_DEFINED
#define __NRSOLVER_DEFINED

#include "ansupport.h"
#include "options.h"
#include "klumatrix.h"
#include "ansolution.h"
#include "status.h"
#include "acct.h"
#include "common.h"

#include "circuit.h"

namespace NAMESPACE {

// // After every topology change or variables change
// rebuild(); 
//
// run(continuePrevious) {
//   initialize()
//   if (not continuePrevious) {
//     zero old solution vector;
//   }
//   do {
//       zero new solution, Jacobian, delta/residual vector;
//       preIteration();
//       buildSystem();
//       load forces;
//       check matrix and rhs;
//       checkResidual();
//       factor and solve;
//       postSolve();
//       check solution;
//       checkDelta() if not in first iteration;
//       check iteration convergence;
//       check convergence;
//       postConvergenceCheck();
//       if (converged) {
//           exit loop;
//       }
//       compute new solution;
//       rotate solutions;
//       postIteration();
//   } while (iteration limit not reached);
// }

class Forces {
public:
    Forces();

    Forces           (const Forces&)  = delete;
    Forces           (      Forces&&) = default;
    Forces& operator=(const Forces&)  = delete;
    Forces& operator=(      Forces&&) = default;

    // Clear forces
    void clear();

    void dump(Circuit& circuit, std::ostream& os) const;
    
    Vector<double> unknownValue_;
    Vector<bool> unknownForced_;
    Vector<double> deltaValue_;
    Vector<std::tuple<UnknownIndex, UnknownIndex>> deltaIndices_;
};


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
        Accounting& acct, 
        KluRealMatrixCore& jac, VectorRepository<double>& solution, 
        NRSettings& settings
    );

    // Clear error
    void clearError() { lastError = Error::OK; }; 

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore, NameResolver* resolver=nullptr) const; 

    // Return value: ok, prevent convergence
    virtual std::tuple<bool, bool> buildSystem(bool continuePrevious) = 0;

    // Return values: ok, residual ok
    virtual std::tuple<bool, bool> checkResidual() = 0;
    
    // Return values: ok, delta ok
    virtual std::tuple<bool, bool> checkDelta() = 0;

    // Rebuild internal structures that depend on topology
    virtual bool rebuild();

    // Initialize run (upsize internal structures)
    // Called once at the beginning of NRSolver::run() 
    // Must set lastError on failure
    virtual bool initialize(bool continuePrevious) = 0;

    // Pre-iteration tasks, called at the beginning of iteration
    // Must set lastError on failure
    virtual bool preIteration(bool continuePrevious) { return true; };

    // Post-solve tasks
    // Must set lastError on failure
    virtual bool postSolve(bool continuePrevious) { return true; };

    // Post convergence check tasks, called after convergence check regardless 
    // of its outcome. 
    // Convergence means that sufficient consecutive iterations converge. 
    // Must set lastError on failure
    virtual bool postConvergenceCheck(bool continuePrevious) { return true; };

    // Post-iteration tasks, called at the end of iteration. 
    // At this point the solution has been rotated and the current solution
    // is the one that will be used for building the next NR system. 
    // Must set lastError on failure
    virtual bool postIteration(bool continuePrevious) { return true; };

    // Post run tasks, called before exit
    virtual bool postRun(bool continuePrevious) { return true; }; 

    // Resize forces repository
    void resizeForces(Int n);

    // Get forces from a forces slot
    Forces& forces(Int ndx);
    const Forces& forces(Int ndx) const;

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

    // Dump solution
    virtual void dumpSolution(std::ostream& os, double* solution, const char* prefix="") {};
    
    // High precision requested
    bool highPrecision;

    // Passed from outside
    KluRealMatrixCore& jac;
    VectorRepository<double>& solution;
    NRSettings& settings;
    Accounting& acct;

    Vector<double*> diagPtrs;
    std::vector<std::vector<std::tuple<double*, double*>>> extraDiags;
    
    std::vector<Forces> forcesList;
    std::vector<bool> forcesEnabled;
    
    Vector<double> rowNorm;

    // Stores residual (before solving), and negative delta (after solving)
    Vector<double> delta;
    
    Int iteration;

    // Convergence check results
    // System build requested no convergence
    bool preventedConvergence;
    // Residual is within tolerances
    bool residualOk;
    // Solution change is within tolerances
    bool deltaOk;
    // Iteration converged (convergence not prevented, residual and delta are ok)
    bool iterationConverged;
    // Sufficient number of consecutive iterations converged (settings.convIter), solver done
    bool converged;

    // Error information
    Error lastError;
    Int errorIteration;
    Int errorForcesIndex;
};

}

#endif
