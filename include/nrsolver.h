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
    Real dampingStep {1.0};
    Int dampingSteps {0};
    bool matrixCheck {};
    bool rhsCheck {};
    bool solutionCheck {};
    Real forceFactor {1e5}; 
} NRSettings;


class NRSolver {
public:
    enum class Error {
        OK, 
        ForcesIndex, 
        EvalAndLoad, 
        LinearSolver, 
        SolutionError, 
        Convergence
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

    // Set extra reference vectors for solution tolerance computation
    void setExtraSolutionToleranceReferences(bool newxref, double* histxref) {
        newxref_ = newxref;
        histxref_ = histxref;
    };

    // Return values: ok, magnitude of residual component with maximal relative magnitude, 
    //                maximal relative magnitude of component, squared L2 norm of relative magnitude vector
    //                corresponding node
    // Relative magnitude is computed wrt. tolerance. 
    std::tuple<bool, double, double, double, Node*> checkResidual(bool* residualOk, bool computeNorms);
    
    // Return values: ok, magnitude of delta component with maximal relative magnitude, 
    //                maximal relative magnitude of component, corresponding node
    // Relative magnitude is computed wrt. tolerance. 
    std::tuple<bool, double, double, Node*> checkDelta(bool* deltaOk, bool computeNorms);

    // Return values: ok, prevent convergence
    virtual std::tuple<bool, bool> computeResidual(bool continuePrevious) = 0;

    // Rebuild internal structures that depend on topology
    virtual bool rebuild();

    // Initialize run (upsize internal structures)
    // Called once at the beginning of NRSolver::run() 
    virtual bool initialize(bool continuePrevious) = 0;

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

    // Passed from outside
    Circuit& circuit;
    KluRealMatrix& jac;
    VectorRepository<double>& states;
    VectorRepository<double>& solution;
    NRSettings& settings;

    // Internal structure for max residual contribution
    Vector<double> maxResidualContribution;
    
    Vector<double*> diagPtrs;
    std::vector<std::vector<std::tuple<double*, double*>>> extraDiags;
    Vector<bool> isFlow;

    std::vector<Forces> forcesList;
    std::vector<bool> forcesEnabled;
    
    Vector<double> rowNorm;

    // Stores residual (before solving), and negative delta (after solving)
    Vector<double> delta;

    bool newxref_;
    double* histxref_;
    
    Int iteration;

    Error lastError;
    Int errorIteration;
};

}

#endif
