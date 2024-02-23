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
    bool infCheck {};
    bool nanCheck {};
    Real forceFactor {1e5}; 
} NRSettings;


class NRSolver {
public:
    NRSolver(
        Circuit& circuit, KluRealMatrix& jac, 
        VectorRepository<double>& states, VectorRepository<double>& solution, 
        NRSettings& settings
    );

    // Return value: ok, prevent convergence
    virtual std::tuple<bool, bool> buildSystem(bool continuePrevious, Status& s=Status::ignore) = 0;

    // Return values: ok, magnitude of residual component with maximal relative magnitude, 
    //                maximal relative magnitude of component, squared L2 norm of relative magnitude vector
    //                corresponding node
    // Relative magnitude is computed wrt. tolerance. 
    virtual std::tuple<bool, double, double, double, Node*> checkResidual(bool* residualOk, bool computeNorms, Status& s=Status::ignore) = 0;
    
    // Return values: ok, magnitude of delta component with maximal relative magnitude, 
    //                maximal relative magnitude of component, corresponding node
    // Relative magnitude is computed wrt. tolerance. 
    virtual std::tuple<bool, double, double, Node*> checkDelta(bool* deltaOk, bool computeNorms, Status& s=Status::ignore) = 0;

    // Return values: ok, prevent convergence
    virtual std::tuple<bool, bool> computeResidual(bool continuePrevious, Status& s=Status::ignore) = 0;

    // Rebuild internal structures that depend on topology
    virtual bool rebuild(Status& s=Status::ignore);

    // Initialize run (upsize internal structures)
    // Called once at the beginning of NRSolver::run() 
    virtual bool initialize(bool continuePrevious, Status& s=Status::ignore);

    // Resize forces repository
    void resizeForces(Int n);

    // Get forces from a forces slot
    Forces& forces(Int ndx);

    // Enable/disable forces slot
    bool enableForces(Int ndx, bool enable, Status& s=Status::ignore);

    // Return values: ok, number of iterations
    bool run(bool continuePrevious, Status& s=Status::ignore);

    // Return number of iterations spent in NR loop
    Int iterations() const { return iteration; }; 

protected:
    // Load forces
    bool loadForces(bool loadJacobian=true, Status& s=Status::ignore); 

    // Passed from outside
    Circuit& circuit;
    KluRealMatrix& jac;
    VectorRepository<double>& states;
    VectorRepository<double>& solution;
    NRSettings& settings;

    Vector<double*> diagPtrs;
    std::vector<std::vector<std::tuple<double*, double*>>> extraDiags;
    Vector<bool> isFlow;

    std::vector<Forces> forcesList;
    std::vector<bool> forcesEnabled;
    
    Vector<double> rowNorm;

    // Stores residual (before solving), and negative delta (after solving)
    Vector<double> delta;
    
    Int iteration;
};

}

#endif
