#ifndef __COREOPNR_DEFINED
#define __COREOPNR_DEFINED

#include "nrsolver.h"
#include "common.h"


namespace NAMESPACE {

struct PreprocessedUserForces {
    std::vector<double> nodeValues;
    std::vector<Node*> nodes;
    std::vector<Id> nodeIds;
    std::vector<double> nodePairValues;
    std::vector<std::tuple<Node*, Node*>> nodePairs;
    std::vector<std::tuple<Id, Id>> nodeIdPairs;

    PreprocessedUserForces() {};

    PreprocessedUserForces           (const PreprocessedUserForces&)  = delete;
    PreprocessedUserForces           (      PreprocessedUserForces&&) = default;
    PreprocessedUserForces& operator=(const PreprocessedUserForces&)  = delete;
    PreprocessedUserForces& operator=(      PreprocessedUserForces&&) = default;

    // Preprocess user specified nodeset/ic parameter values, store node ptrs instead of string names. 
    // If syntax is bad, return error. Otherwise simply store the forced value, 
    // node and id (or node pair and ids). 
    // If node is not found, the corresponding force is ignored. 
    // Checks if all extradiagonal sparsity map entries are present. 
    // This is phase 1 of user forces processing. 
    // Return value: ok, needs to add entries to sparsity map
    std::tuple<bool, bool> set(Circuit& circuit, ValueVector& userForces, Status& s=Status::ignore);

    void clear() {
        nodeValues.clear();
        nodes.clear();
        nodeIds.clear();
        nodePairValues.clear();
        nodePairs.clear();
        nodeIdPairs.clear();
    };
};


class OpNRSolver : public NRSolver {
public:
    // By default OpNRSolver has 2 force slots
    // 0 .. continuation nodesets for sweep and homotopy
    //      cannot contain branch forces
    // 1 .. forces explicitly specified via nodeset analysis parameter
    //      can contain branch forces
    // When created in transient analysis for computing initial point
    // there is an extra slot
    // 2 .. forces specified via ic analysis parameter
    //      can contain branch forces
    // Slots containing branch forces affect the circuit topology. 
    // They need to be set before rebuild() is called. 
    // By default we have 2 slots. Transient analysis requests 3 slots. 
    OpNRSolver(
        Circuit& circuit, KluRealMatrix& jac, 
        VectorRepository<double>& states, VectorRepository<double>& solution, 
        NRSettings& settings, Int forcesSize=2
    ); 

    enum class OpNRSolverError {
        OK, 
        ConflictNode, 
        ConflictDelta, 
        LoadForces, 
    };

    // Clear error
    void clearError() { NRSolver::clearError(); lastOpNRError = OpNRSolverError::OK; }; 

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore, NameResolver* resolver=nullptr) const; 

    // Set forces based on an annotated solution
    bool setForces(Int ndx, const AnnotatedSolution& solution, bool abortOnError);

    // Set forces based on preprocessed user forces
    bool setForces(Int ndx, const PreprocessedUserForces& preprocessed, bool uicMode, bool abortOnError);
    
    virtual bool rebuild();
    virtual bool initialize(bool continuePrevious);
    virtual bool preIteration(bool continuePrevious);
    virtual bool postSolve(bool continuePrevious);
    virtual bool postConvergenceCheck(bool continuePrevious);
    virtual bool postIteration(bool continuePrevious);
    
    virtual std::tuple<bool, bool> buildSystem(bool continuePrevious);
    virtual std::tuple<bool, bool> checkResidual();
    virtual std::tuple<bool, bool> checkDelta();

    // Reset solution maxima and residual maxima
    void resetMaxima();

    // Initialize maxima from another OP NR solver
    void initializeMaxima(OpNRSolver& other);

    // Update historic and global maxima (across history and unknowns)
    void updateMaxima();

    // Return max historic solution vector (needed for LTE convergence check)
    const Vector<double>& historicMaxSolution() const { return historicMaxSolution_; };
    const Vector<double>& historicMaxResidualContribution() const { return historicMaxResidualContribution_; };

    // Return max global historic solution
    const Vector<double>& globalMaxSolution() const { return globalMaxSolution_; };
    const Vector<double>& globalMaxResidualContribution() const { return globalMaxResidualContribution_; };

    // Return point max solution
    const Vector<double>& pointMaxSolution() const { return pointMaxSolution_; };
    const Vector<double>& pointMaxResidualContribution() const { return pointMaxResidualContribution_; };
        
    double* maxResidualContribution() { return maxResidualContribution_.data(); }; 
    
    EvalSetup& evalSetup() { return evalSetup_; };
    LoadSetup& loadSetup() { return loadSetup_; };
    ConvSetup& convSetup() { return convSetup_; };

    virtual void dumpSolution(std::ostream& os, double* solution, const char* prefix="");

protected:
    bool setForceOnUnknown(Forces& f, Node* node, double value);

    // Load forces
    bool loadForces(bool loadJacobian=true); 

    void loadShunts(double gshunt, bool loadJacobian=true);
    bool evalAndLoadWrapper(EvalSetup& evalSetup, LoadSetup& loadSetup);
    
    void setNodesetAndIcFlags(bool continuePrevious);

    Vector<double*> diagPtrs;
    std::vector<std::vector<std::tuple<double*, double*>>> extraDiags;
    
    EvalSetup evalSetup_;
    LoadSetup loadSetup_;
    ConvSetup convSetup_;

    // Passed from outside
    Circuit& circuit;
    VectorRepository<double>& states;
    
    // Internal structures
    Vector<double> dummyStates;
    Vector<double> deviceStates;

    // Internal structure for max residual contribution
    Vector<double> maxResidualContribution_; // maximal residual contributionm at this evaluation for each equation
    
    // What kind of tolerance reference to use
    bool historicSolRef;
    bool globalSolRef;
    bool historicResRef;
    bool globalResRef;
    
    // Historic and global maxima
    Vector<double> historicMaxResidualContribution_; // across produced solutions, maximal value for each equation, updated on external command
    Vector<double> globalMaxResidualContribution_;   // accross produced solutions, maximal value for each nature, updated on external command
    Vector<double> pointMaxResidualContribution_;    // at current solution, maximal value for each nature
    
    Vector<double> historicMaxSolution_; // across produced solutions, maximal value for each unknown, updated on external command
    Vector<double> globalMaxSolution_;   // across produced solutions, maximal value for each nature, updated on external command
    Vector<double> pointMaxSolution_;    // previous solution, maximal value for each nature

    // Flags indicating nodes are flow nodes
    Vector<bool> isFlow;

    // Flag that skipping of device convergence check
    bool skipConvergenceCheck;

    // Solution natures and residual natures are currently limited to 
    //   0 .. voltage
    //   1 .. current
    
    // Convergence check auxiliary results
    double maxResidual; 
    double maxNormResidual; 
    double l2normResidual2;
    Node* maxResidualNode;
    double maxDelta; 
    double maxNormDelta; 
    Node* maxDeltaNode;

    OpNRSolverError lastOpNRError;
    Node* errorNode1;
    Node* errorNode2;
};

}

#endif
