#ifndef __COREOPNR_DEFINED
#define __COREOPNR_DEFINED

#include "nrsolver.h"
#include "common.h"


namespace NAMESPACE {

class OpNRSolver : public NRSolver {
public:
    // By default OpNRSolver has 2 force slots
    // 0 .. continuation nodesets for sweep
    //      cannot contain branch forces
    // 1 .. forces explicitly specified via nodeset analysis parameter
    //      can contain branch forces
    // When created in transient analysis for computing initial point
    // there is an extra slot
    // 2 .. forces specified via ic analysis parameter
    //      can contain branch forces
    // Slots containing branch forces affect the circuit topology. 
    // They need to be set before rebuild() is called. 
    OpNRSolver(
        Circuit& circuit, KluRealMatrix& jac, 
        VectorRepository<double>& states, VectorRepository<double>& solution, 
        NRSettings& settings, Int forcesSize=2
    ); 

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

    // Initialize maxima from another solver
    void initializeMaxima(OpNRSolver& other);

    // Update historic and global maxima (across history and unknowns)
    void updateMaxima();
        
    EvalSetup& evalSetup() { return evalSetup_; };
    LoadSetup& loadSetup() { return loadSetup_; };
    ConvSetup& convSetup() { return convSetup_; };

    // Return max historic solution vector
    const Vector<double>& historicMaxSolution() const { return historicMaxSolution_; };
    const Vector<double>& historicMaxResidualContribution() const { return historicMaxResidualContribution_; };

    // Return max global historic solution
    const Vector<double>& globalMaxSolution() const { return globalMaxSolution_; };
    const Vector<double>& globalMaxResidualContribution() const { return globalMaxResidualContribution_; };

    // Return point max solution
    double pointMaxSolution() const { return pointMaxSolution_; };
    double pointMaxResidualContribution() const { return pointMaxResidualContribution_; };

    double* maxResidualContribution() { return maxResidualContribution_.data(); }; 
    
protected:
    void loadShunts(double gshunt, bool loadJacobian=true);
    bool evalAndLoadWrapper(EvalSetup& evalSetup, LoadSetup& loadSetup);
    
    void setNodesetAndIcFlags(bool continuePrevious);

    virtual void dumpSolution(std::ostream& os, double* solution, const char* prefix="");

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
    Vector<double> maxResidualContribution_; // maximal residual contributionm at this point
    
    // Historic and global maxima
    Vector<double> historicMaxResidualContribution_; // across produced solutions, updated on external command
    Vector<double> globalMaxResidualContribution_; // accross time and all points, updated on external command
                                                   // one component per each residual nature
    double pointMaxResidualContribution_; // at current point solution
    Vector<double> historicMaxSolution_; // across produced solutions, updated on external command
    Vector<double> globalMaxSolution_; // accross time and all points, updated on external command
                                       // one component per each solution nature
    double pointMaxSolution_; // at current point solution

    // Flags indicating nodes are flow nodes
    Vector<bool> isFlow;

    // Flag that forces skipping of device convergence check
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
};

}

#endif
