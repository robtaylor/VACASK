#ifndef __COREHBNR_DEFINED
#define __COREHBNR_DEFINED

#include "nrsolver.h"
#include "klubsmatrix.h"
#include "common.h"


namespace NAMESPACE {

class HBNRSolver : public NRSolver {
public:
    // No forces. 
    HBNRSolver(
        Circuit& circuit, KluBlockSparseRealMatrix& bsjac, 
        VectorRepository<double>& solution, 
        Vector<Real>& spectrum, 
        Vector<Real>& timepoints, 
        DenseMatrix<Real>& XF, 
        DenseMatrix<Real>& XFdot, 
        NRSettings& settings
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
    
    EvalSetup& evalSetup() { return evalSetup_; };
    LoadSetup& loadSetup() { return loadSetup_; };
    
protected:
    bool evalAndLoadWrapper(EvalSetup& evalSetup, LoadSetup& loadSetup);
    
    EvalSetup evalSetup_;
    LoadSetup loadSetup_;

    KluBlockSparseRealMatrix& bsjac;
    Vector<Real>& spectrum; 
    Vector<Real>& timepoints; 
    DenseMatrix<Real>& XF; 
    DenseMatrix<Real>& XFdot; 
    Circuit& circuit;

    // Maximal value in each row of XF, XFdot
    Vector<double> XFrowMax;
    Vector<double> XFdotRowMax;
    
    // Internal structures
    VectorRepository<double> oldSolutionTD;
    VectorRepository<double> oldSolutionTDDot;
    VectorRepository<double> oldSolutionTDtk;
    Vector<double> resistiveResidualTk;
    VectorRepository<double> dummyStatesRepo;

    // Internal structure for max residual contribution
    Vector<double> maxResidualContributionAtTimepoint_; // maximal residual contribution for all equations at a given timepoint
                                                        // filled with maximal resistive contribution in evalAndLoadWrapper()
    Vector<double> maxResidualContribution_; // maximal residual contribution for each equation at each timepoint
    
    // What kind of tolerance reference to use
    // Solution has only global/local reference
    // Point/historic reference distinction is not available 
    // because the solution is in the frequency domain
    bool globalSolRef;
    bool historicResRef;
    bool globalResRef;

    // Historic and global maxima
    Vector<double> historicMaxResidualContribution_; // across produced solutions, maximal value for each euqation
    Vector<double> globalMaxResidualContribution_;   // accross produced solutions, maximal value for each nature
    Vector<double> pointMaxResidualContribution_;    // at current solution, maximal value for each nature
    
    DenseMatrix<double> pointMaxSolution_;  // previous solution, maximal value for each nature
                                            // rows are natures, columns are frequency components

    // Convergence check auxiliary results
    double maxResidual; 
    double maxNormResidual; 
    double l2normResidual2;
    Node* maxResidualNode;
    size_t maxResidualTimepointIndex;
    double maxDelta; 
    double maxNormDelta; 
    Node* maxDeltaNode;
    size_t maxDeltaFreqIndex;
};

}

#endif
