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
        Vector<Complex>& solutionFD,
        Vector<Real>& frequencies, 
        Vector<Real>& timepoints, 
        DenseMatrix<Real>& DDT, 
        DenseMatrix<Real>& DDTcolMajor, 
        DenseMatrix<Real>& APFT, 
        NRSettings& settings
    ); 

    virtual bool rebuild();
    virtual bool initialize(bool continuePrevious);
    virtual bool preIteration(bool continuePrevious);
    virtual bool postSolve(bool continuePrevious);
    virtual bool postConvergenceCheck(bool continuePrevious);
    virtual bool postIteration(bool continuePrevious);
    virtual bool postRun(bool continuePrevious); 
    
    virtual std::tuple<bool, bool> buildSystem(bool continuePrevious);
    virtual std::tuple<bool, bool> checkResidual();
    virtual std::tuple<bool, bool> checkDelta();
    
    EvalSetup& evalSetup() { return evalSetup_; };
    LoadSetup& loadSetup() { return loadSetup_; };
    
protected:
    bool evalAndLoadWrapper(EvalSetup& evalSetup, LoadSetup& loadSetup);
    
    virtual void dumpSolution(std::ostream& os, double* solution, const char* prefix="");
    
    EvalSetup evalSetup_;
    LoadSetup loadSetup_;

    // HB Jacobian
    KluBlockSparseRealMatrix& bsjac;
    
    // Vectors and matrices without a bucket
    std::vector<double>& frequencies;
    Vector<double>& timepoints; 
    DenseMatrix<double>& DDT;
    DenseMatrix<double>& DDTcolMajor;
    DenseMatrix<double>& APFT;
    Vector<Complex>& solutionFD; 
    Circuit& circuit;

    // Internal structures computed at t_k
    // These structures have a bucket because they communicate with 
    // evalAndLoad() which requires vectors to have a bucket. 
    VectorRepository<double> oldSolutionAtTk;
    Vector<double> resistiveResidualAtTk;
    Vector<double> reactiveResidualAtTk;
    Vector<double> maxResidualContributionAtTk_; // maximal residual contribution for all equations at a given timepoint
                                                 // filled with maximal resistive contribution in evalAndLoadWrapper()
    // No bucket, just a dummy
    Vector<double> dummyStates;

    // For all timepoints, no bucket
    Vector<double> resistiveResidual;
    Vector<double> reactiveResidual;
    
    // Internal structure for max residual contribution
    // Has no bucket because it does not communicate with evalAndLoad(). 
    Vector<double> maxResidualContribution_; // maximal residual contribution for each equation at each timepoint
    
    // What kind of tolerance reference to use
    // We support only global/local reference
    // Historic reference is not possible. 
    // We always use point-wise reference. 
    bool globalSolRef;
    bool globalResRef;

    // Global maxima
    DenseMatrix<double> pointMaxResidualContribution_;  // at current solution, maximal value for each nature, each timepoint
                                                        // rows are natures, columns are timepoints
    
    DenseMatrix<double> pointMaxSolution_;  // previous solution, maximal value for each nature, each timepoint
                                            // rows are natures, columns are timepoints

    // Convergence check auxiliary results
    double maxResidual; 
    double maxNormResidual; 
    double l2normResidual2;
    Node* maxResidualNode;
    size_t maxResidualTimepointIndex;
    double maxDelta; 
    double maxNormDelta; 
    Node* maxDeltaNode;
    size_t maxDeltaTimepointIndex;
};

}

#endif
