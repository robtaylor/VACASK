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
    
    // Internal structures
    VectorRepository<double> oldSolutionTD;
    VectorRepository<double> oldSolutionTDDot;
    VectorRepository<double> oldSolutionTDtk;
    Vector<double> resistiveResidualTk;
    VectorRepository<double> dummyStatesRepo;
};

}

#endif
