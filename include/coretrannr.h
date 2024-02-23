#ifndef __CORETRANNR_DEFINED
#define __CORETRANNR_DEFINED

#include "coreopnr.h"
#include "coretrancoef.h"
#include "common.h"


namespace NAMESPACE {

// Transient NR solver is almost identical to OP NR solver
// This one is used for all points, but the first one. 
class TranNRSolver : public OpNRSolver {
public:
    TranNRSolver(
        Circuit& circuit, KluRealMatrix& jac, 
        VectorRepository<double>& states, VectorRepository<double>& solution, 
        NRSettings& settings, IntegratorCoeffs& integCoeffs
    ); 

    virtual bool initialize(bool continuePrevious, Status& s=Status::ignore);

    // No need to override buildSysten() and computeResidual() to set 
    // nodeset and ic flags to false because 
    // nodeset flag is off due to continue mode and 
    // ic flag is off due to forces slot 2 not being present. 
    
private:
    IntegratorCoeffs* integCoeffs;
};

}

#endif
