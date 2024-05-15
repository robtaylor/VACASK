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
    // Slots containing branch forces affect the circuti topology. 
    // They need to be set before rebuild() is called. 
    OpNRSolver(
        Circuit& circuit, KluRealMatrix& jac, 
        VectorRepository<double>& states, VectorRepository<double>& solution, 
        NRSettings& settings, Int forcesSize=2
    ); 

    virtual bool rebuild();
    virtual bool initialize(bool continuePrevious);
    
    virtual std::tuple<bool, bool> buildSystem(bool continuePrevious);
    virtual std::tuple<bool, bool> computeResidual(bool continuePrevious);

    EvalAndLoadSetup& evalSetupSystem() { return elsSystem; };
    EvalAndLoadSetup& evalSetupResidual() { return elsResidual; };

protected:
    void loadShunts(double gshunt, bool loadJacobian=true);
    bool evalAndLoadWrapper(EvalAndLoadSetup& els);

    void setNodesetAndIcFlags(bool continuePrevious);
    
    EvalAndLoadSetup elsSystem;
    EvalAndLoadSetup elsResidual; 

    // Internal structures
    Vector<double> dummyStates;
};

}

#endif
