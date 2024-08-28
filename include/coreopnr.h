#ifndef __COREOPNR_DEFINED
#define __COREOPNR_DEFINED

#include "nrsolver.h"
#include "common.h"


namespace NAMESPACE {

enum class OpNRSolverFlags : uint8_t { 
    Abort = 1,  // Exit analysis immediately, even in the middle of computing a point
    Finish = 2, // Wait until current point is computed to the end, then exit simulation
                 // i.e. for multipoint analyses (sweep, frequency sweep, time sweep) 
                 // wait until current point is computed, then exit
                 // Do not exit sweep. 
    Stop = 4,   // Stop analysis to possibly continue it later
                 // Exit sweep. 
};
DEFINE_FLAG_OPERATORS(OpNRSolverFlags);

class OpNRSolver : public NRSolver, public FlagBase<OpNRSolverFlags> {
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
    virtual bool postSolve(bool continuePrevious);
    
    virtual std::tuple<bool, bool> buildSystem(bool continuePrevious);
    
    EvalSetup& evalSetupSystem() { return esSystem; };
    LoadSetup& loadSetupSystem() { return lsSystem; };
    
protected:
    void loadShunts(double gshunt, bool loadJacobian=true);
    bool evalAndLoadWrapper(EvalSetup& evalSetup, LoadSetup& loadSetup);
    
    void setNodesetAndIcFlags(bool continuePrevious);

    EvalSetup esSystem;
    LoadSetup lsSystem;
    ConvSetup csSystem;
    
    // Internal structures
    Vector<double> dummyStates;
    Vector<double> deviceStates;

    // Flag that forces skipping of device convergence check
    bool skipConvergenceCheck;
};

}

#endif
