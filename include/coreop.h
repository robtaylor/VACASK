#ifndef __ANCOREOP_DEFINED
#define __ANCOREOP_DEFINED

#include "status.h"
#include "circuit.h"
#include "core.h"
#include "klumatrix.h"
#include "output.h"
#include "outrawfile.h"
#include "flags.h"
#include "coreopnr.h"
#include "ansolution.h"
#include "generator.h"
#include "common.h"


namespace NAMESPACE {

// Circuit equations
//              d
//   f(x(t)) + ---- q(x(t)) = 0 
//              dt
// 
//   x(t) .. unknowns
//   f(x) .. resistive residual
//   q(x) .. reactive residual

// Operating point analysis
// Assuming t=0 solves
//   f(x) = 0
//
// Nodesets are specified as a list of values where each nodeset is given
// with 2 or 3 values
// - single node nodesets (...; "<node>"; value; ...)
// - differential nodeset (...; "<node1>"; "<node2>"; value; ...)

typedef struct OperatingPointParameters {
    Value nodeset {Value("")}; // String to specify stored solution slot or
                               // list to specify nodesets
    String store {""};         // Name of stored solution slot to write

    Int writeOutput {1}; // Do we want to write the results to a file
                         // Not exposed as analysis parameter. 

    OperatingPointParameters(); 
} OperatingPointParameters;

enum class OpRunType { 
    OrdinaryOp=0, 
    GshuntStepping=1, 
    GminStepping=2, 
    Spice3GminStepping=3, // Mask = 3, result!=0 for gmin stepping
    SourceStepping=4, // Can be combined with any gmin stepping
    Spice3SourceStepping=8, 
    GSteppingMask=3, 
};
DEFINE_FLAG_OPERATORS(OpRunType);

typedef struct OperatingPointState {
    OperatingPointState() {};
    
    OperatingPointState           (const OperatingPointState&)  = delete;
    OperatingPointState           (      OperatingPointState&&) = default;
    OperatingPointState& operator=(const OperatingPointState&)  = delete;
    OperatingPointState& operator=(      OperatingPointState&&) = default;

    // Annotated solution
    AnnotatedSolution solution;
    // Solution vector with bucket
    Vector<double> stateVector;
    // Is state coherent with current topology
    bool coherent;
    // Is state valid
    bool valid;
} OperatingPointState;

// Operating point core functionality, assumes all circuit parameters and simulator options have been set
// This core uses no other core
class OperatingPointCore : public AnalysisCore {
public:
    using RunType = OpRunType;
    typedef OperatingPointParameters Parameters;
    enum class OpError {
        OK, 
        InitialOp, 
        SteppingSolver, 
        SteppingSteps, 
        NoAlgorithm, 
    };
    
    OperatingPointCore(
        OutputDescriptorResolver& parentResolver, OperatingPointParameters& params, Circuit& circuit, 
        KluRealMatrix& jacobian, VectorRepository<double>& solution, VectorRepository<double>& states
    ); 
    ~OperatingPointCore();
    
    OperatingPointCore           (const OperatingPointCore&)  = delete;
    OperatingPointCore           (      OperatingPointCore&&) = delete;
    OperatingPointCore& operator=(const OperatingPointCore&)  = delete;
    OperatingPointCore& operator=(      OperatingPointCore&&) = delete;

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore) const; 

    bool addDefaultOutputDescriptors();
    bool resolveOutputDescriptors(bool strict, Status& s=Status::ignore);

    std::tuple<bool, bool> preMapping(Status& s=Status::ignore);
    bool populateStructures(Status& s=Status::ignore);

    bool rebuild(Status& s=Status::ignore); 
    bool initializeOutputs(Id name, Status& s=Status::ignore);
    bool run(bool continuePrevious);
    CoreCoroutine coroutine(bool continuePrevious);
    bool finalizeOutputs(Status& s=Status::ignore);
    bool deleteOutputs(Id name, Status& s=Status::ignore);

    size_t stateStorageSize() const;
    void resizeStateStorage(size_t n);
    bool storeState(size_t ndx);
    bool restoreState(size_t ndx);
    void makeStateIncoherent(size_t ndx);

    // Get solver
    OpNRSolver& solver() { return nrSolver; }; 

    void dump(std::ostream& os) const;

protected:
    // Clear error
    void clearError() { AnalysisCore::clearError(); lastOpError = OpError::OK; }; 

    void setError(OpError e) { lastOpError = e; lastError = Error::OK; };
    
    OpError lastOpError;
    RunType errorRunType; 
    Int errorHomotopyIterations;

    bool runSolver(bool continuePrevious);
    
    KluRealMatrix& jac; // Resistive Jacobian
    VectorRepository<double>& solution; // Solution history
    VectorRepository<double>& states; // Circuit states

    OperatingPointState* continueState;

    Forces stateNodesets;
    
    bool continuePrevious;
    bool converged_;

    OutputRawfile* outfile;
    std::vector<OperatingPointState> analysisStateRepository;
    RunType runType;

    PreprocessedUserForces preprocessedNodeset;

private:
    NRSettings nrSettings;
    OpNRSolver nrSolver;

    bool gminStepping(RunType type);
    bool spice3GminStepping();
    bool sourceStepping();
    bool spice3SourceStepping();
    std::string homotopyProgress() const;
    
    OperatingPointParameters& params;
};

}

#endif
