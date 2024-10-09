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


typedef struct OperatingPointState {
    OperatingPointState() {};
    
    OperatingPointState           (const OperatingPointState&)  = delete;
    OperatingPointState           (      OperatingPointState&&) = default;
    OperatingPointState& operator=(const OperatingPointState&)  = delete;
    OperatingPointState& operator=(      OperatingPointState&&) = default;

    // Annotated solution
    AnnotatedSolution solution;
    // States vector
    Vector<double> stateVector;
    // Is state coherent with current topology 
    // Becomes coherent when it is written, 
    // stops being coherent when makeStateIncoherent() is called.
    bool coherent;
    // Is state valid 
    // Becomes valid as soon as something is written into the slot. 
    bool valid;
} OperatingPointState;

// Operating point core functionality, assumes all circuit parameters and simulator options have been set
// This core uses no other core
class OperatingPointCore : public AnalysisCore {
public:
    typedef OperatingPointParameters Parameters;
    enum class OperatingPointError {
        OK, 
        Forces, 
        InitialOp, 
        Homotopy, 
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

    virtual size_t stateStorageSize() const;
    virtual size_t allocateStateStorage(size_t n);
    virtual void deallocateStateStorage(size_t n=0);
    virtual bool storeState(size_t ndx, bool storeDetails=true);
    virtual bool restoreState(size_t ndx);
    virtual void makeStateIncoherent(size_t ndx);

    virtual std::tuple<bool, bool> runSolver(bool continuePrevious);
    virtual Int iterations() const;
    virtual Int iterationLimit(bool continuePrevious) const;

    // Get solver
    OpNRSolver& solver() { return nrSolver; }; 

    void dump(std::ostream& os) const;

protected:
    // Clear error
    void clearError() { AnalysisCore::clearError(); lastOpError = OperatingPointError::OK; }; 

    void setError(OperatingPointError e) { lastOpError = e; lastError = Error::OK; };
    
    OperatingPointError lastOpError;
    Int errorForce;

    KluRealMatrix& jac; // Resistive Jacobian
    VectorRepository<double>& solution; // Solution history
    VectorRepository<double>& states; // Circuit states

    OperatingPointState* continueState;

    Forces stateNodesets;
    
    bool continuePrevious;
    bool converged_;

    OutputRawfile* outfile;
    std::vector<OperatingPointState> analysisStateRepository;
    
    PreprocessedUserForces preprocessedNodeset;

private:
    NRSettings nrSettings;
    OpNRSolver nrSolver;

    OperatingPointParameters& params;
};

}

#endif
