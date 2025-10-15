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
    Value nodeset {Value("")}; // String specifying stored solution slot to read or
                               // list specifying nodesets
    String store {""};         // Name of stored solution slot to write

    Int write {1};             // Write the results to a file

    OperatingPointParameters(); 
} OperatingPointParameters;


// Operating point core functionality, assumes all circuit parameters and simulator options have been set
// This core uses no other core
class OperatingPointCore : public AnalysisCore {
public:
    typedef OperatingPointParameters Parameters;
    enum class OperatingPointError {
        OK, 
        InitialOp, 
        Homotopy, 
        NoAlgorithm, 
    };
    
    OperatingPointCore(
        OutputDescriptorResolver& parentResolver, OperatingPointParameters& params, Circuit& circuit, 
        CommonData& commons, 
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

    virtual bool storeState(size_t ndx, bool storeDetails=true);
    virtual bool restoreState(size_t ndx);
    
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

    CoreStateStorage* continueState;

    Forces stateNodesets;
    
    bool continuePrevious;
    bool converged_;

    OutputRawfile* outfile;
    
    PreprocessedUserForces preprocessedNodeset;

private:
    NRSettings nrSettings;
    OpNRSolver nrSolver;

    OperatingPointParameters& params;
};

}

#endif
