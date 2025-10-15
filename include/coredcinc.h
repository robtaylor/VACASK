#ifndef __ANCOREDCINC_DEFINED
#define __ANCOREDCINC_DEFINED

#include "status.h"
#include "circuit.h"
#include "core.h"
#include "coreop.h"
#include "klumatrix.h"
#include "output.h"
#include "outrawfile.h"
#include "flags.h"
#include "outrawfile.h"
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

// DC incremental small-signal analysis
// Assuming t=0 first solves for operating point (x0)
//   f(x0) = 0
// Then it linearizes the circuit by computing the resistive Jacobian Jr
// (Jacobian of f(x)) at x=x0 and solves
//   Jr dx = du
// where du comprises incremental excitations specified by the mag parameters 
// of independent sources: mag can also be negative. 
// 
// See coreop.h on how to specify nodesets. 

typedef struct DCIncrementalParameters {
    OperatingPointParameters opParams;
    
    Int writeop {0};  // 1 = dump operating point to <analysisname>.op.raw;
    // Nodeset and store parameters of the operating point core 
    // are also exposed. 

    Int write {1};    // Write the results to a file

    DCIncrementalParameters();
} DCIncrementalParameters;


class DCIncrementalCore : public AnalysisCore {
public:
    typedef DCIncrementalParameters Parameters;
    enum class DCIncrementalError {
        OK, 
        EvalAndLoad, 
        MatrixError, 
        SolutionError, 
        OperatingPointError, 
        SingularMatrix, 
    };
    DCIncrementalCore(
        OutputDescriptorResolver& parentResolver, DCIncrementalParameters& params, OperatingPointCore& opCore, Circuit& circuit, 
        CommonData& commons, KluRealMatrix& jacobian, Vector<double>& incrementalSolution
    ); 
    ~DCIncrementalCore();
    
    DCIncrementalCore           (const DCIncrementalCore&)  = delete;
    DCIncrementalCore           (      DCIncrementalCore&&) = delete;
    DCIncrementalCore& operator=(const DCIncrementalCore&)  = delete;
    DCIncrementalCore& operator=(      DCIncrementalCore&&) = delete;

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore) const; 

    bool addDefaultOutputDescriptors();
    bool resolveOutputDescriptors(bool strict);

    bool rebuild(Status& s=Status::ignore); 
    bool initializeOutputs(Id name, Status& s=Status::ignore);
    CoreCoroutine coroutine(bool continuePrevious);
    bool run(bool continuePrevious);
    bool finalizeOutputs(Status& s=Status::ignore);
    bool deleteOutputs(Id name, Status& s=Status::ignore);

    void dump(std::ostream& os) const;

    OperatingPointCore& opCore_;
    OutputRawfile* outfile;

protected:
    // Clear error
    void clearError() { AnalysisCore::clearError(); lastDcIncrError = DCIncrementalError::OK; }; 

    void setError(DCIncrementalError e) { lastDcIncrError = e; lastError = Error::OK; };
    DCIncrementalError lastDcIncrError;
    double errorFreq;
    Status errorStatus;

    KluRealMatrix& jacobian;
    Vector<double>& incrementalSolution;
    DCIncrementalParameters& params;
};

}

#endif
