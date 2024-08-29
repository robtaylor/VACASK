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

typedef struct DcIncrParameters {
    OpParameters opParams;
    
    Int dumpop {0};   // 1 = dump operating point to <analysisname>.op.raw;
    // Nodeset and store parameters of the operating point core 
    // are also exposed. 

    Int writeOutput {1}; // Do we want to write the results to a file
                         // Not exposed as analysis parameter. 

    DcIncrParameters();
} DcIncrParameters;


class DcIncrementalCore : public AnalysisCore {
public:
    typedef DcIncrParameters Parameters;
    enum class DcIncrError {
        OK, 
        EvalAndLoad, 
        MatrixError, 
        SolutionError, 
        OpError, 
        SingularMatrix, 
    };
    DcIncrementalCore(
        Analysis& analysis, DcIncrParameters& params, OperatingPointCore& opCore, Circuit& circuit, 
        KluRealMatrix& jacobian, Vector<double>& incrementalSolution
    ); 
    ~DcIncrementalCore();
    
    DcIncrementalCore           (const DcIncrementalCore&)  = delete;
    DcIncrementalCore           (      DcIncrementalCore&&) = delete;
    DcIncrementalCore& operator=(const DcIncrementalCore&)  = delete;
    DcIncrementalCore& operator=(      DcIncrementalCore&&) = delete;

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore) const; 

    bool addDefaultOutputDescriptors();
    bool resolveOutputDescriptors(bool strict);

    bool rebuild(Status& s=Status::ignore); 
    bool initializeOutputs(Id name);
    CoreCoroutine coroutine(bool continuePrevious);
    bool run(bool continuePrevious);
    bool finalizeOutputs();
    bool deleteOutputs(Id name);

    void dump(std::ostream& os) const;

    OperatingPointCore& opCore_;
    OutputRawfile* outfile;

protected:
    // Clear error
    void clearError() { AnalysisCore::clearError(); lastDcIncrError = DcIncrError::OK; }; 

    void setError(DcIncrError e) { lastDcIncrError = e; lastError = Error::OK; };
    DcIncrError lastDcIncrError;
    double errorFreq;
    Status errorStatus;

    KluRealMatrix& jacobian;
    Vector<double>& incrementalSolution;
    DcIncrParameters& params;
};

}

#endif
