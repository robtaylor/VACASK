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

typedef struct DcIncrParameters {
    OpParameters opParams;
    
    Int dumpop {0};   // // 1 = dump operating point to <analysisname>.op.raw;
    
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
