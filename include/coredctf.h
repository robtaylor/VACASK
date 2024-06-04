#ifndef __ANCOREDCTF_DEFINED
#define __ANCOREDCTF_DEFINED

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

typedef struct DcTfParameters {
    OpParameters opParams;
    
    Value out {""};   // output node or node pair (string vector)
    Int dumpop {0};   // // 1 = dump operating point to <analysisname>.op.raw;

    Int writeOutput {1}; // Do we want to write the results to a file
                         // Not exposed as analysis parameter. 

    DcTfParameters();
} DcTfParameters;


class DcTfCore : public AnalysisCore {
public:
    typedef DcTfParameters Parameters;
    enum class DcTfError {
        OK, 
        NotFound, 
        NotSource, 
        EvalAndLoad, 
        MatrixError, 
        SolutionError, 
        OpError, 
        SingularMatrix, 
    };
    DcTfCore(
        Analysis& analysis, DcTfParameters& params, OperatingPointCore& opCore, std::unordered_map<Id,size_t>& sourceIndex, 
        Circuit& circuit, KluRealMatrix& jacobian, Vector<double>& incrementalSolution, 
        std::vector<Instance*>& sources, Vector<double>& tf, Vector<double>& yin, 
        Vector<double>& zin
    ); 
    ~DcTfCore();
    
    DcTfCore           (const DcTfCore&)  = delete;
    DcTfCore           (      DcTfCore&&) = delete;
    DcTfCore& operator=(const DcTfCore&)  = delete;
    DcTfCore& operator=(      DcTfCore&&) = delete;

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
    void clearError() { AnalysisCore::clearError(); lastDcTfError = DcTfError::OK; }; 

    void setError(DcTfError e) { lastDcTfError = e; lastError = Error::OK; };
    DcTfError lastDcTfError;
    Id errorInstance;
    
    KluRealMatrix& jacobian;
    Vector<double>& incrementalSolution;

    std::unordered_map<Id,size_t>& sourceIndex;
    std::vector<Instance*>& sources;
    Vector<double>& tf;
    Vector<double>& yin;
    Vector<double>& zin;

    DcTfParameters& params;
};

}

#endif
