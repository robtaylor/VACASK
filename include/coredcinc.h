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
    DcIncrementalCore(
        Analysis& analysis, DcIncrParameters& params, OperatingPointCore& opCore, Circuit& circuit, 
        KluRealMatrix& jacobian, Vector<double>& incrementalSolution
    ); 
    ~DcIncrementalCore();
    
    DcIncrementalCore           (const DcIncrementalCore&)  = delete;
    DcIncrementalCore           (      DcIncrementalCore&&) = delete;
    DcIncrementalCore& operator=(const DcIncrementalCore&)  = delete;
    DcIncrementalCore& operator=(      DcIncrementalCore&&) = delete;

    bool addDefaultOutputDescriptors(Status& s=Status::ignore);
    bool resolveOutputDescriptors(bool strict, Status &s);

    bool rebuild(Status& s=Status::ignore); 
    bool initializeOutputs(Id name, Status& s=Status::ignore);
    bool run(bool continuePrevious, Status& s=Status::ignore);
    bool finalizeOutputs(Status& s=Status::ignore);
    bool deleteOutputs(Id name, Status& s=Status::ignore);

    void dump(std::ostream& os) const;

    OperatingPointCore& opCore_;
    OutputRawfile* outfile;

private:
    KluRealMatrix& jacobian;
    Vector<double>& incrementalSolution;
    DcIncrParameters& params;
};

}

#endif
