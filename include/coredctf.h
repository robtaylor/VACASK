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

    std::unordered_map<Id,size_t>& sourceIndex;
    std::vector<Instance*>& sources;
    Vector<double>& tf;
    Vector<double>& yin;
    Vector<double>& zin;

    DcTfParameters& params;
};

}

#endif
