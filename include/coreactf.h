#ifndef __ANCOREACTF_DEFINED
#define __ANCOREACTF_DEFINED

#include "status.h"
#include "circuit.h"
#include "core.h"
#include "coreop.h"
#include "klumatrix.h"
#include "output.h"
#include "flags.h"
#include "outrawfile.h"
#include "common.h"


namespace NAMESPACE {

typedef struct AcTfParameters {
    OpParameters opParams;
    
    Value out {""};   // output node or node pair (string vector)
    Real from {0};    // start frequency for step and dec/oct/lin sweep
    Real to {0};      // stop frequency for step and dec/oct/lin sweep
    Real step {0};    // step size for step sweep
    Id mode {Id()};   // mode for dec/oct/lin sweep
    Int points {0};   // number of points for dec/oct/lin sweep
    Value values {0}; // vector of values for values sweep
    Int dumpop {0};   // // 1 = dump operating point to <analysisname>.op.raw;

    Int writeOutput {1};  // Do we want to write the results to a file
                      // Not exposed as analysis parameter. 

    AcTfParameters();
} AcTfParameters;


class AcTfCore : public AnalysisCore {
public:
    AcTfCore(
        Analysis& analysis, AcTfParameters& params, OperatingPointCore& opCore, std::unordered_map<Id,size_t>& sourceIndex, 
        Circuit& circuit, 
        KluRealMatrix& dcJacobian, VectorRepository<double>& dcSolution, VectorRepository<double>& dcStates, 
        KluComplexMatrix& acMatrix, Vector<Complex>& acSolution, 
        std::vector<Instance*>& sources, Vector<Complex>& tf, Vector<Complex>& yin, Vector<Complex>& zin
    ); 
    ~AcTfCore();
    
    AcTfCore           (const AcTfCore&)  = delete;
    AcTfCore           (      AcTfCore&&) = delete;
    AcTfCore& operator=(const AcTfCore&)  = delete;
    AcTfCore& operator=(      AcTfCore&&) = delete;

    bool addCoreOutputDescriptors(Status& s=Status::ignore);
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
    VectorRepository<double>& dcSolution;
    VectorRepository<double>& dcStates;
    KluRealMatrix& dcJacobian;
    KluComplexMatrix& acMatrix; 
    Vector<Complex>& acSolution;

    std::unordered_map<Id,size_t>& sourceIndex;
    std::vector<Instance*>& sources;
    Vector<Complex>& tf;
    Vector<Complex>& yin;
    Vector<Complex>& zin;

    AcTfParameters& params;
};

}

#endif
