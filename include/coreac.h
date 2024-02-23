#ifndef __ANCOREAC_DEFINED
#define __ANCOREAC_DEFINED

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

typedef struct AcParameters {
    OpParameters opParams;

    Real from {0};    // start frequency for step and dec/oct/lin sweep
    Real to {0};      // stop frequency for step and dec/oct/lin sweep
    Real step {0};    // step size for step sweep
    Id mode {Id()}; // mode for dec/oct/lin sweep
    Int points {0};   // number of points for dec/oct/lin sweep
    Value values {0}; // vector of values for values sweep
    Int dumpop {0};   // // 1 = dump operating point to <analysisname>.op.raw;

    Int writeOutput {1}; // Do we want to write the results to a file
                        // Not exposed as analysis parameter. 

    AcParameters();
} AcParameters;


class AcCore : public AnalysisCore {
public:
    AcCore(
        Analysis& analysis, AcParameters& params, OperatingPointCore& opCore, Circuit& circuit, 
        KluRealMatrix& dcJacobian, VectorRepository<double>& dcSolution, VectorRepository<double>& dcStates, 
        KluComplexMatrix& acMatrix, Vector<Complex>& acSolution
    ); 
    ~AcCore();
    
    AcCore           (const AcCore&)  = delete;
    AcCore           (      AcCore&&) = delete;
    AcCore& operator=(const AcCore&)  = delete;
    AcCore& operator=(      AcCore&&) = delete;

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
    AcParameters& params;
};

}

#endif
