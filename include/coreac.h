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
    typedef AcParameters Parameters;
    enum class AcError {
        OK, 
        Sweeper, 
        SweepCompute, 
        EvalAndLoad, 
        MatrixError, 
        SolutionError, 
        OpError, 
        SingularMatrix, 
        BadFrequency, 
    };
       
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

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore) const; 

    bool addCoreOutputDescriptors();
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
    void clearError() { AnalysisCore::clearError(); lastAcError = AcError::OK; }; 

    void setError(AcError e) { lastAcError = e; lastError = Error::OK; };
    AcError lastAcError;
    double errorFreq;
    Status errorStatus;

    VectorRepository<double>& dcSolution;
    VectorRepository<double>& dcStates;
    KluRealMatrix& dcJacobian;
    KluComplexMatrix& acMatrix;
    Vector<Complex>& acSolution;
    AcParameters& params;
};

}

#endif
