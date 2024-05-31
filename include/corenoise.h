#ifndef __ANCORENOISE_DEFINED
#define __ANCORENOISE_DEFINED

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

typedef struct NoiseParameters {
    OpParameters opParams;
    
    Value out {""};   // output node or node pair (string vector)
    Id in {""};       // input source
    Real from {0};    // start frequency for step and dec/oct/lin sweep
    Real to {0};      // stop frequency for step and dec/oct/lin sweep
    Real step {0};    // step size for step sweep
    Id mode {Id()};   // mode for dec/oct/lin sweep
    Int points {0};   // number of points for dec/oct/lin sweep
    Value values {0}; // vector of values for values sweep
    Int dumpop {0};   // // 1 = dump operating point to <analysisname>.op.raw;

    Int writeOutput {1};  // Do we want to write the results to a file
                      // Not exposed as analysis parameter. 

    NoiseParameters();
} NoiseParameters;


class NoiseCore : public AnalysisCore {
public:
    typedef NoiseParameters Parameters;
    enum class NoiseError {
        OK, 
        NotFound, 
        ContribNotFound, 
        Sweeper, 
        SweepCompute, 
        EvalAndLoad, 
        PsdError, 
        MatrixError, 
        OpError, 
        SingularMatrix, 
        BadFrequency, 
    };
    
    NoiseCore(
        Analysis& analysis, NoiseParameters& params, OperatingPointCore& opCore, 
        std::unordered_map<std::pair<Id, Id>, size_t, IdPairHash>& contributionOffset, 
        Circuit& circuit, 
        KluRealMatrix& dcJacobian, VectorRepository<double>& dcSolution, VectorRepository<double>& dcStates, 
        KluComplexMatrix& acMatrix, Vector<Complex>& acSolution, 
        
        Vector<double>& results, double& powerGain, double& outputNoise
    ); 
    ~NoiseCore();
    
    NoiseCore           (const NoiseCore&)  = delete;
    NoiseCore           (      NoiseCore&&) = delete;
    NoiseCore& operator=(const NoiseCore&)  = delete;
    NoiseCore& operator=(      NoiseCore&&) = delete;

    // Clear error
    void clearError() { AnalysisCore::clearError(); lastNoiseError = NoiseError::OK; }; 

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
    void setError(NoiseError e) { lastNoiseError = e; lastError = Error::OK; };
    NoiseError lastNoiseError;
    double errorFreq;
    Status errorStatus;
    Id errorInstance;
    Id errorContrib;
    
    VectorRepository<double>& dcSolution;
    VectorRepository<double>& dcStates;
    KluRealMatrix& dcJacobian;
    KluComplexMatrix& acMatrix; 
    Vector<Complex>& acSolution;

    // second Id is Id() -> total instance contribution
    std::unordered_map<std::pair<Id, Id>, size_t, IdPairHash>& contributionOffset; 
    // noise contributions
    Vector<double>& results;  
    double& powerGain;
    double& outputNoise;
    
    NoiseParameters& params;
};

}

#endif
