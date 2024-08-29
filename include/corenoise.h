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

// Circuit equations
//              d
//   f(x(t)) + ---- q(x(t)) = 0 
//              dt
// 
//   x(t) .. unknowns
//   f(x) .. resistive residual
//   q(x) .. reactive residual

// Small-signal noise analysis
// Assuming t=0 first solves for operating point (x0)
//   f(x0) = 0
// where f(x) is the resistive residual. 
// Then it linearizes the circuit by computing the resistive Jacobian Jr
// (Jacobian of f(x)) and the reactive Jacobian Jc (Jacobian of q(x)) 
// at x=x0 and solves
//   (Jr + omega Jc) X = U
// for each noise source in the system. U is set to reflect the 
// source's magnitude 1 and phase 0. Additionally this is also done  
// for the input source to compute the power gain from input to output. 
// The obtained X is used for computing 
// - the contribution of that noise source 
//   to the output power spectral density 
// - the contribution of all noise sources within individual instances 
//   to the output power spectral density 
// - output power spectral density
// - gain from the input source to the output
// The output is given as a node or a pair of nodes. 
// Frequency is swept across the given range. 
// See coreac.h for details on the frequency sweep. 
// 
// See coreop.h on how to specify nodesets. 

typedef struct NoiseParameters {
    OpParameters opParams;
    
    Value out {""};   // Output node or node pair (string vector)
    Id in {""};       // Input source
    Real from {0};    // Start frequency for step and dec/oct/lin sweep
    Real to {0};      // Stop frequency for step and dec/oct/lin sweep
    Real step {0};    // Step size for step sweep
    Id mode {Id()};   // Mode for dec/oct/lin sweep
    Int points {0};   // Number of points for dec/oct/lin sweep
    Value values {0}; // Vector of values for values sweep
    Int dumpop {0};   // 1 = dump operating point to <analysisname>.op.raw;
    // Nodeset and store parameters of the operating point core 
    // are also exposed. 

    Int writeOutput {1}; // Do we want to write the results to a file
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
        SolutionError, 
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

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore) const; 

    bool addCoreOutputDescriptors();
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
    void clearError() { AnalysisCore::clearError(); lastNoiseError = NoiseError::OK; }; 

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
