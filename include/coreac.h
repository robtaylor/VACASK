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

// Circuit equations
//              d
//   f(x(t)) + ---- q(x(t)) = 0 
//              dt
// 
//   x(t) .. unknowns
//   f(x) .. resistive residual
//   q(x) .. reactive residual

// AC small-signal analysis
// Assuming t=0 first solves for operating point (x0)
//   f(x0) = 0
// where f(x) is the resistive residual. 
// Then it linearizes the circuit by computing the resistive Jacobian Jr
// (Jacobian of f(x)) and the reactive Jacobian Jc (Jacobian of q(x)) 
// at x=x0 and solves
//   (Jr + omega Jc) X = U
// for 
//   omega = 2 pi f 
// where the range of f is given as
// - given values (values=[...])
// - stepped linear sweep (from, to, step)
// - linear sweep with given number of points 
//   (from, to, points, mode="lin")
// - logarithmic sweep with given number of points per decade 
//   (from, to, points, mode="dec")
// - logarithmic sweep with given number of points per decade 
//   (from, to, points, mode="oct")
// U comprises the AC excitations specified by the mag and phase 
// parameters of independent sources. Phase is given in degrees. 
// Mag can be negative (equivalent to adding 180 degrees to the phase). 
// The resulting X comprises phasors corresponding to sinusoidal responses 
// in the circuit's unknowns. Sinusoidal signal
//   A cos(omega t + phi) 
// corresponds to phasor
//   A exp(j phi)
// where j is the imaginary unit. 
// 
// See coreop.h on how to specify nodesets. 

typedef struct AcParameters {
    OpParameters opParams;

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
    CoreCoroutine coroutine(bool continuePrevious);
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

    double frequency;
};

}

#endif
