#ifndef __ANCORETRAN_DEFINED
#define __ANCORETRAN_DEFINED

#include "status.h"
#include "circuit.h"
#include "core.h"
#include "klumatrix.h"
#include "output.h"
#include "outrawfile.h"
#include "flags.h"
#include "coreop.h"
#include "coretrannr.h"
#include "coretrancoef.h"
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

// Large signal time-domain analysis. 
// icmode="op"  computes the point at t=0 using operating point analysis
//              where initial conditions are forced
// icmode="uic" tries to compute the initial states of reactive components 
//              from the given initial condition and does not perform 
//              operating point analysis. It is equivalent to classic SPICE3 
//              uic transient analysis. 
// 
// See coreop.h on how to specify nodesets. 
// 
// Initial conditions are specified with the same format as nodesets. 

typedef struct TranParameters {
    OperatingPointParameters opParams;
    Real step {0.0};      // Initial timestep
    Real stop {0.0};      // Time up to which the circuit is to be simulated
    Real start {0.0};     // Time at which the results start being recorded
    Real maxstep {0.0};   // Maximal timestep (optional)
    Id icmode {Id()};     // op=op with ic forces, uic=spice uic
    Value ic {Value("")}; // String specifying stored solution slot or
                          // list specifying initial conditions
                          // for transient analysis
    String store {""};    // name of stored solution slot to write transient solution to
    // Nodeset parameter of the operating point core is also exposed. 

    Int writeOutput {1}; // Do we want to write the results to a file
                         // Not exposed as analysis parameter. 
    
    TranParameters();
} TranParameters;

// Operating point core functionality, assumes all circuit parameters and simulator options have been set
// This core uses no other core
class TranCore : public AnalysisCore {
public:
    typedef TranParameters Parameters;
    enum class TranError {
        OK, 
        Tstep, 
        Tstop, 
        Tstart, 
        Method, 
        IcMode, 
        Predictor, 
        Corrector, 
        EvalAndLoad, 
        MatrixError, 
        OpError, 
        UicForces, 
        TimestepTooSmall, 
        BadLteReference, 
        BreakPointPanic, 
    };
    
    TranCore(
        OutputDescriptorResolver& parentResolver, TranParameters& params, OperatingPointCore& opCore, 
        Circuit& circuit, 
        KluRealMatrix& jacobian, VectorRepository<double>& solution, VectorRepository<double>& states
    ); 
    ~TranCore();
    
    TranCore           (const TranCore&)  = delete;
    TranCore           (      TranCore&&) = delete;
    TranCore& operator=(const TranCore&)  = delete;
    TranCore& operator=(      TranCore&&) = delete;

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore) const; 

    bool addCoreOutputDescriptors();
    bool addDefaultOutputDescriptors();
    bool resolveOutputDescriptors(bool strict);

    std::tuple<bool, bool> preMapping(Status& s=Status::ignore);
    bool populateStructures(Status& s=Status::ignore);

    bool rebuild(Status& s=Status::ignore); 
    bool initializeOutputs(Id name, Status& s=Status::ignore);
    void install(ProgressReporter* p);
    CoreCoroutine coroutine(bool continuePrevious);
    bool run(bool continuePrevious);
    bool finalizeOutputs(Status& s=Status::ignore);
    bool deleteOutputs(Id name, Status& s=Status::ignore);

    void dump(std::ostream& os) const;

    static Id icmodeOp;
    static Id icmodeUic;
    static Id methodAM;
    static Id methodBDF;
    static Id methodGear;
    static Id methodEuler;
    static Id methodTrapezoidal;
    static Id methodBDF2;
    static Id methodGear2;

protected:
    // Clear error
    void clearError() { AnalysisCore::clearError(); lastTranError = TranError::OK; }; 

    void setError(TranError e) { lastTranError = e; lastError = Error::OK; };
    TranError lastTranError;
    Id errorId;
    
    KluRealMatrix& jacobian; // Resistive Jacobian
    VectorRepository<double>& solution; // Solution history
    VectorRepository<double>& states; // Circuit states

    Vector<double> predictedSolution;
    Vector<double> scaledLte;

    OutputRawfile* outfile;

    bool finished; 

    PreprocessedUserForces preprocessedIc;
    Forces uicForces;
    
private:
    bool evalAndLoadWrapper(EvalSetup& evalSetup, LoadSetup& loadSetup);
    
    // Update breakpoint, but only if it is after last
    void updateBreakPoint(double& bp, double candidate, double last) { if (candidate<bp && candidate>last) bp = candidate; };

    VectorRepository<double> filteredSolution;

    OperatingPointCore& opCore_;
    NRSettings nrSettings;
    TranNRSolver nrSolver;
    CircularBuffer<double> pastTimesteps;
    IntegratorCoeffs integCoeffs; 
    IntegratorCoeffs predictorCoeffs; 
    CircularBuffer<double> breakPoints;
    double acceptedBoundStep;
    double acceptedHmax;
    
    TranParameters& params;

    size_t nPoints;
    double tk;
};

}

#endif
