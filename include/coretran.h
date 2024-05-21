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

typedef struct TranParameters {
    OpParameters opParams;
    Real step {0.0};
    Real stop {0.0};
    Real start {0.0};
    Real maxStep {0.0};
    Id icMode {"op"};   // op=op with ic forces, uic=spice uic
    Value ic {Value("")};      // string specifying stored solution slot or
                        // list specifying initial conditions
                        // for transient analysis
    String store {""};  // name of stored solution slot to write transient solution to

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
    };
    
    TranCore(
        Analysis& analysis, TranParameters& params, OperatingPointCore& opCore, 
        Circuit& circuit, 
        KluRealMatrix& jacobian, VectorRepository<double>& solution, VectorRepository<double>& states
    ); 
    ~TranCore();
    
    TranCore           (const TranCore&)  = delete;
    TranCore           (      TranCore&&) = delete;
    TranCore& operator=(const TranCore&)  = delete;
    TranCore& operator=(      TranCore&&) = delete;

    // Clear error
    void clearError() { AnalysisCore::clearError(); lastTranError = TranError::OK; }; 

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore) const; 

    bool addCoreOutputDescriptors();
    bool addDefaultOutputDescriptors();
    bool resolveOutputDescriptors(bool strict);

    std::tuple<bool, bool> preMapping(Status& s=Status::ignore);
    bool populateStructures(Status& s=Status::ignore);

    bool rebuild(Status& s=Status::ignore); 
    bool initializeOutputs(Id name);
    bool run(bool continuePrevious);
    bool finalizeOutputs();
    bool deleteOutputs(Id name);

    void dump(std::ostream& os) const;

    static Id icModeOp;
    static Id icModeUic;
    static Id methodAM;
    static Id methodBDF;
    static Id methodGear;
    static Id methodEuler;
    static Id methodTrapezoidal;
    static Id methodBDF2;
    static Id methodGear2;

protected:
    void setError(TranError e) { lastTranError = e; lastError = Error::OK; };
    TranError lastTranError;
    Id errorId;
    
    KluRealMatrix& jacobian; // Resistive Jacobian
    VectorRepository<double>& solution; // Solution history
    VectorRepository<double>& states; // Circuit states

    Vector<double> predictedSolution;
    Vector<double> scaledLte;

    Vector<double> maxSolution;

    OutputRawfile* outfile;

    bool finished; 

    PreprocessedUserForces preprocessedIc;
    Forces uicForces;
    
private:
    bool evalAndLoadWrapper(EvalAndLoadSetup& els);

    // Update breakpoint, but only if it is after last
    void updateBreakPoint(double& bp, double candidate, double last) { if (candidate<bp && candidate>last) bp = candidate; };

    VectorRepository<double> filteredSolution;

    OperatingPointCore& opCore_;
    NRSettings nrSettings;
    TranNRSolver nrSolver;
    CircularBuffer<double> pastTimesteps;
    IntegratorCoeffs integCoeffs; 
    IntegratorCoeffs predictorCoeffs; 
    double reducedOrderTime;
    CircularBuffer<double> breakPoints;
    double acceptedBoundStep;
    double acceptedHmax;
    
    TranParameters& params;

    size_t nPoints;
    double tk;
};

}

#endif
