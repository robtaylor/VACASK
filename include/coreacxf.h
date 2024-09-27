#ifndef __ANCOREACXF_DEFINED
#define __ANCOREACXF_DEFINED

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

// AC small-signal transfer function analysis
// Assuming t=0 first solves for operating point (x0)
//   f(x0) = 0
// where f(x) is the resistive residual. 
// Then it linearizes the circuit by computing the resistive Jacobian Jr
// (Jacobian of f(x)) and the reactive Jacobian Jc (Jacobian of q(x)) 
// at x=x0 and solves
//   (Jr + omega Jc) X = U
// for each independent source in the system. U is set to reflect the 
// source's magnitude 1 and phase 0. 
// The obtained X is used for computing 
// - the small-signal transfer function from the independent source 
//   to the output defined by a node or a pair of nodes, 
// - the small-signal input impedance and admittance felt by a source 
//   at its terminals
// Frequency is swept across the given range. 
// See coreac.h for details on the frequency sweep. 
// 
// See coreop.h on how to specify nodesets. 

typedef struct ACXFParameters {
    OperatingPointParameters opParams;
    
    Value out {""};   // Output node or node pair (string vector)
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

    ACXFParameters();
} ACXFParameters;


class ACXFCore : public AnalysisCore {
public:
    typedef ACXFParameters Parameters;
    enum class ACXFError {
        OK, 
        NotFound, 
        NotSource, 
        Sweeper, 
        SweepCompute, 
        EvalAndLoad, 
        MatrixError, 
        SolutionError, 
        OperatingPointError, 
        SingularMatrix, 
        BadFrequency, 
    };
    
    ACXFCore(
        OutputDescriptorResolver& parentResolver, ACXFParameters& params, OperatingPointCore& opCore, std::unordered_map<Id,size_t>& sourceIndex, 
        Circuit& circuit, 
        KluRealMatrix& dcJacobian, VectorRepository<double>& dcSolution, VectorRepository<double>& dcStates, 
        KluComplexMatrix& acMatrix, Vector<Complex>& acSolution, 
        std::vector<Instance*>& sources, Vector<Complex>& tf, Vector<Complex>& yin, Vector<Complex>& zin
    ); 
    ~ACXFCore();
    
    ACXFCore           (const ACXFCore&)  = delete;
    ACXFCore           (      ACXFCore&&) = delete;
    ACXFCore& operator=(const ACXFCore&)  = delete;
    ACXFCore& operator=(      ACXFCore&&) = delete;

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore) const; 

    bool addCoreOutputDescriptors();
    bool addDefaultOutputDescriptors();
    bool resolveOutputDescriptors(bool strict);

    bool rebuild(Status& s=Status::ignore); 
    bool initializeOutputs(Id name, Status& s=Status::ignore);
    CoreCoroutine coroutine(bool continuePrevious);
    bool run(bool continuePrevious);
    bool finalizeOutputs(Status& s=Status::ignore);
    bool deleteOutputs(Id name, Status& s=Status::ignore);

    void dump(std::ostream& os) const;

    OperatingPointCore& opCore_;
    OutputRawfile* outfile;

protected:
    // Clear error
    void clearError() { AnalysisCore::clearError(); lastAcTfError = ACXFError::OK; }; 

    void setError(ACXFError e) { lastAcTfError = e; lastError = Error::OK; };
    ACXFError lastAcTfError;
    double errorFreq;
    Status errorStatus;
    Id errorInstance;

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

    ACXFParameters& params;

    double frequency;
};

}

#endif
