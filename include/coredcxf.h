#ifndef __ANCOREDCXF_DEFINED
#define __ANCOREDCXF_DEFINED

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

// Circuit equations
//              d
//   f(x(t)) + ---- q(x(t)) = 0 
//              dt
// 
//   x(t) .. unknowns
//   f(x) .. resistive residual
//   q(x) .. reactive residual

// DC small-signal transfer function analysis
// Assuming t=0 first solves for operating point (x0)
//   f(x0) = 0
// Then it linearizes the circuit by computing the resistive Jacobian Jr
// (Jacobian of f(x)) at x=x0 and solves
//   Jr dx = du
// for each independent source in the system. Here du is set to the 
// incremental excitation of 1 originating from a particular independent 
// source for which the equation is being solved. 
// The obtained dx is used for computing 
// - the small-signal transfer function from the independent source 
//   to the output defined by a node or a pair of nodes, 
// - the small-signal input impedance and admittance felt by a source 
//   at its terminals
//
// See coreop.h on how to specify nodesets. 

typedef struct DCXFParameters {
    OperatingPointParameters opParams;
    
    Value out {""};   // Output node or node pair (string vector)
    Int dumpop {0};   // 1 = dump operating point to <analysisname>.op.raw;
    // Nodeset and store parameters of the operating point core 
    // are also exposed. 

    Int writeOutput {1}; // Do we want to write the results to a file
                         // Not exposed as analysis parameter. 

    DCXFParameters();
} DCXFParameters;


class DCXFCore : public AnalysisCore {
public:
    typedef DCXFParameters Parameters;
    enum class DCXFError {
        OK, 
        NotFound, 
        NotSource, 
        EvalAndLoad, 
        MatrixError, 
        SolutionError, 
        OperatingPointError, 
        SingularMatrix, 
    };
    DCXFCore(
        OutputDescriptorResolver& parentResolver, DCXFParameters& params, OperatingPointCore& opCore, std::unordered_map<Id,size_t>& sourceIndex, 
        Circuit& circuit, KluRealMatrix& jacobian, Vector<double>& incrementalSolution, 
        std::vector<Instance*>& sources, Vector<double>& tf, Vector<double>& yin, 
        Vector<double>& zin
    ); 
    ~DCXFCore();
    
    DCXFCore           (const DCXFCore&)  = delete;
    DCXFCore           (      DCXFCore&&) = delete;
    DCXFCore& operator=(const DCXFCore&)  = delete;
    DCXFCore& operator=(      DCXFCore&&) = delete;

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore) const; 

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
    void clearError() { AnalysisCore::clearError(); lastDcTfError = DCXFError::OK; }; 

    void setError(DCXFError e) { lastDcTfError = e; lastError = Error::OK; };
    DCXFError lastDcTfError;
    Id errorInstance;
    
    KluRealMatrix& jacobian;
    Vector<double>& incrementalSolution;

    std::unordered_map<Id,size_t>& sourceIndex;
    std::vector<Instance*>& sources;
    Vector<double>& tf;
    Vector<double>& yin;
    Vector<double>& zin;

    DCXFParameters& params;
};

}

#endif
