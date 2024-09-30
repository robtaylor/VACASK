#ifndef __COREHB_DEFINED
#define __COREHB_DEFINED

#include "densematrix.h"
#include "klubsmatrix.h"
#include "core.h"
#include "corehbnr.h"
#include "outrawfile.h"
#include "corehbnr.h"
#include "value.h"
#include "common.h"

namespace NAMESPACE {

typedef struct HBParameters {
    RealVector freq {};   // Fundamental frequencies (f1, f2, ..., fd) 
    Value nharm {4};      // Number of harmonics for each fundamental frequency 
                          // (H1, H2, ..., Hd) in box truncation
                          // If scalar, applies to all frequencies, 
                          // if vector, components apply to corresponding frequencies. 
    Int immax {0};        // Maximal order of intermodulation products (IM) in diamond truncation. 
                          // If <=0, defaults to largest component of nharm. 
                          // If harmonics is a vector, defaults to its largest component. 
    IntVector imorder {}; // Intermodulation product order 
                          // (reported for each frequency when raw truncation scheme is used)
    Id truncate {Id()};   // Truncation scheme: 
                          // raw     .. values in freq are the frequencies in the spectrum
                          // box     .. box truncation
                          //   kj = 0..Hj, first nonzero kj must be >0
                          // diamond .. diamond truncation (default)
                          //   sum abs(kj) <= IM, first nonzero kj must be >0
    Real samplefac {2};   // Sampling factor in time domain (>=1). 
    Real nper {3};        // Number of periods across which colocation points are selected
    Id sample {Id()};     // Sampling mode (uniform, random), default is random. 

    Int writeOutput {1};  // Do we want to write the results to a file
                          // Not exposed as analysis parameter. 
                             
    HBParameters();
} HBParameters;


class HBCore : public AnalysisCore {
public:
    typedef HBParameters Parameters;
    enum class HBError {
        OK, 
        MatrixError, 
        SolverError, 
    };

    HBCore(
        OutputDescriptorResolver& parentResolver, HBParameters& params, Circuit& circuit, 
        KluBlockSparseRealMatrix& jacobian, VectorRepository<double>& solution
    );
    ~HBCore();
    
    HBCore           (const HBCore&)  = delete;
    HBCore           (      HBCore&&) = delete;
    HBCore& operator=(const HBCore&)  = delete;
    HBCore& operator=(      HBCore&&) = delete;

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore) const; 

    bool addCoreOutputDescriptors();
    bool addDefaultOutputDescriptors();
    bool resolveOutputDescriptors(bool strict, Status& s=Status::ignore);

    std::tuple<bool, bool> requestsRebuild(Status& s = Status::ignore);
    
    bool rebuild(Status& s=Status::ignore); 
    bool initializeOutputs(Id name, Status& s=Status::ignore);
    bool run(bool continuePrevious);
    CoreCoroutine coroutine(bool continuePrevious);
    bool finalizeOutputs(Status& s=Status::ignore);
    bool deleteOutputs(Id name, Status& s=Status::ignore);

    bool buildGrid(Status& s=Status::ignore);
    bool buildColocation(Status& s=Status::ignore);
    bool buildAPFT(Status& s=Status::ignore);
    
    void dump(std::ostream& os) const;

    static Id truncateRaw;
    static Id truncateBox;
    static Id truncateDiamond;
    static Id sampleUniform;
    static Id sampleRandom;

    static bool test();

protected:
    // Clear error
    void clearError() { AnalysisCore::clearError(); lastHbError = HBError::OK; }; 

    void setError(HBError e) { lastHbError = e; lastError = Error::OK; };
    
    bool buildTransformMatrix(DenseMatrix<double>& XF, Status& s=Status::ignore);

    HBError lastHbError;

    KluBlockSparseRealMatrix& bsjac; // Jacobian
    VectorRepository<Real>& solution; // Solution history
    
    OutputRawfile* outfile;

    bool converged_;
    
    struct SpecFreq {
        size_t gridIndex;
        double f;
        bool negated;
        int order;
        bool isHarmonic;
    };

private:
    NRSettings nrSettings;
    HBNRSolver nrSolver;
    VectorRepository<Complex> outputPhasors;
    Complex outputFreq;

    Vector<Complex> solutionFD;

    HBParameters oldParams;
    bool firstBuild;

    HBParameters& params;
    std::vector<SpecFreq> freq;
    std::vector<double> frequencies;
    DenseMatrix<int> grid;
    Vector<double> timepoints;    
    DenseMatrix<double> APFT;
    DenseMatrix<double> IAPFT;
    DenseMatrix<double> DDT;
    DenseMatrix<double> DDTcolMajor;
};

}

#endif
