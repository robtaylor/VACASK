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
    RealVector freq {};    // Fundamental frequencies (f1, f2, ..., fd) 
    Value nharm {4};       // Number of harmonics for each fundamental frequency 
                           // (H1, H2, ..., Hd) in box truncation
                           // If scalar, applies to all frequencies, 
                           // if vector, components apply to corresponding frequencies. 
    Int immax {0};         // Maximal order of intermodulation products (IM) in diamond truncation. 
                           // If <=0, defaults to largest component of nharm. 
                           // If harmonics is a vector, defaults to its largest component. 
    Id truncate {Id()};    // Truncation scheme: 
                           // raw     .. values in freq are the frequencies in the spectrum
                           // box     .. box truncation
                           //   kj = 0..Hj, first nonzero kj must be >0
                           // diamond .. diamond truncation (default)
                           //   sum abs(kj) <= immax, first nonzero kj must be >0
    Real samplefac {2};    // Sampling factor in time domain (>=1). 
    Real nper {3};         // Number of periods across which colocation points are selected
    Id sample {Id()};      // Sampling mode (uniform, random), default is random. 

    // These two are for annotation purpuses only, they do not affect HB simulation. 
    IntVector harmonic {}; // When raw truncation scheme is used this flag indicates 
                           // a frequency in the freq vector is a harmonic. 
                           // If not given, assumes frequencies in freq are all harmonics. 
    IntVector imorder {};  // When raw truncation scheme is used this flag indicates 
                           // the intermodulation product order of each frequency in the freq vector. 
                           // If not set, assumes order is -1 for all frequencies. 
    
 
    Int writeOutput {1};   // Do we want to write the results to a file
                           // Not exposed as analysis parameter. 
                             
    HBParameters();
} HBParameters;


class HBCore : public AnalysisCore {
public:
    typedef HBParameters Parameters;
    enum class HBError {
        OK, 
        NoAlgorithm, 
        MatrixError, 
        SolverError,
        InitialHB, 
        Homotopy,
    };

    HBCore(
        OutputDescriptorResolver& parentResolver, HBParameters& params, Circuit& circuit, 
        KluBlockSparseRealMatrix& jacColoc, KluBlockSparseRealMatrix& jacobian, VectorRepository<double>& solution
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

    virtual bool storeState(size_t ndx, bool storeDetails=true);
    virtual bool restoreState(size_t ndx);
    
    virtual std::tuple<bool, bool> runSolver(bool continuePrevious);
    virtual Int iterations() const;
    virtual Int iterationLimit(bool continuePrevious) const;
        
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

    // Block-sparse matrix with rectangular blocks for 
    // storing Jacobian values at colocation timepoints
    KluBlockSparseRealMatrix& jacColoc;

    // HB Jacobian
    KluBlockSparseRealMatrix& bsjac; 
    VectorRepository<Real>& solution; // Solution history

    CoreStateStorage* continueState;

    OutputRawfile* outfile;
    
    bool converged_;
    
    struct SpecFreq {
        // Index of the frequency grid entry
        size_t gridIndex;
        // Abslute frequency
        double f;
        // If the grid entry results in a negative frequency, this is true
        bool negated;
        // Intermodulation product order
        // For harmonics this is the order of the harmonic. 
        int order;
        // Flag indicating that this frequency is a harmonic 
        // (i.e. all grid coordinates, but one, are 0)
        bool isHarmonic;
    };

private:
    NRSettings nrSettings;
    HBNRSolver nrSolver;

    // Temporary structures for collecting the phasors at a single frequency
    // before they are dumped. This vector has a bucket so that the output 
    // source code is the same as with other analyses. 
    VectorRepository<Complex> outputPhasors;
    Complex outputFreq;

    // Solution in frequency domain, nf phasors for each on of the n unknowns
    // NR solver resizes this vector. This vector has no bucket. 
    Vector<Complex> solutionFD;

    // Previous HB parameters to check if we need to rebuild()
    HBParameters oldParams;
    // Flag indicating rebuild() has not been called yet
    bool firstBuild;

    // HB parameters
    HBParameters& params;

    // Vectors and matrices without a bucket
    // Gridpoints in the frequency grid
    // Rows are gridpoints, columns are fundamental frequency factors
    DenseMatrix<int> grid;
    // Details on each frequency in the spectrum (sorted), references grid entries
    std::vector<SpecFreq> freq;
    // Vector of frequencies including DC (sorted)
    std::vector<double> frequencies;
    // Colocation timepoints (sorted)
    Vector<double> timepoints; 
    // Almost periodic Fourier transform
    DenseMatrix<double> APFT;
    // Inverse almost periodic Fourier transform
    DenseMatrix<double> IAPFT;
    // Derivative wrt time operator for time-domain vectors
    // Results in a time domain vector
    // Computed as IAPFT Omega APFT
    DenseMatrix<double> DDT;
    // DDT operator in column-major order (for cache locality)
    DenseMatrix<double> DDTcolMajor;
};

}

#endif
