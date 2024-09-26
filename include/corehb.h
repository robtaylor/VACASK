#ifndef __COREHB_DEFINED
#define __COREHB_DEFINED

#include "densematrix.h"
#include "core.h"
#include "value.h"
#include "common.h"

namespace NAMESPACE {

typedef struct HbParameters {
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
                          // raw     .. values in freq are the values in the spectrum
                          // box     .. box truncation
                          //   kj = 0..Hj, first nonzero kj must be >0
                          // diamond .. diamond truncation (default)
                          //   sum abs(kj) <= IM, first nonzero kj must be >0
    Real samplefac {2};   // Sampling factor in time domain (>=1). 
    Real nper {3};        // Number of periods across which colocation points are selected
    Id sample {Id()};     // Sampling mode (uniform, random), default is random. 
                             
    HbParameters();
} HbParameters;


class HbCore /*: public AnalysisCore*/ {
public:
    typedef HbParameters Parameters;
    enum class HbError {
        OK, 
        StatusError, 
    };

    HbCore(HbParameters& params) : params(params) {};

    bool buildGrid(Status& s=Status::ignore);
    bool buildColocation(Status& s=Status::ignore);
    bool buildAPFT(Status& s=Status::ignore);
    
    // Recompute the spectrum, check if its length changed
    // Return value: Ok, requesting rebuild
    std::tuple<bool, bool> requestsRebuild(Status& s = Status::ignore);
    
    static Id truncateRaw;
    static Id truncateBox;
    static Id truncateDiamond;
    static Id sampleUniform;
    static Id sampleRandom;

    static bool test();

protected:
    bool buildTransformMatrix(DenseMatrix<double>& XF, Status& s=Status::ignore);
    
    struct SpecFreq {
        size_t gridIndex;
        double f;
        int order;
        bool isHarmonic;
    };

private:
    HbParameters& params;
    Vector<double> spectrum;
    std::vector<SpecFreq> freq;
    DenseMatrix<int> grid;
    Vector<double> timepoints;    
    DenseMatrix<double> APFT;
    DenseMatrix<double> IAPFT;
    DenseMatrix<double> DDT;
    DenseMatrix<double> DDTcolMajor;
};

}

#endif
