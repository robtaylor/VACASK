#include "core.h"
#include "common.h"

namespace NAMESPACE {


typedef struct HbParameters {
    RealVector freq {};   // Fundamental frequencies (f1, f2, ..., fd)
    Value nharm {4};      // Number of harmonics for each fundamental frequency 
                          // (H1, H2, ..., Hd) in box truncation
                          // If scalar, applies to all frequencies, 
                          // if vector, components apply to corresponding frequencies. 
    Int immax {0};        // Maximal order of intermodulation products (IM) in diamond truncation. 
                          // If <=0, defaults to harmonics. 
                          // If harmonics is a vector, defaults to its largest component. 
    IntVector imorder {}; // Intermodulation product order 
                          // (reported for each frequency when raw truncation scheme is used)
    Id truncate {Id()};   // Truncation scheme: 
                          // raw     .. values in freq are the values in the spectrum
                          // box     .. box truncation
                          //   kj = 0..Hj, first nonzero kj must be >0
                          // diamond .. diamond truncation (default)
                          //   sum abs(kj) <= IM, first nonzero kj must be >0
    Real samplfac {2};    // Sampling factor in time domain (>=1). 
    Id samplmode {Id()};  // Sampling mode (uniform, random), default is random. 
                             
    HbParameters();
} HbParameters;


class CoreHb : public AnalysisCore {
public:
    typedef HbParameters Parameters;
    enum class HbError {
        OK, 
    };

protected:
    bool buildSpectrum(Status& s=Status::ignore);
    bool buildCollocation(Status& s=Status::ignore);
    bool buildTransform(Status& s=Status::ignore);

private:
    HbParameters& params;
    Vector<Real> spectrum;
    Vector<Int> intmodOrder;
    Vector<Real> timepoints;
    
};

}