#ifndef __COREHB_DEFINED
#define __COREHB_DEFINED

#include "core.h"
#include "value.h"
#include "common.h"

namespace NAMESPACE {

template<typename T> class DenseRow {
public:
    DenseRow(T* start, size_t nCol) : start(start), nCol(nCol) {};
    
    T& at(size_t col) { return *(start+col); };
    DenseRow<T> advance(size_t n) { return DenseRow<T>(start+(n*nCol), nCol); };

    size_t count() const { return nCol; };

private:
    T* start;
    size_t nCol;
};


template<typename T> class DenseMatrix {
public:
    DenseMatrix() : nRow_(0), nCol_() {};
    DenseMatrix(size_t nRow, size_t nCol) : nRow_(nRow), nCol_(nCol) { data.resize(nRow*nCol); };

    size_t nRow() const { return nRow_; }; 
    size_t nCol() const { return nCol_; }; 

    void resize(size_t nRow, size_t nCol) { nRow_ = nRow; nCol_ = nCol; data.resize(nRow*nCol); }; 

    T& at(size_t row, size_t col) { return data[row*nCol_+col]; }; 
    DenseRow<T> row(size_t i) { return DenseRow<T>(data.data()+nCol_*i, nCol_); }

    void pushRowData(T* rowDataPtr) { for(size_t i=0; i<nCol_; i++) data.push_back(rowDataPtr[i]); };
    DenseRow<T> addRow() { nRow_++; data.resize(nRow_*nCol_); return row(nRow_-1); }

    void zero() { data.assign(data.size(), T()); };

private:
    size_t nRow_;
    size_t nCol_;
    std::vector<T> data;
};


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
    Real samplfac {2};    // Sampling factor in time domain (>=1). 
    Id samplmode {Id()};  // Sampling mode (uniform, random), default is random. 
                             
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

protected:
    
    bool buildCollocation(Status& s=Status::ignore);
    bool buildTransform(Status& s=Status::ignore);

    struct SpecFreq {
        size_t gridIndex;
        double f;
        int order;
        bool isHarmonic;
    };

    static Id truncateRaw;
    static Id truncateBox;
    static Id truncateDiamond;

private:
    HbParameters& params;
    Vector<Real> spectrum;
    std::vector<SpecFreq> freq;
    DenseMatrix<int> grid;
    Vector<Real> timepoints;    
};

}

#endif
