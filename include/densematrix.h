#ifndef __DENSEMATRIX_DEFINED
#define __DENSEMATRIX_DEFINED

#include <stdexcept>
#include <cmath>
#include "common.h"

namespace NAMESPACE {

// A vector view into another vector/matrix
// Due to stride_ it can handle columns, rows, and much more
template<typename T> class VectorView {
public:
    VectorView(T* start, size_t n, size_t stride=1) : start_(start), n_(n), stride_(stride) {};
    
    T& at(size_t col) { return *(start_+col*stride_); };
    
    size_t n() const { return n_; };

    VectorView<T>& operator=(const VectorView<T>& from) { 
        if (n_ != from.n_) {
            throw std::out_of_range("Vector length mismatch.");
        }
        for(decltype(n_) i=0; i<n_; i++) {
            start_[i*stride_] = from.start_[i*stride_];
        }
        return *this;
    };

    VectorView<T>& operator=(const T& from) { 
        for(decltype(n_) i=0; i<n_; i++) {
            start_[i*stride_] = from;
        }
        return *this;
    };

    T& operator[](size_t i) { return *(start_+i*stride_); };

    // Squared norm
    double norm2() const {
        double nrm = 0;
        T* ptr = start_;
        for(size_t i=0; i<n_; i++) {
            auto a = std::abs(*ptr);
            nrm += a*a;
            ptr += stride_;
        }
        return nrm;
    };

    // Norm
    double norm() const { return std::sqrt(norm2()); };

    // Dot product
    // Conjugates other if T is complex, works only for double complex
    T dot(VectorView<T>& other) {
        if (n_ != other.n_) {
            throw std::out_of_range("Vector length mismatch.");
        }
        T sum = 0;
        T* ptr = start_;
        T* ptrOther = other.start_;
        for(size_t i=0; i<n_; i++) {
            if constexpr(std::is_same<T, Complex>::value) {
                sum += *ptr * std::conj(*ptrOther);
            } else {
                sum += *ptr * *ptrOther;
            }
            ptr += stride_;
            ptrOther += other.stride_;
        }
        return sum;
    };

    // Orthogonalize to wrt
    void orthogonalize(VectorView<T>& wrt) {
        // dot() checks vector compatibility
        auto prod = dot(wrt);
        auto nrm2 = wrt.norm2();
        auto fac = prod/nrm2;
        T* ptr = start_;
        T* ptrWrt = wrt.start_;
        for(size_t i=0; i<n_; i++) {
            *ptr -= fac * *ptrWrt;
            ptr += stride_;
            ptrWrt += wrt.stride_;
        }
    };

    // Swap with other, assume vectors have no common elements (e.g. row crossing a column)
    void swap(VectorView<T>& other) {
        if (n_ != other.n_) {
            throw std::out_of_range("Vector length mismatch.");
        }
        T* ptr = start_;
        T* ptrOther = other.start_;
        for(size_t i=0; i<n_; i++) {
            auto tmp = *ptr;
            *ptr = *ptrOther;
            *ptrOther = tmp;
            ptr += stride_;
            ptrOther += other.stride_;
        }
    };

    void swap(VectorView<T>&& other) {
        swap(other);
    }

private:
    T* start_;
    size_t n_;
    size_t stride_;
};


// Dense matrix stored in row-major order
template<typename T> class DenseMatrix {
public:
    DenseMatrix() : nRow_(0), nCol_(0) {};
    DenseMatrix(size_t nRow, size_t nCol) : nRow_(nRow), nCol_(nCol) { data.resize(nRow*nCol); };

    size_t nRow() const { return nRow_; }; 
    size_t nCol() const { return nCol_; }; 

    void resize(size_t nRow, size_t nCol) { nRow_ = nRow; nCol_ = nCol; data.resize(nRow*nCol); }; 

    T& at(size_t row, size_t col) { return data[row*nCol_+col]; }; 
    VectorView<T> row(size_t i) { return VectorView<T>(data.data()+nCol_*i, nCol_, 1); }
    VectorView<T> column(size_t i) { return VectorView<T>(data.data()+i, nRow_, nCol_); }

    void pushRow(T* rowDataPtr) { for(size_t i=0; i<nCol_; i++) data.push_back(rowDataPtr[i]); };
    VectorView<T> addRow() { nRow_++; data.resize(nRow_*nCol_); return row(nRow_-1); }

    void zero() { data.assign(data.size(), T()); };

private:
    size_t nRow_;
    size_t nCol_;
    std::vector<T> data;
};

}

#endif
