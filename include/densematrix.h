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
        T* ptr = start_;
        T* ptrFrom = from.start_;
        for(decltype(n_) i=0; i<n_; i++) {
            *ptr = *ptrFrom;
            ptr += stride_;
            ptrFrom += from.stride_;
        }
        return *this;
    };

    VectorView<T>& operator=(const T& from) { 
        if (n_ != from.n_) {
            throw std::out_of_range("Vector length mismatch.");
        }
        T* ptr = start_;
        T* ptrFrom = from.start_;
        for(decltype(n_) i=0; i<n_; i++) {
            *ptr = *ptrFrom;
            ptr += stride_;
            ptrFrom += from.stride_;
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


template<typename T> class DenseMatrixView {
public:
    DenseMatrixView() 
        : start_(nullptr), nRow_(0), nCol_(0), rowStride_(0), colStride_(0) {};
    DenseMatrixView(T* start, size_t nRow, size_t nCol, size_t rowStride, size_t colStride) 
        : start_(start), nRow_(nRow), nCol_(nCol), rowStride_(rowStride_), colStride_(colStride) {}; 

    size_t nRow() const { return nRow_; }; 
    size_t nCol() const { return nCol_; }; 

    T& at(size_t row, size_t col) { return *(start_ + row*rowStride_ + col*colStride_); };
    VectorView<T> row(size_t i) { return VectorView<T>(start_ + i*rowStride_, nCol_, colStride_); };
    VectorView<T> column(size_t i) { return VectorView<T>(start_ + i*colStride_, nRow_, rowStride_); };

    void zero() { 
        for(size_t i=0; i<nRow_; i++) {
            for(size_t j=0; j<nCol_; j++) {
                at(i, j) = 0;
            }
        }
    };

    // TODO: optimize this further
    void vecMul(VectorView<T>& other, VectorView<T>& result) {
        for(size_t i=0; i<nRow_; i++) {
            auto x = row(i).dot(other);
            result[i] = x;
        }
    };
    
protected:
    T* start_;
    size_t nRow_;
    size_t nCol_;
    size_t rowStride_;
    size_t colStride_;
};


// Dense matrix stored in row-major order
template<typename T> class DenseMatrix : public DenseMatrixView<T> {
public:
    enum class Major { Row=0, Column=1 }; 

    using DenseMatrixView<T>::start_;
    using DenseMatrixView<T>::nRow_;
    using DenseMatrixView<T>::nCol_;
    using DenseMatrixView<T>::rowStride_;
    using DenseMatrixView<T>::colStride_;
    
    DenseMatrix() 
        : DenseMatrixView<T>(nullptr, 0, 0, 1, 1), major_(Major::Row) {};

    DenseMatrix(size_t nRow, size_t nCol, Major major=Major::Row) 
        : DenseMatrixView<T>(nullptr, nRow_, nCol_, 1, 1), major_(major) { 
        data.resize(nRow*nCol); 
        start_ = data.data();
        setStride();
    };

    // Does not reorder existing data, data moves in mysterious ways index-wise
    void resize(size_t nRow, size_t nCol) { 
        nRow_ = nRow; 
        nCol_ = nCol; 
        data.resize(nRow*nCol);
        setStride();
    }; 
    
    // Specializations for DenseMatrix
    T& at(size_t row, size_t col) { 
        switch (major_) {
            case Major::Row:
                return data[row*nCol_+col]; 
            case Major::Column:
            default:
                return data[row+nRow_*col]; 
        }
    }; 
    
    // Specializations for DenseMatrix
    VectorView<T> row(size_t i) { 
        switch (major_) {
            case Major::Row:
                return VectorView<T>(data.data()+nCol_*i, nCol_, 1); 
            case Major::Column:
            default:
                return VectorView<T>(data.data()+i, nCol_, nRow_); 
        }   
    };
    
    // Specializations for DenseMatrix
    VectorView<T> column(size_t i) { 
        switch (major_) {
            case Major::Row:
                return VectorView<T>(data.data()+i, nRow_, nCol_); 
            case Major::Column:
            default:
                return VectorView<T>(data.data()+i*nRow_, nRow_, 1); 
        }
    };

    // Specializations for DenseMatrix
    void zero() { data.assign(data.size(), T()); };
    
    // Row major only
    VectorView<T> addRow() { 
        if (major_==Major::Column) {
            throw std::out_of_range("Rows cannot be added to column major matrices.");
        }
        nRow_++; 
        data.resize(nRow_*nCol_); 
        return row(nRow_-1); 
    };
    
private:
    void setStride() {
        switch (major_) {
            case Major::Row:
                rowStride_ = nCol_;
                colStride_ = 1;
                break;
            case Major::Column:
                rowStride_ = 1;
                colStride_ = nRow_;
                break;
        }
    };

    Major major_;
    std::vector<T> data;
};

}

#endif
