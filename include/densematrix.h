#ifndef __DENSEMATRIX_DEFINED
#define __DENSEMATRIX_DEFINED

#include <stdexcept>
#include <vector>
#include <cmath>
#include <optional>
#include <type_traits>
#include <iostream>
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
    };

    void swap(size_t i, size_t j) {
        T tmp = at(j);
        at(j) = at(i);
        at(i) = tmp;
    };

    void scale(T factor) {
        T* ptr = start_;
        for(size_t i=0; i<n_; i++) {
            *ptr *= factor;
            ptr += stride_;
        }
    };

    // Add scaled vector, assume vectors have no common components
    void addScaled(VectorView<T>& other, T factor) {
        if (n_ != other.n_) {
            throw std::out_of_range("Vector length mismatch.");
        }
        T* ptr = start_;
        T* ptrOther = other.start_;
        for(size_t i=0; i<n_; i++) {
            *ptr += *ptrOther * factor;
            ptr += stride_;
            ptrOther += other.stride_;
        }
    };

    void dump(std::ostream& os)  {
        T* ptr = start_;
        for(size_t i=0; i<n_; i++) {
            os << *ptr << " ";
            ptr += stride_;
        }
        os << "\n";
    };
    
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
        : start_(start), nRow_(nRow), nCol_(nCol), rowStride_(rowStride), colStride_(colStride) {}; 

    size_t nRows() const { return nRow_; }; 
    size_t nCols() const { return nCol_; }; 

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

    void identity() {
        for(size_t i=0; i<nRow_; i++) {
            for(size_t j=0; j<nCol_; j++) {
                if (i==j) {
                    at(i, j) = 1;
                } else {
                    at(i, j) = 0;
                }
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

    bool destructiveSolve(VectorView<T>& rhs) {
        return solveCore(&rhs, nullptr);
    };

    bool destructiveSolve(DenseMatrixView<T>& rhs) {
        return solveCore(&rhs, nullptr);
    };

    bool factor(VectorView<size_t>& rowPerm) {
        return solveCore(static_cast<VectorView<T>*>(nullptr), &rowPerm);
    };

    // Vector must be distinct from result
    void multiply(VectorView<T>& vector, VectorView<T>& result) {
        if (nCol_!=vector.n()) {
            throw std::out_of_range("Matrix is not compatible with vector.");
        }
        if (nRow_!=result.n()) {
            throw std::out_of_range("Result is not compatible with product.");
        }
        for(size_t i=0; i<nRow_; i++) {
            result[i] = row(i).dot(vector);
        }
    };

    // Result must be distinct from this and other
    void multiply(DenseMatrixView<T>& other, DenseMatrixView<T>& result) {
        if (nCol_!=other.nRow_) {
            throw std::out_of_range("Matrices are not compatible.");
        }
        if (nRow_!=result.nRow_ || other.nCol_!=result.nCol_) {
            throw std::out_of_range("Result is not compatible with product.");
        }
        for(size_t i=0; i<nRow_; i++) {
            for(size_t j=0; j<other.nCol_; j++) {
                auto otherCol = other.column(j);
                result.at(i, j) = row(i).dot(otherCol);
            }
        }
    };

    // Destructive invert, result must be distinct from this
    bool destructiveInvert(DenseMatrixView<T>& result) {
        if (nRow_!=result.nRow_ || nCol_!=result.nCol_) {
            throw std::out_of_range("Matrices are not compatible.");
        }
        if (nRow_!=nCol_) {
            throw std::out_of_range("Matrix is not square.");
        }
        result.identity();
        return solveCore(&result);
    };

    void dump(std::ostream& os)  {
        for(size_t i=0; i<nRow_; i++) {
            for(size_t j=0; j<nCol_; j++) {
                os << at(i, j) << " ";
            }
            os << "\n";
        }
    };
    
protected:
    T* start_;
    size_t nRow_;
    size_t nCol_;
    size_t rowStride_;
    size_t colStride_;

private:
    // Destructive solve/factor core 
    // Replaces matrix content with LU decomposition. 
    // Solution is placed in rhs. 
    // Stores row permutation vector. 
    // Returns true on success. 
    template<typename RhsType> bool solveCore(RhsType* rhs, VectorView<size_t>* rowPerm = nullptr) {
        auto n = nCol_;
        if constexpr(std::is_same<RhsType, VectorView<T>>::value) {
            if (rhs && rhs->n()!=n) {
                throw std::out_of_range("Vector length does not match matrix size.");
            }
        } else if constexpr(std::is_same<RhsType, DenseMatrixView<T>>::value) {
            if (rhs && rhs->nRows()!=n) {
                throw std::out_of_range("Vector length does not match matrix size.");
            }
        } else {
            if (rhs) {
                throw std::out_of_range("Bad rhs type.");
            }
        }
        if (rowPerm && rowPerm->n()!=n) {
            throw std::out_of_range("Row permutation vector length does not match matrix size.");
        }
        if (nRow_!=n) {
            throw std::out_of_range("Matrix is not square.");
        }
        if (rowPerm) {
            // Initialize row permutation
            for(size_t i=0; i<n; i++) {
                (*rowPerm)[i] = i;
            }
        }
        
        // Eliminate
        for(size_t i=0; i<n-1; i++) {
            // Find pivot
            auto pivCol = column(i);
            size_t pivI = i;
            T pivot = pivCol[pivI];
            for(size_t j=i; j<n; j++) {
                auto cand = std::abs(pivCol[j]);
                if (cand>pivot) {
                    pivot = cand;
                    pivI = j;
                }
            }
            pivot = pivCol[pivI];
            if (pivot==0) {
                return false;
            }
            
            // Swap with pivot
            if (pivI!=i) {
                auto pivRow = row(pivI);
                auto pivDest = row(i);
                pivRow.swap(pivDest);
                if (rowPerm) {
                    rowPerm->swap(i, pivI);
                }
                if (rhs) {
                    if constexpr(std::is_same<RhsType, VectorView<T>>::value) {
                        rhs->swap(i, pivI);
                    } else {
                        rhs->row(i).swap(rhs->row(pivI));
                    }
                }
            }

            // Eliminate in column i, rows i+1..n-1
            for(size_t j=i+1; j<n; j++) {
                auto pivRow = row(i);
                auto targetRow = row(j);
                T fac = -targetRow[i]/pivRow[i];
                for(size_t k=i+1; k<n; k++) {
                    targetRow[k] += fac*pivRow[k];
                }
                targetRow[i] = -fac;
                if (rhs) {
                    if constexpr(std::is_same<RhsType, VectorView<T>>::value) {
                        rhs->at(j) += fac*rhs->at(i);
                    } else {
                        auto rhsRow = rhs->row(i);
                        rhs->row(j).addScaled(rhsRow, fac);
                    }
                }
            }
        }

        // Back-substitute
        if (rhs) {
            for(size_t i=0; i<n; i++) {
                auto irow = n-1-i;
                auto matRow = row(irow);
                if constexpr(std::is_same<RhsType, VectorView<T>>::value) {
                    for(size_t j=irow+1; j<n; j++) {
                        rhs->at(irow) -= rhs->at(j)*matRow[j];
                    }
                    rhs->at(irow) /= matRow[irow];
                } else {
                    auto rhsRowi = rhs->row(i);
                    for(size_t j=irow+1; j<n; j++) {
                        auto rhsRowj = rhs->row(j);
                        auto fac = -matRow[j];
                        rhsRowi.addScaled(rhsRowj, fac);
                    }
                    rhsRowi.scale(1.0/matRow[irow]);
                }
            }
        }
        return true;
    };
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
        : DenseMatrixView<T>(nullptr, nRow, nCol, 1, 1), major_(major) { 
        data_.resize(nRow*nCol); 
        start_ = data_.data();
        setStride();
    };

    DenseMatrix(std::vector<T>&& from, size_t nRow, size_t nCol, Major major=Major::Row) { 
        if (nRow*nCol != from.size()) {
            throw std::out_of_range("Matrix size inconsistent with data.");
        }
        data_ = std::move(from);
        start_ = data_.data();
        nRow_ = nRow;
        nCol_ = nCol;
        setStride();
    };
    
    // Does not reorder existing data, data moves in mysterious ways index-wise
    void resize(size_t nRow, size_t nCol) { 
        data_.resize(nRow*nCol);
        nRow_ = nRow; 
        nCol_ = nCol; 
        start_ = data_.data();
        setStride();
    }; 
    
    // Specializations for DenseMatrix
    T& at(size_t row, size_t col) { 
        switch (major_) {
            case Major::Row:
                return data_[row*nCol_+col]; 
            case Major::Column:
            default:
                return data_[row+nRow_*col]; 
        }
    }; 
    
    // Specializations for DenseMatrix
    VectorView<T> row(size_t i) { 
        switch (major_) {
            case Major::Row:
                return VectorView<T>(data_.data()+nCol_*i, nCol_, 1); 
            case Major::Column:
            default:
                return VectorView<T>(data_.data()+i, nCol_, nRow_); 
        }   
    };
    
    // Specializations for DenseMatrix
    VectorView<T> column(size_t i) { 
        switch (major_) {
            case Major::Row:
                return VectorView<T>(data_.data()+i, nRow_, nCol_); 
            case Major::Column:
            default:
                return VectorView<T>(data_.data()+i*nRow_, nRow_, 1); 
        }
    };

    // Specializations for DenseMatrix
    void zero() { data_.assign(data_.size(), T()); };
    
    // Row major only
    VectorView<T> addRow() { 
        if (major_==Major::Column) {
            throw std::out_of_range("Rows cannot be added to column major matrices.");
        }
        nRow_++; 
        data_.resize(nRow_*nCol_); 
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
    std::vector<T> data_;
};

}

#endif
