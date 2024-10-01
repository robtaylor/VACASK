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
    VectorView(std::vector<T>& v) : start_(v.data()), n_(v.size()), stride_(1) {};
    VectorView(std::vector<T>& v, size_t offset, size_t length, size_t stride) 
        : start_(v.data()+offset), n_(length), stride_(stride) {};
    VectorView(T* start, size_t n, size_t stride=1) : start_(start), n_(n), stride_(stride) {};
    
    // Access to members
    const T& at(size_t col) const { return *(start_+col*stride_); };
    T& at(size_t col) { return *(start_+col*stride_); };

    const T& operator[](size_t i) const { return *(start_+i*stride_); };
    T& operator[](size_t i) { return *(start_+i*stride_); };

    // Length
    size_t n() const { return n_; };

    // Assign elements from another VectorView
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

    // Assign same value to all elements
    VectorView<T>& operator=(const T& from) { 
        T* ptr = start_;
        for(decltype(n_) i=0; i<n_; i++) {
            *ptr = from;
            ptr += stride_;
        }
        return *this;
    };

    // Apply function to each element, put result in result
    void apply(T (*func)(T), VectorView<T>& result) const {
        for(size_t i=0; i<n_; i++) {
            result.at(i) = func(at(i));
        }
    };

    // Apply function to each element in place
    void apply(T (*func)(T)) {
        for(size_t i=0; i<n_; i++) {
            at(i) = func(at(i));
        }
    };

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
    T dot(const VectorView<T>& other) const {
        if (n_ != other.n_) {
            throw std::out_of_range("Vector length mismatch.");
        }
        T sum = 0;
        const T* ptr = start_;
        const T* ptrOther = other.start_;
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
    void orthogonalize(const VectorView<T>& wrt) {
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

    // Swap elements with other
    // Assume vectors have no common elements (e.g. row crossing a column)
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

    // Swap elements with other (other is a rvalue reference)
    // Assume vectors have no common elements (e.g. row crossing a column)
    void swap(VectorView<T>&& other) {
        swap(other);
    };

    // Swap i-th and j-th elements
    void swap(size_t i, size_t j) {
        T tmp = at(j);
        at(j) = at(i);
        at(i) = tmp;
    };

    // Scale by a factor
    void scale(T factor) {
        T* ptr = start_;
        for(size_t i=0; i<n_; i++) {
            *ptr *= factor;
            ptr += stride_;
        }
    };

    // Add scaled vector
    void addScaled(const VectorView<T>& other, T factor) {
        if (n_ != other.n_) {
            throw std::out_of_range("Vector length mismatch.");
        }
        T* ptr = start_;
        const T* ptrOther = other.start_;
        for(size_t i=0; i<n_; i++) {
            *ptr += *ptrOther * factor;
            ptr += stride_;
            ptrOther += other.stride_;
        }
    };

    // Add scaled vector
    void add(const VectorView<T>& other) {
        if (n_ != other.n_) {
            throw std::out_of_range("Vector length mismatch.");
        }
        T* ptr = start_;
        const T* ptrOther = other.start_;
        for(size_t i=0; i<n_; i++) {
            *ptr += *ptrOther;
            ptr += stride_;
            ptrOther += other.stride_;
        }
    };

    // Write scaled vector
    void writeScaled(const VectorView<T>& other, T factor) {
        if (n_ != other.n_) {
            throw std::out_of_range("Vector length mismatch.");
        }
        T* ptr = start_;
        const T* ptrOther = other.start_;
        for(size_t i=0; i<n_; i++) {
            *ptr = *ptrOther * factor;
            ptr += stride_;
            ptrOther += other.stride_;
        }
    };

    // Add scaled vector, store in result
    // Result can be *self
    // Assume result is not other
    void addScaled(const VectorView<T>& other, T factor, VectorView<T>& result) {
        if (n_ != other.n_) {
            throw std::out_of_range("Vector lengths do not match.");
        }
        if (n_ != result.n_) {
            throw std::out_of_range("Result length does not match vector.");
        }
        T* ptr = start_;
        const T* ptrOther = other.start_;
        T* ptrResult = result.start_;
        for(size_t i=0; i<n_; i++) {
            *ptrResult = *ptr + *ptrOther * factor;
            ptr += stride_;
            ptrOther += other.stride_;
            ptrResult += result.stride_;
        }
    };

    // Maximal absolute element
    double maxAbs() const {
        double m = 0;
        T* ptr = start_;
        for(size_t i=0; i<n_; i++) {
            auto c = std::abs(*ptr);
            if (c>m) {
                m = c;
            }
            ptr += stride_;
        }
        return m;
    };

    void dump(std::ostream& os) const {
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

template<typename T> std::ostream& operator<<(std::ostream& os, const VectorView<T>& obj) {
    for(decltype(obj.n()) i=0; i<obj.n(); i++) {
        if (i>0) {
            os << " ";
        }
        os << obj[i];
    }
    return os;
}

// Deduction guides for VectorView
template<typename T> VectorView(std::vector<T>& v) -> VectorView<T>;
template<typename T> VectorView(std::vector<T>& v, size_t offset, size_t length, size_t stride) -> VectorView<T>;
template<typename T> VectorView(T* start, size_t n, size_t stride=1) -> VectorView<T>;

template<typename T> class DenseMatrixView {
public:
    // Default constructor, uninitialized view
    DenseMatrixView() 
        : start_(nullptr), nRow_(0), nCol_(0), rowStride_(0), colStride_(0) {};

    // Construct from array
    DenseMatrixView(T* start, size_t nRow, size_t nCol, size_t rowStride, size_t colStride) 
        : start_(start), nRow_(nRow), nCol_(nCol), rowStride_(rowStride), colStride_(colStride) {}; 

    // Size
    size_t nRows() const { return nRow_; }; 
    size_t nCols() const { return nCol_; }; 

    // Element access
    const T& at(size_t row, size_t col) const { return *(start_ + row*rowStride_ + col*colStride_); };
    T& at(size_t row, size_t col) { return *(start_ + row*rowStride_ + col*colStride_); };
    
    // Row access
    const VectorView<T> row(size_t i) const { return VectorView<T>(start_ + i*rowStride_, nCol_, colStride_); };
    VectorView<T> row(size_t i) { return VectorView<T>(start_ + i*rowStride_, nCol_, colStride_); };
    
    // Column access
    const VectorView<T> column(size_t i) const { return VectorView<T>(start_ + i*colStride_, nRow_, rowStride_); };
    VectorView<T> column(size_t i) { return VectorView<T>(start_ + i*colStride_, nRow_, rowStride_); };

    // Assign elements from another MatrixView
    DenseMatrixView<T>& operator=(const DenseMatrixView<T>& other) {
        if (nRow_!=other.nRow_ || nCol_!=other.nCol_) {
            throw std::out_of_range("Matrices do not match.");
        }
        for(size_t i=0; i<nRow_; i++) {
            row(i) = other.row(i);
        }
        return *this;
    };

    // Assign value to all elements
    DenseMatrixView<T>& operator=(const T& val) {
        for(size_t i=0; i<nRow_; i++) {
            row(i) = val;
        }
        return *this;
    };

    // Set to zero
    void zero() { 
        *this = 0;
    };

    // Set to identity
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

    // Solve Ax = rhs, destroy A, result in rhs
    // Use partial pivoting
    bool destructiveSolve(VectorView<T>& rhs) {
        return solveCore(&rhs, nullptr);
    };

    // Solve Ax = Rhs, destroy A, result in Rhs
    // Use partial pivoting
    bool destructiveSolve(DenseMatrixView<T>& rhs) {
        return solveCore(&rhs, nullptr);
    };

    // Perform LU decomposition in place, return row permutation vector
    // Use partial pivoting
    bool factor(VectorView<size_t>& rowPerm) {
        return solveCore(static_cast<VectorView<T>*>(nullptr), &rowPerm);
    };

    // Multiply with vector, store result in result
    // Vector must be distinct from result
    void multiply(const VectorView<T>& vector, VectorView<T>& result) const {
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

    // Multiply with matrix, store result in result
    // Result must be distinct from this and other
    void multiply(const DenseMatrixView<T>& other, DenseMatrixView<T>& result) const {
        if (nCol_!=other.nRow_) {
            throw std::out_of_range("Matrices are not compatible.");
        }
        if (nRow_!=result.nRow_ || other.nCol_!=result.nCol_) {
            throw std::out_of_range("Result is not compatible with product.");
        }
        for(size_t i=0; i<nRow_; i++) {
            for(size_t j=0; j<other.nCol_; j++) {
                result.at(i, j) = row(i).dot(other.column(j));
            }
        }
    };

    // Add scaled other matrix, put result in result
    // Result must be distinct from this and other
    void addScaled(const DenseMatrixView<T>& other, T factor, DenseMatrixView<T>& result) {
        if (nRow_!=other.nRow_ || nCol_!=other.nCol_) {
            throw std::out_of_range("Matrices are not compatible.");
        }
        if (nRow_!=result.nRow_ || nCol_!=result.nCol_) {
            throw std::out_of_range("Result is not compatible with matrix.");
        }
        for(size_t i=0; i<nRow_; i++) {
            auto rrow = result.row(i);
            row(i).addScaled(other.row(i), factor, rrow);
        }
    };

    // Add other matrix, put result in result
    // Result must be distinct from this and other
    void add(DenseMatrixView<T>& other, DenseMatrixView<T>& result) {
        addScaled(other, 1.0, result);
    };

    // Subtract other matrix, put result in result
    // Result must be distinct from this and other
    void subtract(DenseMatrixView<T>& other, DenseMatrixView<T>& result) {
        addScaled(other, -1.0, result);
    };

    // Apply function to each element, put result in result
    void apply(T (*func)(T), DenseMatrixView<T>& result) const {
        for(size_t i=0; i<nRow_; i++) {
            auto rrow = result.row(i);
            row(i).apply(func, rrow);
        }
    };

    // Apply function to each element in place
    void apply(T (*func)(T), DenseMatrixView<T>& result) {
        for(size_t i=0; i<nRow_; i++) {
            row(i).apply(func);
        }
    };

    // Absolute maximal element
    double maxAbs() const {
        double m = 0;
        for(size_t i=0; i<nRow_; i++) {
            auto c = row(i).maxAbs();
            if (c>m) {
                m = c;
            }
        }
        return m;
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

    void dump(std::ostream& os) const {
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
            auto p = std::abs(pivCol[pivI]);
            for(size_t j=i; j<n; j++) {
                auto cand = std::abs(pivCol[j]);
                if (cand>p) {
                    p = cand;
                    pivI = j;
                }
            }
            auto pivot = pivCol[pivI];
            if (p==0) {
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
            for(size_t cnt=0; cnt<n; cnt++) {
                auto i = n-1-cnt;
                auto matRow = row(i);
                if constexpr(std::is_same<RhsType, VectorView<T>>::value) {
                    for(size_t j=i+1; j<n; j++) {
                        rhs->at(i) -= rhs->at(j)*matRow[j];
                    }
                    rhs->at(i) /= matRow[i];
                } else {
                    auto rhsRowi = rhs->row(i);
                    for(size_t j=i+1; j<n; j++) {
                        auto rhsRowj = rhs->row(j);
                        auto fac = -matRow[j];
                        rhsRowi.addScaled(rhsRowj, fac);
                    }
                    rhsRowi.scale(1.0/matRow[i]);
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
    
    // Default constructor, uninitialized matrix
    DenseMatrix() 
        : DenseMatrixView<T>(nullptr, 0, 0, 1, 1), major_(Major::Row) {};

    // Copy constructor
    DenseMatrix(const DenseMatrix<T>& A) {
        major_ = A.major_;
        data_ = A.data_;
        start_ = data_.data();
        nRow_ = A.nRow_;
        nCol_ = A.nCol_;
        setStride();
    };

    // Move constructor
    DenseMatrix(DenseMatrix<T>&& A) {
        major_ = A.major_;
        data_ = std::move(A.data_);
        start_ = data_.data();
        nRow_ = A.nRow_;
        nCol_ = A.nCol_;
        setStride();
    };

    // Size-based constructor
    DenseMatrix(size_t nRow, size_t nCol, Major major=Major::Row) 
        : DenseMatrixView<T>(nullptr, nRow, nCol, 1, 1), major_(major) { 
        data_.resize(nRow*nCol); 
        start_ = data_.data();
        setStride();
    };

    // Move-construct from vector and size
    DenseMatrix(std::vector<T>&& from, size_t nRow, size_t nCol, Major major=Major::Row) { 
        if (nRow*nCol != from.size()) {
            throw std::out_of_range("Matrix size inconsistent with data.");
        }
        major_ = major;
        data_ = std::move(from);
        start_ = data_.data();
        nRow_ = nRow;
        nCol_ = nCol;
        setStride();
    };

    // Copy-construct from vector and size
    DenseMatrix(const std::vector<T>& from, size_t nRow, size_t nCol, Major major=Major::Row) { 
        if (nRow*nCol != from.size()) {
            throw std::out_of_range("Matrix size inconsistent with data.");
        }
        major_ = major;
        data_ = from;
        start_ = data_.data();
        nRow_ = nRow;
        nCol_ = nCol;
        setStride();
    };

    // Copy assignment
    DenseMatrix<T>& operator=(const DenseMatrix<T>& other) {
        major_ = other.major_;
        data_ = other.data_;
        start_ = data_.data();
        nRow_ = other.nRow_;
        nCol_ = other.nCol_;
        setStride();
        return *this;
    };

    // Move assignment
    DenseMatrix<T>& operator=(DenseMatrix<T>&& other) {
        major_ = other.major_;
        data_ = other.data_;
        start_ = std::move(data_.data());
        nRow_ = other.nRow_;
        nCol_ = other.nCol_;
        setStride();
        return *this;
    };
    
    // Resize, does not reorder elements (content is invalidated)
    void resize(size_t nRow, size_t nCol, Major major=Major::Row) { 
        major_ = major;
        data_.resize(nRow*nCol);
        nRow_ = nRow; 
        nCol_ = nCol; 
        start_ = data_.data();
        setStride();
    }; 
    
    // Override for DenseMatrix
    T& at(size_t row, size_t col) { 
        switch (major_) {
            case Major::Row:
                return data_[row*nCol_+col]; 
            case Major::Column:
            default:
                return data_[row+nRow_*col]; 
        }
    }; 
    
    // Override for DenseMatrix
    VectorView<T> row(size_t i) { 
        switch (major_) {
            case Major::Row:
                return VectorView<T>(data_.data()+nCol_*i, nCol_, 1); 
            case Major::Column:
            default:
                return VectorView<T>(data_.data()+i, nCol_, nRow_); 
        }   
    };
    
    // Override for DenseMatrix
    VectorView<T> column(size_t i) { 
        switch (major_) {
            case Major::Row:
                return VectorView<T>(data_.data()+i, nRow_, nCol_); 
            case Major::Column:
            default:
                return VectorView<T>(data_.data()+i*nRow_, nRow_, 1); 
        }
    };

    // Override for DenseMatrix
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

    static bool test();

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
