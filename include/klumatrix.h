#ifndef __KLUMATRIX_DEFINED
#define __KLUMATRIX_DEFINED

#include <suitesparse/klu.h>
#include <unordered_map>
#include <complex>
#include <type_traits>
#include <optional>
#include "status.h"
#include "identifier.h"
#include "flags.h"
#include "hash.h"
#include "acct.h"
#include "common.h"


namespace NAMESPACE {

// TODO: make multiple matrices share the same sparsity pattern without allocating a new copy
//       of KLU sparsity pattern

// Abstract class that resolves a row/column index into a name
// Should be defined by the user of the sparse matrix
class NameResolver {
public:
    NameResolver() {};

    // Indices in matrix are 0-based, operator() must take this into account
    virtual Id operator()(MatrixEntryIndex u) = 0;
};


// Index pair holding row and column position
typedef std::pair<EquationIndex,UnknownIndex> MatrixEntryPosition;


// Hash function for MatrixEntryPosition
 typedef struct MatrixEntryPositionHash {
    auto operator()(const MatrixEntryPosition& p) const -> size_t {
        if constexpr(sizeof(MatrixEntryPosition)<=sizeof(size_t)) {
            // Faster hash if entry coordinates are 32+32 bits on a machine with a 64-bit size_t
            return std::hash<size_t>{}(
                (static_cast<size_t>(p.first) << (sizeof(size_t)*8/2)) +
                p.second
            );
        } else {
            // Standard approach
            return hash_val(p.first, p.second);
        }
    }
} MatrixEntryPositionHash;


// Sparsity map - maps MatrixEntyPosition to an index in a linear array
class SparsityMap {
public:
    SparsityMap() {};

    SparsityMap           (const SparsityMap&)  = delete;
    SparsityMap           (      SparsityMap&&) = delete;
    SparsityMap& operator=(const SparsityMap&)  = delete;
    SparsityMap& operator=(      SparsityMap&&) = delete;

    // Type of map from MatrixEntryPosition into a linear array index
    typedef std::unordered_map<MatrixEntryPosition, MatrixEntryIndex, MatrixEntryPositionHash> Map;

    // Clear
    void clear();
    
    // Size
    MatrixEntryIndex size() const { return smap.size(); };

    // Insert, the returned pointer points to an integer that 
    // will containt the entry index after enumerate() is called
    // bool value indicates if a pre-existing entry was found
    std::tuple<MatrixEntryIndex*, bool> insert(EquationIndex e, UnknownIndex u) {
        auto [it, inserted] = smap.insert({std::make_pair(e, u), 0});
        return std::make_tuple(&(it->second), true);
    };

    // Find, bool value indicates if an entry was found
    std::tuple<MatrixEntryIndex, bool> find(const MatrixEntryPosition& mep) const {
        auto it = smap.find(mep);
        if (it==smap.end()) {
            return std::make_tuple(0, false);
        }
        return std::make_tuple(it->second, true);
    };
    
    // Get vector of sorted matrix entry positions
    std::vector<MatrixEntryPosition>& positions() { return ordering; };
    const std::vector<MatrixEntryPosition>& positions() const { return ordering; };

    // Enumerate entries
    void enumerate();

    // Dump map
    void dump(int indent, std::ostream& os) const;

private:
    Map smap;
    std::vector<MatrixEntryPosition> ordering;
};


// Element component when element is complex
enum class Component { Real=1, Imaginary=2 };
DEFINE_FLAG_OPERATORS(Component);


// Matrix binding interface for accessing element indices and pointers
// Because Circuit::bind() should by matrix type agnostic we need this interface. 
// Analyses can use different types of matrices, but instances must be able to 
// handle them all in the same way via this interface. 
// Assumes the underlying type of element is either double or std::complex<double> (Complex)
// This interface is used by the Device::bind() method to bind instances 
// to matrix elements and their components. 
template<typename IndexType> class MatrixAccess {
public:
    // Return array holding matrix nonzero elements
    // For complex matrices returns nullptr
    virtual double* valueArray() = 0;

    // Return array holding matrix nonzero elements
    // For real matrices returns nullptr
    virtual Complex* cxValueArray() = 0;

    // Return index into element array corresponding to row, column given by mep (1-based)
    // For block matrices, mep is the block position (1-based) and 
    // blockMep is the position of the element within the block (1-based). 
    // If blockMep is not given, (1,1) is assumed. 
    // Return value: index, found
    virtual std::tuple<IndexType, bool> valueIndex(const MatrixEntryPosition& mep, const std::optional<MatrixEntryPosition>& blockMep=std::nullopt) const = 0;

    // Return pointer to element's component
    // For block matrices, mep is the block position (1-based) and 
    // blockMep is the position of the element within the block (1-based). 
    // If blockMep is not given, (1,1) is assumed. 
    // Returns bucket if element is not found 
    // Returns nullptr if imaginary part is requested from a real matrix
    virtual double* valuePtr(const MatrixEntryPosition& mep, Component comp=Component::Real, const std::optional<MatrixEntryPosition>& blockMep=std::nullopt) = 0;

    // Return pointer to element's component (complex matrix)
    // For block matrices, mep is the block position (1-based) and 
    // blockMep is the position of the element within the block (1-based). 
    // If blockMep is not given, (1,1) is assumed. 
    // Returns complex bucket if element is not found 
    // Returns nullptr if matrix is real 
    virtual Complex* cxValuePtr(const MatrixEntryPosition& mep, const std::optional<MatrixEntryPosition>& blockMep=std::nullopt) = 0;
};


template<typename IndexType, typename ValueType> class KluMatrixCore : public MatrixAccess<IndexType> {
public: 
    using Common = typename std::conditional<std::is_same<int32_t,IndexType>::value, klu_common, klu_l_common>::type;
    using Symbolic = typename std::conditional<std::is_same<int32_t,IndexType>::value, klu_symbolic, klu_l_symbolic>::type;
    using Numeric = typename std::conditional<std::is_same<int32_t,IndexType>::value, klu_numeric, klu_l_numeric>::type;

    enum class Error {
        OK, 
        Memory, 
        Defaults, 
        Analysis, 
        Factorization, 
        Refactorization, 
        ReciprocalPivotGrowth, 
        ReciprocalCondEstimate, 
        MatrixInfNan, 
        VectorInfNan, 
        Solve
    };
    
    KluMatrixCore();
    
    KluMatrixCore           (const KluMatrixCore&)  = delete;
    KluMatrixCore           (      KluMatrixCore&&) = delete;
    KluMatrixCore& operator=(const KluMatrixCore&)  = delete;
    KluMatrixCore& operator=(      KluMatrixCore&&) = delete;

    virtual ~KluMatrixCore();

    // Matrix binding interface
    // Block element position is ignored
    virtual double* valueArray();
    virtual Complex* cxValueArray();
    virtual std::tuple<IndexType, bool> valueIndex(const MatrixEntryPosition& mep, const std::optional<MatrixEntryPosition>& blockMep=std::nullopt) const;
    virtual double* valuePtr(const MatrixEntryPosition& mep, Component comp=Component::Real, const std::optional<MatrixEntryPosition>& blockMep=std::nullopt) ;
    virtual Complex* cxValuePtr(const MatrixEntryPosition& mep, const std::optional<MatrixEntryPosition>& blockMep=std::nullopt);

    // Set accounting structure
    void setAccounting(Accounting& accounting) { acct = &accounting; }; 

    // Turn off accounting
    void noAccounting() { acct = nullptr; };

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore, NameResolver* resolver=nullptr) const; 

    // Error element
    std::tuple<IndexType, IndexType> errorElement() const {
        for(IndexType col=0; col<AN; col++) {
            if (AP[col]<=errorIndex && errorIndex<AP[col+1]) {
                return std::make_tuple(AI[errorIndex], col); 
            }
        }
        return std::make_tuple(-1, -1);
    }

    // Error row
    IndexType errorRow() const { return errorIndex; }; 

    // Error row
    IndexType errorRank() const { return errorRank_; };

    // Rebuild it based on the given sparsity map, set to zero, clear error
    bool rebuild(SparsityMap& m, EquationIndex n);

    // Checks if matrix is valid (rebuild completed successfully)
    bool valid() const { return symbolic; };

    // Returns a pointer to element (component), if element is not found returns pointer to bucket
    // Assumes the undelying type is double or std::complex<double> (Complex)
    // This method is used when the type of the matrix is known. 
    double* elementPtr(const MatrixEntryPosition& mep, Component comp=Component::Real) {
        auto [entryOffset, found] = smap->find(mep);
        if (!found) {
            return reinterpret_cast<double*>(&bucket_);
        }
        if constexpr(std::is_same<ValueType, Complex>::value) {
            return (comp==Component::Imaginary) ? 
                reinterpret_cast<double*>(Ax+entryOffset)+1 : 
                reinterpret_cast<double*>(Ax+entryOffset);
        } else {
            return (comp==Component::Imaginary) ? nullptr : (Ax+entryOffset);
        }
    };

    // Returns internal data array or a real matrix
    ValueType* data() { return reinterpret_cast<ValueType*>(Ax); };

    // Return number of nonzeros
    IndexType nnz() const { return AP[AN]; };

    // Set entries to 0, clear error
    void zero(Component what=Component::Real|Component::Imaginary);

    // Factorization
    bool factor();
    bool refactor();
    bool isFactored() const { return numeric; };

    // Reciprocal pivot growth
    bool rgrowth(double& rgrowth);

    // Cheap reciprocal condition number estimation
    bool rcond(double& rcond);

    // Check matrix for inf/nan
    bool isFinite(bool infCheck=false, bool nanCheck=true);

    // Check vector for inf/nan
    bool isFinite(ValueType* vec, bool infCheck=false, bool nanCheck=true);

    // Maximal element in row
    bool rowMaxNorm(double* maxNorm);
    
    // Solve after factorization, result is stored in rhs
    bool solve(ValueType* b);
    
    // Matrix-vector product, result is stored in res
    bool product(ValueType* vec, ValueType* res);
    
    // Residual (Ax-b), stored in res
    bool residual(ValueType* x, ValueType* b, ValueType* res);
    
    // Structural rank
    IndexType structuralRank() const { return common.structural_rank; };

    // Numerical rank
    IndexType numericalRank() const { return common.numerical_rank; };

    // Singular column
    IndexType singularColumn() const { return common.singular_col; };

    // Computes offset of a nonzero element
    std::tuple<IndexType, bool> nonzeroOffset(EquationIndex row, UnknownIndex col);

    // Dump nonzero pattern
    void dumpSparsity(std::ostream& os);

    // Dump nonzero pattern vectors
    void dumpSparsityTables(std::ostream& os);

    // Dump entry values
    void dumpEntries(std::ostream& os);

    // Dump matrix and optional rhs with given column width and precision
    void dump(std::ostream& os, ValueType* rhs=nullptr, int colw=12, int prec=2);
    
    // Dump vector
    void dumpVector(std::ostream& os, ValueType* v, int colw=12, int prec=2);
    
protected:
    Accounting* acct;
    bool isComplex_;
    IndexType AN;
    IndexType* AP;
    IndexType* AI;
    ValueType* Ax;
    Symbolic* symbolic;
    Numeric* numeric;
    Common common;
    SparsityMap* smap;
    
    ValueType bucket_;
    
    // Clear error
    void clearError() { lastError = Error::OK; }; 

    Error lastError;
    IndexType errorIndex;
    IndexType errorRank_;
    bool errorNan;
};

// Default KLU matrix flavor
typedef KluMatrixCore<MatrixEntryIndex, double> KluRealMatrix;
typedef KluMatrixCore<MatrixEntryIndex, Complex> KluComplexMatrix;
typedef MatrixAccess<MatrixEntryIndex> KluMatrixAccess;
}

#endif
