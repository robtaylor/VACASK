#ifndef __KLUMATRIX_DEFINED
#define __KLUMATRIX_DEFINED

#include <suitesparse/klu.h>
#include <unordered_map>
#include <complex>
#include <type_traits>
#include "status.h"
#include "identifier.h"
#include "flags.h"
#include "common.h"


namespace NAMESPACE {

// TODO: make multiple matrices share the same sparsity pattern without allocating a new copy
//       of KLU sparsity pattern

class NameResolver {
public:
    NameResolver() {};

    // Indices in matrix are 0-based, operator() must take this into account
    virtual Id operator()(MatrixEntryIndex u) = 0;
};

class SparsityMap {
public:
    SparsityMap() {};

    SparsityMap           (const SparsityMap&)  = delete;
    SparsityMap           (      SparsityMap&&) = delete;
    SparsityMap& operator=(const SparsityMap&)  = delete;
    SparsityMap& operator=(      SparsityMap&&) = delete;

    // Internal map type
    typedef std::pair<EquationIndex,UnknownIndex> MatrixEntryPosition;
    
    typedef struct MatrixEntryPositionHash {
        auto operator()(const MatrixEntryPosition& p) const -> size_t {
            return std::hash<size_t>{}(
                (static_cast<size_t>(p.first) << (sizeof(size_t)/2)) ^ 
                p.second
            );
        }
    } MatrixEntryPositionHash;

    typedef std::unordered_map<MatrixEntryPosition, MatrixEntryIndex, MatrixEntryPositionHash> Map;

    // Clear
    void clear();
    
    // Size
    MatrixEntryIndex size() const { return smap.size(); };

    // Insert, the returned pointer points to an integer that 
    // will containt the entry index after enumerate() is called
    // bool value indicates if a pre-existing entry was found
    std::tuple<MatrixEntryIndex*, bool> insert(EquationIndex e, UnknownIndex u);

    // Find, bool value indicates if an entry was found
    std::tuple<MatrixEntryIndex, bool> find(EquationIndex e, UnknownIndex u) const;

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


typedef enum Component { RealPart=1, ImagPart=2 } Component;
DEFINE_FLAG_OPERATORS(Component);

template<typename IndexType, typename ValueType> class KluMatrixCore {
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

    ~KluMatrixCore();

    // Clear error
    void clearError() { lastError = Error::OK; }; 

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
    double* elementPtr(EquationIndex e, UnknownIndex u, Component comp=Component::RealPart) {
        auto [entryOffset, found] = smap->find(e, u);
        if (found) {
            if constexpr(std::is_same<ValueType, Complex>::value) {
                return (comp==Component::ImagPart) ? 
                    reinterpret_cast<double*>(Ax+entryOffset)+1 : 
                    reinterpret_cast<double*>(Ax+entryOffset);
            } else {
                return (comp==Component::ImagPart) ? nullptr : (Ax+entryOffset);
            }
        } else {
            return &bucket_;
        }
    };

    // Returns internal data array or a real matrix
    ValueType* data() { return reinterpret_cast<ValueType*>(Ax); };

    // Return number of nonzeros
    IndexType nnz() const { return AP[AN]; };

    // Set entries to 0, clear error
    void zero(Component what=Component::RealPart|Component::ImagPart);

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
    
private:
    bool isComplex_;
    IndexType AN;
    IndexType* AP;
    IndexType* AI;
    ValueType* Ax;
    Symbolic* symbolic;
    Numeric* numeric;
    Common common;
    SparsityMap* smap;
    
    double bucket_;
    
    Error lastError;
    IndexType errorIndex;
    IndexType errorRank_;
    bool errorNan;
};

// Default KLU matrix flavor
typedef KluMatrixCore<MatrixEntryIndex, double> KluRealMatrix;
typedef KluMatrixCore<MatrixEntryIndex, Complex> KluComplexMatrix;

}

#endif
