#ifndef __KLUBSMATRIX_DEFINED
#define __KLUBSMATRIX_DEFINED

#include <unordered_map>
#include <complex>
#include <type_traits>
#include <optional>
#include "status.h"
#include "densematrix.h"
#include "klumatrix.h"
#include "identifier.h"
#include "flags.h"
#include "hash.h"
#include "acct.h"
#include "common.h"


namespace NAMESPACE {

// Blocks of size nb x nb
// mep          .. coordinates of the block
// blockMep     .. element coordinates within the block
// block origin .. element with blockMep=(1,1)
template<typename IndexType, typename ValueType> 
class KluBlockSparseMatrixCore : public KluMatrixCore<IndexType, ValueType> {
public: 
    KluBlockSparseMatrixCore();
    
    KluBlockSparseMatrixCore           (const KluBlockSparseMatrixCore&)  = delete;
    KluBlockSparseMatrixCore           (      KluBlockSparseMatrixCore&&) = delete;
    KluBlockSparseMatrixCore& operator=(const KluBlockSparseMatrixCore&)  = delete;
    KluBlockSparseMatrixCore& operator=(      KluBlockSparseMatrixCore&&) = delete;

    virtual ~KluBlockSparseMatrixCore();

    // Matrix binding interface
    // If blockMep is not given the block origin is returned
    virtual double* valueArray();
    virtual Complex* cxValueArray();
    virtual std::tuple<IndexType, bool> valueIndex(const MatrixEntryPosition& mep, const std::optional<MatrixEntryPosition>& blockMep=std::nullopt) const;
    virtual double* valuePtr(const MatrixEntryPosition& mep, Component comp=Component::Real, const std::optional<MatrixEntryPosition>& blockMep=std::nullopt) ;
    virtual Complex* cxValuePtr(const MatrixEntryPosition& mep, const std::optional<MatrixEntryPosition>& blockMep=std::nullopt);

    // BlockSparseMatrixCore specific interface
    // Returns a dense matrix view of a block. 
    // Storage is column major due to KLU. 
    // A block column occupies a consecutive block of memory. 
    // Consecutive columns of the same bloc do not generally occupy a continuous block of memory. 
    // Column stride depends on the number of dense blocks in a column of dense blocks. 
    //   column stride = number of dense blocks in the column x nb
    // Returns DenseMatrixView of block, found flag
    // If the block is not found the dense matrix view has row and column stride 0. 
    // All elements refer to the bucket. 
    std::tuple<DenseMatrixView<ValueType>, bool> block(const MatrixEntryPosition& mep) {
        auto [nzPosition, found] = elementIndex(mep);
        if (!found) {
            return std::make_tuple(DenseMatrixView<ValueType>(&bucket_, nb_, nb_, 0, 0), false);
        }
        // KLU organizes elements in column major order
        // row stride is 1, column stride depends on the column of dense blocks
        auto [row, col] = mep;
        return std::make_tuple(
            DenseMatrixView<ValueType>(Ax+nzPosition, nb_, nb_, 1, blockColumnStride[col]), 
            true
        );
    };

    // Rebuild it based on the given sparsity map of dense blocks, 
    // dense block size nb x nb, set elements to zero, clear error
    bool rebuild(SparsityMap& m, EquationIndex n, EquationIndex nb);

    // Returns the linear nonzero element index coresponding to dense block
    // at block position mep (1-based), block element position blockMep (1-based). 
    // If blockMep is not given assumes (1, 1), i.e. block origin. 
    // Returns index, found. found=true if element exists. 
    std::tuple<IndexType, bool> elementIndex(const MatrixEntryPosition& mep, const std::optional<MatrixEntryPosition>& blockMep=std::nullopt) const {
        auto [entryOffset, found] = smap->find(mep);
        if (!found) {
            return std::make_tuple(0, false);
        }
        
        // Get 0-based block position
        auto [row, col] = mep;
        row--;
        col--;
        // Get index of the first dense block in column
        auto firstBlockInColumn = denseColumnBegin[col];
        // Which block in column is this
        auto blockInColumn = entryOffset - firstBlockInColumn;
        // Compute element index of block origin
        auto nzPosition = blockColumnOrigin[col] + nb_*blockInColumn;

        // Do we have a blockMep
        if (blockMep.has_value()) {
            // Get 0-based dense block element position
            auto [brow, bcol] = blockMep.value();
            brow--;
            bcol--;
            // Compute element position
            nzPosition += bcol * blockColumnStride[col] + brow;
        }
        return std::make_tuple(nzPosition, true); 
    };

    // Returns a pointer to element (component), if element is not found returns pointer to bucket. 
    // Assumes the undelying type is double or std::complex<double> (Complex)
    // This method is used when the type of the matrix is known. 
    // If blockMap is not given returns the element at the origin of a dense block. 
    double* elementPtr(const MatrixEntryPosition& mep, Component comp=Component::Real, const std::optional<MatrixEntryPosition>& blockMep=std::nullopt) {
        auto [nzPosition, found] = elementIndex(mep, blockMep);
        if (!found) {
            return reinterpret_cast<double*>(&bucket_);
        }
        // Return pointer
        if constexpr(std::is_same<ValueType, Complex>::value) {
            return (comp==Component::Imaginary) ? 
                reinterpret_cast<double*>(Ax+nzPosition)+1 : 
                reinterpret_cast<double*>(Ax+nzPosition);
        } else {
            return (comp==Component::Imaginary) ? nullptr : (Ax+nzPosition);
        }
    };

    IndexType blocksInColumn() const { return n_; };
    IndexType blockRows() const { return nb_; };
    
    void dumpBlockSparsity(std::ostream& os);

protected:
    using Error = KluMatrixCore<IndexType, ValueType>::Error;
    using KluMatrixCore<IndexType, ValueType>::smap;
    using KluMatrixCore<IndexType, ValueType>::AN;
    using KluMatrixCore<IndexType, ValueType>::AP;
    using KluMatrixCore<IndexType, ValueType>::AI;
    using KluMatrixCore<IndexType, ValueType>::Ax;
    using KluMatrixCore<IndexType, ValueType>::bucket_;
    using KluMatrixCore<IndexType, ValueType>::lastError;
    using KluMatrixCore<IndexType, ValueType>::common;
    using KluMatrixCore<IndexType, ValueType>::symbolic;

    // Number of blocks in row/column
    IndexType n_;

    // Number of rows/columns in a dense block
    IndexType nb_;

    // Origin of dense block column (origin of topmost dense block).
    // This is the linear index of the nonzero element at the topmost block's origin. 
    // Array has n entries. 
    IndexType* blockColumnOrigin;
    
    // Number of nnz elements to skip to reach the element 
    // in the next column of the same row of a dense block. 
    // Depends on the number of dense blocks in a column. 
    //   number of dense blocks in the column x nb
    // Array has n elements. 
    IndexType* blockColumnStride;

    // Linear index of the first dense block in a column of dense blocks
    // when the dense blocks are ordered in a column major ordering. 
    // has n+1 elements where the n+1-th element is the number of dense blocks. 
    IndexType* denseColumnBegin;
};

// Default KLU matrix flavor
typedef KluBlockSparseMatrixCore<MatrixEntryIndex, double> KluBlockSparseRealMatrix;
typedef KluBlockSparseMatrixCore<MatrixEntryIndex, Complex> KluBlockSparseComplexMatrix;

}

#endif
