#include "klubsmatrix.h"
#include "common.h"
#include <iomanip>
#include <algorithm>
#include <type_traits>

namespace NAMESPACE {

template<typename IndexType, typename ValueType> KluBlockSparseMatrixCore<IndexType, ValueType>::KluBlockSparseMatrixCore(bool largeBucket) 
    : denseColumnBegin(nullptr), blockColumnOrigin(nullptr), blockColumnStride(nullptr), blockBucket_(nullptr), largeBucket_(largeBucket) {
}

template<typename IndexType, typename ValueType> KluBlockSparseMatrixCore<IndexType, ValueType>::~KluBlockSparseMatrixCore() {
    if (denseColumnBegin) {
        delete [] denseColumnBegin;
        denseColumnBegin = nullptr;
    }
    if (blockColumnOrigin) {
        delete [] blockColumnOrigin;
        blockColumnOrigin = nullptr;
    }
    if (blockColumnStride) {
        delete [] blockColumnStride;
        blockColumnStride = nullptr;
    }
    if (blockBucket_ && largeBucket_) {
        delete [] blockBucket_;
        blockBucket_ = nullptr;
    }
}

template<typename IndexType, typename ValueType> 
double* KluBlockSparseMatrixCore<IndexType, ValueType>::valueArray() {
    if constexpr(std::is_same<ValueType, Complex>::value) {
        return nullptr;
    } else {
        return KluMatrixCore<IndexType, ValueType>::data();
    }
} 

template<typename IndexType, typename ValueType> 
Complex* KluBlockSparseMatrixCore<IndexType, ValueType>::cxValueArray() {
    if constexpr(std::is_same<ValueType, Complex>::value) {
        return KluMatrixCore<IndexType, ValueType>::data();
    } else {
        return nullptr;
    }
} 

template<typename IndexType, typename ValueType> 
std::tuple<IndexType, bool> KluBlockSparseMatrixCore<IndexType, ValueType>::valueIndex(
    const MatrixEntryPosition& mep, const std::optional<MatrixEntryPosition>& blockMep
) const {
    return elementIndex(mep, blockMep);
}

template<typename IndexType, typename ValueType> 
double* KluBlockSparseMatrixCore<IndexType, ValueType>::valuePtr(
    const MatrixEntryPosition& mep, Component comp, const std::optional<MatrixEntryPosition>& blockMep
) {
    return elementPtr(mep, comp, blockMep);
}

template<typename IndexType, typename ValueType> 
Complex* KluBlockSparseMatrixCore<IndexType, ValueType>::cxValuePtr(
    const MatrixEntryPosition& mep, const std::optional<MatrixEntryPosition>& blockMep
) {
    if constexpr(std::is_same<ValueType, Complex>::value) {
        auto [nzPosition, found] = elementIndex(mep, blockMep);
        if (found) {
            return Ax+nzPosition;
        } else {
            return blockBucket_;
        }
    } else {
        return nullptr;
    }
}

template<typename IndexType, typename ValueType> 
bool KluBlockSparseMatrixCore<IndexType, ValueType>::rebuild(SparsityMap& m, EquationIndex n, EquationIndex nbRow, UnknownIndex nbCol) {
    KluMatrixCore<IndexType, ValueType>::clearError();
    
    this->~KluBlockSparseMatrixCore();
    this->KluMatrixCore<IndexType, ValueType>::~KluMatrixCore();

    n_ = n;
    nbRow_ = nbRow;
    nbCol_ = nbCol;

    smap = &m;
    
    // Set number of columns of elements
    AN = n_*nbCol_;

    // Number of nonzeros (for now, assume all dense blocks are fully dense)
    auto nnz_ = m.size()*nbRow_*nbCol_;

    // Allocate arrays
    AP = new IndexType[AN+1];
    AI = new IndexType[nnz_];
    denseColumnBegin = new IndexType[n+1];
    blockColumnOrigin = new IndexType[n];
    blockColumnStride = new IndexType[n];
    if (!AP || !AI || !denseColumnBegin || !blockColumnOrigin || !blockColumnStride) {
        this->~KluBlockSparseMatrixCore();
        this->KluMatrixCore<IndexType, ValueType>::~KluMatrixCore();
        lastError = Error:: Memory;
        return false;
    }

    // Element column index
    decltype(nnz_) atCol = 0;

    // Nonzero element index
    decltype(nnz_) atNz = 0;
    
    // Collect indices of the beginnings of columns in the sparsity map positions vector. 
    // Entries are already sorted by column first, then row
    // so they are in the same order as the dense blocks appear in column-major ordering. 
    denseColumnBegin[0] = 0;
    atCol = 0;
    auto& positions = m.positions();
    for(size_t posNdx=0; posNdx<positions.size(); posNdx++) {
        auto [row, col] = positions[posNdx];
        // Make column index 0-based
        col--;
        // Reached next column
        if (col!=atCol) {
            // Make sure empty columns of blocks are also handled
            for(; atCol<col;) {
                atCol++;
                denseColumnBegin[atCol] = posNdx;
            }
        }
    }
    // Trailing empty columns
    // n-1<atCol is an internal error!
    for(; atCol<n-1;) {
        atCol++;
        denseColumnBegin[atCol] = positions.size();
    }
    // Number of dense blocks
    denseColumnBegin[n] = positions.size();

    // Index of a column of elements
    atCol = 0;
    // Index of nonzero element
    atNz = 0;
    // Number of dense blocks
    auto blockCount = denseColumnBegin[n];
    // Iterate through columns of dense blocks
    // This also iterates throuh columns with no dense blocks
    for(decltype(n) blockColNdx=0; blockColNdx<n; blockColNdx++) {
        // Beginning and end of a column of blocks
        auto colBeginNdx = denseColumnBegin[blockColNdx];
        auto colEndNdx = denseColumnBegin[blockColNdx+1];

        // How many dense blocks do we have in this column
        auto blocksInColumn = colEndNdx-colBeginNdx;

        // Add column origin and stride
        blockColumnOrigin[blockColNdx] = atNz;
        blockColumnStride[blockColNdx] = blocksInColumn*nbRow_;

        // No blocks in this block column, nothing to do
        if (blocksInColumn==0) {
            continue;
        }

        // Iterate through subcolumns of each dense block
        // This loop runs even if there are no dense blocks in this column of dense blocks
        // Therefore AP is filled with indices correctly even in this case
        for(decltype(nbCol_) subColNdx=0; subColNdx<nbCol_; subColNdx++) {
            // Add index of first nonzero element in column
            AP[atCol] = atNz;

            // For each dense block in column of dense blocks
            for(auto blkPos = colBeginNdx; blkPos < colEndNdx ; blkPos++) {
                // Get block row and column index
                // We do not need blkCol - it is useful for debugging
                auto [blkRow, blkCol] = positions[blkPos];
                // These indices are 1-based, make them 0-based
                blkRow--;
                blkCol--;

                // For each row in dense block
                IndexType rowIndex = blkRow * nbRow_;
                for(decltype(nbRow_) subRowNdx=0; subRowNdx<nbRow_; subRowNdx++) {
                    // Write row index
                    AI[atNz] = rowIndex;

                    // Advance row index
                    rowIndex++;

                    // Advance nonzero element index
                    atNz++;
                }
            }

            // Advance column index
            atCol++;
        }
    }

    // Add final AP entry
    AP[atCol] = atNz;
    
    // Allocate array for nozero element values
    Ax = new ValueType[nnz_];
    if (largeBucket_) {
        blockBucket_ = new ValueType[nbRow_*nbCol_];
    } else {
        blockBucket_ = &bucket_;
    }
    if (!Ax || (largeBucket_ && !blockBucket_)) {
        this->~KluBlockSparseMatrixCore();
        this->KluMatrixCore<IndexType, ValueType>::~KluMatrixCore();
        lastError = Error::Memory;
        return false;
    }

    // Zero array
    KluMatrixCore<IndexType, ValueType>::zero();
    
    // Set up KLU structures
    int st;
    if constexpr(std::is_same<int32_t, IndexType>::value) {
        st = klu_defaults(&common);
    } else {
        st = klu_l_defaults(&common);
    }
    if (!st) {
        lastError = Error::Defaults;
        return false;
    }

    if constexpr(std::is_same<int32_t, IndexType>::value) {
        symbolic = klu_analyze(AN, AP, AI, &common);
    } else {
        symbolic = klu_l_analyze(AN, AP, AI, &common);
    }
    if (!symbolic) {
        lastError = Error::Analysis;
        return false;
    }
    
    return true;
}

template<typename IndexType, typename ValueType> 
void KluBlockSparseMatrixCore<IndexType, ValueType>::dumpBlockSparsity(std::ostream& os) {
   for(IndexType row=0; row<n_; row++) {
        for(IndexType col=0; col<n_; col++) {
            auto [_, found] = block(MatrixEntryPosition(row+1, col+1));
            if (found) {
                std::cout << "x";
            } else {
                std::cout << ".";
            }
        }
        std::cout << "\n";
    } 
}


// Instantiate template class for int32 and int64 indices, double and Complex values
template class KluBlockSparseMatrixCore<int32_t, double>;
template class KluBlockSparseMatrixCore<int32_t, Complex>;
template class KluBlockSparseMatrixCore<int64_t, double>;
template class KluBlockSparseMatrixCore<int64_t, Complex>;

}
