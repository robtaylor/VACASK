#include "klumatrix.h"
#include "common.h"
#include <iomanip>
#include <algorithm>
#include <type_traits>

namespace NAMESPACE {

void SparsityMap::clear() {
    smap.clear();
    ordering.clear();
}

/*
std::tuple<MatrixEntryIndex*, bool> SparsityMap::insert(EquationIndex e, UnknownIndex u) {
    auto [it, inserted] = smap.insert({std::make_pair(e, u), 0});
    return std::make_tuple(&(it->second), true);
}

std::tuple<MatrixEntryIndex, bool> SparsityMap::find(EquationIndex e, UnknownIndex u) const {
    auto it = smap.find(std::make_pair(e, u));
    if (it==smap.end()) {
        return std::make_tuple(0, false);
    }
    return std::make_tuple(it->second, true);
}
*/

// We could make ordering more efficient - now we must order up to n^2 elements 
// which is on average n^2 log(n^2) operations with std::sort(). 
// If we collect them by columns and then sort each column separately, 
// the average complexity with std::sort() is n n log(n) operations which is half 
// of the previous. 
void SparsityMap::enumerate() {
    // Prepare a vector of map keys
    ordering.clear();
    for(auto it=smap.begin(); it!=smap.end(); ++it) {
        ordering.push_back(it->first);
    }
    
    // Order them
    struct {
        bool operator()(const MatrixEntryPosition& lhs, const MatrixEntryPosition& rhs) const {
            // Compare first by column (unknown), then by row (equation)
            return (lhs.second < rhs.second) || ((lhs.second == rhs.second) && (lhs.first < rhs.first));
        }
    } comparison;
    
    std::sort(ordering.begin(), ordering.end(), comparison);

    // Traverse keys, enumerate entries
    MatrixEntryIndex num = 0;
    for(auto it=ordering.begin(); it!=ordering.end(); ++it) {
        auto e = it->first;
        auto u = it->second;
        smap[*it] = num;
        num++;
    }
}

void SparsityMap::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    for(auto it=ordering.cbegin(); it!=ordering.end(); ++it) {
        auto [offs, found] = find(it->first, it->second);
        os << pfx << "(" << it->first << ", " << it->second << ") : ";
        if (found) {
            os << offs;
        } else {
            os << "?";
        }
        os << "\n";
    }
}


template<typename IndexType, typename ValueType> KluMatrixCore<IndexType, ValueType>::KluMatrixCore() 
    : AN(0), AP(nullptr), AI(nullptr), numeric(nullptr), symbolic(nullptr),  
      Ax(nullptr), smap(nullptr), lastError(Error::OK), acct(nullptr) {
    // Sanity check: IndexType can only be int32_t or int64_t
    static_assert(
        std::is_same<IndexType, int>::value || std::is_same<IndexType, int64_t>::value, 
        "Klu matrix index type is neither int32_t nor int64_t."
    );
    // Sanity check: ValueType can only be double or std::complex<double>
    static_assert(
        std::is_same<ValueType, double>::value || std::is_same<ValueType, std::complex<double>>::value, 
        "Klu matrix value type is neither double nor std::complex<double>."
    );
}

template<typename IndexType, typename ValueType> KluMatrixCore<IndexType, ValueType>::~KluMatrixCore() {
    if (numeric) {
        if constexpr(std::is_same<ValueType, Complex>::value) {
            if constexpr(std::is_same<int32_t,IndexType>::value) {
                klu_z_free_numeric(&numeric, &common);
            } else {
                klu_zl_free_numeric(&numeric, &common);
            }
        } else {
            if constexpr(std::is_same<int32_t,IndexType>::value) {
                klu_free_numeric(&numeric, &common);
            } else {
                klu_l_free_numeric(&numeric, &common);
            }
        }
        numeric = nullptr;
    }
    if (symbolic) {
        if constexpr(std::is_same<int32_t,IndexType>::value) {
            klu_free_symbolic(&symbolic, &common);
        } else {
            klu_l_free_symbolic(&symbolic, &common);
        }
        symbolic = nullptr;
    }
    if (AI) {
        delete [] AI;
        AI = nullptr;
    }
    if (AP) {
        delete [] AP;
        AP = nullptr;
    }
    if (Ax) {
        delete [] Ax;
        Ax = nullptr;
    }
}

template<typename IndexType, typename ValueType> double* KluMatrixCore<IndexType, ValueType>::valueArray() {
    if constexpr(std::is_same<ValueType, Complex>::value) {
        return nullptr;
    } else {
        return data();
    }
} 

template<typename IndexType, typename ValueType> 
Complex* KluMatrixCore<IndexType, ValueType>::cxValueArray() {
    if constexpr(std::is_same<ValueType, Complex>::value) {
        return data();
    } else {
        return nullptr;
    }
} 

template<typename IndexType, typename ValueType> 
std::tuple<IndexType, bool> KluMatrixCore<IndexType, ValueType>::valueIndex(IndexType row, IndexType col) {
    return smap->find(row, col);
}

template<typename IndexType, typename ValueType> 
double* KluMatrixCore<IndexType, ValueType>::valuePtr(IndexType row, IndexType col, Component comp) {
    return elementPtr(row, col, comp);
}

template<typename IndexType, typename ValueType> 
Complex* KluMatrixCore<IndexType, ValueType>::cxValuePtr(IndexType row, IndexType col) {
    if constexpr(std::is_same<ValueType, Complex>::value) {
        auto [entryOffset, found] = smap->find(row, col);
        if (found) {
            return Ax+entryOffset;
        } else {
            return &bucket_;
        }
    } else {
        return nullptr;
    }
}

template<typename IndexType, typename ValueType> bool KluMatrixCore<IndexType, ValueType>::rebuild(SparsityMap& m, EquationIndex n) {
    clearError();
    
    this->~KluMatrixCore();

    smap = &m;
    
    AN = n;
    auto nnz_ = m.size();
    AP = new IndexType[n+1];
    AI = new IndexType[nnz_];
    if (!AP || !AI) {
        this->~KluMatrixCore();
        lastError = Error::Memory;
        return false;
    }

    decltype(nnz_) atCol = 0;
    decltype(nnz_) atNz = 0;
    
    // Go through allocated Jacobian entries
    // These entries are already sorted by column first, then row
    // so they are in the same order as they appear in jacI
    // Note that Jacobian indices are zero-based, so we subtract 1 from row and column index
    bool start = true;
    for(auto it=m.positions().begin(); it!=m.positions().end(); ++it) {
        auto row = it->first;
        auto col = it->second;

        // Skip entries that have zero index (they correspond to ground)
        if (!row || !col) {
            continue;
        }

        // KLU sparsity pattern indices are 0-based
        row -= 1;
        col -= 1;

        // First entry, or new column
        if (start || atCol!=col) {
            // Starting new / first column
            AP[col] = atNz;
            atCol = col;
            start = false;
        }
        AI[atNz] = row;
        atNz++;
    }
    AP[atCol+1] = nnz_;
    
    // New does not call default constructor for builtin types
    // i.e. doubles are not initialized to 0. 
    // If, however we do
    //   new double[...]() 
    // then initialization takes place. 
    Ax = new ValueType[nnz_];
    
    if (!Ax) {
        this->~KluMatrixCore();
        lastError = Error::Memory;
        return false;
    }
    zero();
    
    // Maybe move this to clear
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

template<typename IndexType, typename ValueType> void KluMatrixCore<IndexType, ValueType>::zero(Component what) {
    auto nnz_ = AP[AN];
    if constexpr(std::is_same<ValueType, Complex>::value) {
        if (what==(Component::RealPart|Component::ImagPart)) {
            for(IndexType i=0; i<nnz_; i++) {
                Ax[i] = 0.0;
            }
        } else if (what==Component::RealPart) {
            for(IndexType i=0; i<nnz_; i++) {
                Ax[i].real(0.0);
            }
        } else if (what==Component::ImagPart) {
            for(IndexType i=0; i<nnz_; i++) {
                Ax[i].imag(0.0);
            }
        }
    } else {
        if ((what&Component::RealPart)==Component::RealPart) {
            for(IndexType i=0; i<nnz_; i++) {
                Ax[i] = 0.0;
            }
        }
    }
    // Clear error
    clearError();
}

template<typename IndexType, typename ValueType> bool KluMatrixCore<IndexType, ValueType>::factor() {
    auto t0 = Accounting::wclk();
    if (acct) {
        if constexpr(std::is_same<ValueType, Complex>::value) {
            acct->acctNew.cxfactor++;
        } else {
            acct->acctNew.factor++;
        }
    }

    clearError();

    if (numeric) {
        if constexpr(std::is_same<int32_t, IndexType>::value) {
            klu_free_numeric(&numeric, &common);
        } else {
            klu_l_free_numeric(&numeric, &common);
        }
    }
    if constexpr(std::is_same<ValueType, Complex>::value) {
        if constexpr(std::is_same<int32_t, IndexType>::value) {
            numeric = klu_z_factor(AP, AI, reinterpret_cast<double*>(Ax), symbolic, &common);
        } else {
            numeric = klu_zl_factor(AP, AI, reinterpret_cast<double*>(Ax), symbolic, &common);
        }
    } else {
        if constexpr(std::is_same<int32_t, IndexType>::value) {
            numeric = klu_factor(AP, AI, Ax, symbolic, &common);
        } else {
            numeric = klu_l_factor(AP, AI, Ax, symbolic, &common);
        }
    }
    bool isSingular = common.status==KLU_SINGULAR;
    auto nr = numericalRank();
    if (acct) {
        if constexpr(std::is_same<ValueType, Complex>::value) {
            acct->acctNew.tcxfactor += Accounting::wclkDelta(t0);
        } else {
            acct->acctNew.tfactor += Accounting::wclkDelta(t0);
        }
    }
    // Check status and numerical rank if it was computed
    if (!numeric || isSingular || (nr>=0 && nr!=AN)) {
        lastError = Error::Factorization;
        errorIndex = singularColumn();
        errorRank_ = numericalRank();
        return false;
    }
    return true;
}

template<typename IndexType, typename ValueType> bool KluMatrixCore<IndexType, ValueType>::refactor() {
    auto t0 = Accounting::wclk();
    if (acct) {
        if constexpr(std::is_same<ValueType, Complex>::value) {
            acct->acctNew.cxrefactor++;
        } else {
            acct->acctNew.refactor++;
        }
    }

    clearError();

    if (!numeric) {
        return factor();
    }
    int st;
    if constexpr(std::is_same<ValueType, Complex>::value) {
        if constexpr(std::is_same<int32_t, IndexType>::value) {
            st = klu_z_refactor(AP, AI, reinterpret_cast<double*>(Ax), symbolic, numeric, &common);
        } else {
            st = klu_zl_refactor(AP, AI, reinterpret_cast<double*>(Ax), symbolic, numeric, &common);
        }
    } else {
        if constexpr(std::is_same<int32_t, IndexType>::value) {
            st = klu_refactor(AP, AI, Ax, symbolic, numeric, &common);
        } else {
            st = klu_l_refactor(AP, AI, Ax, symbolic, numeric, &common);
        }
    }
    bool isSingular = common.status==KLU_SINGULAR;
    auto nr = numericalRank();
    if (acct) {
        if constexpr(std::is_same<ValueType, Complex>::value) {
            acct->acctNew.tcxrefactor += Accounting::wclkDelta(t0);
        } else {
            acct->acctNew.trefactor += Accounting::wclkDelta(t0);
        }
    }
    // Check status and numerical rank if it was computed
    if (!st || isSingular || (nr>=0 && nr!=AN)) {
        lastError = Error::Refactorization;
        errorRank_ = numericalRank();
        return false;
    }
    return true;
}

template<typename IndexType, typename ValueType> bool KluMatrixCore<IndexType, ValueType>::rgrowth(double& rgrowth) {
    clearError();

    int st;
    if constexpr(std::is_same<ValueType, Complex>::value) {
        if constexpr(std::is_same<int32_t, IndexType>::value) {
            st = klu_z_rgrowth(AP, AI, reinterpret_cast<double*>(Ax), symbolic, numeric, &common);
        } else {
            st = klu_zl_rgrowth(AP, AI, reinterpret_cast<double*>(Ax), symbolic, numeric, &common);
        }
    } else {
        if constexpr(std::is_same<int32_t, IndexType>::value) {
            st = klu_rgrowth(AP, AI, Ax, symbolic, numeric, &common);
        } else {
            st = klu_l_rgrowth(AP, AI, Ax, symbolic, numeric, &common);
        }
    }
    if (!st) {
        lastError = Error::ReciprocalPivotGrowth;
        return false;
    }
    rgrowth = common.rgrowth;
    return true;
}

template<typename IndexType, typename ValueType> bool KluMatrixCore<IndexType, ValueType>::rcond(double& rcond) {
    clearError();

    int st;
    if constexpr(std::is_same<ValueType, Complex>::value) {
        if constexpr(std::is_same<int32_t, IndexType>::value) {
            st = klu_z_rcond(symbolic, numeric, &common);
        } else {
            st = klu_zl_rcond(symbolic, numeric, &common);
        }
    } else {
        if constexpr(std::is_same<int32_t, IndexType>::value) {
            st = klu_rcond(symbolic, numeric, &common);
        } else {
            st = klu_l_rcond(symbolic, numeric, &common);
        }
    }
    if (!st) {
        lastError = Error::ReciprocalCondEstimate;
        return false;
    }
    rcond = common.rcond;
    return true;
}

template<typename IndexType, typename ValueType> bool KluMatrixCore<IndexType, ValueType>::isFinite(bool infCheck, bool nanCheck) {
    clearError();

    if (!infCheck && !nanCheck) {
        return true;
    }
    // Check matrix
    bool gotInf = false;
    bool gotNan = false;
    auto nnz = AP[AN];
    IndexType i;
    for(i=0; i<nnz; i++) {
        if constexpr(std::is_same<ValueType, Complex>::value) {
            if (nanCheck) {
                gotNan = gotNan || (std::isnan(Ax[i].real()) || std::isnan(Ax[i].imag()));
            } 
            if (infCheck) {
                gotInf = gotInf || (std::isinf(Ax[i].real()) || std::isinf(Ax[i].imag()));
            }
        } else {
            if (nanCheck) {
                gotNan = gotNan || std::isnan(Ax[i]);
            }
            if (infCheck) {
                gotInf = gotInf || std::isinf(Ax[i]);
            }
        }
        if (gotInf || gotNan) {
            break;
        }
    }
           
    if (gotInf || gotNan) {
        lastError = Error::MatrixInfNan;
        errorIndex = i;
        errorNan = gotNan;
        return false;
    }
    return true;
}

template<typename IndexType, typename ValueType> bool KluMatrixCore<IndexType, ValueType>::isFinite(ValueType* vec, bool infCheck, bool nanCheck) {
    clearError();

    if (!infCheck && !nanCheck) {
        return true;
    }
    for(IndexType i=0; i<AN; i++) {
        bool gotInf = false;
        bool gotNan = false;
        if constexpr(std::is_same<ValueType, Complex>::value) {
            if (nanCheck && (std::isnan(vec[i].real()) || std::isnan(vec[i].imag()))) {
                // NaN found
                gotNan = true;
            } else if (infCheck && (std::isinf(vec[i].real()) || std::isinf(vec[i].imag()))) {
                // Inf found
                gotInf = true;
            }
        } else {
            if (nanCheck && std::isnan(vec[i])) {
                // NaN found
                gotNan = true;
            } else if (infCheck && std::isinf(vec[i])) {
                // Inf found
                gotInf = true;
            }
        }
        if (gotInf || gotNan) {
            lastError = Error::VectorInfNan;
            errorIndex = i;
            errorNan = gotNan;
            return false;
        }
    }  
    return true;
}

template<typename IndexType, typename ValueType> bool KluMatrixCore<IndexType, ValueType>::rowMaxNorm(double* maxNorm) {
    // Zero out result
    for(IndexType i=0; i<AN; i++) {
        maxNorm[i] = 0.0;
    }

    // Go through entries
    IndexType col1, col2;
    auto nnz = AP[AN];
    for(IndexType i=0; i<nnz; i++) {
        auto row = AI[i];
        double nrm;
        if constexpr(std::is_same<double, ValueType>::value) {
            nrm = std::abs(Ax[i]);
        } else {
            nrm = std::sqrt(Ax[i].real()*Ax[i].real() + Ax[i].imag()*Ax[i].imag());
        }
        if (nrm>maxNorm[row]) {
            maxNorm[row] = nrm;
        }
    }

    return true;
}

template<typename IndexType, typename ValueType> bool KluMatrixCore<IndexType, ValueType>::solve(ValueType* b) {
    auto t0 = Accounting::wclk();
    if (acct) {
        if constexpr(std::is_same<ValueType, Complex>::value) {
            acct->acctNew.cxsolve++;
        } else {
            acct->acctNew.solve++;
        }
    }

    clearError();
    
    int st;
    if constexpr(std::is_same<ValueType, Complex>::value) {
        if constexpr(std::is_same<int32_t, IndexType>::value) {
            st = klu_z_solve(symbolic, numeric, AN, 1, reinterpret_cast<double*>(b), &common);
        } else {
            st = klu_zl_solve(symbolic, numeric, AN, 1, reinterpret_cast<double*>(b), &common);
        }
    } else {
        if constexpr(std::is_same<int32_t, IndexType>::value) {
            st = klu_solve(symbolic, numeric, AN, 1, b, &common);
        } else {
            st = klu_l_solve(symbolic, numeric, AN, 1, b, &common);
        }
    }

    if (acct) {
        if constexpr(std::is_same<ValueType, Complex>::value) {
            acct->acctNew.tcxsolve += Accounting::wclkDelta(t0);
        } else {
            acct->acctNew.tsolve += Accounting::wclkDelta(t0);
        }
    }
    
    if (!st) {
        lastError = Error::Solve;
        return false;
    }
    return true;
}

// Both vectors must be distinct
template<typename IndexType, typename ValueType> bool KluMatrixCore<IndexType, ValueType>::product(ValueType* vec, ValueType* res) {
    // Zero out result
    for(IndexType i=0; i<AN; i++) {
        res[i] = 0.0;
    }

    // Go through entries
    IndexType col1, col2;
    for(IndexType col=0; col<AN; col++) {
        col1 = AP[col];
        col2 = AP[col+1];
        for(IndexType i=col1; i<col2; i++) {
            auto row = AI[i];
            res[row] += Ax[i]*vec[col];
        }
    }

    return true;
}

// All 3 vectors must be distinct
template<typename IndexType, typename ValueType> bool KluMatrixCore<IndexType, ValueType>::residual(ValueType* x, ValueType* b, ValueType* res) {
    product(x, res);
    for(IndexType i=0; i<AN; i++) {
        res[i] -= b[i];
    }
    return true;
}

template<typename IndexType, typename ValueType> std::tuple<IndexType, bool> KluMatrixCore<IndexType, ValueType>::nonzeroOffset(EquationIndex row, UnknownIndex col) {
    IndexType i1 = AP[col];
    IndexType i2 = AP[col+1];

    // Last entry for this column is at i2-1
    i2 = i2-1;
    // i1>i2 ... column is empty
    if (i1>i2) {
        return std::make_tuple(0, false);
    }
    // Check endpoints
    if (AI[i1]==row) {
        return std::make_tuple(i1, true);
    }
    if (AI[i2]==row) {
        return std::make_tuple(i2, true);
    }

    // Bisect for row between i1 and i2
    // At the beginning of the loop body both endpoints i1 and i2 are already checked
    // Therefore if i2-i1==1 we are done.
    while (i2-i1>1) {
        IndexType ic = (i1+i2)/2;
        if (AI[ic]==row) {
            return std::make_tuple(ic, true);
        } else if (AI[ic]>row) {
            i2 = ic;
        } else {
            i1 = ic;
        }
    }
    
    return std::make_tuple(0, false);
}

template<typename IndexType, typename ValueType> void KluMatrixCore<IndexType, ValueType>::dumpSparsity(std::ostream& os) {
   for(IndexType row=0; row<AN; row++) {
        if (row>0) {
            std::cout << std::endl;
        }
        for(IndexType col=0; col<AN; col++) {
            auto [offs, found] = nonzeroOffset(row, col);
            if (found) {
                std::cout << "x";
            } else {
                std::cout << ".";
            }
        }
    } 
}

template<typename IndexType, typename ValueType> void KluMatrixCore<IndexType, ValueType>::dumpSparsityTables(std::ostream& os) {
    os << "Ap: ";
    for(IndexType i=0; i<=AN; i++) {
        os << AP[i] << " ";
    }
    os << std::endl;
    os << "Ai: ";
    for(IndexType i=0; i<AP[AN]; i++) {
        os << AI[i] << " ";
    }
}

template<typename IndexType, typename ValueType> void KluMatrixCore<IndexType, ValueType>::dumpEntries(std::ostream& os) {
    os << "Ax: ";
    for(IndexType i=0; i<AP[AN]; i++) {
        os << Ax[i] << " ";
    }
}

template<typename IndexType, typename ValueType> void KluMatrixCore<IndexType, ValueType>::dump(std::ostream& os, ValueType* rhs, int colw, int prec) {
    std::ios oldState(nullptr);
    oldState.copyfmt(std::cout);
    os << std::scientific << std::setprecision(prec);

    os << "     ";
    for(IndexType col=0; col<AN; col++) {
        if constexpr(std::is_same<ValueType, Complex>::value) {
            os << std::setw(colw*2+1) << col;
        } else {
            os << std::setw(colw) << col;
        }
    }
    if (rhs) {
        if constexpr(std::is_same<ValueType, Complex>::value) {
            os << std::setw(colw*2+1) << "rhs";
        } else {
            os << std::setw(colw) << "rhs";
        }
    }
    std::cout << std::endl;
    for(IndexType row=0; row<AN; row++) {
        if (row>0) {
            os << std::endl;
        }
        os << std::setw(4) << row << ":";
        for(IndexType col=0; col<AN; col++) {
            double* ptr;
            auto [offs, found] = nonzeroOffset(row, col);
            if (found) {
                if constexpr(std::is_same<ValueType, Complex>::value) {
                    os << std::setw(colw) << (Ax+offs)->real();
                    if ((Ax+offs)->imag()>=0) {
                        os << "+" << std::setw(colw-1) << (Ax+offs)->imag();
                    } else {
                        os << std::setw(colw) << (Ax+offs)->imag();
                    }
                    os << "i";
                } else {
                    os << std::setw(colw) << *(Ax+offs); 
                }
            } else {
                if constexpr(std::is_same<ValueType, Complex>::value) {
                    os << std::setw(colw*2+1) << "0+0i";
                } else {
                    os << std::setw(colw) << 0;
                }
            }
        }
        if (rhs) {
            if constexpr(std::is_same<ValueType, Complex>::value) {
                os << std::setw(colw) << rhs[row].real();
                if (rhs[row].imag()>=0) {
                    os << "+" << std::setw(colw-1) << rhs[row].imag();
                } else {
                    os << std::setw(colw) << rhs[row].imag();
                }
                os << "i";
            } else {
                os << std::setw(colw) << rhs[row];
            }
        }
    }

    std::cout.copyfmt(oldState);
}

template<typename IndexType, typename ValueType> void KluMatrixCore<IndexType, ValueType>::dumpVector(std::ostream& os, ValueType* v, int colw, int prec) {
    std::ios oldState(nullptr);
    oldState.copyfmt(std::cout);
    os << std::scientific << std::setprecision(prec);

    for(IndexType i=0; i<AN; i++) {
        if (i>0) {
            os << std::endl;
        }
        if constexpr(std::is_same<ValueType, Complex>::value) {
            os << std::setw(colw) << v[i].real();
            if (v[i].imag()>=0) {
                os << "+" << std::setw(colw-1) << v[i].imag();
            } else {
                os << std::setw(colw) << v[i].imag();
            }
            os << "i";
        } else {    
            os << std::setw(colw) << v[i];
        }
    }
    
    std::cout.copyfmt(oldState);
}

template<typename IndexType, typename ValueType> bool KluMatrixCore<IndexType, ValueType>::formatError(Status& s, NameResolver* resolver) const {
    std::string txt;
    IndexType row, col;
    switch (lastError) {
        case Error::Memory:
            s.set(Status::LinearSolver, "Out of memory.");
            return false;
        case Error::Defaults:
            s.set(Status::LinearSolver, "Cannot set up KLU defaults.");
            return false;
        case Error::Analysis:
            s.set(Status::LinearSolver, "KLU matrix analysis failed. Probably the matrix is singular.");
            return false;
        case Error::ReciprocalPivotGrowth:
            s.set(Status::LinearSolver, "Failed to compute reciprocal pivot growth.");
            return false;
        case Error::ReciprocalCondEstimate:
            s.set(Status::LinearSolver, "Failed to compute reciprocal condition number estimate.");
            return false;
        case Error::Solve:
            s.set(Status::LinearSolver, "Failed to solve factorized system.");
            return false;
        case Error::Factorization:
            txt = "Factorization failed, size="+std::to_string(AN);
            if (errorRank_>=0) {
                txt += ", rank="+std::to_string(errorRank_);
            }
            if (resolver) {
                txt += std::string(", zero pivot @ node '")+std::string((*resolver)(errorIndex))+"'" ;
            } else {
                txt += std::string(", zero pivot @ column ")+std::to_string(errorIndex+1);
            }
            txt += ".";
            s.set(Status::LinearSolver, txt);
            return false;
        case Error::Refactorization:
            txt = "Refactorization failed, size="+std::to_string(AN);
            if (errorRank_>=0) {
                txt += ", rank="+std::to_string(errorRank_);
            }
            txt += ".";
            s.set(Status::LinearSolver, txt);
            return false;
        case Error::MatrixInfNan:
            if (errorNan) {
                txt = "NaN found in matrix";
            } else {
                txt = "Inf found in matrix";
            }
            std::tie(row, col) = errorElement();
            if (resolver) {
                txt +=   ", row node '"+std::string((*resolver)(row))+"'"
                       + ", column node '"+std::string((*resolver)(col))+"'";
            } else {
                txt += ", row "+std::to_string(row+1)+", column "+std::to_string(col+1);
            }
            txt +=".";
            s.set(Status::LinearSolver, txt);
            return false;
        case Error::VectorInfNan:
            if (errorNan) {
                txt = "NaN found in vector";
            } else {
                txt = "Inf found in vector";
            }
            if (resolver) {
                txt += ", row node '"+std::string((*resolver)(errorIndex))+"'";
            } else {
                txt += ", row "+std::to_string(errorIndex+1);
            }
            txt += ".";
            return false;
    }
    return true;
}

// Instantiate template class for int32 and int64 indices, double and Complex values
template class KluMatrixCore<int32_t, double>;
template class KluMatrixCore<int32_t, Complex>;
template class KluMatrixCore<int64_t, double>;
template class KluMatrixCore<int64_t, Complex>;

}
