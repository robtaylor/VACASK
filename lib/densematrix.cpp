#include "densematrix.h"
#include "common.h"

namespace NAMESPACE {

template<typename T> bool DenseMatrix<T>::test() {
    // Test status
    bool ok = true;

    // Small matrix inverse
    DenseMatrix<T> S({1.0, 2, 3, 4}, 2, 2);
    std::cout << "Small matrix\n";
    S.dump(std::cout);
    std::cout << "\n";
    DenseMatrix<T> S1 = S;
    DenseMatrix<T> SI(2, 2);
    S1.destructiveInvert(SI);
    DenseMatrix<T> SIexact({-2, 1, 1.5, -0.5}, 2, 2);
    std::cout << "Small matrix inverse\n";
    SI.dump(std::cout);
    std::cout << "\n";
    for(size_t i=0; i<2; i++) {
        bool exit = false;
        for(size_t j=0; j<2; j++) {
            if (std::abs(SI.at(i,j)-SIexact.at(i,j))>1e-12) {
                ok = false;
                exit = true;
                std::cout << "Small matrix inverse failed\n";
                break;
            }
        }
        if (exit) {
            break;
        }
    }

    // Small matrix right-multiplied by inverse
    DenseMatrix<T> B(2, 2);
    S.multiply(SI, B);
    std::cout << "Small matrix right-multiplied by inverse\n";
    B.dump(std::cout);
    std::cout << "\n";
    for(size_t i=0; i<2; i++) {
        bool exit = false;
        for(size_t j=0; j<2; j++) {
            if (
                i==j && std::abs(B.at(i,j)-1.0)>1e-12 ||
                i!=j && std::abs(B.at(i,j))>1e-12
            ) {
                ok = false;
                exit = true;
                std::cout << "Small matrix right-multiplied by inverse failed\n";
                break;
            }
        }
        if (exit) {
            break;
        }
    }

    size_t n = 4;

    // Invert a scaled identity
    DenseMatrix<T> D1(n, n);
    D1.zero();
    for(size_t i=0; i<n; i++) {
        D1.at(i, i) = i+1;
    }
    DenseMatrix<T> D = D1;
    std::cout << "Diagonal matrix\n";
    D1.dump(std::cout);
    std::cout << "\n";
    
    DenseMatrix<T> ID(n, n);
    D1.destructiveInvert(ID);
    std::cout << "Inverted diagonal matrix\n";
    ID.dump(std::cout);
    std::cout << "\n";
    for(size_t i=i; i<n; i++) {
        if (std::abs(ID.at(i, i)-1.0/(i+1))>1e-12) {
            ok = false;
            std::cout << "Inverted diagonal matrix failed\n";
            break;
        }
    }

    // Matrix multiply
    DenseMatrix<T> PR(n, n);
    D.multiply(ID, PR);
    std::cout << "Diagonal matrix right-multiplied by inverse\n";
    PR.dump(std::cout);
    std::cout << "\n";
    for(size_t i=i; i<n; i++) {
        if (std::abs(PR.at(i, i)-1.0)>1e-12) {
            ok = false;
            std::cout << "Diagonal matrix right-multiplied by inverse failed\n";
            break;
        }
    }

    // Generic matrix multiplication
    DenseMatrix<T> A(n, n);
    for(size_t i=0; i<n; i++) {
        for(size_t j=0; j<n; j++) {
            A.at(i,j) = 1.0/(i+j+1);
        }
    }
    std::cout << "Generic matrix\n";
    A.dump(std::cout);
    std::cout << "\n";

    // Right-multiply by diagonal
    A.multiply(D, PR);
    std::cout << "Generic matrix right-multiplied by diagonal matrix\n";
    PR.dump(std::cout);
    std::cout << "\n";
    for(size_t i=0; i<n; i++) {
        for(size_t j=0; j<n; j++) {
            if (std::abs(PR.at(i,j)-A.at(i,j)*D.at(j,j))>1e-12) {
                ok = false;
                std::cout << "Generic matrix right-multiplied by diagonal matrix failed\n";
                break;
            }
        }
    }

    // Left-multiply by diagonal
    D.multiply(A, PR);
    std::cout << "Generic matrix left-multiplied by diagonal matrix\n";
    PR.dump(std::cout);
    std::cout << "\n";
    for(size_t i=0; i<n; i++) {
        for(size_t j=0; j<n; j++) {
            if (std::abs(PR.at(i,j)-A.at(i,j)*D.at(i,i))>1e-12) {
                ok = false;
                std::cout << "Generic matrix left-multiplied by diagonal matrix failed\n";
                break;
            }
        }
    }

    // Invert
    DenseMatrix<T> A1 = A;
    DenseMatrix<T> IA(n, n);
    A1.destructiveInvert(IA);
    A.multiply(IA, PR);
    std::cout << "Generic matrix right-multiplied by inverse matrix\n";
    PR.dump(std::cout);
    std::cout << "\n";
    for(size_t i=0; i<n; i++) {
        for(size_t j=0; j<n; j++) {
            if (
                i==j && std::abs(PR.at(i,j)-1.0)>1e-12 ||
                i!=j && std::abs(PR.at(i,j))>1e-12    
            ) {
                ok = false;
                std::cout << "Generic matrix right-multiplied by inverse matrix failed\n";
                break;
            }
        }
    }


    std::cout << "Dense matrix test " << (ok ? "OK" : "FAILED") << "\n";
    return ok;
}

template class DenseMatrix<double>;
template class DenseMatrix<Complex>;

}
