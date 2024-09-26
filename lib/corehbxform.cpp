#include <vector>
#include <algorithm>
#include <random>
#include "corehb.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

// Amplitude vector (2m-1 components), first component is DC, 
// the remaining components are real and imaginary parts of 
// the phasors corresponding to m-1 nonzero frequencies. 
// 
//   [A0 A1r A1i A2r A2i ... A(m-1)r A(m-1)i]

// Transformation of an amplitude vector to a vector of signal values 
// at n colocation points.
// n rows, 2m-1 columns
// Should satisfy n=2m-1, but during colocation points selection n>=2m-1.

// See chapter 2 in:
//   Kundert, White, Sangiovanni-Vincentelli: 
//   Steady-state methods for simulating analog and microwave circuits, 
//   Springer, 1990. 

// Notable difference: the book works with magnitudes of cosine and sine 
// contributions. We, however, work with real and imaginary part of a phasor. 
// Because 
//   Re(X exp(jwt)) = Xr cos(jwt) - Xi sin(jwt) 
// the sin(jwt)) contribution is replaced with -sin(jwt). 

// We use the same subtraction of multiples of 2 pi from the phase of base 
// frequencies as in the book. If trigonometric functions have a decent 
// implementation this should not be neccessary. 

bool HbCore::buildTransformMatrix(DenseMatrix<double>& XF, Status& s) {
    auto n = timepoints.size();
    auto m = freq.size();
    auto ncoef = 2*m-1;
    
    XF.resize(n, ncoef);
    
    // Storage of (2 pi f t) factors for base frequencies
    auto nBase = params.freq.size();
    double baseFac[nBase];
    
    // Loop through timepoints
    for(decltype(n) i=0; i<n; i++) {
        auto row = XF.row(i);
        auto t = timepoints[i];
        
        // For base frequencies compute phase at t (2 pi f t)
        // Subtract integer multiple of 2 pi to keep phase small
        for(decltype(nBase) k=0; k<nBase; k++) {
            auto prod = params.freq[k]*t; 
            auto frac = prod - std::trunc(prod);
            baseFac[k] = 2 * std::numbers::pi * frac;
        }

        // DC
        row.at(0) = 1.0;

        // Nonzero frequencies
        for(decltype(m) j=1; j<m; j++) {
            // Grid entry
            auto weights = grid.row(freq[j].gridIndex);
            // Assemble phase from base frequency contributions
            double phase = 0;
            for(decltype(nBase) k=0; k<nBase; k++) {
                phase += weights.at(k)*baseFac[k];
            }
            // Compute cosine and sine component
            row.at(j*2-1) =  std::cos(phase);
            row.at(j*2)   = -std::sin(phase);
        }
    }

    return true;
}

bool HbCore::buildAPFT(Status& s) {
    auto n = timepoints.size();
    auto m = freq.size();
    auto ncoef = 2*m-1;

    if (!buildTransformMatrix(IAPFT, s)) {
        return false;
    }
    
    // Make a copy that will be destroyed during matrix inversion
    DenseMatrix<double> coeffs = IAPFT;
    
    // Invert to obtain APFT
    APFT.resize(n, n);
    if (!coeffs.destructiveInvert(APFT)) {
        s.set(Status::Analysis, "Failed to compute forward transform matrix.");
        return false;
    }

    // Construct Omega matrix (for computing the derivative wrt. time on a spectrum)
    // This matrix is the time derivative matrix that operates on frequency domain vectors.  
    // Assumes the first component is DC magnitude and the remaining ones are (cosine, -sine) 
    // magnitudes, i,e, (Re, Im) parts of a phasor
    // First 6 rows and columns are
    //   0 0   0   0    0  .
    //   0 0  -w_1 0    0  .
    //   0 w_1 0   0    0  .
    //   0 0   0   0   -w_2 .
    //   0 0   0   w_2  0  .
    //   . .  .  .   .  .
    // where wi is (2 pi fi). 
    DenseMatrix<double> Omega(n, n);
    Omega.zero();
    for(decltype(m) i=1; i<m; i++) {
        auto base = 1+(i-1)*2;
        auto omega = 2*std::numbers::pi*freq[i].f;
        Omega.at(base, base+1) = -omega;
        Omega.at(base+1, base) = omega;
    }

    // Compute DDT matrix as IAPFT * Omega * APFT
    DDT.resize(n, n);
    DenseMatrix<double> tmp(n, n);
    IAPFT.multiply(Omega, tmp);
    tmp.multiply(APFT, DDT);

    // Reorganize DDT in column major form for better cache locality in solver
    DDTcolMajor.resize(n, n, DenseMatrix<double>::Major::Column);
    for(decltype(n) i=0; i<n; i++) {
        DDTcolMajor.row(i) = DDT.row(i);
    }

    return true;
}

bool HbCore::test() {
    Status s;
    bool ok = true;

    HbParameters p;
    p.freq = {1000, 100000};
    p.nharm = 4;
    p.truncate = "diamond";
    p.sample = "random";
    p.samplefac = 4;

    HbCore hb(p);

    if (ok && !hb.buildGrid(s)) {
        ok = false;
        std::cout << "Failed to build grid: " << s.message() << "\n";
    }
    
    if (ok && !hb.buildColocation(s)) {
        ok = false;
        std::cout << "Failed to select colocation points: " << s.message() << "\n";
    } 
    
    if (!hb.buildAPFT(s)) {
        ok = false;
        std::cout << "Failed to build APFT: " << s.message() << "\n";
    }

    if (ok) {
        double delta;
        auto n = hb.timepoints.size();
        DenseMatrix<double> result(n, n);

        // APFT*IAPFT
        hb.APFT.multiply(hb.IAPFT, result);
        DenseMatrix<double> I(n, n);
        I.identity();
        result.subtract(I, result);
        delta = result.maxAbs();
        std::cout << "APFT * IAPFT - I :: delta = " << delta << "\n";
        std::cout << "\n";
        if (delta>1e-12) {
            ok = false;
            std::cout << "IAPFT inverse failed\n";
        }

        // APFT of first non zero frequency cosine
        std::vector<double> v(n, 0.0);
        std::vector<double> vres(n, 0.0);
        auto f = hb.freq[1].f;
        auto mag = 10;
        for(size_t i=0; i<n; i++) {
            auto t = hb.timepoints[i];
            v[i] = mag*std::cos(2*std::numbers::pi*f*t);
        }

        auto vv = VectorView<double>(v);
        auto vvres = VectorView<double>(vres);
        hb.APFT.multiply(vv, vvres);
        auto norm = vvres.maxAbs();
        std::cout << "APFT of cosine at f1\n";
        vvres.dump(std::cout);
        std::cout << "\n";
        for(size_t i=0; i<n; i++) {
            if (
                i==1 && std::abs(vvres[i]-mag)/norm>1e-12 ||
                i!=1 && std::abs(vvres[i])/norm>1e-12
            ) {
                ok = false;
                std::cout << "APFT failed\n";
                break;
            }
        }

        // Derivative of first nonzero frequency cosine wrt time
        hb.DDT.multiply(vv, vvres);
        norm = vvres.maxAbs();
        std::cout << "APFT of DDT of cosine at f1\n";
        auto spec = std::vector<double>(n, 0);
        auto vspec = VectorView<double>(spec);
        hb.APFT.multiply(vvres, vspec);
        vspec.dump(std::cout);
        std::cout << "\n";
        for(size_t i=0; i<n; i++) {
            auto t = hb.timepoints[i];
            auto exact = -mag*2*std::numbers::pi*f*std::sin(2*std::numbers::pi*f*t);
            if (std::abs(vvres[i] - exact)/norm>1e-12) {
                ok = false;
                std::cout << "DDT failed\n";
                break;
            }
        }

    }

    if (!ok) {
        std::cout << "HB core test failed.\n";
    }

    return ok;
}

}
