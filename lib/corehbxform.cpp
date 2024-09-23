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

// Notable difference: the book works with sums of cosine and sine contributions. 
// We, however, work with real and imaginary part of a phasor. Therefore the 
// sine contribution is negative. 

// We use the same subtraction of multiples of 2 pi from the phase of base 
// frequencies as in the book. If trigonometric functions have a decent 
// implementation this should not be neccessary. 

bool HbCore::buildTransformMatrix(Status& s) {
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


// Tranformation of an amplitude vector to a vector of 
// signal derivative values at n colocation points. 
// n rows, 2m-1 columns
// Should satisfy n=2m-1. 
// We avoid recomputing sines and cosines. 
// Must call buildTransformMatrix() first. 

bool HbCore::buildDdtTransformMatrix(Status& s) {
    auto nt = XF.nRows();
    auto nc = XF.nCols();
    auto nf = freq.size();

    if (nt!=nc) {
        s.set(Status::Analysis, "Number of rows and number of columns in transform matrix do not match.");
        return false;
    }
    
    if (nf*2-1 != nc) {
        s.set(Status::Analysis, "Number of frequencies and number FD coefficients do not match.");
        return false;
    }

    XFdot.resize(nt, nc);
    for(decltype(nt) i=0; i<nt; i++) {
        auto rowXF = XF.row(i);
        auto rowXFdot = XFdot.row(i);
        
        // DC 
        rowXFdot[0] = 0.0;
        
        // Nonzero frequencies
        for(decltype(nc) j=1; j<nf; j++) {
            auto omega = 2 * std::numbers::pi * freq[j].f;
            rowXFdot[2*j-1] =   omega * rowXF[2*j];   // -sin
            rowXFdot[2*j]   = - omega * rowXF[2*j-1]; // -cos
        }
    }   

    return true;
}

}
