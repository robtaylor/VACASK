#include <vector>
#include <algorithm>
#include <random>
#include "corehb.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

// Selection of colocation points for harmonic balance
// See Appendix B, chapter 2.1 in: 
//   Kundert, White, Sangiovanni-Vincentelli: 
//   Steady-state methods for simulating analog and microwave circuits, 
//   Springer, 1990. 

bool HbCore::buildColocation(Status& s) {
    auto hb_debug = 2;

    // Includes DC
    auto n = freq.size();

    // Must have 2 timepoints for each nonzero frequency and one for DC
    auto nt = 2*n-1;

    // Maximal frequency
    auto fmax = freq.back().f;
    auto fmin = freq[1].f;

    if (params.samplefac<1.0) {
        s.set(Status::BadArguments, "samplefac mus be >=1.");
    }
    
    // Number of samples
    size_t nsam = std::ceil(params.samplefac*nt);

    // Build initial sample pool
    auto range = params.nper/fmin;
    timepoints.clear();
    if (params.sample==HbCore::sampleUniform) {
        if (hb_debug>1) {
            Simulator::dbg() << "Generating pool of " << nsam << " uniformly distributed points, tmax=" << range << "\n" ;
        }
        for(decltype(nsam) i=0; i<nsam; i++) {
            timepoints.push_back(i/fmin);
        }
    } else if (params.sample==HbCore::sampleRandom) {
        std::mt19937_64 gen;
        gen.seed(1);
        std::uniform_real_distribution dist(0.0, 1.0);
        // Select across 3 periods of fmin
        if (hb_debug>1) {
            Simulator::dbg() << "Generating pool of " << nsam << " random points, tmax=" << range << "\n" ;
        }
        for(decltype(nsam) i=0; i<nsam; i++) {
            timepoints.push_back(dist(gen)*range);
        }
    } else {
        s.set(Status::BadArguments, "Unknown samplmode.");
        return false;
    }

    // Build fd->td transform matrix from timepoints
    buildTransformMatrix();
    
    // Number of cadidate rows
    auto ncand = XF.nRow();

    if (hb_debug>1) {
        Simulator::dbg() << "Colocation points (" << nt << "):\n" ;
    }

    // Number of kept rows, first row is always kept
    decltype(ncand) nkeep = 1;

    if (hb_debug>1) {
        Simulator::dbg() << "  t=" << timepoints[0] << " xform norm2=" << XF.row(0).norm2() << "\n" ;
    }

    // Keep nt best rows
    for(decltype(ncand) i=0; i<ncand-1; i++) {
        // Orthogonalize rows i+1 .. ncand-1, find the one with largest norm
        double maxNorm2 = 0;
        decltype(ncand) ndx = 0;
        auto wrt = XF.row(i);
        for(decltype(ncand) j=i+1; j<ncand; j++) {
            auto row = XF.row(j);
            row.orthogonalize(wrt);
            auto nrm2 = row.norm2();
            if (nrm2>maxNorm2) {
                maxNorm2 = nrm2;
                ndx = j;
            }
        }
        
        if (maxNorm2<=0) {
            s.set(Status::Analysis, "Zero norm encountered while computing colocation points.");
            return false;
        }

        // Swap row ndx with row i+1
        if (hb_debug>1) {
            Simulator::dbg() << "  t=" << timepoints[ndx] << " xform norm2=" << maxNorm2 << "\n" ;
        }
        XF.row(i+1).swap(XF.row(ndx));
        nkeep++;

        // Also swap timepoints
        auto tmp = timepoints[i+1];
        timepoints[i+1] = timepoints[ndx];
        timepoints[ndx] = tmp;
        
        // Stop if we have nt timepoints
        if (nkeep>=nt) {
            break;
        }
    }

    // Keep first nt timepoints
    timepoints.resize(nt);

    // Sort them
    std::sort(
        timepoints.begin(), timepoints.end(), 
        [](const double& t1, const double& t2) { return t1<t2; }
    );

    return true;
}

}
