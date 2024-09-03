#include <vector>
#include <algorithm>
#include "corehb.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

// Generation of harmonics and intermodulation products
// See Chapter 1.1 in: 
//   Kundert, White, Sangiovanni-Vincentelli: 
//   Steady-state methods for simulating analog and microwave circuits, 
//   Springer, 1990. 

bool HbCore::buildGrid(Status& s) {
    auto n = params.freq.size();

    Int hb_debug = 2;
    auto hb_freqtol = 1e-14;
    
    if (params.freq.size()<1) {
        s.set(Status::BadArguments, "freq must have at least one component.");
        return false;
    }
    
    auto fundamentals = params.freq.data();
    
    // Check freq
    for(decltype(n) i=0; i<n; i++) {
        if (params.freq[i]==0.0) {
            s.set(Status::BadArguments, "Zero frequency should not be specified explicitly.");
            return false;
        }
    }

    // Check nharm
    Int* nharm;
    Int nharmCount;
    Int nharmMax = 0;
    Int nharmScalar;
    bool nharmCommon = false;
    if (params.nharm.type()==ValueType::Int) {
        if (params.nharm.val<Int>()<=0) {
            s.set(Status::BadArguments, "nharm must be >0.");
            return false;
        }
        nharmScalar = params.nharm.val<Int>();
        nharmCommon = true;
        nharmMax = nharmScalar;
    } else if (params.nharm.type()==ValueType::IntVec) {
        for(auto nh : params.nharm.val<IntVector>()) {
            if (nh<=0) {
                s.set(Status::BadArguments, "nharm components must be >0.");
                return false;
            }
            if (nh>nharmMax) {
                nharmMax = nh;
            }
        }
        nharm = params.nharm.val<IntVector>().data();
        nharmCount = params.nharm.val<IntVector>().size();
        if (nharmCount!=n) {
            s.set(Status::BadArguments, "Number of nharm components must match number of freq components.");
            return false;
        }
    } else {
        s.set(Status::BadArguments, "nharm must be an integer or an integer vector.");
        return false;
    }

    // Build grid
    if (params.truncate==HbCore::truncateRaw) {
        // Raw
        // Check imorder
        if (params.imorder.size()!=0 && params.imorder.size()!=n) {
            s.set(Status::BadArguments, "If given imorder must match freq in size.");
            return false;
        }
        grid.resize(n+1, n); // first entry is DC (0Hz)
        grid.zero();
        freq.push_back({
            .gridIndex = 0, 
            .f = 0, 
            .order = 0, 
            .isHarmonic = true
        });
        for(decltype(n) i=0; i<n; i++) {
            grid.at(i+1, i) = 1.0;
            Int order = -1;
            if (params.imorder.size()>0) {
                order = params.imorder[i];
                if (order<0) {
                    s.set(Status::BadArguments, "imorder components must be >=0.");
                    return false;
                }
            }
            freq.push_back({
                .gridIndex = i+1, 
                .f = fundamentals[i], 
                .order = order, 
                .isHarmonic = false, 
            });
        }
    } else if (params.truncate==HbCore::truncateBox || params.truncate==HbCore::truncateDiamond) {
        // Box and diamond
        grid.resize(0, n); // Empty table
        std::vector<Int> cnt(n);
        std::vector<Int> end(n);
        Int lastChanged = 0;
        cnt[0] = 0;
        end[0] = nharmCommon ? nharmScalar+1 : nharm[0]+1;

        // Compute immax
        auto immax = params.immax>0 ? params.immax : nharmMax;
        
        while (true) {
            // Do we need to build ranges
            if (lastChanged<n-1) {
                // Check if all previous coordinates are 0
                bool allZero = true;
                for(decltype(lastChanged) i=0; i<lastChanged+1; i++) {
                    if (cnt[i] != 0) {
                        allZero = false;
                        break;
                    }
                }
                // From lastChanged+1 to n-1, set up ranges
                for(decltype(lastChanged) i=lastChanged+1; i<n; i++) {
                    if (allZero) {
                        cnt[i] = 0;
                    } else {
                        cnt[i] = nharmCommon ? -nharmScalar : -nharm[i];
                    }
                    end[i] = nharmCommon ? nharmScalar+1 : nharm[i]+1;
                }
            }

            // Compute properties
            Int order = 0;
            double f = 0;
            Int nnz = 0;
            for(decltype(n) i=0; i<n; i++) {
                order += std::abs(cnt[i]);
                f += cnt[i]*fundamentals[i];
            }
            f = std::abs(f);

            // Check immax
            // Not optimal for diamond truncation because we traverse the whole box and 
            // leave out frequencies with order above imax. 
            // But then again, HB spends a lot more time solving the problem. 
            if (!(params.truncate==HbCore::truncateDiamond && order>immax)) {
                // Construct component
                auto row = grid.addRow();
                for(decltype(n) i=0; i<n; i++) {
                    row.at(i) = cnt[i];
                    if (cnt[i]!=0) {
                        nnz++;
                    }
                }
                freq.push_back({
                    .gridIndex = grid.nRow()-1, 
                    .f = f, 
                    .order = order, 
                    .isHarmonic = nnz<=1, 
                });
            }

            // Advance, count up because size_t is unsigned
            for(decltype(n) i=0; i<n; i++) {
                lastChanged = n-1-i;
                cnt[lastChanged]++;
                if (cnt[lastChanged]<end[lastChanged]) {
                    break;
                }
            }

            // Check if done
            if (cnt[0]>=end[0]) {
                // Done
                break;
            }
        }
    } else {
        s.set(Status::BadArguments, "Unknown spectrum truncation method.");
        return false;
    }

    if (hb_debug>1) {
        Simulator::out() << "Raw HB frequency grid\n";
        auto nn = grid.nRow();
        for(decltype(nn) i=0; i<nn; i++) {
            auto row = grid.row(i);
            
            std::cout << "  #" << i << " [";
            auto nel = row.n();
            for(decltype(nel) j=0; j<nel; j++) {
                std::cout << row.at(j) << " ";
            }
            std::cout << "]";
            
            std::cout << " f=" << freq[i].f;
            std::cout << " order=" << freq[i].order;
            if (freq[i].isHarmonic) {
                std::cout << " harmonic";
            }
            std::cout << "\n";
        }
    }


    // Remove duplicate frequencies
    // Lower order im products are kept over higher order ones
    // Harmonics are kept over im products
    // Lower index is kept over higher index
    auto nf = freq.size();
    std::vector<bool> removed(nf, false);
    for(decltype(nf) i=0; i<nf-1; i++) {
        // Is i removed
        if (removed[i]) {
            // Go to next i
            continue;
        }
        for(decltype(nf) j=i+1; j<nf; j++) {
            // Is j removed
            if (removed[j]) {
                // Go to next j
                continue;
            }
            // std::cout << i << " " << j << " " << freq[i].f << " " << freq[j].f << "\n";
            auto df = std::abs(freq[i].f - freq[j].f);
            auto tol = std::max(freq[i].f, freq[j].f);
            if (df <= tol*hb_freqtol) {
                // Frequencies match, compare order
                if (freq[i].order < freq[j].order) {
                    // j has higher order, keep i, mark j as removed, continue
                    removed[j] = true;
                    if (hb_debug>1) {
                        Simulator::out() << "Removing #" << j << " (higher order)\n";
                    }
                 } else if (freq[i].order > freq[j].order) {
                    // i has higher order, keep j, mark i as removed, exit inner loop
                    removed[i] = true;
                    if (hb_debug>1) {
                        Simulator::out() << "Removing #" << i << " (higher order)\n";
                    }
                    break;
                } else {
                    // Same order
                    // Check harmonic
                    if (freq[i].isHarmonic && !freq[j].isHarmonic) {
                        // i is harmonic, j is not, keep i, mark j as removed, continue
                        removed[j] = true;
                        if (hb_debug>1) {
                            Simulator::out() << "Removing #" << j << " (not harmonic)\n";
                        }
                    } else if (!freq[i].isHarmonic && freq[j].isHarmonic) {
                        // j is harmonic, i is not, keep j, mark i as removed, exit inner loop
                        removed[i] = true;
                        if (hb_debug>1) {
                            Simulator::out() << "Removing #" << i << " (not harmonic)\n";
                        }
                        break;
                    } else {
                        // Harmonic satus is the same, keep the one with lower index (i), mark j as removed
                        removed[j] = true;
                        if (hb_debug>1) {
                            Simulator::out() << "Removing #" << j << " (higher index)\n";
                        }
                    }
                }
            }
        }
    }

    // Rebuild freq vector
    decltype(nf) dest = 0;
    for(decltype(nf) i=0; i<nf; i++) {
        if (!removed[i]) {
            if (i!=dest) {
                freq[dest] = freq[i];
            }
            dest++;
        }
    }
    freq.resize(dest);

    if (dest<2) {
        s.set(Status::BadArguments, "Too few frequencies in spectrum.");
        return false;
    }

    // Sort freq vector by frequency
    std::sort(
        freq.begin(), freq.end(), 
        [](const SpecFreq& a, const SpecFreq& b) { return a.f < b.f; }
    );

    if (hb_debug>0) {
        Simulator::out() << "HB spectrum, " << freq.size() << " frequencies\n";
        auto nn = grid.nRow();
        for(auto& fd : freq) {
            std::cout << "  #" << fd.gridIndex << " [";
            auto row = grid.row(fd.gridIndex);
            auto nel = row.n();
            for(decltype(nel) j=0; j<nel; j++) {
                std::cout << row.at(j) << " ";
            }
            std::cout << "]";
            
            std::cout << " f=" << fd.f;
            std::cout << " order=" << fd.order;
            if (fd.isHarmonic) {
                std::cout << " harmonic";
            }
            std::cout << "\n";
        }
    }
    
    return true;
}

}
