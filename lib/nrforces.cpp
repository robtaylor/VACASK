#include "nrsolver.h"
#include "common.h"

namespace NAMESPACE {

Forces::Forces() {
}

void Forces::clear() {
    unknownValue_.clear();
    unknownForced_.clear();
    deltaValue_.clear();
    deltaIndices_.clear();
}

void Forces::dump(Circuit& circuit, std::ostream& os) const {
    auto n = unknownForced_.size();
    for(decltype(n) i=1; i<n; i++) {
        if (!unknownForced_[i]) {
            continue;
        }
        os << i << " : " << unknownValue_[i] << "\n";
    }
    auto nd = deltaIndices_.size();
    for(decltype(nd) i=0; i<nd; i++) {
        auto [u1, u2] = deltaIndices_[i];
        os << circuit.reprNode(u1)->name() << ", " << circuit.reprNode(u2)->name()
            << " : " << deltaValue_[i] << "\n"; 
    }
}

}
