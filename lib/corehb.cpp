#include <vector>
#include <algorithm>
#include "core.h"
#include "corehb.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

Id HbCore::truncateRaw = Id::createStatic("raw");
Id HbCore::truncateBox = Id::createStatic("box");
Id HbCore::truncateDiamond = Id::createStatic("diamond");

Id HbCore::sampleUniform = Id::createStatic("uniform");
Id HbCore::sampleRandom = Id::createStatic("random");

HbParameters::HbParameters() {
    truncate = HbCore::truncateDiamond;
    sample = HbCore::sampleRandom;
}

// Analysis asks cores if they request a rebuild. 
// HB core replies that it does if the spectrum changes. 
// Along with changed spectrum this function recomputes
// - colocation points
// - transform matrices
std::tuple<bool, bool> HbCore::requestsRebuild(Status& s) {
    auto oldSpectrum = spectrum;
    
    // Recompute spectrum
    if (!buildGrid(s)) {
        return std::make_tuple(false, false);
    }

    // See if it changed
    bool changed = oldSpectrum.size() != spectrum.size();
    if (!changed) {
        auto nf = spectrum.size();
        for(decltype(nf) i=0; i<nf; i++) {
            if (oldSpectrum[i]!=spectrum[i]) {
                changed = true;
                break;
            }
        }
    }

    if (changed) {
        // Recompute colocation
        if (!buildColocation(s)) {
            return std::make_tuple(false, changed);
        }

        // Recompute transforms
        if (!buildTransformMatrix(s)) {
            return std::make_tuple(false, changed);
        }

        if (!buildDdtTransformMatrix(s)) {
            return std::make_tuple(false, changed);
        }
    }

    return std::make_tuple(true, changed);
}

}
