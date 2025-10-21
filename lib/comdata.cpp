#include "comdata.h"
#include "common.h"


namespace NAMESPACE {

CommonData::CommonData() {
    reset();
}

void CommonData::reset() {
    // Default obtained from SimulatorOptions, also set by homotopy
    gmin = 0.0;
    gshunt = 0.0;

    // Set by homotopy (by default set so that they have no effect)
    gdev = 0.0;
    sourcescalefactor = 1.0;
    
    // Set by analysis
    analysis_name = "";
    analysis_type = "";
    allowContinueStateBypass = false; // Allow the analysis to force bypass in the first iteration of NR
                                      // with continuation when continueState is used without forcing it. 
                                      // This is set to true for all but the first point of the innermost sweep
                                      // if nr_contbypass is enabled and the innermost sweep allows continuation.  
    
    // Set by analysis core
    iteration = 0;
    requestForcedBypass = false;      // Request forced bypass in next NR iteration for all bypassable devices 
                                      // regardless of their converged state. 
}

void CommonData::fromOptions(const SimulatorOptions& options) {
    reset();
    gmin = options.gmin;
    gshunt = options.gshunt;
}

void CommonData::resetTolerances(UnknownIndex n) {
    unknown_abstol.resize(n+1);
    unknown_idt_abstol.resize(n+1);
    unknown_natureId.resize(n+1);
    unknown_idt_natureId.resize(n+1);
    residual_abstol.resize(n+1);
    residual_idt_abstol.resize(n+1);
    residual_natureId.resize(n+1);
    residual_idt_natureId.resize(n+1);
    for(UnknownIndex i=0; i<=n; i++) {
        setTolerances(i,
            NatureTolerance(
                NatureRegistry::noNature, 
                std::numeric_limits<double>::infinity()
            ), 
            NatureTolerance(
                NatureRegistry::noNature, 
                std::numeric_limits<double>::infinity()
            ), 
            NatureTolerance(
                NatureRegistry::noNature, 
                std::numeric_limits<double>::infinity()
            ), 
            NatureTolerance(
                NatureRegistry::noNature, 
                std::numeric_limits<double>::infinity()
            ) 
        );
    }
}

void CommonData::setTolerances(UnknownIndex i, NatureTolerance u, NatureTolerance idt_u, NatureTolerance r, NatureTolerance idt_r) {
    updateTolerances(i, u, idt_u, r, idt_r, true);
}

void CommonData::updateTolerances(UnknownIndex i, NatureTolerance u, NatureTolerance idt_u, NatureTolerance r, NatureTolerance idt_r, bool force) {
    if (force || unknown_abstol[i]>u.abstol) {
        unknown_abstol[i] = u.abstol;
        unknown_natureId[i] = u.id;
    }
    if (force || unknown_idt_abstol[i]>idt_u.abstol) {
        unknown_idt_abstol[i] = idt_u.abstol;
        unknown_idt_natureId[i] = idt_u.id;
    }
    if (force || residual_abstol[i]>r.abstol) {
        residual_abstol[i] = r.abstol;
        residual_natureId[i] = r.id;
    }
    if (force || residual_idt_abstol[i]>idt_r.abstol) {
        residual_idt_abstol[i] = idt_r.abstol;
        residual_idt_natureId[i] = idt_r.id;
    }
}

void CommonData::defaultTolerances(UnknownIndex i, NatureTolerance u, NatureTolerance idt_u, NatureTolerance r, NatureTolerance idt_r) {
    if (unknown_natureId[i]==NatureRegistry::noNature) {
        unknown_abstol[i] = u.abstol;
        unknown_natureId[i] = u.id;
    }
    if (unknown_idt_natureId[i]==NatureRegistry::noNature) {
        unknown_idt_abstol[i] = idt_u.abstol;
        unknown_idt_natureId[i] = idt_u.id;
    }
    if (residual_natureId[i]==NatureRegistry::noNature) {
        residual_abstol[i] = r.abstol;
        residual_natureId[i] = r.id;
    }
    if (residual_idt_natureId[i]==NatureRegistry::noNature) {
        residual_idt_abstol[i] = idt_r.abstol;
        residual_idt_natureId[i] = idt_r.id;
    }
}

void CommonData::scaleTolerances(double scl) {
    for(UnknownIndex i=0; i<unknown_abstol.size(); i++) {
        unknown_abstol[i] *= scl;
        unknown_idt_abstol[i] *= scl;
        residual_abstol[i] *= scl;
        residual_idt_abstol[i] *= scl;
    }
}

bool CommonData::enumerateNatures(Status& s) {
    auto n = unknown_natureId.size();
    unknown_natureIndex.resize(n+1);
    unknown_idt_natureIndex.resize(n+1);
    residual_natureIndex.resize(n+1);
    residual_idt_natureIndex.resize(n+1);
    for(decltype(n) i=0; i<n; i++) {
        unknown_natureIndex[i] = natures.natureIndex(unknown_natureId[i]);
        unknown_idt_natureIndex[i] = natures.natureIndex(unknown_idt_natureId[i]);
        residual_natureIndex[i] = natures.natureIndex(residual_natureId[i]);
        residual_idt_natureIndex[i] = natures.natureIndex(residual_idt_natureId[i]);
    }
    return true;
}

std::tuple<NatureTolerance, NatureTolerance, NatureTolerance, NatureTolerance> CommonData::getTolerances(UnknownIndex i) {
    return std::make_tuple(
        NatureTolerance(
            unknown_natureId[i], 
            unknown_abstol[i]
        ), 
        NatureTolerance(
            unknown_idt_natureId[i], 
            unknown_idt_abstol[i]
        ), 
        NatureTolerance(
            residual_natureId[i], 
            residual_abstol[i]
        ), 
        NatureTolerance(
            residual_idt_natureId[i], 
            residual_idt_abstol[i]
        )
    );
}

}

