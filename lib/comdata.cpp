#include "comdata.h"
#include "common.h"


namespace NAMESPACE {

CommonData::CommonData() {
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
    gmin = options.gmin;
    gshunt = options.gshunt;
}

void CommonData::resetTolerances(UnknownIndex n) {
    unknown_abstol.resize(n+1);
    unknown_idt_abstol.resize(n+1);
    residual_abstol.resize(n+1);
    residual_idt_abstol.resize(n+1);
    for(UnknownIndex i=0; i<=n; i++) {
        setTolerances(i, 
            std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 
            std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()
        );
    }
}

void CommonData::setTolerances(UnknownIndex i, double abstol, double idt_abstol, double res_abstol, double res_idt_abstol) {
    unknown_abstol[i] = abstol;
    unknown_idt_abstol[i] = idt_abstol;
    residual_abstol[i] = res_abstol;
    residual_idt_abstol[i] = res_idt_abstol;
}

void CommonData::updateTolerances(UnknownIndex i, double abstol, double idt_abstol, double res_abstol, double res_idt_abstol) {
    if (unknown_abstol[i]>abstol) {
        unknown_abstol[i] = abstol;
    }
    if (unknown_idt_abstol[i]>idt_abstol) {
        unknown_idt_abstol[i] = idt_abstol;
    }
    if (residual_abstol[i]>res_abstol) {
        residual_abstol[i] = res_abstol;
    }
    if (residual_idt_abstol[i]>res_idt_abstol) {
        residual_idt_abstol[i] = res_idt_abstol;
    }
}

void CommonData::defaultTolerances(UnknownIndex i, double abstol, double idt_abstol, double res_abstol, double res_idt_abstol) {
    if (std::isinf(unknown_abstol[i])) {
        unknown_abstol[i] = abstol;
    }
    if (std::isinf(unknown_idt_abstol[i])) {
        unknown_idt_abstol[i] = idt_abstol;
    }
    if (std::isinf(residual_abstol[i])) {
        residual_abstol[i] = res_abstol;
    }
    if (std::isinf(residual_idt_abstol[i])) {
        residual_idt_abstol[i] = res_idt_abstol;
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

std::tuple<double, double, double, double> CommonData::getTolerances(UnknownIndex i) {
    return std::make_tuple(
        unknown_abstol[i], 
        unknown_idt_abstol[i], 
        residual_abstol[i], 
        residual_idt_abstol[i]
    );
}

// Nature registry, one per circuit
// Issues ids for nature information (currently nature name)
// Discipline nature name is <discipline>.potential, <discipline>.flow
// Registry maps nature (id) -> number
// Each device registers its natures at load, receives nature id
// When asked for tolerances, it returns tolerance and nature id
// At tolerance setup id is stored in the corresponding vector
// If a tolerance is updated, the node's nature changes to the latest nature 
// After tolerances are collected ids are converted to numbers via registry
//
// Spice natures are always registered before natures from devices
// name              number  tol     idt nature
// ::spice_voltage   0       vntol   ::spice_flux
// ::spice_current   1       abstol  ::spice_charge
// ::spice_flux      2       fluxtol -
// ::spice_charge    3       chgtol  -
//
// Nature numbers are stored in CommonData vectors
// They are used for setting up global maxima (point and historic)

}

