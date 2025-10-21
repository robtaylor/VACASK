#ifndef __COMDATA_DEFINED
#define __COMDATA_DEFINED

#include "value.h"
#include "options.h"
#include "natures.h"
#include "common.h"


namespace NAMESPACE {

// Common data used during analysis and elaboration
typedef struct CommonData {
    CommonData           (const CommonData&)  = delete;
    CommonData           (      CommonData&&) = delete;
    CommonData& operator=(const CommonData&)  = delete;
    CommonData& operator=(      CommonData&&) = delete;
    
    // Homotopy support
    Real sourcescalefactor;
    Real gmin;     // gmin applied in parallel to nonlinear branches
    Real gdev;     // extra gmin applied during homotopy, 
                   // Computed and passed by the simulator. 
                   // Usually models don't use it
                   // therefore our homotopy algorithms modify only gmin. 
    Real gshunt;   // Conductance connected between potential nodes and the ground
    
    // Iteration counter
    Int iteration;

    // Analysis information
    String analysis_name;
    String analysis_type;
    
    // Bypass control
    bool allowContinueStateBypass;
    bool requestForcedBypass; 

    // Natures subset
    // This is a subset of all natures. These natures are used by this circuit. 
    NaturesSubset natures;

    // Static abstols for unknowns, component 0 is for gound
    std::vector<double> unknown_abstol;
    std::vector<double> unknown_idt_abstol;
    // Nature identifiers for unknown and its idt
    std::vector<NatureId> unknown_natureId;
    std::vector<NatureId> unknown_idt_natureId;
    // Nature indices for unknown and its idt
    std::vector<NaturesSubset::NatureIndex> unknown_natureIndex;
    std::vector<NaturesSubset::NatureIndex> unknown_idt_natureIndex;

    // Static abstols for residuals, component 0 is for gound
    std::vector<double> residual_abstol;
    std::vector<double> residual_idt_abstol;
    // Nature identifiers for residual and its idt
    std::vector<NatureId> residual_natureId;
    std::vector<NatureId> residual_idt_natureId;
    // Nature indices for residual and its idt
    std::vector<NaturesSubset::NatureIndex> residual_natureIndex;
    std::vector<NaturesSubset::NatureIndex> residual_idt_natureIndex;

    CommonData();
    void reset();
    void fromOptions(const SimulatorOptions& options);
    void resetTolerances(UnknownIndex n);
    void updateTolerances(UnknownIndex i, NatureTolerance u, NatureTolerance idt_u, NatureTolerance r, NatureTolerance idt_r, bool force=false);
    void defaultTolerances(UnknownIndex i, NatureTolerance u, NatureTolerance idt_u, NatureTolerance r, NatureTolerance idt_r);
    void setTolerances(UnknownIndex i, NatureTolerance u, NatureTolerance idt_u, NatureTolerance r, NatureTolerance idt_r);
    void scaleTolerances(double scl);
    bool enumerateNatures(Status& s=Status::ignore);
    std::tuple<NatureTolerance, NatureTolerance, NatureTolerance, NatureTolerance> getTolerances(UnknownIndex i);
} CommonData;

}

#endif
