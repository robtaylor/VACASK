#ifndef __COMDATA_DEFINED
#define __COMDATA_DEFINED

#include "value.h"
#include "options.h"
#include "common.h"


namespace NAMESPACE {

    // Common data used during analysis and elaboration
typedef struct CommonData {
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

    // Static abstols for unknowns, component 0 is for gound
    std::vector<double> unknown_abstol;
    std::vector<double> unknown_idt_abstol;

    // Static abstols for residuals, component 0 is for gound
    std::vector<double> residual_abstol;
    std::vector<double> residual_idt_abstol;

    CommonData();
    void fromOptions(const SimulatorOptions& options);
    void resetToleranceVectors(UnknownIndex n);
    void updateTolerances(UnknownIndex i, double abstol, double idt_abstol, double res_abstol, double res_idt_abstol);
    void defaultTolerances(UnknownIndex i, double abstol, double idt_abstol, double res_abstol, double res_idt_abstol);
    void setTolerances(UnknownIndex i, double abstol, double idt_abstol, double res_abstol, double res_idt_abstol);
} CommonData;

}

#endif
