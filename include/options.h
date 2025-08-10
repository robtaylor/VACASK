#ifndef __OPTIONS_DEFINED
#define __OPTIONS_DEFINED

#include "parameterized.h"
#include "value.h"
#include "common.h"


namespace NAMESPACE {

// Simulator options
typedef struct SimulatorOptions  {
    Real temp;
    Real tnom;
    Real gmin; 
    Real gshunt; 
    Real minr;
    Real scale;
    Real reltol; 
    Real abstol;
    Real vntol;
    Real chgtol;
    Real fluxtol;
    Id relrefsol;
    Id relrefres;
    Id relreflte;
    Id relref;
    Real restol;
    Real vnrestol;
    Int matrixcheck;
    Int rhscheck;
    Int solutioncheck;
    Real rcondcheck;
    Int sweep_pointmarker;
    Int sweep_debug;
    int nr_debug;
    Int nr_bypass; 
    Real nr_convtol;
    Real nr_bypasstol;
    Int nr_conviter;
    Int nr_residualcheck;
    Real nr_damping;
    Real nr_force;
    Int nr_contbypass;
    Int homotopy_debug;
    Int homotopy_gminsteps;
    Int homotopy_srcsteps;
    Real homotopy_gminfactor;
    Real homotopy_startgmin;
    Real homotopy_maxgmin;
    Real homotopy_mingmin;
    Real homotopy_maxgminfactor;
    Real homotopy_mingminfactor;
    Real homotopy_srcstep;
    Real homotopy_srcscale;
    Real homotopy_minsrcstep;
    Real homotopy_sourcefactor;
    Int op_debug; 
    Int op_itl;
    Int op_itlcont;
    Int op_skipinitial;
    std::vector<Id> op_homotopy;
    std::vector<Id> op_srchomotopy;
    Int op_nsiter;
    Int smsig_debug;
    Int tran_debug;
    Id tran_method;
    Int tran_maxord;
    Real tran_fs;
    Real tran_ffmax;
    Real tran_fbr;
    Real tran_rmax;
    Int tran_minpts;
    Int tran_itl;
    Real tran_ft;
    Int tran_predictor;
    Real tran_redofactor;
    Real tran_lteratio;
    Int tran_spicelte;
    Real tran_xmu;
    Int tran_trapltefilter;
    Int hb_debug;
    Int hb_itl;
    Int hb_itlcont;
    Int hb_skipinitial;
    std::vector<Id> hb_homotopy;
    Id rawfile;
    Int strictoutput;
    Int strictsave;
    Int strictforce;
    Int accounting;
    
    SimulatorOptions();

    // Isn't C++20 great?
    bool operator==(const SimulatorOptions& other) const = default;
    bool operator!=(const SimulatorOptions& other) const = default; 

    static Id relrefPointLocal;
    static Id relrefLocal;
    static Id relrefPointGlobal;
    static Id relrefGlobal;
    static Id relrefRelref;
    
    static Id relrefAlllocal;
    static Id relrefSigglobal;
    static Id relrefAllglobal;

    static Id rawfileAscii;
    static Id rawfileBinary; 

    // Options that affect mapping (node collapsing)
    static std::unordered_map<Id, ParameterIndex> mappingAffectingOptions;
    // Options that affect hierarchical parameters 
    static std::unordered_map<Id, ParameterIndex> parametrizationAffectingOptions;
    // Options that affect topology
    static std::unordered_map<Id, ParameterIndex> hierarchyAffectingOptions;

    static bool staticInitialize();

    bool optionsDiffer(std::unordered_map<Id, ParameterIndex>& optionsList, SimulatorOptions& opt);
} SimulatorOptions;


// Simulator internals
// TODO: make this a local temporary structure that is created in analysis, 
//       passed on to the core, and filled out by analysis and core.
//       Problem is in setup() that needs this structure. Higher up 
//       setup is called by elaborate() and elaborateChanges(), which in 
//       turn are also called from the command interpreter. 
typedef struct SimulatorInternals {
    Real sourcescalefactor;
    Real gmin;     // gmin applied in parallel to nonlinear branches
    Real gdev;     // extra gmin applied during homotopy, 
                   // Computed and passed by the simulator. 
                   // Usually models don't use it
                   // therefore our homotopy algorithms modify only gmin. 
    Real gshunt;   // Conductance connected between potential nodes and the ground
    Int iteration;
    String analysis_name;
    String analysis_type;
    // String cwd;

    bool allowContinueStateBypass;
    bool requestForcedBypass; 
    
    SimulatorInternals();
    void fromOptions(const SimulatorOptions& options);
} SimulatorInternals;

}

#endif
