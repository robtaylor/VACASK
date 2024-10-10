#include "options.h"
#include "introspection.h"
#include "simulator.h"
#include "homotopy.h"
#include "common.h"


namespace NAMESPACE {

Id SimulatorOptions::relrefPointLocal = Id::createStatic("pointlocal"); // at given time for each unknown separately
Id SimulatorOptions::relrefLocal = Id::createStatic("local"); // maximum over past time for each unknown separately
Id SimulatorOptions::relrefPointGlobal = Id::createStatic("pointglobal"); // at given time, maximum over all unknowns
Id SimulatorOptions::relrefGlobal = Id::createStatic("global"); // maximum over past time, maximum over all unknowns
Id SimulatorOptions::relrefRelref = Id::createStatic("relref"); // let relref option decide on tolerance reference

Id SimulatorOptions::relrefAlllocal = Id::createStatic("alllocal"); // relref value alllocal
Id SimulatorOptions::relrefSigglobal = Id::createStatic("sigglobal"); // relref value sigglobal
Id SimulatorOptions::relrefAllglobal = Id::createStatic("allglobal"); // relref value allglobal

Id SimulatorOptions::rawfileAscii = Id::createStatic("ascii");
Id SimulatorOptions::rawfileBinary = Id::createStatic("binary");

// Default options
// TODO: options validation
SimulatorOptions::SimulatorOptions() { 
    temp = 27; // ambient temperature in C
    tnom = 27; // default device parameter measurement temperature in C 
    gmin = 1e-12; // >=0, device gmin
    gshunt = 0.0; // >=0, 0 = off, >0 shunt conductance from potential nodes to ground
    minr = 0.0; // >=0
    scale = 1.0; // >0
    reltol = 1e-3; // 0<x<1, Relative tolerance 
    abstol = 1e-12; // >0, absolute current tolerance in A
    vntol = 1e-6; // >0, absolute voltage tolerance in V
    chgtol = 1e-15; // >0, charge tolerance in As, default is 1mV across 1pF
    fluxtol = 1e-15; // >0, flux tolerance in Vs, default is 1uA across 1nH
    relrefsol = relrefRelref; // reference value for solution delta reltol
                             // pointlocal = separate for each unknown, each timepoint
                             // local = separate for each unknown, maximum over past timepoints
                             // pointglobal = maximum over all unknowns, separate for each timepoint
                             // global = maximum over all unknowns, maximum over past timepoints
                             // relref = let recoreopnrlref option decide 
    relrefres = relrefRelref; // reference value for residual reltol
    relreflte = relrefRelref; // reference value for lte reltol
    relref = relrefAlllocal;
    restol = 1e-12; // >0, residual tolerance (A, applied to potential nodes)
    vnrestol = 1e-6; // >0, residual tolerance (V, applied to flow nodes)

    matrixcheck = 0; // check matrix for inf and nan
    rhscheck = 1; // check rhs vector for inf and nan
    solutioncheck = 1; // check solution vector for inf and nan
    
    sweep_pointmarker=0; // 1 = yelds analysis execution in stepped mode before each sweep point
                         // Should be enabled in cosimulation when the analog simulator is the slave. 
                         // Before each sweep point Analysis::resume() will stop and return 
                         // AnalysisState::SweepPoint. At that point the digital simulator should reset
                         // its state to the initial state at t=0. 
    sweep_debug = 0; // 1 = debug sweep, >=2 print details
    
    nr_debug = 0;  // >0  enables nonlinear solver debugging
                   // >=1 print messages
                   // >=2 print linear system
                   // >=3 print new solution
                   // >=4 print old solution
    nr_bypass = 0; // 1 = disable core evaluation for bypassed instances 
                   // To be bypassed an instance must converge and 
                   // the instance inputs change between iterations must be within tolerances. 
                   // As soon as the inputs change is outside tolerances instance is no longer bypassed. 
    nr_convtol = 0.01;   // Tolerance factor applied to residuals for instance convergence check. Should be <1.
                         // 1.0 corresponds to set tolerances (abstol, vntol, ...), <1.0 makes them more strict. 
    nr_bypasstol = 0.01; // Tolerance factor applied to instance inputs for instance bypass check. Should be <1.
                         // 1.0 corresponds to set tolerances (abstol, vntol, ...), <1.0 makes them more strict. 
    nr_conviter = 1; // >0, number of consecutive convergent iterations before convergence is confirmed
    nr_residualcheck = 1; // check residual beside unknowns change to establish convergence 
    nr_damping = 1.0; // 0<x<=1, Newton-Raphson damping factor (<=1)
    nr_force = 1e5; // x>0, forcing factor for nodesets and initial conditions
    nr_contbypass = 1; // allow forced bypass of instance evaluation 
                       // in the first NR iteration when continuation mode is enabled

    homotopy_debug = 0; // >0 enables homotopy debugging
    homotopy_gminsteps = 100; // >1, <=0 disables gmin stepping
    homotopy_srcsteps = 100; // >1, <=0 disables source stepping
    homotopy_gminfactor = 10.0; // >0, initial gmin stepping factor for dynamic gmin stepping
    homotopy_maxgminfactor = 10.0; // >=op_gminfactor, maximal gmin stepping factor for dynamic gmin stepping
    homotopy_mingminfactor = 1.00005; // 1<x<op_gminfactor, maximal gmin stepping factor for dynamic gmin stepping
                                      // give up when gmin step fails and factor falls below this value
    homotopy_startgmin = 1e-3; // >mingmin, value at which dynamic gmin stepping starts
    homotopy_maxgmin = 1e2;   // >mingmin, if op dynamic gmin stepping failes to solve the circuit above this value of gmin, it fails
    homotopy_mingmin = 1e-15; // >0, value where dynamic gmin stepping stops if gmin/gshunt are set to 0
    homotopy_srcstep = 0.001; // >0, initial source step for dynamic source stepping
    homotopy_minsrcstep = 1e-7; // 0<x<srcstep, source step at which dynamic source stepping gives up
    homotopy_sourcefactor = 1.0; // x = homotopy_sourcefactor * sourcescalefactor where sourcescalefactor is set by 
                                 // the source stepping homotopy algorithm is the scaling factor for all independent sources. 
                                 // This option makes it possible to set up manual homotopy via dc sweep. 
                                 // For normal simulation this should be 1.0. 
    
    op_debug = 0; // >0 enables op analysis debugging
                  // >=1 print iteration type, homotopy information, convergence report
                  // >=2 print continuation mode information
    op_itl = 100;  // >0, maximal number of iterations in non-continuation mode
    op_itlcont = 50; // >0, maximal number of iterations in continuation mode

    op_skipinitial = 0; // 1 = no initial op, go straight to homotopy
    op_homotopy = { "gdev", "gshunt", "src" };
                        // list of homotopy algorithms to apply in operating point analysis
    op_srchomotopy = { "gdev", "gshunt" };
                        // list of homotopy algorithms to apply in operating point analysis
                        // when source stepping fails at sourcefactor=0
    
    op_nsiter = 1; // Number of iterations during which nodesets are applied


    smsig_debug = 0; // 1=debug small signal analyses (dcinc, dcxf, ac, acxf, noise), 
                     // >=100 print linear system
                     // >=101 print matrix before matrix checks are performed
    
    tran_debug = 0; // 1 = debug steps, 2 = debug solver
    tran_method = "trap"; // am, bdf, gear - Adams-Moulton, backward differentiation (Gear)
                          //   (use maxord to limit maximal integration order)
                          // euler (ignore maxord, use am maxord=1), 
                          // trap, am2 (ignore maxord, use am maxord=2)
                          // gear2, bdf2 (ignore maxord, use bdf maxord=2)
    tran_maxord = 2;
    tran_fs = 0.25; // 0<x<=0.5, fraction of timestep for first timepoint and 
                    // for maximal timestep between two breakpoints
                    // Used for multiplying step to get first timestep. 
    tran_ffmax = 0.25; // >=0, fraction of maximal excitation frequency's period
                       // to which the initial time step is limited. 
                       // Setting it to 0 turns off step limiting with maximal frequency. 
    tran_fbr = 0.2501; // 0<x<=1/3, fraction of distance between consecutive breakpoints
                     // to which the step is limited from above. 
                     // To guarantee at least 3 points between two breakpoints it 
                     // must be below 1/3. 0.2501 guarantees 4 points with the last 
                     // timestep slightly shorter to avoid step halving before breakpoint. 
                     // We need 2 past points for LTE estimation of a first order algorithm
                     // and also 2 past points for first order predictor. 
                     // Because the point at t=0 with UIC is not good, we cannot use it
                     // and we need a minimum of 2 intermediate points. This means we need 
                     // a step that is at most 1/3 of the distance between breakpoints. 
                     // Therefore 0<tran_fbr<=1/3. 
    tran_rmax = 0.0; // ratio of step to timestep upper limit
                     // <1 means that the timestep has no upper limit based on step
    tran_minpts = 50; // Minimal number of timepoints from start to stop (excluding start). 
                      // <1 means that timestep is not limited by simulation interval
                      // If rmax and minpts are both disabled, the simulation interval is 
                      // the upper limit to timestep. 
                      // The upper limit to timestep can be further reduced by specifying maxstep. 
    tran_itl = 10;  // >0, maximal number of iterations per timepoint
    tran_ft = 0.25; // 0<x<1, factor for cutting the timestep when iterations>tran_itl
    tran_predictor = 0; // 1 .. use predictor for computing initial guess 
                        // 0 .. use previous solution as initial guess
    tran_redofactor = 2.5; // If timestep/lte_computed_timestep > redofactor reject timepoint. 
                           // No LTE-based rejections take place if redofactor=0
    tran_lteratio = 3.5; // x>1, LTE overestimation factor (greater values mean more loose LTE tolerance)
    tran_spicelte = 0; // 0 .. correct LTE handling
                       // 1 .. (incorrect) SPICE-like LTE handling
    tran_xmu = 0.5; // for trapezoidal algorithm (Adams-Moulton of order 2) chooses a mixture between trapezoidal and Euler
                    // 0 = pure Euler, 0.5 = pure trapezoidal
    tran_trapltefilter = 1; // enable trap ringing filter for predictor and LTE computation, 
                            // applied only when Adams-Moulton algorithm of order 2 is used 
    
    hb_debug = 0; // 0 = none, 1 = solver and homotopy runs, 2 = homotopy and continuation internals
    hb_itl = 100; // >0, maximal number of iterations in non-continuation mode
    hb_itlcont = 50; // >0, maximal number of iterations in continuation mode
    hb_skipinitial = 0; // 1 = no initial hb, go straight to homotopy
    hb_homotopy = { "src" }; // list of homotopy algorithms to apply in hb analysis
    
    rawfile = "binary"; // ascii or binary
    strictoutput = 2; // 0 = leave output files in place after error, 
                      // 1 = delete output files on error
                      // 2 = delete output files before analysis
    strictsave = 1; // 0 = failure to bind an output binds to constant 0
                    // 1 = if cannot bind at first binding (ordinary analysis, first sweep point), 
                    //     signal error, later failed binding attempts bind to constant 0
                    // 2 = if cannot bind always signal error
    strictforce = 1; // 0 = ignore conflicting nodeset/ic and print warning
                     // 1 = abort on nodeset/ic conflicts
    accounting = 0; // 0 = accounting enabled outside analyses only
                    // 1 = accounting enabled outside and inside analyses
                    
}

SimulatorInternals::SimulatorInternals() {
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
    initalizeLimiting = false;
    requestForcedBypass = false;      // Request forced bypass in next NR iteration for all bypassable devices 
                                      // regardless of their converged state. 
}

void SimulatorInternals::fromOptions(const SimulatorOptions& options) {
    gmin = options.gmin;
    gshunt = options.gshunt;
}


// Introspection for options structure
template<> int Introspection<SimulatorOptions>::setup() {
    registerMember(temp);
    registerMember(tnom);
    
    registerMember(gmin);
    registerMember(gshunt);
    registerMember(minr);
    
    registerMember(scale);
    
    registerMember(reltol);
    registerMember(abstol);
    registerMember(vntol);
    registerMember(chgtol);
    registerMember(fluxtol);
    registerMember(relrefsol);
    registerMember(relrefres);
    registerMember(relreflte);
    registerMember(relref);
    registerMember(restol);
    registerMember(vnrestol);

    registerMember(matrixcheck);
    registerMember(rhscheck);
    registerMember(solutioncheck);
    
    registerMember(sweep_pointmarker);
    registerMember(sweep_debug);
    
    registerMember(nr_debug);
    registerMember(nr_bypass);
    registerMember(nr_convtol);
    registerMember(nr_bypasstol);
    registerMember(nr_conviter);
    registerMember(nr_residualcheck);
    registerMember(nr_damping);
    registerMember(nr_force);
    registerMember(nr_contbypass);

    registerMember(homotopy_debug);
    registerMember(homotopy_gminsteps);
    registerMember(homotopy_srcsteps);
    registerMember(homotopy_gminfactor);
    registerMember(homotopy_maxgminfactor);
    registerMember(homotopy_mingminfactor);
    registerMember(homotopy_startgmin);
    registerMember(homotopy_maxgmin);
    registerMember(homotopy_mingmin);
    registerMember(homotopy_srcstep);
    registerMember(homotopy_minsrcstep);
    registerMember(homotopy_sourcefactor);
    
    registerMember(op_debug);
    registerMember(op_itl);
    registerMember(op_itlcont);
    registerMember(op_skipinitial);
    registerMember(op_homotopy);
    registerMember(op_srchomotopy);
    registerMember(op_nsiter);

    registerMember(smsig_debug);

    registerMember(tran_debug);
    registerMember(tran_method);
    registerMember(tran_maxord);
    registerMember(tran_fs);
    registerMember(tran_ffmax);
    registerMember(tran_fbr);
    registerMember(tran_rmax);
    registerMember(tran_minpts);
    registerMember(tran_itl);
    registerMember(tran_ft);
    registerMember(tran_predictor);
    registerMember(tran_redofactor);
    registerMember(tran_lteratio);
    registerMember(tran_spicelte);
    registerMember(tran_xmu);
    registerMember(tran_trapltefilter);

    registerMember(hb_debug);
    registerMember(hb_itl);
    registerMember(hb_itlcont);
    registerMember(hb_skipinitial);
    registerMember(hb_homotopy);
    
    registerMember(rawfile);
    registerMember(strictoutput);
    registerMember(strictsave);
    registerMember(strictforce);
    registerMember(accounting);

    return 0;
}
instantiateIntrospection(SimulatorOptions);

std::unordered_map<Id, ParameterIndex> SimulatorOptions::mappingAffectingOptions;
std::unordered_map<Id, ParameterIndex> SimulatorOptions::hierarchyAffectingOptions;

bool SimulatorOptions::staticInitialize() {
    for(auto it : std::initializer_list<Id>{
        Id::createStatic("tnom"),
        Id::createStatic("temp"),
        Id::createStatic("scale"),
        Id::createStatic("minr"), 
    } ) {
        auto [ndx, found] = Introspection<SimulatorOptions>::index(it);
        mappingAffectingOptions.insert({it, static_cast<ParameterIndex>(ndx)});
    }

    for(auto it : std::initializer_list<Id>{}) {
        auto [ndx, found] = Introspection<SimulatorOptions>::index(it);
        hierarchyAffectingOptions.insert({it, static_cast<ParameterIndex>(ndx)});
    }

    return true;
}

static bool dummy = SimulatorOptions::staticInitialize();

bool SimulatorOptions::optionsDiffer(std::unordered_map<Id, ParameterIndex>& optionsList, SimulatorOptions& opt) {
    for(auto& it : optionsList) {
        if (!Introspection<SimulatorOptions>::memberEqual(*this, opt, it.second)) {
            return true;
        }
    }
    return false;
}


}
