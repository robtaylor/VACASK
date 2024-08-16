#include "options.h"
#include "introspection.h"
#include "simulator.h"
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
    
    sweep_debug = 0; // 1 = debug sweep, >=2 print details

    op_debug = 0; // 0 = none, 1 = NRSolver and homotopy runs, 2 = homotopy and continuation internals, 
                  // 100  = nrdebug=1, print NRSolver internals (progress)
                  // 101  = nrdebug=2, print linear system
                  // 102  = nrdebug=3, print solutions
                  // >102 = nrdebug>=4, print old solution before building system

    nr_bypass = 0; // 1 = disable core evaluation for bypassed instances 
                   // To be bypassed an instance must converge and 
                   // the instance inputs change between iterations must be within tolerances. 
                   // As soon as the inputs change is outside tolerances instance is no longer bypassed. 
    nr_convtol = 1.0;   // Tolerance factor applied to residuals for instance convergence check. Should be <1.
                        // 1.0 corresponds to set tolerances (abstol, vntol, ...), <1.0 makes them more strict. 
    nr_bypasstol = 1.0; // Tolerance factor applied to instance inputs for instance bypass check. Should be <1.
                        // 1.0 corresponds to set tolerances (abstol, vntol, ...), <1.0 makes them more strict. 
    nr_conviter = 1; // >0, number of consecutive convergent iterations before convergence is confirmed
    nr_residualcheck = 1; // check residual beside unknowns change to establish convergence 
    nr_damping = 1.0; // 0<x<=1, Newton-Raphson damping factor (<=1)
    nr_force = 1e5; // x>0, forcing factor for nodesets and initial conditions
    
    op_itl = 100;  // >0, maximal number of iterations in non-continuation mode
    op_itlcont = 50; // >0, maximal number of iterations in continuation mode

    op_skipinitial = 0; // 1 = no initial op, go straight to gmin stepping
    op_skipgmin = 0; // 1 = skip gmin stepping and go straight to source stepping
    op_skipsrc = 0;  // 1 = skip source stepping
    op_skiphomotopy = 0; // 1 = skip homotopy methods (gmin and source stepping)
                         // at least one of the methods must be enabled (initial, gmin stepping, source stepping)
    op_spice3gmin = 0; // 1 = use spice3-style gmin stepping
    op_spice3src = 0; // 1 = use spice3-style source stepping
    op_sourcefactor = 1.0; // x = op_sourcefactor * sourcescalefactor where sourcescalefactor is set by 
                           // the source stepping homotopy algorithm is the scaling factor for all independent sources. 
                           // The product (x) is available in the SimulatorInternals data structure as sourcescalefactor 
                           // and should be used for computing the resistive residual of all independent sources. 
                           // This option makes it possible to set up manual homotopy via dc sweep. 
                           // For normal simulation this should be 1.0. 
    op_gshuntalg = 1; // 1= in non-spice3 mode do device gmin stepping followed by gshunt stepping, 
                   // 0= device gmin stepping only
    op_gminfactor = 10.0; // >0, initial gmin stepping factor for dynamic gmin stepping
    op_maxgminfactor = 10.0; // >=op_gminfactor, maximal gmin stepping factor for dynamic gmin stepping
    op_mingminfactor = 1.00005; // 1<x<op_gminfactor, maximal gmin stepping factor for dynamic gmin stepping
                                // give up when gmin step fails and factor falls below this value
    op_startgmin = 1e-3; // >mingmin, value at which dynamic gmin stepping starts
    op_maxgmin = 1e2;   // >mingmin, if op dynamic gmin stepping failes to solve the circuit above this value of gmin, it fails
    op_mingmin = 1e-15; // >0, value where dynamic gmin stepping stops if gmin/gshunt are set to 0
    op_gminsteps = 100; // >1, <=0 is equivalent to op_skipgmin=1
    
    op_srcstep = 0.001; // >0, initial source step for dynamic source stepping
    op_minsrcstep = 1e-7; // 0<x<srcstep, source step at which dynamic source stepping gives up
    op_gminsteps = 100; // >1, <=0 is equivalent to op_skipgmin=1
    op_srcsteps = 100; // >1, <=0 is equvalent to op_skipsrc=1

    op_nsiter = 1; // Number of iterations during which nodesets are applied


    smsig_debug = 0; // 1=debug small signal analyses (dcinc, dctf, ac, actf, noise), 
                     // >=100 print linear system
    
    tran_debug = 0; // 1 = debug steps, 2 = debug solver
                    // 100 = print linear system, 
                    // >100 print solutions
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
    sourcescalefactor = 1.0;
    gmin = 0.0;
    gdev = 0.0;
    gshunt = 0.0;
    iteration = 0;
    analysis_name = "";
    analysis_type = "";
    cwd = Simulator::startupPath();
    initalizeLimiting = false;
    highPrecision = false; // request high precision from simulator (prevents bypass for all bypassable devices)
    forceBypass = false;   // force bypass in next NR iteration for all bypassable devices regardless of their 
                           // converged state, has lower precedence than highPrecision
    frequency = 0.0;
    time = 0.0;
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
    
    registerMember(sweep_debug);
    
    registerMember(nr_bypass);
    registerMember(nr_convtol);
    registerMember(nr_bypasstol);
    registerMember(nr_conviter);
    registerMember(nr_residualcheck);
    registerMember(nr_damping);
    registerMember(nr_force);
    
    registerMember(op_debug);
    registerMember(op_itl);
    registerMember(op_itlcont);
    registerMember(op_skipinitial);
    registerMember(op_skipgmin);
    registerMember(op_skipsrc);
    registerMember(op_skiphomotopy);
    registerMember(op_spice3gmin);
    registerMember(op_spice3src);
    registerMember(op_sourcefactor);
    registerMember(op_gshuntalg);
    registerMember(op_gminfactor);
    registerMember(op_maxgminfactor);
    registerMember(op_mingminfactor);
    registerMember(op_startgmin);
    registerMember(op_maxgmin);
    registerMember(op_mingmin);
    registerMember(op_srcstep);
    registerMember(op_minsrcstep);
    registerMember(op_gminsteps);
    registerMember(op_srcsteps);
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
