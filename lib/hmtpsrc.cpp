#include "hmtpsrc.h"
#include "hmtpgmin.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

SourceStepping::SourceStepping(Circuit& circuit, AnalysisCore& core, const std::vector<Id>& alternativeHomotopy)
    : Homotopy(circuit, core), alternativeHomotopy(alternativeHomotopy) {
}

static Homotopy* createHomotopy(Circuit& circuit, AnalysisCore& core, Homotopy& parent, Int debug, Id name) {
    if (name==Homotopy::src || name==Homotopy::spice3Src) {
        if (debug>0) {
            Simulator::dbg() << parent.formatProgress() << ", source stepping cannot be used.\n";
        }
        return nullptr;
    } else if (name==Homotopy::gdev) {
        if (debug>0) {
            Simulator::dbg() << parent.formatProgress() << ", trying gdev stepping.\n";
        }
        return new GminStepping(circuit, core, true);
    } else if (name==Homotopy::gshunt) {
        if (debug>0) {
            Simulator::dbg() << parent.formatProgress() << ", trying gshunt stepping.\n";
        }
        return new GminStepping(circuit, core, false);
    } else if (name==Homotopy::spice3Gmin) {
        if (debug>0) {
            Simulator::dbg() << parent.formatProgress() << ", trying SPICE3 gmin stepping.\n";
        }
        return new Spice3GminStepping(circuit, core);
    } else {
        if (debug>0) {
            Simulator::dbg() << parent.formatProgress() << ", unknown homotopy '"+std::string(name)+"'.\n";
        }
        return nullptr;
    }
}

std::tuple<bool, bool> SourceStepping::run() {
    auto& options = circuit.simulatorOptions().core();
    auto debug = options.homotopy_debug;
    double goodFactor = 0.0;
    bool continuation = false;
    bool leave = false;
    double raise = options.homotopy_srcstep;
    bool converged = false;
    auto srcsteps = options.homotopy_srcsteps;
    
    if (srcsteps<1) {
        // Skip it
        return std::make_tuple(false, false);
    }

    // Allocate storage for analysis state
    auto stateNdx = core.allocateStateStorage(1);

    // Set source factor to 0
    circuit.simulatorInternals().sourcescalefactor = goodFactor;
    itCount=0;
    
    // Try solving without gmin stepping at source scale factor = 0.0
    std::tie(converged, leave) = core.runSolver(continuation);
    if (debug>0) {
        Simulator::dbg() << formatProgress() << ", initial OP " 
            << (converged ? "converged in " : "failed to converge in ")
            << core.iterations() << " solver iteration(s)" << ".\n";
    }

    if (!leave && !converged) {
        // Failed at sourcefactor=0, now go through all alternative homotopies
        Homotopy* homotopy;
        for(auto it : alternativeHomotopy) {
            homotopy = createHomotopy(circuit, core, *this, debug, it);
            if (!homotopy) {
                continue;
            }
            std::tie(converged, leave) = homotopy->run();
            if (debug>0) {
                Simulator::dbg() << formatProgress() << ", " << (converged ? "converged" : "failed to converge") << ".\n";
            }
            delete homotopy;
            if (leave || converged) {
                break;
            }
        }
    } 

    if (!leave && converged) {
        // Now try source stepping
        
        // Set continuation mode
        continuation = true;

        // Store old solution
        core.storeState(stateNdx, false);
        
        // Source stepping loop
        while (true) {
            // New source factor, limit to 1.0
            auto srcfactor = goodFactor+raise;
            if (srcfactor>1.0) {
                srcfactor = 1.0;
            }
            circuit.simulatorInternals().sourcescalefactor = srcfactor;

            // If previous runSolver() converged we can force bypass in first iteration
            // If not, we do not force it. 
            circuit.simulatorInternals().requestForcedBypass = converged && options.nr_contbypass;
            std::tie(converged, leave) = core.runSolver(continuation);
            itCount++;
            if (debug>0) {
                Simulator::dbg() << formatProgress() << ", step " << itCount << " " 
                    << (converged ? "converged in " : "failed to converge in ")
                    << core.iterations() << " solver iteration(s)" << ".\n";
            }

            if (leave) {
                converged = false;
                break;
            }

            if (converged) {
                // Converged
                // Store point
                goodFactor = srcfactor;
                core.storeState(stateNdx, false);
                
                if (goodFactor >= 1.0) {
                    // Success
                    break;
                }
                
                if (core.iterations() <= core.iterationLimit(true)/4) {
                    raise = raise*options.homotopy_srcscale;
                }

                if (core.iterations() > core.iterationLimit(true)*3/4) {
                    raise = raise/options.homotopy_srcscale;
                    if (raise < options.homotopy_minsrcstep) {
                        raise = options.homotopy_minsrcstep;
                    }
                }
            } else {
                // Not converged
                raise = raise*0.5;

                // Restore old solution
                core.restoreState(stateNdx);
                
                if (raise<options.homotopy_minsrcstep) {
                    // Give up
                    break;
                }
            }

            if (itCount>=srcsteps) {
                if (debug>0) {
                    Simulator::dbg() << "Number of available source steps exhausted, stopping.\n";
                }
                break;
            }
        }
    }

    converged = goodFactor>=1.0;

    if (leave) {
        // Leaving early, did not converge
        // Add a status message at a higher level
        converged = false;
    }

    // Restore source scale factor (in case we exited the loop early)
    circuit.simulatorInternals().sourcescalefactor = 1.0;

    // Deallocate state storage
    core.deallocateStateStorage(stateNdx);

    return std::make_tuple(converged, leave);
}

std::string SourceStepping::formatProgress() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(2);
    std::string txt = std::string("Homotopy: srcfact=");
    ss.str(""); ss << circuit.simulatorInternals().sourcescalefactor;
    txt += ss.str();
    return txt;
}


Spice3SourceStepping::Spice3SourceStepping(Circuit& circuit, AnalysisCore& core)
    : Homotopy(circuit, core) {
}

std::tuple<bool, bool> Spice3SourceStepping::run() {
    auto& options = circuit.simulatorOptions().core();
    auto debug = options.homotopy_debug;
    bool converged = false;

    auto srcsteps = options.homotopy_srcsteps;

    if (srcsteps<1) {
        // Skip it
        return std::make_tuple(false, false);
    }

    // First point without continuation
    bool continuation = false;

    bool leave = false;

    itCount = 0;
    decltype(srcsteps) i;
    for(i=0; i<=srcsteps; i++) {
        // Set source scale factor
        if (i<srcsteps) {
            circuit.simulatorInternals().sourcescalefactor = (double)i/srcsteps;
        } else {
            // In last step restore source scale factor to 1.0 to avoid roundoff errors
            circuit.simulatorInternals().sourcescalefactor = 1.0; 
        }

        // If previous runSolver() converged we can force bypass in first iteration
        // If not, we do not force it. 
        circuit.simulatorInternals().requestForcedBypass = converged && options.nr_contbypass;
        std::tie(converged, leave) = core.runSolver(continuation);
        itCount++;
        if (debug>0) {
            Simulator::dbg() << formatProgress() << ", step " << itCount << " " 
                << (converged ? "converged in " : "failed to converge in ") 
                << core.iterations() << " solver iteration(s)" << ".\n";
        }
        if (leave || !converged) {
            // Leave
            converged = false;
            break;
        }

        continuation = true;
    }

    if (leave) {
        // Leaving early, did not converge
        // Add a status message at a higher level
        converged = false;
    }

    // Restore source scale factor (in case we exited the loop early)
    circuit.simulatorInternals().sourcescalefactor = 1.0;

    return std::make_tuple(converged, leave);
}

std::string Spice3SourceStepping::formatProgress() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(2);
    std::string txt = std::string("Homotopy: SPICE3 srcfact=");
    ss.str(""); ss << circuit.simulatorInternals().sourcescalefactor;
    txt += ss.str();
    return txt;
}

}

