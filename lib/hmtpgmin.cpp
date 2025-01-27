#include "hmtpgmin.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

GminStepping::GminStepping(Circuit& circuit, AnalysisCore& core, bool gdev)
    : Homotopy(circuit, core), gdev(gdev) {
}

std::tuple<bool, bool> GminStepping::run() {
    auto& options = circuit.simulatorOptions().core();
    auto debug = options.homotopy_debug;

    if (options.homotopy_gminsteps<1) {
        // Skip it
        return std::make_tuple(false, false);
    }

    // Allocate storage for analysis state
    auto stateNdx = core.allocateStateStorage(1);

    bool leave = false;
    
    // Stepping starts with zero initial solution and state
    double goodGmin;
    
    // Initially solve without previous solution (from all-zeros)
    bool continuation = false;
    
    // Store initial value of variable swept during homotopy
    double originalGmin;
    if (gdev) {
        originalGmin = circuit.simulatorInternals().gdev;
    } else {
        originalGmin = circuit.simulatorInternals().gshunt;
    }
    
    auto factor = options.homotopy_gminfactor;
    auto atGmin = options.homotopy_startgmin;
    auto maxGmin = options.homotopy_maxgmin;
    bool converged = false;

    // Target gmin value
    double targetGmin;
    if (gdev) {
        // Stepping gdev which is added to gmin option, step it down to gmin/10
        // We stop when effective gmin is 11/10 gmin option
        targetGmin = circuit.simulatorInternals().gmin / 10;
    } else {
        // Stepping gshunt, go to original gshunt value
        targetGmin = circuit.simulatorInternals().gshunt;
    }
    if (targetGmin==0.0) {
        targetGmin = circuit.simulatorOptions().core().homotopy_mingmin;
    }

    auto gminsteps = options.homotopy_gminsteps;
    itCount=0;
    while (true) {
        if (gdev) {
            circuit.simulatorInternals().gdev = atGmin;
        } else {
            circuit.simulatorInternals().gshunt = atGmin;
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
        if (leave) {
            converged = false;
            break;
        }
        
        if (converged) {
            // Converged
            continuation = true;
            if (atGmin <= targetGmin) {
                // Success
                break;
            }

            // Store point
            core.storeState(stateNdx, false);
            goodGmin = atGmin;
            
            // Solved with a few iterations
            if (core.iterations() <= core.iterationLimit(true)/4) {
                factor = factor * std::sqrt(factor);
                if (factor > options.homotopy_maxgminfactor) { 
                    factor = options.homotopy_maxgminfactor;
                }
            }

            // Solved with many iterations
            if (core.iterations() > core.iterationLimit(true)*3/4) {
                factor = std::sqrt(factor);
                if (factor < options.homotopy_mingminfactor) { 
                    factor = options.homotopy_mingminfactor;
                }
            }

            // Update gmin
            if (atGmin/factor < targetGmin) {
                factor = atGmin/targetGmin;
                atGmin = targetGmin;
            } else {
                atGmin = atGmin/factor;
            }
        } else {
            // Not converged
            if (!continuation) {
                // Not in continuation mode (no good solution yet)
                // Increase gmin
                atGmin = atGmin*factor;

                // Give up if gmin too large, append last analysis status
                if (atGmin>maxGmin) {
                    break;
                }
            } else {
                // Continuation mode, (have a good solution)
                // Decrease factor
                factor = std::pow(factor, 0.25);

                if (factor < options.homotopy_mingminfactor) { 
                    // Give up
                    break;
                }

                // Restore last successful point
                core.restoreState(stateNdx);
                atGmin = goodGmin;
            }
        }
        if (itCount>=gminsteps) {
            if (debug>0) {
                Simulator::dbg() << "Number of available gmin steps exhausted, stopping.\n";
            }
            break;
        }
    }

    // Restore original gmin
    if (gdev) {
        circuit.simulatorInternals().gdev = originalGmin;
    } else {
        circuit.simulatorInternals().gshunt = originalGmin;
    }

    // One final 
    if (!leave) {
        // Did not leave early
        if (converged) {
            // Converged, we can force bypass in first iteration
            circuit.simulatorInternals().requestForcedBypass = true;

            // Try with original gmin, but this time with continuation
            std::tie(converged, leave) = core.runSolver(continuation);
            itCount++;
            if (debug>0) {
                Simulator::dbg() << formatProgress() << ", final step " 
                    << (converged ? "converged in " : "failed to converge in ")
                    << core.iterations() << " solver iteration(s)"<< ".\n";
            }
        }
    } else {
        // Leaving early, did not converge
        // Add a status message at a higher level
        converged = false;
    }

    // Deallocate state storage
    core.deallocateStateStorage(1);

    return std::make_tuple(converged, leave);
}
    
std::string GminStepping::formatProgress() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(2);
    std::string txt = std::string("Homotopy: ") + (gdev ? "gdev=" : "gshunt=");
    if (gdev) {
        ss.str(""); ss << circuit.simulatorInternals().gmin + circuit.simulatorInternals().gdev;
    } else {
        ss.str(""); ss << circuit.simulatorInternals().gshunt;
    }
    txt += ss.str();
    return txt;
}


Spice3GminStepping::Spice3GminStepping(Circuit& circuit, AnalysisCore& core)
    : Homotopy(circuit, core) {
}

std::tuple<bool, bool> Spice3GminStepping::run() {
    auto& options = circuit.simulatorOptions().core();
    auto debug = options.homotopy_debug;
    bool converged = false;

    if (options.homotopy_gminsteps<1) {
        // Skip it
        return std::make_tuple(false, false);
    }

    // End value
    auto gshuntMin = circuit.simulatorInternals().gshunt;
    if (gshuntMin==0) {
        gshuntMin = options.homotopy_mingmin;
    }
    // Step
    auto logFactor = std::log(options.homotopy_gminfactor);
    
    // Number of steps
    auto gminsteps = options.homotopy_gminsteps;
    
    // Start value
    auto gshuntMax = gshuntMin * std::exp(logFactor*gminsteps);
    if (gshuntMax>options.homotopy_maxgmin) {
        gshuntMax = options.homotopy_maxgmin;
        logFactor = std::log(gshuntMax/gshuntMin)/gminsteps;
    }

    // First point without continuation
    bool continuation = false;

    // Store gshunt
    auto originalGshunt = circuit.simulatorInternals().gshunt;

    bool leave = false;

    itCount = 0;
    for(decltype(gminsteps) i=0; i<=gminsteps; i++) {
        // Set gshunt
        circuit.simulatorInternals().gshunt = gshuntMax / std::exp(logFactor*i);

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
            // Give up
            converged = false;
            break;
        }

        continuation = true;
    }

    // Restore gshunt
    circuit.simulatorInternals().gshunt = originalGshunt;

    // One final try with original gshunt, but this time with continuation
    if (!leave) {
        if (converged) {
            // Converged, we can force bypass in first iteration
            circuit.simulatorInternals().requestForcedBypass = true;

            // Try with original gmin, but this time with continuation
            std::tie(converged, leave) = core.runSolver(continuation);
            itCount++;
            if (debug>0) {
                Simulator::dbg() << formatProgress() << ", final step  " 
                    << (converged ? "converged in " : "failed to converge in ")
                    << core.iterations() << " solver iteration(s)" << ".\n";
            }
        } 
    } else {
        // Leaving early, did not converge
        // Add a status message at a higher level
        converged = false;
    }
    
    return std::make_tuple(converged, leave);
}

std::string Spice3GminStepping::formatProgress() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(2);
    std::string txt = std::string("Homotopy: SPICE3 gshunt=");
    ss.str(""); ss << circuit.simulatorInternals().gshunt;
    txt += ss.str();
    return txt;
}

}
