#include <cmath>
#include <iomanip>
#include "coreop.h"
#include "devbase.h"
#include "simulator.h"
#include "common.h"


namespace NAMESPACE {

// TODO: Bypass can cause problems with homotopy because small steps taken by homotopy
//       can be mistaken for instances converging and consequently being bypassed. 
//       This can lead to homotopy failure. 

// Each homotopy algorithm returns the error message of last failed nrSolver run. 
// The message reporting which homotopy algorithm failed is added one level higher

// Based on Ngspice dynamic_gmin() and new_gmin()
bool OperatingPointCore::gminStepping(RunType type) {
    auto debug = circuit.simulatorOptions().core().op_debug;
    bool leave = false;
    runType = runType | type;
    // Stepping starts with zero initial solution and state
    Vector<double> goodSolution;
    Vector<double> goodState;
    double goodGmin;
    
    // Initially solve without previous solution (from all-zeros)
    bool continuation = false;
    
    // Store initial value of variable swept during homotopy
    double originalGmin;
    if (type == RunType::GshuntStepping) {
        originalGmin = circuit.simulatorInternals().gshunt;
    } else {
        originalGmin = circuit.simulatorInternals().gdev;
    }
    
    auto& options = circuit.simulatorOptions().core();
    auto factor = options.op_gminfactor;
    auto atGmin = options.op_startgmin;
    auto maxGmin = options.op_maxgmin;
    bool converged = false;

    // Target gmin value
    double targetGmin;
    if (type == RunType::GshuntStepping) {
        // STepping gshunt, go to original gshunt value
        targetGmin = circuit.simulatorInternals().gshunt;
    } else {
        // Stepping gdev which is added to gmin option, step it down to gmin/10
        // We stop when effective gmin is 11/10 gmin option
        targetGmin = circuit.simulatorInternals().gmin / 10;
    }
    if (targetGmin==0.0) {
        targetGmin = circuit.simulatorOptions().core().op_mingmin;
    }

    auto gminsteps = options.op_gminsteps;
    decltype(gminsteps) itCount=0;
    while (true) {
        if (type == RunType::GshuntStepping) {
            circuit.simulatorInternals().gshunt = atGmin;
        } else {
            circuit.simulatorInternals().gdev = atGmin;
        }
        // If previous runSolver() converged we can force bypass in first NR iteration
        // If not, we do not force it. 
        circuit.simulatorInternals().forceBypass = converged && options.nr_contbypass;
        converged = runSolver(continuation);
        if (!converged) {
            setError(OpError::SteppingSolver);
            errorRunType = runType;
            errorHomotopyIterations = itCount;
        }
        leave = nrSolver.checkFlags(OpNRSolver::Flags::Abort);
        itCount++;
        if (debug>1) {
            Simulator::dbg() << homotopyProgress() << ", step " << itCount << " " 
                << (converged ? "converged" : "failed to converge") << ", " 
                << nrSolver.iterations() << " NR iteration(s)" << ".\n";
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
            goodSolution = solution.vector();
            goodState = states.vector();
            goodGmin = atGmin;
            
            // Solved with a few iterations
            if (nrSolver.iterations() <= options.op_itlcont/4) {
                factor = factor * std::sqrt(factor);
                if (factor > options.op_maxgminfactor) { 
                    factor = options.op_maxgminfactor;
                }
            }

            // Solved with many iterations
            if (nrSolver.iterations() > options.op_itlcont*3/4) {
                factor = std::sqrt(factor);
                if (factor < options.op_mingminfactor) { 
                    factor = options.op_mingminfactor;
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

                if (factor < options.op_mingminfactor) { 
                    // Give up
                    break;
                }

                // Restore last successful point
                solution.vector() = goodSolution;
                states.vector() = goodState;
                atGmin = goodGmin;
            }
        }
        if (itCount>=gminsteps) {
            if (debug>1) {
                Simulator::dbg() << "Number of available gmin steps exhausted, stopping.\n";
            }
            break;
        }
    }

    // Restore original gmin
    if (type == RunType::GshuntStepping) {
        circuit.simulatorInternals().gshunt = originalGmin;
    } else {
        circuit.simulatorInternals().gdev = originalGmin;
    }

    // One final 
    if (!leave) {
        // Did not leave early
        if (converged) {
            // Converged, we can force bypass in first NR iteration
            circuit.simulatorInternals().forceBypass = true;

            // Try with original gmin, but this time with continuation
            converged = runSolver(continuation);
            if (!converged) {
                setError(OpError::SteppingSolver);
                errorRunType = runType;
                errorHomotopyIterations = itCount;
            }
            leave = nrSolver.checkFlags(OpNRSolver::Flags::Abort);
            itCount++;
            if (debug>1) {
                Simulator::dbg() << homotopyProgress() << ", final step " 
                    << (converged ? "converged" : "failed to converge") << ", " 
                    << nrSolver.iterations() << " NR iteration(s)"<< ".\n";
            }
        } else {
            // Failed 
            setError(OpError::SteppingSteps);
            errorRunType = type;
            errorHomotopyIterations = itCount;
        }
    } else {
        // Leaving early, did not converge
        // Add a status message at a higher level
        converged = false;
    }

    runType = runType & ~type;
    return converged;
}

// Based on Ngspice spice3_gmin()
bool OperatingPointCore::spice3GminStepping() {
    auto& options = circuit.simulatorOptions().core();
    auto debug = options.op_debug;
    runType = runType | RunType::Spice3GminStepping;
    bool converged = false;

    // End value
    auto gshuntMin = circuit.simulatorInternals().gshunt;
    if (gshuntMin==0) {
        gshuntMin = options.op_mingmin;
    }
    // Step
    auto logFactor = std::log(options.op_gminfactor);
    
    // Number of steps
    auto gminsteps = options.op_gminsteps;
    
    // Start value
    auto gshuntMax = gshuntMin * std::exp(logFactor*gminsteps);
    if (gshuntMax>options.op_maxgmin) {
        gshuntMax = options.op_maxgmin;
        logFactor = std::log(gshuntMax/gshuntMin)/gminsteps;
    }

    // First point without continuation
    bool continuation = false;

    // Store gshunt
    auto originalGshunt = circuit.simulatorInternals().gshunt;

    bool leave = false;

    decltype(gminsteps) itCount = 0;
    for(decltype(gminsteps) i=0; i<=gminsteps; i++) {
        // Set gshunt
        circuit.simulatorInternals().gshunt = gshuntMax / std::exp(logFactor*i);

        // If previous runSolver() converged we can force bypass in first NR iteration
        // If not, we do not force it. 
        circuit.simulatorInternals().forceBypass = converged && options.nr_contbypass;
        converged = runSolver(continuation);
        if (!converged) {
            setError(OpError::SteppingSolver);
            errorRunType = runType;
            errorHomotopyIterations = itCount;
        }
        leave = nrSolver.checkFlags(OpNRSolver::Flags::Abort);
        itCount++;
        if (debug>1) {
            Simulator::dbg() << homotopyProgress() << ", step " << itCount << " " 
                << (converged ? "converged" : "failed to converge") << ", " 
                 << nrSolver.iterations() << " NR iteration(s)" << ".\n";
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
            // Converged, we can force bypass in first NR iteration
            circuit.simulatorInternals().forceBypass = true;

            // Try with original gmin, but this time with continuation
            converged = runSolver(continuation);
            if (!converged) {
                setError(OpError::SteppingSolver);
                errorRunType = runType;
                errorHomotopyIterations = itCount;
            }
            leave = nrSolver.checkFlags(OpNRSolver::Flags::Abort);
            itCount++;
            if (debug>1) {
                Simulator::dbg() << homotopyProgress() << ", final step  " 
                    << (converged ? "converged" : "failed to converge")  << ", " 
                    << nrSolver.iterations() << " NR iteration(s)" << ".\n";
            }
        } else {
            // Failed
            setError(OpError::SteppingSteps);
            errorRunType = RunType::Spice3GminStepping;
            errorHomotopyIterations = itCount;
        }
    } else {
        // Leaving early, did not converge
        // Add a status message at a higher level
        converged = false;
    }
    
    runType = runType & ~RunType::Spice3GminStepping;
    return converged;
}

// Based on Ngspice gillespie_src()
bool OperatingPointCore::sourceStepping() {
    auto& options = circuit.simulatorOptions().core();
    auto debug = options.op_debug;
    runType = runType | RunType::SourceStepping;
    double goodFactor = 0.0;
    bool continuation = false;
    bool leave = false;
    double raise = options.op_srcstep;
    bool converged = false;
    
    // Set source factor to 0
    circuit.simulatorInternals().sourcescalefactor = goodFactor;

    auto srcsteps = options.op_srcsteps;
    decltype(srcsteps) itCount=0;
    
    // Try solving without gmin stepping at source scale factor = 0.0
    converged = runSolver(continuation);
    if (!converged) {
        setError(OpError::SteppingSolver);
        errorRunType = runType;
        errorHomotopyIterations = itCount;
    }
    leave = nrSolver.checkFlags(OpNRSolver::Flags::Abort);
    if (debug>1) {
        Simulator::dbg() << homotopyProgress() << ", initial OP " 
            << (converged ? "converged" : "failed to converge") << ", " 
            << nrSolver.iterations() << " NR iteration(s)" << ".\n";
    }

    // Do gmin stepping even if op_skipgmin is 1
    if (!leave && !converged && options.op_gminsteps>1) {
        // Failed, now try with gmin stepping
        if (!options.op_spice3gmin) {
            // New algorithms
            // Device gmin first
            converged = gminStepping(RunType::GminStepping);
            leave = nrSolver.checkFlags(OpNRSolver::Flags::Abort);
            if (debug>0) {
                Simulator::dbg() << homotopyProgress() << ", initial gmin stepping " << (converged ? "converged" : "failed to converge") << ".\n";
            }
            if (!converged && !leave && circuit.simulatorOptions().core().op_gshuntalg) {
                // Diagonal gshunt second
                converged = gminStepping(RunType::GshuntStepping);
                leave = nrSolver.checkFlags(OpNRSolver::Flags::Abort);
                if (debug>0) {
                    Simulator::dbg() << homotopyProgress() << ", initial gshunt stepping " << (converged ? "converged" : "failed to converge") << ".\n";
                }
            } 
        } else {
            // Spice3 gmin stepping
            converged = spice3GminStepping();
            leave = nrSolver.checkFlags(OpNRSolver::Flags::Abort);
            if (debug>0) {
                Simulator::dbg() << homotopyProgress() << ", initial spice3 gmin stepping " << (converged ? "converged" : "failed to converge") << ".\n";
            }
        }
    } 

    if (!leave && converged) {
        // Now try source stepping
        
        // Set continuation mode
        continuation = true;

        // Store old solution
        Vector<double> goodSolution;
        Vector<double> goodState;

        goodSolution = solution.vector();
        goodState = states.vector();
        
        // Source stepping loop
        while (true) {
            // New source factor, limit to 1.0
            auto srcfactor = goodFactor+raise;
            if (srcfactor>1.0) {
                srcfactor = 1.0;
            }
            circuit.simulatorInternals().sourcescalefactor = srcfactor*options.op_sourcefactor;

            // If previous runSolver() converged we can force bypass in first NR iteration
            // If not, we do not force it. 
            circuit.simulatorInternals().forceBypass = converged && options.nr_contbypass;
            converged = runSolver(continuation);
            if (!converged) {
                setError(OpError::SteppingSolver);
                errorRunType = runType;
                errorHomotopyIterations = itCount;
            }
            leave = nrSolver.checkFlags(OpNRSolver::Flags::Abort);
            itCount++;
            if (debug>1) {
                Simulator::dbg() << homotopyProgress() << ", step " << itCount << " " 
                    << (converged ? "converged" : "failed to converge") << ", " 
                    << nrSolver.iterations() << " NR iteration(s)" << ".\n";
            }

            if (leave) {
                converged = false;
                break;
            }

            if (converged) {
                // Converged
                // Store point
                goodFactor = srcfactor;
                goodSolution = solution.vector();
                goodState = states.vector();

                if (goodFactor >= 1.0) {
                    // Success
                    break;
                }
                
                if (nrSolver.iterations() <= options.op_itlcont/4) {
                    raise = raise*1.5;
                }

                if (nrSolver.iterations() > options.op_itlcont*3/4) {
                    raise = raise*0.5;
                    if (raise < options.op_minsrcstep) {
                        raise = options.op_minsrcstep;
                    }
                }
            } else {
                // Not converged
                raise = raise*0.5;

                // Restore old solution
                solution.vector() = goodSolution;
                states.vector() = goodState;

                if (raise<options.op_minsrcstep) {
                    // Give up
                    break;
                }
            }

            if (itCount>=srcsteps) {
                if (debug>1) {
                    Simulator::dbg() << "Number of available source steps exhausted, stopping.\n";
                }
                break;
            }
        }
    }

    converged = goodFactor>=1.0;

    if (!leave) {
        if (itCount>=srcsteps) {
            // Failed
            setError(OpError::SteppingSteps);
            errorRunType = RunType::SourceStepping;
            errorHomotopyIterations = itCount;
        }
    } else {
        // Leaving early, did not converge
        // Add a status message at a higher level
        converged = false;
    }

    // Restore source scale factor (in case we exited the loop early)
    circuit.simulatorInternals().sourcescalefactor = 1.0*options.op_sourcefactor;

    runType = runType & ~RunType::SourceStepping;
    return converged;
}

// Based on Ngspice spice3_src()
bool OperatingPointCore::spice3SourceStepping() {
    auto& options = circuit.simulatorOptions().core();
    auto debug = options.op_debug;
    runType = runType | RunType::Spice3SourceStepping;
    bool converged = false;

    auto srcsteps = options.op_srcsteps;

    // First point without continuation
    bool continuation = false;

    bool leave = false;

    decltype(srcsteps) itCount = 0;
    decltype(srcsteps) i;
    for(i=0; i<=srcsteps; i++) {
        // Set source scale factor
        if (i<srcsteps) {
            circuit.simulatorInternals().sourcescalefactor = (double)i/srcsteps*options.op_sourcefactor;
        } else {
            // In last step restore source scale factor to 1.0 to avoid roundoff errors
            circuit.simulatorInternals().sourcescalefactor = 1.0*options.op_sourcefactor; 
        }

        // If previous runSolver() converged we can force bypass in first NR iteration
        // If not, we do not force it. 
        circuit.simulatorInternals().forceBypass = converged && options.nr_contbypass;
        converged = runSolver(continuation);
        if (!converged) {
            setError(OpError::SteppingSolver);
            errorRunType = runType;
            errorHomotopyIterations = itCount;
        }
        leave = nrSolver.checkFlags(OpNRSolver::Flags::Abort);
        itCount++;
        if (debug>1) {
            Simulator::dbg() << homotopyProgress() << ", step " << itCount << " " 
                << (converged ? "converged" : "failed to converge") << ", " 
                << nrSolver.iterations() << " NR iteration(s)" << ".\n";
        }
        if (leave || !converged) {
            // Leave
            converged = false;
            break;
        }

        continuation = true;
    }

    if (!leave) {
        if (itCount>srcsteps) {
            // Failed
            setError(OpError::SteppingSteps);
            errorRunType = RunType::Spice3SourceStepping;
            errorHomotopyIterations = itCount;
}
    } else {
        // Leaving early, did not converge
        // Add a status message at a higher level
        converged = false;
    }

    // Restore source scale factor (in case we exited the loop early)
    circuit.simulatorInternals().sourcescalefactor = 1.0;

    runType = runType & ~RunType::Spice3SourceStepping;
    return converged;
}

std::string OperatingPointCore::homotopyProgress() const {
    // Formatting is a pain, after gcc libc++ supports format from c++20 this will be much simpler
    std::stringstream ss;
    ss << std::scientific << std::setprecision(2);
    std::string txt = "Homotopy, ";
    if (static_cast<bool>(runType & RunType::SourceStepping)) {
        ss.str(""); ss << circuit.simulatorInternals().sourcescalefactor;
        txt += "srcfact="+ss.str();
        if (static_cast<bool>(runType & RunType::GSteppingMask)) {
            txt += ", ";
        }
    } else if  (static_cast<bool>(runType & RunType::Spice3SourceStepping)) {
        ss.str(""); ss << circuit.simulatorInternals().sourcescalefactor;
        txt += "spice3 srcfact="+ss.str();
    }
    switch (runType & RunType::GSteppingMask) {
        case RunType::GminStepping:
            ss.str(""); ss << circuit.simulatorInternals().gmin + circuit.simulatorInternals().gdev;
            txt += "gmin="+ss.str();
            break;
        case RunType::GshuntStepping:
            ss.str(""); ss << circuit.simulatorInternals().gshunt;
            txt += "gshunt="+ss.str();
            break;
        case RunType::Spice3GminStepping:
            ss.str(""); ss << circuit.simulatorInternals().gshunt;
            txt += "spice3 gmin="+ss.str();
            break;
    }
    return txt;
}

}
