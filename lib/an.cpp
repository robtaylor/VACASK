#include "an.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

Analysis::Analysis(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : name_(name), circuit(circuit), sweeper(nullptr), 
      ptAnalysis(ptAnalysis), resolver(circuit), commonSaves(nullptr) {
}

Analysis::~Analysis() {
}

void Analysis::setSaves(PTSavesVector* saves) {
    // Store common saves
    commonSaves = saves;
}

void Analysis::setParametrization(const PTParameterMap* optionsMap) {
    // Store only options that are defined with expressions
    parameterizedOptions.clear();
    for(auto& it : *optionsMap) {
        if (std::holds_alternative<const PTParameterExpression*>(it.second)) {
            parameterizedOptions.insert(it);
        }
    }
}

bool Analysis::registerFactory(Id name, Analysis::AnalysisFactory factory) {
    auto [it, inserted] = getRegistry().insert(std::make_pair(name, factory));
    return inserted;
}

Analysis* Analysis::create(
    PTAnalysis& ptAnalysis, 
    PTSavesVector* commonSaves, 
    const PTParameterMap* optionsMap, 
    Circuit& circuit, Status& s
) {
    auto it = getRegistry().find(ptAnalysis.typeName());
    if (it==getRegistry().end()) {
        s.set(Status::NotFound, "Analysis type '"+std::string(ptAnalysis.typeName())+"' not found.");
        s.extend(ptAnalysis.location());
        return nullptr;
    }

    // Create analysis
    auto factory = it->second;
    auto* an = factory(ptAnalysis, circuit, s);
    if (!an) {
        return nullptr;
    }

    // Store common saves
    an->setSaves(commonSaves);

    // Set up parametrization
    an->setParametrization(optionsMap);

    // Set analysis parameters based on global circuit context
    auto [ok, changed] = an->parameters().setParameters(ptAnalysis.parameters(), circuit.variableEvaluator(), s);
    if (!ok) {
        delete an;
        return nullptr;
    }

    // Store simulator options state
    an->originalSimOptions.core() = circuit.simulatorOptions().core();

    // Set up initial simulator options for analysis
    an->simOptions.core() = an->originalSimOptions.core();
    
    // Add sweeps
    // Sweep settings can depend on global circuit parameters
    // But if during sweep a global parameter changed and the sweep settings depend on it
    // they will not be affected by the change. Basically sweep parameters are taken 
    // as they are before the sweep is initiated. 
    // To perform an analysis with sweep settings changed due to changed global parameters 
    // an analysis has to be recreated. 
    auto& anSweeps = ptAnalysis.sweeps();
    for(auto it=anSweeps.data().cbegin(); it!=anSweeps.data().cend(); ++it) {
        IStruct<SweepSettings> sw;
        sw.core().name = it->name();
        sw.core().location = it->location();
        auto [ok, changed] = sw.setParameters(it->parameters(), circuit.variableEvaluator(), s);
        if (!ok) {
            s.extend("Error in sweep setup for analysis '"+std::string(ptAnalysis.name())+"'.");
            s.extend(it->location());
            delete an;
            return nullptr;
        }
        an->addSweep(std::move(sw.core()));
    }

    // Build analysis state repository
    an->resizeAnalysisStateStorage(anSweeps.data().size());

    return an;
}

bool Analysis::addOutputDescriptors(Status& s) {
    // Clear old descriptors in all cores
    clearOutputDescriptors();

    // Add sweep variables
    auto nSweeps = sweeps_.size();
    for(decltype(nSweeps) i=0; i<nSweeps; i++) {
        addCommonOutputDescriptor(OutputDescriptor(OutdSweepvar, sweeper->sweepName(i), i)); 
    }

    // Add core output descriptors
    addCoreOutputDescriptors();

    // Add saves
    bool strict;
    
    // Add common saves from circuit, unknown saves are ignored, known ones are checked for syntax
    strict = false;
    if (commonSaves) {
        for(auto saves : *commonSaves) {
            for (auto it = saves->saves().cbegin(); it != saves->saves().cend(); ++it) {
                // Ignore failures, unless strict is true
                Status saveSt;
                auto ok = resolveSave(*it, strict, saveSt);
                // Semantic error in save directive results in ok=false, it is always an error
                // Failure to add a descriptor in strict mode results in ok=false and thus an error
                // Unsupported descriptor type in strict mode results in ok=false and thus error
                // Unsupported save directive is not an error unless in strict mode
                if (strict && !ok) {
                    s.set(saveSt);
                    return false;
                }
            }
        }
    }

    // Add default saves if needed
    addDefaultOutputDescriptors(s);

    return true;
}

bool Analysis::resolveOutputDescriptor(const OutputDescriptor& descr, Output::SourcesList& srcs, bool strict, Status& s) {
    // Abstract analysis handles only sweep variables
    switch (descr.type) {
        case OutdSweepvar: 
            srcs.emplace_back(sweeper.get(), descr.ndx);
            break;
        default:
            s.set(Status::Internal, "Unknown output descriptor type.");
            return false;
    }
    return true;
}

bool Analysis::addSweep(const SweepSettings& sw, Status& s) {
    sweeps_.push_back(sw);
    return true;
}

bool Analysis::addSweep(SweepSettings&& sw, Status& s) {
    sweeps_.push_back(std::move(sw));
    return true;
}

bool Analysis::run(Status& s) {
    // Output descriptors are created with analysis
    // Binding must be done here, just before the core analysis is run

    // Set simulator internals
    auto& options = circuit.simulatorOptions().core(); 
    SimulatorInternals internals;
    internals.fromOptions(options);
    internals.analysis_name = std::string(name_);
    internals.analysis_type = std::string(ptAnalysis.typeName());
    internals.initalizeLimiting = false;
    circuit.simulatorInternals() = internals;

    // Clear Abort, Finish, and Stop flags
    // Responders 
    //   Abort (is an error)
    //     innermost (NR) loop of analysis
    //     homotopy loops
    //     frequency/time sweep loop
    //     sweep loop
    //   Finish, Stop
    //     frequency/time sweep loop
    //     sweep loop
    circuit.clearFlags(Circuit::Flags::Abort);
    circuit.clearFlags(Circuit::Flags::Finish);
    circuit.clearFlags(Circuit::Flags::Stop);

    auto sweepDebug = circuit.simulatorOptions().core().sweep_debug;

    bool runOk = false;

    // Do we have sweep(s)
    if (sweeps_.size()>0) {
        // Sweep required
        if (sweepDebug>0) {
            Simulator::dbg() << "Starting sweep, analysis '" << std::string(name_) << "'.\n";
        }

        // Create sweeper
        sweeper = std::make_unique<ParameterSweeper>(sweeps_, s);
        if (!sweeper->isValid()) {
            s.extend("Failed to initialize sweep.");
            return false;
        }

        // Bind sweeper to actual parameters and options
        if (!sweeper->bind(circuit, simOptions, s)) {
            s.extend("Failed to bind sweep parameters.");
            return false;
        }

        // Collect current values so we can restore them later
        if (!sweeper->storeState(s)) {
            s.extend("Failed to store initial circuit state.");
            return false;
        }

        // Reset sweeper
        sweeper->reset();

        // Variable indicating that the outputs are bound
        bool outputsBound = false;

        // Varible indicating that the output is initializes
        bool outputInitialized = false;

        // Main sweep loop
        bool haveStoredState = false;
        // Initially the outermost sweep index increases
        Int advancedSweepIndex = 0;
        do {
            // Assume analysis failed
            runOk = false;

            if (sweepDebug>0) {
                Simulator::dbg() << "Sweep point: " << sweeper->progress() << ".\n";
            }

            // Set current sweep point
            auto [ok, hierarchyChanged, needsCoreRebuild] = circuit.elaborateChanges(
                sweeper.get(), ParameterSweeper::WriteValues::Sweep, 
                this, &simOptions, 
                &parameterizedOptions, 
                s
            );
            if (!ok) {
                s.extend("Failed to set circuit state.");
                runOk = false;
                break;
            }

            // Delete output files if strictoutput>=2 and this is the first sweep point
            auto strictoutput = circuit.simulatorOptions().core().strictoutput;
            if (strictoutput>=2 && !outputInitialized) {
                deleteOutputs();
            }

            // Need to re-bind sweeper if hierarchy changed
            if (hierarchyChanged) {
                if (!sweeper->bind(circuit, simOptions, s)) {
                    s.extend("Failed to re-bind sweep parameters.");
                    runOk = false;
                    break;
                }
            }

            // If outputs not bound yet or core needs rebuilding, bind them to actual quantities
            bool coreRebuilt = false;
            if (!outputsBound || needsCoreRebuild) {
                if (sweepDebug>1) {
                    Simulator::dbg() << "Rebuilding analysis internals and invalidating stored analysis states.\n";
                }

                // Only the first time (while outputsBound=false)
                if (!outputsBound) {
                    if (!addOutputDescriptors(s)) {
                        s.extend("Failed to construct list of outputs.");
                        runOk = false;
                        break;
                    }
                }

                // Every time core needs rebuild, rebind outputs
                bool strict = (!outputsBound && circuit.simulatorOptions().core().strictsave>0) ||
                              (outputsBound && circuit.simulatorOptions().core().strictsave>1);
                if (!resolveOutputDescriptors(strict, s)) {
                    s.extend("Failed to bind analysis outputs.");
                    runOk = false;
                    break;
                }

                // Rebuild core
                if (!rebuildCores(s)) {
                    s.extend("Failed to rebuild analysis structures.");
                    runOk = false;
                    break;
                }
                outputsBound = true;
                coreRebuilt = true;

                // Core rebuild makes all stored states incoherent with current topology
                for(decltype(sweepCount()) i=0; i<sweepCount(); i++) {
                    makeStateIncoherent(i);
                }
            }
            
            // Initialize outputs (core needs to be rebuilt for this)
            if (!outputInitialized) {
                if (!initializeOutputs(s)) {
                    s.extend("Failed to initialize results output.");
                    runOk = false;
                    break;
                }
                outputInitialized = true;
            }

            // Restore analysis state of the sweep whose index increased 
            // during the last call to advance()
            bool restoredState = false;
            restoredState = restoreState(advancedSweepIndex);
            
            if (sweepDebug>1) {
                Simulator::dbg() << "Invoking analysis" << (restoredState ? " in continuation mode" : "") << ".\n"; 
            }

            // Do pre-analysis computations
            if (!circuit.preAnalysis(s)) {
                s.extend("Pre-analysis computations failed.");
                runOk = false;
                break;
            }

            // Invoke analysis, dump results
            // Continue mode on if state has been restored
            runOk = runCores(restoredState, s);
            if (circuit.checkFlags(Circuit::Flags::Abort)) {
                // This is an error, stop sweep, return error
                s.set(Status::AbortRequested, std::string("Analysis '")+std::string(name_)+"' - aborted.");
                runOk = false;
                break;
            } else if (circuit.checkFlags(Circuit::Flags::Finish)) {
                // This is not an error, write debug notification and continue sweep
                if (sweepDebug>1) {
                    Simulator::dbg() << "Analysis '"+std::string(name_)+"' - finish requested.\n";
                }
            } else if (circuit.checkFlags(Circuit::Flags::Stop)) {
                // This is not an error, write debug notification and stop sweep
                // In future we will have to implement some kind of continue analysis option
                if (sweepDebug>1) {
                    Simulator::dbg() << "Analysis '"+std::string(name_)+"' - stop requested.\n";
                    Simulator::dbg() << "Exiting sweep @ "+sweeper->progress()+".\n";
                }
                break;
            } else if (!runOk) {
                s.set(Status::AnalysisFailed, std::string("Analysis '")+std::string(name_)+"' failed.");
                break;
            }

            // Advance sweeper
            bool finished;
            std::tie(finished, advancedSweepIndex) = sweeper->advance();
            if (finished) {
                if (sweepDebug>0) {
                    Simulator::dbg() << "Sweep finished.\n";
                }
                break;
            }

            // If we don't have stored state, initialize states of all sweeps (start with sweep 0)
            // otherwise store state of the sweep whose index has been advanced and all its inner sweeps
            auto updateStatesFrom = advancedSweepIndex;
            if (!haveStoredState) {
                updateStatesFrom = 0;
                haveStoredState = true;
            }

            // Update stored analysis state for the sweep that advanced its index and all its inner sweeps
            // Store state only if continuation is used for that particular sweep
            for(auto i=updateStatesFrom; i<sweepCount(); i++) {
                if (sweeps_.at(i).continuation) {
                    if (sweepDebug>1) {
                        Simulator::dbg() << "Updating analysis state for sweep level " << (i+1) << ".\n";
                    }

                    storeState(i);
                }
            }

            // Check Finish and Stop flags, for now they do the same as we have no means to contnue analysis
            if (circuit.checkFlags(Circuit::Flags::Finish) || circuit.checkFlags(Circuit::Flags::Stop)) {
                break;
            }
        } while (true);

        if (!runOk) { 
            s.extend("Sweep aborted @ "+sweeper->progress()+".");
        } 

        // Finalize outputs
        if (outputInitialized) {
            if (!finalizeOutputs(s)) {
                s.extend("Failed to finalize results output.");
                runOk = false;
            }
        }

        // Delete output files on failure if strictoutput>=1 and analysis failed
        auto strictoutput = circuit.simulatorOptions().core().strictoutput;
        if (!runOk && strictoutput>=1) {
            deleteOutputs();
        }

        // Restore original state
        auto [ok, hierarchyChanged, needsCoreRebuild] = circuit.elaborateChanges(
            sweeper.get(), ParameterSweeper::WriteValues::StoredState, 
            this, &originalSimOptions, 
            &parameterizedOptions, 
            s
        );
        if (!ok) {
            s.extend("Failed to restore initial circuit state.");
            runOk = false;
        }
    } else {
        // No sweep
        bool outputInitialized = false;

        // Set analysis options
        auto [ok, hierarchyChanged, needsCoreRebuild] = circuit.elaborateChanges(
            nullptr, ParameterSweeper::WriteValues::Sweep, 
            this, &simOptions, 
            &parameterizedOptions, 
            s
        );
        if (!ok) {
            s.extend("Failed to set circuit state.");
            runOk = false;
        }

        // Delete output files if strictoutput>=2
        auto strictoutput = circuit.simulatorOptions().core().strictoutput;
        if (strictoutput>=2) {
            deleteOutputs();
        }
        
        // Create list of output descriptors
        if (ok) {
            ok = addOutputDescriptors(s);
            if (!ok) {
                s.extend("Failed to construct list of outputs.");
                runOk = false;
            }
        }
    
        // Bind outputs
        if (ok) {
            // This happens only once per analysis run so treat it as first binding
            bool strict = circuit.simulatorOptions().core().strictsave>0;
            ok = resolveOutputDescriptors(strict, s);
            if (!ok) {
                s.extend("Failed to bind analysis outputs.");
                runOk = false;
            }
        }

        // Rebuild core structures for simulation (always, because there is only one core run)
        if (ok) {
            ok = rebuildCores(s);
            if (!ok) {
                s.extend("Failed to rebuild analysis structures.");
                runOk = false;
            }    
        }

        // Initialize outputs
        if (ok) {
            outputInitialized = initializeOutputs(s);
            if (!outputInitialized) {
                s.extend("Failed to initialize results output.");
                runOk = false;
            }

            // Do pre-analysis computations
            bool preAnalysisOk = false;
            if (outputInitialized) {
                preAnalysisOk =  circuit.preAnalysis(s);
                if (!preAnalysisOk) {
                    s.extend("Pre-analysis computations failed.");
                    runOk = false;
                }
            }
            
            // Rune core simulation, dump results
            if (outputInitialized && preAnalysisOk) {
                // Continue mode off
                runOk = runCores(false, s);
                if (circuit.checkFlags(Circuit::Flags::Abort)) {
                    s.set(Status::AbortRequested, std::string("Abort requsted."));
                    runOk = false;
                } else if (circuit.checkFlags(Circuit::Flags::Finish)) {
                    // This is not an error. 
                    if (sweepDebug>1) {
                        Simulator::dbg() << "Analysis '"+std::string(name_)+"' - finish requested.\n";
                    }
                } else if (circuit.checkFlags(Circuit::Flags::Stop)) {
                    // This is not an error. 
                    // In future we will have to implement some kind of continue analysis option
                    if (sweepDebug>1) {
                        Simulator::dbg() << "Analysis '"+std::string(name_)+"' - stop requested.\n";
                    }
                } else if (!runOk) {
                    s.set(Status::AnalysisFailed, std::string("Analysis '")+std::string(name_)+"' failed.");
                }
            }

            // Finalize outputs
            if (outputInitialized) {
                ok = finalizeOutputs(s);
                if (!ok) {
                    s.extend("Failed to finalize results output.");
                    runOk = false;
                }
            }

            // Delete output files on failure if strictoutput>=1 and analysis failed
            if (!ok && strictoutput>=1) {
                deleteOutputs();
            }
        }

        
        // Restore original state (interactive simulator mode)
        std::tie(ok, hierarchyChanged, needsCoreRebuild) = circuit.elaborateChanges(
            nullptr, ParameterSweeper::WriteValues::StoredState, 
            this, &originalSimOptions, 
            &parameterizedOptions, 
            s
        );
        if (!ok) {
            s.extend("Failed to restore initial circuit state.");
            runOk = false;
        }
    }

    return runOk;
}

std::tuple<bool, bool> Analysis::updateParameterExpressions(Status& s) {
    return parameters().setParameters(ptAnalysis.parameters().expressions(), circuit.variableEvaluator(), s);
}

void Analysis::dump(std::ostream& os) const {
    os << "Analysis " << std::string(name());
    os << std::endl; 

    os << "  Simulator options:" << std::endl;
    simOptions.dump(os, "    ");
    auto& params = parameters();
    if (params.parameterCount()>0) {
        os << std::endl;
        os << "  Analysis parameters:" << std::endl;
        params.dump(os, "    ");
    }
    os << std::endl;
}

}
