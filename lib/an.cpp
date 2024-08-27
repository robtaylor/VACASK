#include "an.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

Analysis::Analysis(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : name_(name), circuit(circuit), sweeper(circuit, ptAnalysis.sweeps()), 
      ptAnalysis(ptAnalysis), commonSaves(nullptr), progressReporter(nullptr) {
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
    
    // Build analysis state repository
    an->resizeAnalysisStateStorage(an->sweepCount());

    return an;
}

bool Analysis::addOutputDescriptors(Status& s) {
    // Clear old descriptors in all cores
    clearOutputDescriptors();

    // Add sweep variables
    auto nSweeps = sweepCount();
    for(decltype(nSweeps) i=0; i<nSweeps; i++) {
        if (!addCommonOutputDescriptor(OutputDescriptor(OutdSweepvar, ptAnalysis.sweeps().data()[i].name(), i))) {
            s.set(Status::Analysis, "Failed to add output descriptor for sweep variable '"+std::string(ptAnalysis.sweeps().data()[i].name())+"'.");
            return false;
        }
    }

    // Add core output descriptors
    if (!addCoreOutputDescriptors(s)) {
        return false;
    }

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

    // Add default saves if needed (no saves specified)
    // Ignore errors and conflicts
    addDefaultOutputDescriptors();

    return true;
}

bool Analysis::resolveOutputDescriptor(const OutputDescriptor& descr, Output::SourcesList& srcs, bool strict) {
    // Abstract analysis handles only sweep variables
    switch (descr.type) {
        case OutdSweepvar: 
            srcs.emplace_back(&sweeper, descr.ndx);
            break;
        default:
            DBGCHECK(true, "Unknown output descriptor type.");
            return false;
    }
    return true;
}

bool Analysis::updateSweeper(Status& s) {
    return sweeper.update(advancedSweepIndex, s);
}

AnalysisCoroutine Analysis::coroutine(Status& s) {
    // Output descriptors are created with analysis
    // Binding must be done here, just before the core analysis is run

    // Clear instance flags: Converged, Bypassed, HasDeviceHistory
    circuit.applyInstanceFlags(
        Instance::Flags::HasDeviceHistory |
        Instance::Flags::Converged |
        Instance::Flags::Bypassed, 
        Instance::NoFlags
    );

    // Set simulator internals
    auto& options = circuit.simulatorOptions().core(); 
    SimulatorInternals internals;
    internals.fromOptions(options);
    internals.analysis_name = std::string(name_);
    internals.analysis_type = std::string(ptAnalysis.typeName());
    internals.initalizeLimiting = false;
    circuit.simulatorInternals() = internals;
    // By default high precision is not requested (bypass possible)
    internals.highPrecision = false;

    // Are we in debug mode
    auto debugMode = options.sweep_debug || options.op_debug || options.smsig_debug || options.tran_debug;
    if (debugMode && progressReporter) {
        progressReporter->disable();
    }

    // Stop and Finish do not interrupt the NR solver, 
    // nor do they interrupt the OP homotopy algorithms. 
    // They are left to be handled at a higher level (transient analysis). 
    // Abort interrupts immediately the NR solver, or the homotopy algorithm.

    auto sweepDebug = circuit.simulatorOptions().core().sweep_debug;

    // Skip this step to avoid wasting time
    // Mark as started
    // co_yield AnalysisState::Ready;

    // Do we have sweep(s)
    if (ptAnalysis.sweeps().data().size()>0) {
        // Sweep required

        if (progressReporter && !debugMode) {
            progressReporter->setValueFormat(ProgressReporter::ValueFormat::Fixed, 0);
            progressReporter->setValueDecoration("point# ", "");
            sweeper.install(progressReporter);
        }
        
        if (sweepDebug>0) {
            Simulator::dbg() << "Starting sweep, analysis '" << std::string(name_) << "'.\n";
        }

        // Setup sweeper
        if (!sweeper.setup(s)) {
            s.extend("Failed to set up sweep for analysis '"+std::string(name_)+"'.");
            co_yield AnalysisState::Aborted;
        }
        preSweepValuesStored = true;
        
        // Bind sweeper to actual parameters and options
        if (!sweeper.bind(circuit, simOptions, s)) {
            s.extend("Failed to bind sweep parameters.");
            co_yield AnalysisState::Aborted;
        }

        // Collect current values so we can restore them later
        if (!sweeper.storeState(s)) {
            s.extend("Failed to store initial circuit state.");
            co_yield AnalysisState::Aborted;
        }

        // Reset sweeper
        sweeper.reset();

        // Variable indicating that the outputs are bound
        bool outputsBound = false;

        // Main sweep loop
        bool haveStoredState = false;
        // Initially the outermost sweep index increases
        advancedSweepIndex = 0;
        do {
            if (sweepDebug>0) {
                Simulator::dbg() << "Sweep point: " << sweeper.progress() << ".\n";
            }

            // Indicate we have a new sweep point
            if (options.sweep_pointmarker) {
                co_yield AnalysisState::SweepPoint;
            }

            // Set current sweep point
            auto [ok, hierarchyChanged, needsCoreRebuild] = circuit.elaborateChanges(
                &sweeper, ParameterSweeper::WriteValues::Sweep, 
                this, &simOptions, 
                &parameterizedOptions, 
                // TODO: for now ignore devReq and Abort, Finish, Stop
                nullptr, 
                s
            );
            if (!ok) {
                s.extend("Failed to set circuit state.");
                co_yield AnalysisState::Aborted;
            }

            // Need to re-bind sweeper if hierarchy changed
            if (hierarchyChanged) {
                if (!sweeper.bind(circuit, simOptions, s)) {
                    s.extend("Failed to re-bind sweep parameters.");
                    co_yield AnalysisState::Aborted;
                }
            }

            // Do not allow analysis to use forced bypass by default
            circuit.simulatorInternals().allowForcedBypass = false;

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
                        co_yield AnalysisState::Aborted;
                    }
                }

                // Every time core needs rebuild, rebind outputs
                bool strict = (!outputsBound && circuit.simulatorOptions().core().strictsave>0) ||
                              (outputsBound && circuit.simulatorOptions().core().strictsave>1);
                if (!resolveOutputDescriptors(strict, s)) {
                    s.extend("Failed to bind analysis outputs.");
                    co_yield AnalysisState::Aborted;
                }

                // Rebuild core
                if (!rebuildCores(s)) {
                    s.extend("Failed to rebuild analysis structures.");
                    co_yield AnalysisState::Aborted;
                }
                outputsBound = true;
                coreRebuilt = true;

                // After core rebuild device history is invalidated
                circuit.applyInstanceFlags(
                    Instance::Flags::HasDeviceHistory |
                    Instance::Flags::Converged |
                    Instance::Flags::Bypassed, 
                    Instance::NoFlags
                );

                // Core rebuild makes all stored states incoherent with current topology
                for(decltype(sweepCount()) i=0; i<sweepCount(); i++) {
                    makeStateIncoherent(i);
                }
            } else {
                // If there was no topology change or core rebuild
                // and this is not the first point of the innermost sweep 
                // we allow the analysis to force the bypass of the first NR iteration
                // because the circuti state is the same as after the last computed 
                // NR iteration. The analysis decides on its own 
                // if it is going to force the bypass or not. Operating point 
                // analysis and all small-signal analyses allow the bypass
                // if ordinary continue mode is used. 
                if (sweeper.innermostSweepPosition()>0 && sweeper.continuation(sweeper.count()-1)) {
                    circuit.simulatorInternals().allowForcedBypass = (options.sweep_innerbypass!=0);
                }
            }
            
            // Initialize outputs (core needs to be rebuilt for this)
            if (!outputInitialized) {
                if (!initializeOutputs(s)) {
                    s.extend("Failed to initialize results output.");
                    co_yield AnalysisState::Aborted;
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
                co_yield AnalysisState::Aborted;
            }

            // Create coroutine, continue mode on if state has been restored
            auto cc = coreCoroutine(restoredState);
            while (cc) {
                bool exitLoop;
                auto state = cc.resume();
                switch (state) {
                    case CoreState::Aborted:
                        // Aborting a core aborts the analysis (and sweep)
                        formatCoreError(s);
                        s.extend("Analysis '"+std::string(name_)+"' aborted.");
                        s.extend("Sweep aborted @ "+sweeper.progress()+".");
                        co_yield AnalysisState::Aborted;
                    case CoreState::Finished:
                        // If a core finishes, the analysis is not finished 
                        // because we may still have some sweep points to process. 
                        if (sweepDebug>1) {
                            Simulator::dbg() << "Analysis '"+std::string(name_)+"' - finish requested.\n";
                            Simulator::dbg() << "Finishing sweep @ "+sweeper.progress()+".\n";
                        }
                        // We just exit the inner loop (stop running the core)
                        exitLoop = true;
                        break;
                    case CoreState::Stopped:
                        // Stopping a core stops the analysis (and sweep)
                        if (sweepDebug>1) {
                            Simulator::dbg() << "Analysis '"+std::string(name_)+"' - stop requested.\n";
                            Simulator::dbg() << "Stopping sweep @ "+sweeper.progress()+".\n";
                        }
                        co_yield AnalysisState::Stopped;
                    default:
                        // Any other state we don't know aborts the analysis
                        s.set(Status::AbortRequested, "Yielded value not supported. Aborting analysis '"+std::string(name_)+"'.");
                        co_yield AnalysisState::Aborted;
                }
                if (exitLoop) {
                    break;
                }
            }

            // Advance sweeper
            bool finished;
            std::tie(finished, advancedSweepIndex) = sweeper.advance();
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
                if (sweeper.continuation(i)) {
                    if (sweepDebug>1) {
                        Simulator::dbg() << "Updating analysis state for sweep level " << (i+1) << ".\n";
                    }
                    storeState(i);
                }
            }
        } while (true);
    } else {
        // No sweep
        advancedSweepIndex = 0;

        if (progressReporter && !debugMode) {
            analysisCore().install(progressReporter);
        }
        
        // Indicate we have a new sweep point
        if (options.sweep_pointmarker) {
            co_yield AnalysisState::SweepPoint;
        }
        
        // Set analysis options
        auto [ok, hierarchyChanged, needsCoreRebuild] = circuit.elaborateChanges(
            nullptr, ParameterSweeper::WriteValues::Sweep, 
            this, &simOptions, 
            &parameterizedOptions, 
            // TODO: for now ignore devReq and Abort, Finish, Stop
            nullptr, 
            s
        );
        if (!ok) {
            s.extend("Failed to set circuit state.");
            co_yield AnalysisState::Aborted;
        }

        // Create list of output descriptors
        ok = addOutputDescriptors(s);
        if (!ok) {
            s.extend("Failed to construct list of outputs.");
            co_yield AnalysisState::Aborted;
        }
        
        // Bind outputs
        // This happens only once per analysis run so treat it as first binding
        bool strict = circuit.simulatorOptions().core().strictsave>0;
        ok = resolveOutputDescriptors(strict, s);
        if (!ok) {
            s.extend("Failed to bind analysis outputs.");
            co_yield AnalysisState::Aborted;
        }

        // Rebuild core structures for simulation (always, because there is only one core run)
        ok = rebuildCores(s);
        if (!ok) {
            s.extend("Failed to rebuild analysis structures.");
            co_yield AnalysisState::Aborted;
        }
        // After core rebuild device history is invalidated
        circuit.applyInstanceFlags(
            Instance::Flags::HasDeviceHistory |
            Instance::Flags::Converged |
            Instance::Flags::Bypassed, 
            Instance::NoFlags
        );
        
        // Initialize outputs
        outputInitialized = initializeOutputs(s);
        if (!outputInitialized) {
            s.extend("Failed to initialize results output.");
            co_yield AnalysisState::Aborted;
        }

        // Do pre-analysis computations
        bool preAnalysisOk = false;
        if (outputInitialized) {
            preAnalysisOk =  circuit.preAnalysis(s);
            if (!preAnalysisOk) {
                s.extend("Pre-analysis computations failed.");
                co_yield AnalysisState::Aborted;
            }
        }
            
        // Run core simulation, dump results
        if (outputInitialized && preAnalysisOk) {
            // Create coroutine, continue mode off
            auto cc = coreCoroutine(false);
            while (cc) {
                auto state = cc.resume();
                switch (state) {
                    case CoreState::Aborted:
                        formatCoreError(s);
                        s.extend("Analysis '"+std::string(name_)+"' aborted.");
                        co_yield AnalysisState::Aborted;
                    case CoreState::Finished:
                        if (sweepDebug>1) {
                            Simulator::dbg() << "Analysis '"+std::string(name_)+"' - finish requested.\n";
                        }
                        // THis is not a swept analysis - if core finishes, so does the analysis
                        co_yield AnalysisState::Finished;
                    case CoreState::Stopped:
                        if (sweepDebug>1) {
                            Simulator::dbg() << "Analysis '"+std::string(name_)+"' - stop requested.\n";
                        }
                        co_yield AnalysisState::Stopped;
                    default:
                        s.set(Status::AbortRequested, "Yielded value not supported. Aborting analysis.");
                        co_yield AnalysisState::Aborted;
                }
            }
        }
    }

    // If we reach this point the analysis is finished
    co_yield AnalysisState::Finished;
}

bool Analysis::start(Status& s) {
    // Do we have a coroutine that is running? 
    // A finished coroutine can be restarted. 
    if (coroutine_.isValid()) {
        // Cannot restart an already started analysis
        s.set(Status::AbortRequested, "Analysis is running. Cannot restart it.");
        return false;
    }
    coroutine_ = coroutine(s);
    
    // No pre-sweep values stored
    preSweepValuesStored = false;

    // Output is not initialized
    outputInitialized = false;

    // Last coroutine state
    lastCoroutineState = AnalysisState::Uninitilized;

    return true;
}

AnalysisState Analysis::resume() {
    // Cannot do this on a coroutine that is not running
    if (!coroutine_.isValid() || coroutine_.done()) {
        return AnalysisState::Aborted;
    }
    lastCoroutineState = coroutine_.resume();
    return lastCoroutineState;
}

bool Analysis::finish(Status& s) {
    bool ok = true;
    if (!coroutine_.isValid()) {
        s.set(Status::AbortRequested, "Analysis is not running.");
        return false;
    }

    // Cleanup old files in case analysis was aborted
    if (lastCoroutineState==AnalysisState::Aborted) {
        // Delete output files, even if outputs were not initialized
        auto strictoutput = circuit.simulatorOptions().core().strictoutput;
        if (strictoutput>=1) {
            deleteOutputs();
        }
        ok = false;
    } else {
        // Finalize output files
        if (outputInitialized) {
            ok = finalizeOutputs(s);
            if (!ok) {
                s.extend("Failed to finalize results output.");
            }
        }
    }
    
    // Restore original parameter values (if stored), ignore Abort, Finish, Stop requests
    if (preSweepValuesStored) {
        Status tmps;
        auto sptr = ptAnalysis.sweeps().data().size()>0 ? &sweeper : nullptr;
        auto [ok, hierarchyChanged, needsCoreRebuild] = circuit.elaborateChanges(
            sptr, ParameterSweeper::WriteValues::StoredState, 
            this, &originalSimOptions, 
            &parameterizedOptions, 
            nullptr, 
            tmps
        );
        if (!ok) {
            s.extend(tmps.message());
            s.extend("Failed to restore initial circuit state.");
            ok = false;
        }
    }

    // Clear coroutine
    coroutine_ = AnalysisCoroutine();
    return ok;
}

bool Analysis::run(Status& s) {
    if (!start(s)) {
        return false;
    }
    bool ok = true;
    while (true) {
        if (coroutine_.done()) {
            // Coroutine reached its end
            // This is an internal error because before exit the state must be eiter Aborted or Finished
            s.set(Status::AbortRequested, "Internal error. Coroutine exited without yielding Aborted or Finished.");
            ok = false;
            break;
        }
        auto state = resume();
        if (state==AnalysisState::Aborted) {
            // Aborting analysis
            // Do not collect status to avoid overwriting what resume() has put there
            finish();
            ok = false;
            break;
        } else if (state==AnalysisState::Finished) {
            // Finish requested
            ok = finish(s);
            break;
        }
    }
    return ok;
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
