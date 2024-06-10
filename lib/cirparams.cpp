#include "circuit.h"
#include "an.h"
#include "common.h"


namespace NAMESPACE {

// Propagate parameter values down the hierarchy
// Check for hierarchy changes and rebuild parts of hierarchy if needed
std::tuple<bool, bool> Circuit::propagateDownHierarchy(Status& s) {
    // Normally propagation is applied only to instances whose parameters changed
    // force overrides this and applies propagation to all instances
    
    // Did hierarchy change?
    bool hierarchyChanged = false;

    // Do we need to check for hierarchy change in all instances? 
    bool forceHierarchyChangeCheck = checkFlags(Flags::HierarchyAffectingOptionsChanged);

    // Do we need to propagate parameters on all instances
    bool forceParameterPropagation = checkFlags(Flags::VariablesChanged); 

    // If hierarchy changed, at which toplevel instance did it change
    // A value greater or equal to the number of toplevel instances means 
    // there has been no hierarchy change. 
    auto n = toplevelInstances_.size(); 
    size_t changeAt = n;

    // Get context marker to revert stack in case of abort
    auto initialContextMarker = paramEvaluator_.contextMarker();

    // Propagate variables to toplevel circuits and check for hierarchy change in toplevel instances
    for(decltype(n) i=0; i<n; i++) {
        // Propagate variables
        
        // Propagate to the global context of the toplevel instance. 
        // Nothing to propagate to parameters because toplevel instances are
        // instantiated without any given instance parameters (defaults are used). 
        // The context is added to the path and recomputed if propagation is forced (variables changed). 
        auto [ok, contextMarker] = toplevelInstances_[i]->enterContext(*this, &(toplevelContext_[i]), true, forceParameterPropagation, s);
        if (!ok) {
            Instance::revertContext(*this, initialContextMarker);
            return std::make_tuple(false, hierarchyChanged);
        }
        
        // Check for hierarchy change
        if (forceParameterPropagation || forceHierarchyChangeCheck) {
            // Check for hierarchy change
            auto [ok, changed] = toplevelInstances_[i]->subhierarchyChanged(*this, s);
            hierarchyChanged |= changed;
            if (!ok) {
                // Error, revert to initial stack state and abort
                Instance::revertContext(*this, initialContextMarker);
                return std::make_tuple(false, hierarchyChanged);
            }
        }
        
        // Leave context only if this is not the default toplevel instance
        // The latter must remain in the stack and the path because its parameters are accessible 
        // in the whole hierarchy, even inside other toplevel instances. 
        if (i>0) {
            Instance::revertContext(*this, contextMarker);
        }

        // If hierarchy changed
        if (hierarchyChanged) {
            changeAt = i;
            break;
        }
    }

    // If there was a hierarchy change this loop will iterate because n>changeAt. 
    // Delete everything in reverse order for i=n-1 down to and including changeAt. 
    for(decltype(n) i=n-1; i-->changeAt;) {
        // Load the context, add to the path, but do not recompute. 
        // No need to check for error because loading a context cannot produce an error. 
        // Do this only for toplevel instances other than the default instance (0). 
        size_t contextMarker;
        if (i>0) {
            bool ok;
            std::tie(ok, contextMarker) = toplevelInstances_[i]->enterContext(*this, &(toplevelContext_[i]), true, false, s);
        }
        // Delete hierarchy
        if (!toplevelInstances_[i]->deleteHierarchy(*this, s))  {
            // Error, revert to initial stack state and abort
            Instance::revertContext(*this, initialContextMarker);
            return std::make_tuple(false, hierarchyChanged);
        }
        // Leave context, but only for toplevel instances other than defalt (0)
        if (i>0) {
            Instance::revertContext(*this, contextMarker);
        }
    }
    
    // Force parameter propagation 
    // - if variables changed
    // - if hierarchy affecting options changed
    //   because they may cause hierarchy change anywhere in the hierarchy
    bool force = checkFlags(Flags::VariablesChanged) || checkFlags(Flags::HierarchyAffectingOptionsChanged);

    // Now propagate parameters for i=0 to changeAt-1
    // In this part toplevel instances did not experience topology change
    for(decltype(n) i=0; i<changeAt; i++) {
        auto toplevelInstance = toplevelInstances_[i];
        // Parameter change in default toplevel instance must be propagated across all instances
        if (i==0) {
            force |= toplevelInstance->checkFlags(Instance::Flags::ParamsChanged);
        }
        // Load toplevel instance context and add it to path. 
        // No need to recompute it, we already did that. 
        // Do this only for i>0. 
        // Context of the default toplevel instance (i=0) is already there. 
        if (i>0) {
            paramEvaluator_.contextStack().enter(&(toplevelContext_[i]), true);
        }
        for(auto it=toplevelInstance->hierarchyBegin(); it!=toplevelInstance->hierarchyEnd(); ++it) {
            // Is the instance we are processing the toplevel instance? 
            bool isToplevel = (&*it == toplevelInstance);
            // Do we need to propagate? 
            if (force || it->checkFlags(Instance::Flags::ParamsChanged)) {
                size_t contextMarker;
                if (!isToplevel) {
                    // For toplevel instance the context is already there
                    // For others enter context, do not add it to the path, and recompute it. 
                    bool ok;
                    std::tie(ok, contextMarker) = it->enterContext(*this, nullptr, false, true, s);
                    if (!ok) {
                        Instance::revertContext(*this, initialContextMarker);
                        return std::make_tuple(false, hierarchyChanged);
                    }
                }
                // std::cout << "Context of " << std::string(it->name()) << std::endl;
                // std::cout << paramEvaluator_.contextStack() << std::endl;

                // Check if there was a subhierarchy change
                auto [ok, instanceSubhierarchyChanged] = it->subhierarchyChanged(*this, s);
                if (!ok) {
                    Instance::revertContext(*this, initialContextMarker);
                    return std::make_tuple(false, hierarchyChanged);
                }
                if (hierarchyChanged) {
                    // Rebuild instance's subhierarchy
                    // Delete subhierarchy 
                    if (!it->deleteHierarchy(*this, s)) {
                        s.extend("Failed to delete subhierarchy.");
                        Instance::revertContext(*this, initialContextMarker);
                        return std::make_tuple(false, hierarchyChanged);
                    }
                    // Build subhierarchy
                    InstantiationData idata(&*it);
                    if (!it->buildHierarchy(*this, paramEvaluator_, idata, s)) {
                        s.extend("Failed to rebuild subhierarchy.");
                        Instance::revertContext(*this, initialContextMarker);
                        return std::make_tuple(false, hierarchyChanged);
                    }
                } else {
                    // No need to rebuild hierarchy, just propagate parameters to immediate descendents
                    auto ok = it->propagateParameters(*this, paramEvaluator_, s);
                    if (!ok) {
                        Instance::revertContext(*this, initialContextMarker);
                        return std::make_tuple(false, hierarchyChanged);
                    }
                }
                // Remember subhierarchy change
                hierarchyChanged |= instanceSubhierarchyChanged;
                // If subhierarchy changed, stop descending along this instance
                if (instanceSubhierarchyChanged) {
                    it.stopDescent();
                }

                // Leave context, but only if this is not the toplevel instance
                if (!isToplevel) {
                    Instance::revertContext(*this, contextMarker);
                }
            }
        }
        // Leave context of toplevel instance, but not for default toplevel instance
        if (i>0) {
            paramEvaluator_.contextStack().exit();
        }
    }

    // Rebuild for i=changeAt to n-1
    for(decltype(n) i=changeAt; i<n; i++) {
        // Load the context, add to the path, but do not recompute. 
        // No need to check for error because loading a context cannot produce an error. 
        // Do this only for toplevel instances other than the default instance (0). 
        size_t contextMarker;
        if (i>0) {
            bool ok;
            std::tie(ok, contextMarker) = toplevelInstances_[i]->enterContext(*this, &(toplevelContext_[i]), true, false, s);
        }
        // Build hierarchy
        InstantiationData idata(toplevelInstances_[i]);
        if (!toplevelInstances_[i]->buildHierarchy(*this, paramEvaluator_, idata, s))  {
            // Error, revert to initial stack state and abort
            Instance::revertContext(*this, initialContextMarker);
            return std::make_tuple(false, hierarchyChanged);
        }
        // Leave context
        if (i>0) {
            paramEvaluator_.contextStack().exit();
        }
    }

    // Revert to initial stack state
    Instance::revertContext(*this, initialContextMarker);

    return std::make_tuple(true, hierarchyChanged);
}

// If sweeper is not nullptr
//   - if what is ParameterSweeper::WriteValues::Sweep
//       set everything to current sweeper state
//    - else 
//       set everything to stored sweeper state
// Writes to instances/models/options to which the sweeper has been bound. 
// 
// If an is not nullptr
//   - update analysis parameters according to new variable values
// 
// opt specifies simulator options that will be updated with
// variables and options expressions from optionsMap. 
// If opt is nullptr circuit's simulator options are used. 
// 
// Simulator options pointed to by opt are written to circuit's simulator options. 
// 
// Add sparsity map entries and state vector slots requested by analysis. 
// 
// Make circuit consistent. 

std::tuple<bool, bool, bool> Circuit::elaborateChanges(
    ParameterSweeper* sweeper, int advancedSweepIndex, ParameterSweeper::WriteValues what, 
    Analysis* an, IStruct<SimulatorOptions>* options, 
    PTParameterMap* optionsMap, 
    Status& s
) {
    auto t0 = Accounting::wclk();
    tables_.accounting().acctNew.sweep.chgelab++; 
    
    // Mark circuit as unelaborated
    clearFlags(Flags::Elaborated);

    // First, write variables if we have a sweeper
    if (sweeper) {
        auto [ok, changed] = sweeper->write(ParameterSweeper::ParameterFamily::Variable, what, s);
        if (!ok) {
            s.extend("Failed to write swept variables.");
            tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
            return std::make_tuple(false, false, false);
        }
    }

    // Did circuit variables change (also includes changes before swept global parameters are written)
    bool variablesChanged = checkFlags(Circuit::Flags::VariablesChanged);

    // Propagate variable changes to inner sweeps
    if (sweeper && variablesChanged) {
        if (!an->updateSweeper(advancedSweepIndex, s)) {
            return std::make_tuple(false, false, false);
        }
    }
    
    // Simulator options structure that will be updated and written to the circuit
    auto optPtr = options ? options : &simOptions;
    
    // If variables changed, rebuild global context, update options and analysis parameters
    if (variablesChanged) { 
        if (!updateGlobalContext(s)) {
            s.extend("Failed to update global context.");
            tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
            return std::make_tuple(false, false, false);
        }
        
        // Options expressions -> options
        if (optionsMap && optionsMap->size()>0) {
            if (auto [ok, changed] = optPtr->setParameters(*optionsMap, variableEvaluator_, s); !ok) {
                s.extend("Failed to apply options expressions.");
                tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
                return std::make_tuple(false, false, false);
            }
        }
        
        // Do we have an analysis
        if (an) {
            // Update analysis parameters with expressions
            if (auto [ok, changed] = an->updateParameterExpressions(s); !ok) {
                s.extend("Failed to update analysis parameters.");
                tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
                return std::make_tuple(false, false, false);
            }
        }
    }

    // Write swept options to places where scalar sweeps are bound
    if (sweeper) { 
        if (auto [ok, _] = sweeper->write(ParameterSweeper::ParameterFamily::Option, what, s); !ok) {
            s.extend("Failed to write swept options.");
            tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
            return std::make_tuple(false, false, false);
        }
    }

    // Write options to the circuit (if needed)
    if (optPtr!=&simOptions) {
        setOptions(*optPtr);
    }
    bool hierarchyAffectingOptionsChanged = checkFlags(Circuit::Flags::HierarchyAffectingOptionsChanged);
    bool mappingAffectingOptionsChanged = checkFlags(Circuit::Flags::MappingAffectingOptionsChanged);

    // Write instance and model parameters if we have a sweeper
    if (sweeper) {
        auto [ok, changed] = sweeper->write(ParameterSweeper::ParameterFamily::Instance | ParameterSweeper::ParameterFamily::Model, what, s); 
        if (!ok) {
            s.extend("Failed to write swept instance and model parameters.");
            tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
            return std::make_tuple(false, false, false);
        }
        if (changed) {
            setFlags(Flags::HierarchyParametersChanged);
        }
    }

    // Add the state of HierarchyParametersChanged flag
    auto hierarchyParametersChanged = checkFlags(Flags::HierarchyParametersChanged);

    // Propagate parameters and check for hierarchy changes if one of the following is true
    // - variables changed (force propagation on all instances)
    // - hierarchy affecting options changed (propagate only on instances whose parameters changed)
    // - swept instance/model parameters changed (propagate only on instances whose parameters changed)
    bool hierarchyChanged = false;
    if (variablesChanged || hierarchyParametersChanged || hierarchyAffectingOptionsChanged) {
        bool ok;
        std::tie(ok, hierarchyChanged) = propagateDownHierarchy(s);
        if (!ok) {
            s.extend("Failed to propagate changes down hierarchy.");
            tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
            return std::make_tuple(false, false, false);
        }
    }

    // If hierarchy changed, rebuild entity lists and node ordering, rebind sweep
    if (hierarchyChanged) {
        if (!buildEntityLists(s)) {
            s.extend("Failed to rebuild entity lists.");
            tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
            return std::make_tuple(false, false, false);
        }

        // Order nodes
        if (!nodeOrdering(s)) {
            s.extend("Failed to rebuild node ordering.");
            tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
            return std::make_tuple(false, false, false);
        }
    }

    // Now we are ready for calling setup

    // Do a setup if
    // - hierarchy changed
    // - variables changed
    // - hierarchical parameters changed
    // - mappingAffectingOptionsChanged (force full setup)
    bool unknownsChanged = false;
    bool sparsityChanged = false;
    if (hierarchyChanged || variablesChanged || hierarchyParametersChanged || mappingAffectingOptionsChanged) {
        bool ok;
        std::tie(ok, unknownsChanged, sparsityChanged) = setup(mappingAffectingOptionsChanged, s);
        if (!ok) {
            s.extend("Circuit setup failed.");
            tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
            return std::make_tuple(false, false, false);
        }
    }

    // Check if analysis needs to add entries to sparsity map and states vector
    bool analysisNeedsSparsity = false;
    if (an) {
        bool preMappingOk;
        std::tie(preMappingOk, analysisNeedsSparsity) = an->preMapping(s);
        if (!preMappingOk) {
            s.extend("Failed to pre-map analysis.");
            tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
            return std::make_tuple(false, false, false);
        }
    }

    // Perform mapping of unknowns to nodes if
    // - hierarchy changed
    // - node collapsing changed
    // - analysis needs to add sparsity map entries
    bool mapUnknownsNeeded = hierarchyChanged || unknownsChanged;
    if (mapUnknownsNeeded) {
        if (!mapUnknowns(s)) {
            s.extend("Failed to map unknowns to nodes.");
            tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
            return std::make_tuple(false, false, false);
        }
    }

    // Perform sparsity building and states allocation 
    // - unknowns were remapped
    // - analysis requires sparsity map entries
    bool buildNeeded = mapUnknownsNeeded || analysisNeedsSparsity; 
    if (buildNeeded) {
        if (!buildSparsityAndStates(s)) {
            s.extend("Failed to create sparsity pattern and allocate states.");
            tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
            return std::make_tuple(false, false, false);
        }
    }

    // Populate structures with parts needed by the analysis
    if (an && analysisNeedsSparsity && !an->populateStructures(s)) { 
        s.extend("Failed to map analysis.");
        tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
        return std::make_tuple(false, false, false);
    }

    // Enumerate Jacobian entries
    if (buildNeeded && !enumerateSystem(s)) {
        s.extend("System enumeration failed.");
        tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
        return std::make_tuple(false, false, false);
    }

    // Circuit is now in a consistent state
    clearFlags(Circuit::Flags::VariablesChanged);
    clearFlags(Circuit::Flags::HierarchyAffectingOptionsChanged);
    clearFlags(Circuit::Flags::MappingAffectingOptionsChanged);
    clearFlags(Circuit::Flags::HierarchyParametersChanged);

    // Mark circuit as elaborated
    setFlags(Flags::Elaborated);
    
    tables_.accounting().acctNew.sweep.tchgelab += Accounting::wclkDelta(t0);
    return std::make_tuple(true, hierarchyChanged, mapUnknownsNeeded || buildNeeded);
}

}
