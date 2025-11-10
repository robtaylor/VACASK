#include <iomanip>
#include <cmath>
#include <filesystem>
#include "coreop.h"
#include "simulator.h"
#include "hmtpgmin.h"
#include "hmtpsrc.h"
#include "common.h"

namespace NAMESPACE {

// Default parameters
OperatingPointParameters::OperatingPointParameters() {
}

// Introspection for parameters structure
template<> int Introspection<OperatingPointParameters>::setup() {
    registerMember(nodeset);
    registerMember(store); 
    registerMember(write); 
    return 0;
}
instantiateIntrospection(OperatingPointParameters);


OperatingPointCore::OperatingPointCore(
    OutputDescriptorResolver& parentResolver, OperatingPointParameters& params, Circuit& circuit, 
    CommonData& commons, 
    KluRealMatrix& jacobian, VectorRepository<double>& solution, VectorRepository<double>& states
) : AnalysisCore(parentResolver, circuit, commons), params(params), outfile(nullptr), 
      nrSolver(circuit, commons, jacobian, states, solution, nrSettings), 
      jac(jacobian), solution(solution), states(states), 
      converged_(false), continueState(nullptr) {
}

OperatingPointCore::~OperatingPointCore() {
    delete outfile;
}

bool OperatingPointCore::resolveOutputDescriptors(bool strict, Status& s) {
    // Clear output sources
    outputSources.clear();
    // Resolve output descriptors
    bool ok = true;
    for (auto it = outputDescriptors.cbegin(); it != outputDescriptors.cend(); ++it) {
        Node *node;
        Instance *inst;
        switch (it->type) {
        case OutdSolComponent:
            ok = addRealVarOutputSource(strict, it->id, solution);
            break;
        case OutdOutvar:
            ok = addOutvarOutputSource(strict, it->idId.id1, it->idId.id2);
            break;
        default:
            // Delegate to parent
            ok = parentResolver.resolveOutputDescriptor(*it, outputSources, strict);
            break;
        }
        if (!ok) {
            break;
        }
    }
    return ok;
}

bool OperatingPointCore::addDefaultOutputDescriptors() {
    // If output is suppressed, skip all this work
    if (!params.write || Simulator::noOutput()) {
        return true;
    }
    if (savesCount==0) {
        return addAllUnknowns(PTSave("default", Id(), Id()));
    }
    return true;
}

bool OperatingPointCore::initializeOutputs(Id name, Status& s) {
    if (!params.write || Simulator::noOutput()) {
        return true;
    }
    // Create output file if not created yet
    if (!outfile) {
        outfile = new OutputRawfile(
            name, outputDescriptors, outputSources,
            (circuit.simulatorOptions().core().rawfile==SimulatorOptions::rawfileBinary ? OutputRawfile::Flags::Binary : OutputRawfile::Flags::None) |
                OutputRawfile::Flags::Padded);
        outfile->setTitle(circuit.title());
        outfile->setPlotname("Operating Point");
    }
    outfile->prologue();

    return true;
}

bool OperatingPointCore::finalizeOutputs(Status &s) {
    if (outfile) {
        outfile->epilogue();
        delete outfile;
        outfile = nullptr;
    }

    // Write DC solution to repository if analysis is OK
    if (converged_ && params.store.length()>0) {
        auto sol = circuit.newStoredSolution("dc", params.store);
        sol->setNames(circuit);
        sol->setValues(solution.vector());
    }

    return true;
}

bool OperatingPointCore::deleteOutputs(Id name, Status &s) {
    if (!params.write || Simulator::noOutput()) {
        return true;
    }

    // Cannot assume outfile is available
    auto fname = std::string(name)+".raw";
    if (std::filesystem::exists(fname)) {
        std::filesystem::remove(fname);
    }
    return true;
}
    
bool OperatingPointCore::storeState(size_t ndx, bool storeDetails) {
    auto& repo = coreStates.at(ndx);
    // Store current solution as annotated solution
    if (storeDetails) {
        repo.solution.setNames(circuit);
    } else {
        repo.solution.clearNames();
    }
    repo.solution.setValues(solution.vector());
    // Store current state
    repo.solution.setAuxData(states.vector());
    // Stored state is coherent and valid
    repo.coherent = true;
    repo.valid = true;
    return true;
}

bool OperatingPointCore::restoreState(size_t ndx) {
    auto& state = coreStates.at(ndx);
    if (state.valid) {
        // State is valid
        continueState = &state;
        return true;
    } else { 
        // Nothing to restore, do not use continuation mode
        return false;
    }
}

std::tuple<bool, bool> OperatingPointCore::preMapping(Status& s) {
    // Go through nodesets. Decode, set, and check delta part. 
    // No need to check unknowns part because diagonal entries are
    // always allocated by Circuit. 
    auto& nsParam = params.nodeset;
    if (nsParam.type()==Value::Type::String) {
        // It is a string. 
        // Will retrieve DC solution and set it as nodeset. 
        // Nothing to do here because DC solutions do not introduce delta foces. 
        return std::make_tuple(true, false);
    } else if (nsParam.type()==Value::Type::ValueVec) {
        // It is a list
        // Retrieve it and check if we need extra matrix entries
        auto& valVec = nsParam.val<ValueVector>();
        auto [ok, needsMapping] = preprocessedNodeset.set(circuit, valVec, s);
        if (!ok) {
            s.extend("Failed to preprocess nodesets.");
        }
        return std::make_tuple(ok, needsMapping);
    } else {
        // Unsupported type, signal an error.
        s.set(Status::BadArguments, "Nodeset must be a list or a string.");
        return std::make_tuple(false, false);
    }
}

bool OperatingPointCore::populateStructures(Status& s) {
    // Go through node pairs and add entries to sparsity map
    for(auto& pair : preprocessedNodeset.nodePairs) {
        auto [node1, node2] = pair;
        if (!node1 | !node2) {
            // One of the two nodes was not found, ignore pair
            continue;
        }
        
        // Forcing entries are resistive
        if (auto [_, ok] = circuit.createJacobianEntry(node1, node2, EntryFlags::Resistive, s); !ok) {
            return false;
        }
        
        if (auto [_, ok] = circuit.createJacobianEntry(node2, node1, EntryFlags::Resistive, s); !ok) {
            return false;
        }
    }

    return true;
}


bool OperatingPointCore::rebuild(Status& s) {
    // Bind Jacobian entries
    // Resistive parts bound to entries of jac, reactive parts not bound
    if (!circuit.bind(&jac, Component::Real, std::nullopt, nullptr, Component::Real, std::nullopt, s)) {
        return false;
    }

    // Prepare NR solver settings
    auto& options = circuit.simulatorOptions().core();
    nrSettings = NRSettings {
        .debug = options.nr_debug, 
        .itlim = options.op_itl, 
        .itlimCont = options.op_itlcont, 
        .convIter = options.nr_conviter, 
        .residualCheck = bool(options.nr_residualcheck),  
        .dampingFactor = options.nr_damping, 
        .matrixCheck = bool(options.matrixcheck), 
        .rhsCheck = bool(options.rhscheck), 
        .solutionCheck = bool(options.solutioncheck), 
        .forceFactor = options.nr_force, 
    };

    // Create Forces from preprocessed nodesets in solver force slot 1. 
    // This has to be done now for forces that include delta forces
    // (e.g. user nodesets and ics) because nrSolver.rebuild() 
    // (which is called before nrSolver.run()) requires 
    // nrSolver.forcesList[].deltaIndices() to be set up. 
    // Setting up is done by nrSolver.forces().set(). 
    // We also handle ordinary nodesets from solution repository. 
    auto strictforce = circuit.simulatorOptions().core().strictforce; 
    if (params.nodeset.type()==Value::Type::String) {
        String& solutionName = params.nodeset.val<String>();
        if (solutionName.length()>0) {
            // Get solution from repository
            auto solPtr = circuit.storedSolution("dc", solutionName);
            if (!solPtr) {
                // No nodesets
                nrSolver.forces(1).clear();
                Simulator::wrn() << "Warning, solution '"+solutionName+"' not found. No user nodesets applied.\n";
            } else {
                // Nodesets from solution repository
                if (!nrSolver.setForces(1, *solPtr, strictforce)) {
                    // Abort if strictforce is set
                    if (strictforce) {
                        nrSolver.formatError(s);
                        return false;
                    }
                }
            }
        } else {
            // No nodesets, clear slot
            nrSolver.forces(1).clear();
        }
    } else if (params.nodeset.type()==Value::Type::ValueVec) {
        // A list with possibly delta forces
        if (!nrSolver.setForces(1, preprocessedNodeset, false, strictforce)) {
            // Abort on error if strictforce is set
            if (strictforce) {
                nrSolver.formatError(s);
                return false;
            }
        }
    } else {
        // Error
        s.set(Status::BadArguments, "Nodeset must be a list or a string.");
        return false;
    }
    
    // Rebuild NR solver structures
    if (!nrSolver.rebuild()) {
        s.set(Status::NonlinearSolver, "Failed to rebuild internal structures of nonlinear solver.");
        return false;
    }

    return true;
}

std::tuple<bool, bool> OperatingPointCore::runSolver(bool continuePrevious) {
    auto& options = circuit.simulatorOptions().core();
    auto strictforce = options.strictforce; 
    // Assume no initial state given, start with standard initial point. 
    // Coherence information is set by an.cpp and homotopy. 
    // They are responsible for detecting topology changes/rebuilds. 
    // Homotopy sweeps parameters that do not cause topology changes/rebuilds. 
    bool runInContinueMode = false;
    // Handle continuation
    if (continuePrevious) {
        // Continue mode
        if (continueState &&
            continueState->valid && continueState->coherent &&
            continueState->solution.values().size()==circuit.unknownCount()+1 &&
            continueState->solution.auxData().size()==circuit.statesCount() 
        ) {
            // Continue a state
            // State is valid, coherent, and its lengths match those of the solver vectors
            // Restore current state
            solution.vector() = continueState->solution.values();
            states.vector() = continueState->solution.auxData();
            runInContinueMode = true;
            // No forces applied
            nrSolver.enableForces(0, false);
            nrSolver.enableForces(1, false);
            if (options.op_debug>1) {
                Simulator::dbg() << "OP using ordinary continue mode with stored analysis state.\n";
            }
            // Use forced bypass if allowed
            commons.requestForcedBypass = commons.allowContinueStateBypass;
        } else if (continueState && continueState->valid) {
            // Stored analysis state is valid, but not coherent with current circuit, 
            // its lengths may not match those of the solver vectors. 
            // Use forces to continue, but set no initial states vector nor initial solution. 
            // Ignore forces conflicts arising from stored solution. Slot 0 is for sweep continuation and homotopy. 
            // There should be no such conflicts as we are applying forces to nodes only, not node deltas. 
            strictforce = false;
            if (!nrSolver.setForces(0, continueState->solution, strictforce)) {
                if (strictforce) {
                    // Failed, strictforce set
                    return std::make_tuple(false, false);
                }
            }
            nrSolver.enableForces(0, true);
            // Disable user-specified forces
            nrSolver.enableForces(1, false);
            if (options.op_debug>1) {
                Simulator::dbg() << "OP using forced continue mode with stored analysis state.\n";
            }
            // Forced bypass is not allowed
            commons.requestForcedBypass = false;
        } else {
            // No valid state, continue with whatever is in solution and states vector
            runInContinueMode = true;
            // No forces applied
            nrSolver.enableForces(0, false);
            nrSolver.enableForces(1, false);
            if (options.op_debug>1) {
                Simulator::dbg() << "OP using ordinary continue mode with previous solution.\n";
            }
            // Forced bypass is allowed if set outside
        }
        // Continue state is spent after first use
        continueState = nullptr;
    } else {
        // Continue mode not requested
        // Disable continuation forces in slot 0
        nrSolver.enableForces(0, false);

        // Apply forces specified by user in slot 1
        nrSolver.enableForces(1, true); 
        
        if (options.op_debug>1) {
            Simulator::dbg() << "OP using standard initial solution with forced nodesets.\n";
        }
    }

    auto converged = nrSolver.run(runInContinueMode);
    auto abort = nrSolver.checkFlags(OpNRSolver::Flags::Abort);
    return std::make_tuple(converged, abort);
}

Int OperatingPointCore::iterations() const {
    return nrSolver.iterations();
}

Int OperatingPointCore::iterationLimit(bool continuePrevious) const {
    return continuePrevious ? nrSettings.itlimCont : nrSettings.itlim;
}

// System of equations is 
//   g(x) = 0
//   |      
//   nonlinear vector-valued function at new solution
//
// Linearized system
//   g(xold) + G(xold) (xnew - xold) = 0
//   |         |
//   |         resistive Jacobian
//   function value (residual) at xold
//
// xold .. previous NR solution
// xnew .. new NR solution
//
// The system can be solved as 
//   xnew - xold = - inverse(G(xold)) g(xold)
//   xnew = xold - inverse(G(xold)) g(xold)
//        = xold + inverse(G(xold)) (-g(xold))
//
// We load the Jacobian and residuals, and solve
//   G(xold) deltax = -g(xold)
// 
// Then we compute new approximate solution as
//   xnew = xold + deltax
// TODO: verify reactive residual and limiting

// Abort is handled and reported here
// Finish and Stop are ignored in OP analysis as there is nothing to interrupt
// Error messages accumulate across
//   - initial operating point
//   - dynamic gmin stepping
//   - dynamic gshunt steppint
//   - spice3 gmin stepping
//   - source stepping
//   - spice3 source stepping
// Each message starts with a nrSolver error message
CoreCoroutine OperatingPointCore::coroutine(bool continuePrevious) {
    initProgress(1, 0);

    clearError();

    auto& options = circuit.simulatorOptions().core();
    converged_ = false;
    auto debug = options.op_debug;
    bool leave = false;
    bool tried = false;

    // Set time to 0
    nrSolver.evalSetup().time = 0.0;

    // Make sure repositories are large enough
    auto n = circuit.unknownCount();
    solution.upsize(2, n+1);
    states.upsize(2, circuit.statesCount());
    
    // Initial plain op
    auto skipinitial = options.op_skipinitial;
    if (!skipinitial) {    
        tried = true;
        std::tie(converged_, leave) = runSolver(continuePrevious);
        if (!converged_) {
            setError(OperatingPointError::InitialOp);
        }
        if (debug>0) {
            if (converged_) {
                Simulator::dbg() << "OP core algorithm converged in " << std::to_string(nrSolver.iterations()) << " NR iteration(s).\n";
            } else {
                Simulator::dbg() << "OP core algorithm failed to converge in " << std::to_string(nrSolver.iterations()) << " NR iteration(s).\n";
            }
        }
    }

    // Try homotopy
    if (!converged_ && !leave && options.op_homotopy.size()>0) {
        Homotopy* homotopy;
        for(auto it : options.op_homotopy) {
            if (it==Homotopy::gdev) {
                if (debug>0) {
                    Simulator::dbg() << "Trying gdev stepping.\n";
                }
                homotopy = new GminStepping(circuit, *this, true);
            } else if (it==Homotopy::gshunt) {
                if (debug>0) {
                    Simulator::dbg() << "Trying gshunt stepping.\n";
                }
                homotopy = new GminStepping(circuit, *this, false);
            } else if (it==Homotopy::spice3Gmin) {
                if (debug>0) {
                    Simulator::dbg() << "Trying SPICE3 gmin stepping.\n";
                }
                homotopy = new Spice3GminStepping(circuit, *this);
            } else if (it==Homotopy::src) {
                if (debug>0) {
                    Simulator::dbg() << "Trying source stepping.\n";
                }
                homotopy = new SourceStepping(circuit, *this, options.op_srchomotopy);
            } else if (it==Homotopy::spice3Src) {
                if (debug>0) {
                    Simulator::dbg() << "Trying SPICE3 source stepping.\n";
                }
                homotopy = new Spice3SourceStepping(circuit, *this);
            } else {
                if (debug>0) {
                    Simulator::dbg() << "Unknown homotopy '"+std::string(it)+"'.\n";
                }
                homotopy = nullptr;
            }
            if (!homotopy) {
                continue;
            }
            // Run
            tried = true;
            std::tie(converged_, leave) = homotopy->run();
            if (debug>0) {
                if (converged_) {
                    Simulator::dbg() << "Homotopy converged in " << std::to_string(homotopy->stepCount()) << " step(s).\n";
                } else {
                    Simulator::dbg() << "Homotopy failed to converge in " << std::to_string(homotopy->stepCount()) << " step(s).\n";
                }
            }
            delete homotopy;
            if (leave || converged_) {
                break;
            }
        }
        if (!converged_) {
            setError(OperatingPointError::Homotopy);
        }
    }
    
    if (!leave) {
        // Did not leave early
        if (!tried) {
            // No algorithm tried
            setError(OperatingPointError::NoAlgorithm);
        } else if (converged_) {
            // Tried and converged, write results
            if (outfile && params.write) {
                outfile->addPoint();
            }
        }
    } else {
        // Leaving early, did not converge
        // Add a status message one level higher
        converged_ = false;
    }

    setProgress(1);

    // OP analysis can only Abort or Finish
    if (converged_) {
        co_yield CoreState::Finished;
    } else {
        co_yield CoreState::Aborted;
    }
}

bool OperatingPointCore::run(bool continuePrevious) {
    auto c = coroutine(continuePrevious);
    bool ok = true;
    while (!c.done()) {
        if (c.resume()==CoreState::Aborted) {
            ok = false;
            break;
        };
    }
    return ok;
}

bool OperatingPointCore::formatError(Status& s) const {
    auto nr = UnknownNameResolver(circuit);
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4);

    // Delegate to NRSolver (which in turn delegates to KluMatrix)
    auto solverError = nrSolver.formatError(s, &nr);
    
    // First, handle AnalysisCore errors
    if (lastError!=Error::OK) {
        AnalysisCore::formatError(s);
        return false;
    }
    
    // Then handle OperatingPointCore errors
    switch (lastOpError) {
        case OperatingPointError::InitialOp:
            s.extend("Initial OP analysis failed.");
            return false;
        case OperatingPointError::Homotopy:
            s.set(Status::Analysis, "Homotopy failed.");
            return false;
        case OperatingPointError::NoAlgorithm:
            s.set(Status::Analysis, "No operating point algorithm tried."); 
            return false;
    }
    return solverError;
}

void OperatingPointCore::dump(std::ostream& os) const {
    AnalysisCore::dump(os);
    os << "  Solver ";
    if (converged_) {
        os << "converged";
    } else {
        os << "not converged";
    }
    os << ", last NR run completed with " << std::to_string(nrSolver.iterations()) << " iterations";
    os << std::endl;
    os << "  Results: " << "\n";
    auto n = circuit.unknownCount();
    for(decltype(n) i=1; i<=n; i++) {
        auto rn = circuit.reprNode(i);
        os << "    " << rn->name() << " : " << solution.data()[i] << "\n";
    }
}

}


/*
  op analysis
  - start mode
    start with all zero solution
    use nodeset parameter for nodesets
  - continue mode
    ignore nodeset parameter
    - consistent topology
      start with stored solution
    - inconsistent topology
      use nodesets based on stored solution
      start with all zero solution
  tran analysis
  - icmode = uic
    - start mode and continue mode
      ignore nodeset parameter
      use ic parameter to build initial solution vector
  - icmode = op
    - start mode
      start with all zero solution
      use nodesets in first iterations
      use ic in all iterations (force)
    - continue mode
      ignore nodeset parameter
      - consistent topology
        start with stored solution
        use ic in all iterations (force)
      - inconsistent topology
        start with all zero solution
        use nodesets based on stored solution in first iterations 
        use ic in all iterations (force)
*/
