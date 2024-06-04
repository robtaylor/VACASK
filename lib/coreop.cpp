#include <iomanip>
#include <cmath>
#include <filesystem>
#include "coreop.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

// Default parameters
OpParameters::OpParameters() {
}

// Introspection for parameters structure
template<> int Introspection<OpParameters>::setup() {
    registerMember(nodeset);
    registerMember(store); 
    return 0;
}
instantiateIntrospection(OpParameters);


OperatingPointCore::OperatingPointCore(
    Analysis& analysis, OpParameters& params, Circuit& circuit, 
    KluRealMatrix& jacobian, VectorRepository<double>& solution, VectorRepository<double>& states
) : AnalysisCore(analysis, circuit), params(params), outfile(nullptr), 
      nrSolver(circuit, jac, states, solution, nrSettings), 
      jac(jacobian), solution(solution), states(states), 
      converged(false), continueState(nullptr) {
}

OperatingPointCore::~OperatingPointCore() {
    delete outfile;
}

// Implement this in every derived class so that calls to 
// resolveOutputDescriptor() will be inlined. 
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
        case OutdOpvar:
            ok = addOpvarOutputSource(strict, it->idId.id1, it->idId.id2);
            break;
        default:
            // Delegate to parent
            ok = analysis.resolveOutputDescriptor(*it, outputSources, strict);
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
    if (!params.writeOutput) {
        return true;
    }
    if (savesCount==0) {
        return addAllUnknowns(PTSave(Loc::bad, "default", Id(), Id()));
    }
    return true;
}

bool OperatingPointCore::initializeOutputs(Id name, Status& s) {
    if (!params.writeOutput) {
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
    if (!params.writeOutput) {
        return true;
    }
    outfile->epilogue();
    delete outfile;
    outfile = nullptr;

    // Write DC solution to repository if analysis is OK
    if (converged && params.store.length()>0) {
        Id label = params.store;
        circuit.storeDcSolution(params.store, solution.vector());
    }

    return true;
}

bool OperatingPointCore::deleteOutputs(Id name, Status &s) {
    if (!params.writeOutput) {
        return true;
    }

    // Cannot assume outfile is available
    auto fname = std::string(name)+".raw";
    if (std::filesystem::exists(fname)) {
        std::filesystem::remove(fname);
    }
    return true;
}
    
size_t OperatingPointCore::stateStorageSize() const { 
    return analysisStateRepository.size(); 
}

void OperatingPointCore::resizeStateStorage(size_t n) { 
    analysisStateRepository.resize(n); 
    // Initially no state is coherent nor valid
    for(auto& it : analysisStateRepository) {
        it.coherent = false;
        it.valid = false;
    }
}

bool OperatingPointCore::storeState(size_t ndx) {
    auto& repo = analysisStateRepository.at(ndx);
    // Store current solution as annotated solution
    repo.solution.set(circuit, solution.vector());
    // Store current state
    repo.stateVector = states.vector();
    // Stored state is coherent and valid
    repo.coherent = true;
    repo.valid = true;
    return true;
}

bool OperatingPointCore::restoreState(size_t ndx) {
    auto& state = analysisStateRepository.at(ndx);
    if (state.valid) {
        // State is valid
        continueState = &state;
        return true;
    } else { 
        // Nothing to restore, do not use continuation mode
        return false;
    }
}

void OperatingPointCore::makeStateIncoherent(size_t ndx) {
    analysisStateRepository.at(ndx).coherent = false;
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
        auto [ptr1, ok1] = circuit.createJacobianEntry(node1, node2, s);
        if (!ok1) {
            return false;
        }
        auto [ptr2, ok2] = circuit.createJacobianEntry(node2, node1, s);
        if (!ok2) {
            return false;
        }
    }

    return true;
}


bool OperatingPointCore::rebuild(Status& s) {
    // Bind Jacobian entries
    // Resistive parts bound to entries of jac, reactive parts not bound
    if (!circuit.bind(&jac, nullptr, Component::RealPart, nullptr, nullptr, Component::RealPart, s)) {
        return false;
    }

    // Prepare NR solver settings
    auto& options = circuit.simulatorOptions().core();
    nrSettings = NRSettings {
        .debug = options.op_debug >= 100 ? options.op_debug-100+1 : 0, 
        .itlim = options.op_itl, 
        .itlimCont = options.op_itlcont, 
        .convIter = options.nr_conviter, 
        .residualCheck = bool(options.nr_residualcheck),  
        .dampingFactor = options.nr_damping, 
        .dampingStep = options.nr_dampingstep, 
        .dampingSteps = options.nr_dampingsteps, 
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
            auto solPtr = circuit.retrieveDcSolution(solutionName);
            if (!solPtr) {
                // No nodesets
                nrSolver.forces(1).clear();
                Simulator::wrn() << "Warning, solution '"+solutionName+"' not found. No user nodesets applied.\n";
            } else {
                // Nodesets from solution repository
                if (!nrSolver.forces(1).set(circuit, *solPtr, strictforce)) {
                    // Abort if strictforce is set
                    if (strictforce) {
                        nrSolver.forces(1).formatError(s);
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
        if (!nrSolver.forces(1).set(circuit, preprocessedNodeset, false, strictforce)) {
            // Abort on error if strictforce is set
            if (strictforce) {
                nrSolver.forces(1).formatError(s);
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

bool OperatingPointCore::runSolver(bool continuePrevious) {
    auto& options = circuit.simulatorOptions().core();
    bool hasInitialStates;
    // Handle continuation
    if (continuePrevious) {
        // Continue mode
        if (continueState &&
            continueState->valid && continueState->coherent &&
            continueState->solution.values().size()==circuit.unknownCount()+1 &&
            continueState->stateVector.size()==circuit.statesCount() 
        ) {
            // Continue a state
            // State is valid, coherent, and its lengths match those of the solver vectors
            // Restore current state
            solution.vector() = continueState->solution.values();
            states.vector() = continueState->stateVector;
            hasInitialStates = true;
            // No nodesets applied
            nrSolver.enableForces(0, false);
            nrSolver.enableForces(1, false);
            if (options.op_debug>1) {
                Simulator::dbg() << "OP using ordinary continue mode with stored state.\n";
            }
        } else if (continueState && continueState->valid) {
            // Stored analysis state is not coherent with current circuit, 
            // its lengths may not match those of the solver vectors, 
            // but it is valid. 
            // Use nodesets to continue, but set no initial states vector nor initial solution. 
            // Ignore nodeset conflicts arising from stored solution. 
            // There should be no such conflicts as we are applying nodesets to nodes only, not node deltas. 
            nrSolver.forces(0).set(circuit, continueState->solution, false);
            nrSolver.enableForces(0, true);
            // Disable user-specified nodesets
            nrSolver.enableForces(1, false);
            if (options.op_debug>1) {
                Simulator::dbg() << "OP using nodeset continue mode with stored state.\n";
            }
        } else {
            // No valid state, continue with whatever is in solution and states vector
            hasInitialStates = true;
            // No nodesets applied
            nrSolver.enableForces(0, false);
            nrSolver.enableForces(1, false);
            if (options.op_debug>2) {
                Simulator::dbg() << "OP using ordinary continue mode with previous solution.\n";
            }
        }

        // State is spent after first use
        continueState = nullptr;
    } else {
        // No initial state given, start with standard initial point
        hasInitialStates = false;
        
        // Disable continuation nodesets in slot 0
        nrSolver.enableForces(0, false);

        // Apply nodesets specified by user in slot 1
        nrSolver.enableForces(1, true); 
        
        if (options.op_debug>1) {
            Simulator::dbg() << "OP using standard initial solution and state.\n";
        }
    }

    return nrSolver.run(hasInitialStates);
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
bool OperatingPointCore::run(bool continuePrevious, Status& s) {
    clearError();

    auto& options = circuit.simulatorOptions().core();
    auto& internals = circuit.simulatorInternals();
    converged = false;
    runType = RunType::OrdinaryOp;
    auto debug = options.op_debug;
    bool leave = false;
    bool tried = false;

    // Set time to 0
    internals.time = 0.0;

    // Make sure repositories are large enough
    auto n = circuit.unknownCount();
    solution.upsize(2, n+1);
    states.upsize(2, circuit.statesCount());
    
    // Initial plain op
    auto skipinitial = options.op_skipinitial;
    if (!skipinitial) {    
        tried = true;
        converged = runSolver(continuePrevious);
        if (!converged) {
            setError(OpError::InitialOp);
        }
        leave = circuit.checkFlags(Circuit::Flags::Abort);
        if (debug>0) {
            if (converged) {
                Simulator::dbg() << "OP core algorithm converged in " << std::to_string(nrSolver.iterations()) << " NR iteration(s).\n";
            } else {
                Simulator::dbg() << "OP core algorithm failed to converge in " << std::to_string(nrSolver.iterations()) << " NR iteration(s).\n";
            }
        }
    }

    auto skiphomotopy = options.op_skiphomotopy;

    // Try gmin stepping
    auto skipgmin = options.op_skipgmin;
    auto gminsteps = options.op_gminsteps;
    auto spice3gmin = options.op_spice3gmin;
    if (!converged && !leave && !skipgmin && !skiphomotopy && gminsteps>1) {
        tried = true;
        if (!spice3gmin) {
            // New algorithms
            // Device gmin first
            converged = gminStepping(RunType::GminStepping);
            leave = circuit.checkFlags(Circuit::Flags::Abort);
            if (debug>0) {
                Simulator::dbg() << "Gmin stepping " << (converged ? "succeeded" : "failed") << ".\n";
            }
            if (!converged && !leave && options.op_gshuntalg) {
                // Diagonal gshunt second
                converged = gminStepping(RunType::GshuntStepping);
                leave = circuit.checkFlags(Circuit::Flags::Abort);
                if (debug>0) {
                    Simulator::dbg() << "Gshunt stepping " << (converged ? "succeeded" : "failed") << ".\n";
                }
            }
        } else {
            // Spice3 gmin stepping
            converged = spice3GminStepping();
            leave = circuit.checkFlags(Circuit::Flags::Abort);
            if (debug>0) {
                Simulator::dbg() << "Spice3 Gmin stepping " << (converged ? "succeeded" : "failed") << ".\n";
            }
        }
    }

    // Try source stepping
    auto skipsrc = options.op_skipsrc;
    auto srcsteps = options.op_srcsteps;
    auto spice3src = options.op_spice3gmin;
    if (!converged && !leave && !skipsrc && !skiphomotopy && srcsteps>1) { 
        tried = true;
        if (!spice3src) {
            converged = sourceStepping();
            leave = circuit.checkFlags(Circuit::Flags::Abort);
            if (debug>0) {
                Simulator::dbg() << "Source stepping " << (converged ? "succeeded" : "failed") << ".\n";
            }
        } else {
            converged = spice3SourceStepping();
            leave = circuit.checkFlags(Circuit::Flags::Abort);
            if (debug>0) {
                Simulator::dbg() << "Spice3 source stepping " << (converged ? "succeeded" : "failed") << ".\n";
            }
        }
    }

    if (!leave) {
        // Did not leave early
        if (!tried) {
            // No algorithm tried
            setError(OpError::NoAlgorithm);
        } else if (converged) {
            // Tried and converged, write results
            if (outfile && params.writeOutput) {
                outfile->addPoint();
            }
        }
    } else {
        // Leaving early, did not converge
        // Add a status message one level higher
        converged = false;
    }

    return converged;
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
        case OpError::InitialOp:
            s.extend("Initial OP analysis failed.");
            return false;
        case OpError::SteppingSolver:
        case OpError::SteppingSteps:
            if (lastOpError==OpError::SteppingSteps) {
                s.set(Status::Analysis, "Homotopy reached step limit.");
            }
            switch (errorRunType) {
                case RunType::GminStepping:
                    s.extend(
                        homotopyProgress()+", dynamic "+(errorRunType==RunType::GminStepping ? "gmin" : "gsunt")+
                        " stepping failed after "+std::to_string(errorHomotopyIterations)+" step(s)."
                    );
                    break;
                case RunType::Spice3GminStepping:
                    s.extend(
                        homotopyProgress()+", SPICE3 gmin stepping failed after "+
                        std::to_string(errorHomotopyIterations)+" step(s)."
                    );
                    break;
                case RunType::SourceStepping:
                    s.extend(
                        homotopyProgress()+", dynamic source stepping failed after "+
                        std::to_string(errorHomotopyIterations)+" step(s)."
                    );
                    break;
                case RunType::Spice3SourceStepping:
                    s.extend(
                        homotopyProgress()+", SPICE3 source stepping failed after "+
                        std::to_string(errorHomotopyIterations)+" step(s)."
                    );
                    break;
            }
            return false;
        case OpError::NoAlgorithm:
            s.set(Status::Analysis, "No operating point algorithm tried."); 
            return false;
    }
    return solverError;
}

void OperatingPointCore::dump(std::ostream& os) const {
    AnalysisCore::dump(os);
    os << "  Solver ";
    if (converged) {
        os << "converged";
    } else {
        os << "not converged";
    }
    os << ", last NR run completed with " << std::to_string(nrSolver.iterations()) << " iterations";
    os << std::endl;
    os << "  Unknowns vector: " << std::endl;
    circuit.dumpSolution(os, solution.data(), "    ");
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
  - icMode = uic
    - start mode and continue mode
      ignore nodeset parameter
      use ic parameter to build initial solution vector
  - icMode = op
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
