#include "coreopnr.h"
#include "simulator.h"
#include "common.h"
#include <iomanip>

namespace NAMESPACE {

// Equations and limiting
//
// SPICE formulation
//   G    .. Jacobian
//   x    .. unknowns
//   f(x) .. residual
// 
//   G_i x_{i+1} = - f(x_i) + G_i x_i
// 
// Delta formulation
//   G_i (x_{i+1} - x_i) = - f(x_i)
// 
// SPICE formulation with limiting
//   xl   .. limited unknowns
//   Gl   .. Jacobian at limited unknowns (G(xl))
//   
//   Gl_i x_{i+1} = -f(xl_i) + Gl_i xl_i
//
//   Gl_i x_{i+1} = -f(xl_i) + Gl_i xl_i + Gl_i x_i - Gl_i x_i
//
//   Gl_i (x_{i+1} - x_i) = -f(xl_i) + Gl_i (xl_i - x_i) 
// 
// In the spirit of nrsolver.cpp we rewrite this as
//
//   Gl_i (x_i - x_{i+1}) = f(xl_i) - Gl_i (xl_i - x_i) 
//                          -------   -----------------
//                          residual  RHS linearized residual contribution (stored by OSDI models)
// 
// load_residual_resist()  ... loads residual (computed at the point of limiting)
// load_limit_rhs_resist() ... loads negated RHS linearized residual contribution
// 
// We first load residuals by calling load_residual_resist(), 
// then subtract RHS limited residual by calling load_limit_rhs_resist(). 
//
// The above equation can also be interpreted as
//   Gl_i (x_i - x_{i+1}) = f(xl_i) + Gl_i (x_i - xl_i)
//                          -----------------------------
//                          residual linearized at xl_i and computed at x_i

// Slots 0 (current) and -1 (future) are used for NR solver
// Slots 1, 2, ... correspond to past values (at t_{k}, t_{k-1}, ...)
// Therefore historyOffset needs to be set to 1

std::tuple<bool, bool> PreprocessedUserForces::set(Circuit& circuit, ValueVector& userForces, Status& s) {
    clear();

    // 0 -> 1 -> 2
    // 
    // 0 = have nothing, 
    // 1 = have node1, 
    // 2 = have node2, 
    // 3 = have number for single node (end)
    // 4 = have number for node pair (end)
    int state = 0; 
    Node* node1;
    Node* node2;
    Id id1, id2;
    double value; 
    size_t nsNdx = 0;
    bool haveAllEntries = true;
    for(auto& it : userForces) {
        switch (state) {
            case 0:
                node1 = node2 = nullptr;
                if (it.type()!=Value::Type::String) {
                    s.set(Status::BadArguments, "Expecting a string at position "+std::to_string(nsNdx)+".");
                    return std::make_tuple(false, false);
                }
                id1 = it.val<String>();
                node1 = circuit.findNode(id1);
                // Node not found is an error
                // if (!node1) {
                //     s.set(Status::BadArguments, "Cannot find node '"+std::string(id1)+"' during force preprocessing.");
                //     return std::make_tuple(false, false);
                // }
                state = 1;
                break;
            case 1:
                switch (it.type()) {
                    case Value::Type::Int:
                        value = it.val<Int>();
                        state = 3;
                        break;
                    case Value::Type::Real:
                        value = it.val<Real>();
                        state = 3;
                        break;
                    case Value::Type::String:
                        id2 = it.val<String>();
                        node2 = circuit.findNode(id2);
                        // Node not found is an error
                        // if (!node2) {
                        //     s.set(Status::BadArguments, "Cannot find node '"+std::string(id2)+"' during force preprocessing.");
                        //     return std::make_tuple(false, false);
                        // }
                        state = 2;
                        break;
                    default:
                        s.set(Status::BadArguments, "Expecting a string, an integer, or a real at position "+std::to_string(nsNdx)+".");
                        return std::make_tuple(false, false);
                }
                break;
            case 2:
                switch (it.type()) {
                    case Value::Type::Int:
                        value = it.val<Int>();
                        state = 4;
                        break;
                    case Value::Type::Real:
                        value = it.val<Real>();
                        state = 4;
                        break;
                    default:
                        s.set(Status::BadArguments, "Expecting an integer or a real at position "+std::to_string(nsNdx)+".");
                        return std::make_tuple(false, false);
                }
                break;
        }
        if (state==3) {
            // Have single node, ignore force if node is not found
            if (node1) {
                nodes.push_back(node1);
                nodeIds.push_back(id1);
                nodeValues.push_back(value);
            }
            state = 0;
        } else if (state==4) {
            // Have node pair, ignore force if node is not found
            if (node1 && node2) {
                nodePairs.push_back(std::make_tuple(node1, node2)); 
                nodeIdPairs.push_back(std::make_tuple(id1, id2));
                nodePairValues.push_back(value); 
            }
            
            // Check existence of extradiagonal entries, but only if both nodes were found. 
            // No need to check if haveAllEntries is already false. 
            if (node1 && node2 && haveAllEntries) {
                auto u1 = node1->unknownIndex();
                auto u2 = node2->unknownIndex();
                auto entry12 = circuit.sparsityMap().find(MatrixEntryPosition(u1, u2));
                auto entry21 = circuit.sparsityMap().find(MatrixEntryPosition(u2, u1));
                haveAllEntries = haveAllEntries && entry12 && entry21;
            }
            state = 0;
        }
        nsNdx++;
    }
    
    return std::make_tuple(true, !haveAllEntries);
}

    
OpNRSolver::OpNRSolver(
    Circuit& circuit, CommonData& commons, KluRealMatrix& jac, 
    VectorRepository<double>& states, VectorRepository<double>& solution, 
    NRSettings& settings, Int forcesSize
) : circuit(circuit), commons(commons), states(states), 
    NRSolver(circuit.tables().accounting(), jac, solution, settings) {
    // Slot 0 is for sweep continuation and homotopy (set via CoreStateStorage object)
    // Slot 1 is 
    resizeForces(forcesSize);

    // For constructing the linearized system in NR loop
    evalSetup_ = EvalSetup {
        // Inputs
        .solution = &solution, 
        .states = &states, 

        // Signal this is static and DC analysis
        .staticAnalysis = true, 
        .dcAnalysis = true, 

        // Evaluation 
        .enableLimiting = true, 
        .evaluateResistiveJacobian = true, 
        .evaluateResistiveResidual = true, 
        .evaluateLinearizedResistiveRhsResidual = true, 
        .evaluateOutvars = true, 
    };

    loadSetup_ = LoadSetup {
        .states = &states, 
        .loadResistiveJacobian = true, 
    };
}

bool OpNRSolver::setForces(Int ndx, const AnnotatedSolution& solution, bool abortOnError) {
    // Get forces
    auto& f = forces(ndx);

    // Clear forced values
    f.clear();
    
    // Number of unknowns
    auto n = circuit.unknownCount();

    // Make space for variable forces (also include bucket)
    f.unknownValue_.resize(n+1);
    f.unknownForced_.resize(n+1);
    
    bool error = false;

    // Go through all solution components, excluding ground
    auto nSol = solution.values().size();
    // Ignore components that do not have a name
    nSol = std::min(nSol, solution.names().size());
    for(decltype(nSol) i=1; i<nSol; i++) {
        // Node
        auto name = solution.names()[i];
        auto value = solution.values()[i];
        Node* node = circuit.findNode(name);
        if (!node) {
            // Node not found
            continue;
        }

        if (!setForceOnUnknown(f, node, value)) {
            error = true;
            break;
        }
    }

    return !error;
}

bool OpNRSolver::setForces(Int ndx, const PreprocessedUserForces& preprocessed, bool uicMode, bool abortOnError) {
    // Get forces
    Forces& f = forces(ndx);

    // Clear forced values
    f.unknownValue_.clear();
    f.unknownForced_.clear();
    f.deltaValue_.clear();
    f.deltaIndices_.clear();

    // Number of unknowns
    auto n = circuit.unknownCount();

    // Make space for forces on unknowns, set them by default to 0
    f.unknownValue_.resize(n+1, 0.0);
    f.unknownForced_.resize(n+1, false);

    bool error = false;
    
    // Set forces on unknowns
    auto nNodeForces = preprocessed.nodes.size();
    for(decltype(nNodeForces) i=0; i<nNodeForces; i++) {
        // Check if node was found
        auto node = preprocessed.nodes[i];
        auto value = preprocessed.nodeValues[i];
        if (!setForceOnUnknown(f, node, value)) {
            error = true;
            if (abortOnError) {
                return false;
            }
        }
    }
    
    // Set delta forces of the form v(x,0) or v(0,x), check node pairs
    auto nDeltaForces = preprocessed.nodePairs.size();  
    for(decltype(nNodeForces) i=0; i<nDeltaForces; i++) {
        // Check if both nodes were found
        auto [node1, node2] = preprocessed.nodePairs[i];
        auto [id1, id2] = preprocessed.nodeIdPairs[i]; 
        
        // Get unknowns and value
        auto u1 = node1->unknownIndex();
        auto u2 = node2->unknownIndex();
        auto value = preprocessed.nodePairValues[i];

        // Check if both nodes are ground? 
        if (u1==0 && u2==0) {
            // If yes, ignore force
            continue;
        } else if (u1==0) {
            // Check if first node is ground, convert it to a force on an unknown
            // v(0,x) = value -> v(x)=-value
            if (!setForceOnUnknown(f, node2, -value)) {
                error = true;
                if (abortOnError) {
                    return false;
                }
            }
        } else if (u2==0) {
            // v(x,0) = value -> v(x)=value
            if (!setForceOnUnknown(f, node1, value)) {
                error = true;
                if (abortOnError) {
                    return false;
                }
            }
        }
    }

    // Set real delta forces
    for(decltype(nNodeForces) i=0; i<nDeltaForces; i++) {
        // Check if both nodes were found
        auto [node1, node2] = preprocessed.nodePairs[i];
        auto [id1, id2] = preprocessed.nodeIdPairs[i]; 
        
        // Get unknowns and value
        auto u1 = node1->unknownIndex();
        auto u2 = node2->unknownIndex();
        auto value = preprocessed.nodePairValues[i];

        // Check if both nodes are ground
        if (u1==0 && u2==0) {
            // Skip this node pair
            continue;
        } else if (u1==0) {
            // v(0,x) = value -> v(x)=-value, already handled
            continue;
        } else if (u2==0) {
            // v(x,0) = value -> v(x)=value, already handled
            continue;
        } else {
            // Actual delta force 
            // Is force already set on both nodes
            if (f.unknownForced_[u1] && f.unknownForced_[u2]) {
                // Both nodes are forced
                // Does delta force conflict with node forces
                if (f.unknownValue_[u1]-f.unknownValue_[u2]!=value) {
                    lastOpNRError = OpNRSolverError::ConflictDelta;
                    errorNode1 = node1;
                    errorNode2 = node2;
                    error = true;
                    if (abortOnError) {
                        return false;
                    }
                } else {
                    // Matches node forces, no need to add it, skip
                    continue;
                }
            } else {
                // At least one node is not forced yet
                if (uicMode) {
                    // As UIC forces, apply to nodes
                    // One of the nodes is not forced
                    if (f.unknownForced_[u1]) {
                        // Force node2 to node1-value
                        f.unknownValue_[u2] = f.unknownValue_[u1] - value;
                        f.unknownForced_[u2] = true;
                    } else if (f.unknownForced_[u2]) {
                        // Force node1 to node2+value
                        f.unknownValue_[u1] = f.unknownValue_[u2] + value;
                        f.unknownForced_[u1] = true;
                    } else {
                        // Force u2 to 0 and u1 to value
                        f.unknownValue_[u1] = value;
                        f.unknownValue_[u2] = 0;
                        f.unknownForced_[u1] = true;
                        f.unknownForced_[u2] = true;
                    }
                } else {
                    // As delta forces
                    f.deltaValue_.push_back(value);
                    f.deltaIndices_.push_back(std::make_tuple(u1, u2));
                }        
            }
        }
    }

    return true;
}

bool OpNRSolver::setForceOnUnknown(Forces& f, Node* node, double value) {
    // Unknown
    auto u = node->unknownIndex();
    // Is it a ground node? 
    if (u==0) {
        // If yes, ignore the force. 
        return true;
    }
    // Is it conflicting with a previous nodeset
    if (f.unknownForced_[u] && f.unknownValue_[u]!=value) {
        lastOpNRError = OpNRSolverError::ConflictNode;
        errorNode1 = node;
        return false;
    }
    f.unknownValue_[u] = value;
    f.unknownForced_[u] = true;

    return true;
}


bool OpNRSolver::rebuild() {
    // Call parent's rebuild
    if (!NRSolver::rebuild()) {
        // Assume parent has set the error flag
        return false;
    }

    // Allocate space in vetors
    auto n = circuit.unknownCount();
    dummyStates.resize(circuit.statesCount());
    maxResidualContribution_.resize(n+1);       
    historicMaxSolution_.resize(n+1);
    globalMaxSolution_.resize(commons.natures.count()); // Number of solution natures
    pointMaxSolution_.resize(commons.natures.count());  // Number of solution natures
    historicMaxResidualContribution_.resize(n+1);
    globalMaxResidualContribution_.resize(commons.natures.count()); // Number of residual natures
    pointMaxResidualContribution_.resize(commons.natures.count());  // Number of residual natures
    resetMaxima();

    // Get diagonal and extradiagonal pointers for forces
    diagPtrs.resize(n+1);
    
    // Bind diagonal matrix elements
    // Needed for forcing unknown values and setting gshunts
    for(decltype(n) i=0; i<n; i++) {
        // We know the matrix type so we can use the elementPtr() non-virtual function
        diagPtrs[i+1] = jac.elementPtr(MatrixEntryPosition(i+1, i+1), Component::Real);
    }

    // Bind extradiagonal matrix entries for forced deltas
    extraDiags.resize(forcesList.size());
    auto nForces = forcesList.size();
    for(decltype(nForces) iForce=0; iForce<nForces; iForce++) {
        auto& deltaIndices = forcesList[iForce].deltaIndices_; 
        auto nDelta = deltaIndices.size();
        auto& ptrs = extraDiags[iForce];
        ptrs.clear();
        for(decltype(nDelta) i=0; i<nDelta; i++) {
            auto [u1, u2] = deltaIndices[i];
            // We know the matrix type so we can use the elementPtr() non-virtual function
            ptrs.push_back(
                std::make_tuple(
                    jac.elementPtr(MatrixEntryPosition(u1, u2), Component::Real),
                    jac.elementPtr(MatrixEntryPosition(u2, u1), Component::Real)
                )
            );
        }
    }

    // Build flow node flags
    // TODO: do something with this... use natures instead of PotentialNode flag of a Node
    isFlow.resize(n+1);
    for(decltype(n) i=1; i<=n; i++) {
        auto* node = circuit.reprNode(i);
        isFlow[i] = node->maskedFlags(Node::Flags::NodeTypeMask)==Node::Flags::PotentialNode;
    }
    
    return true;
}

bool OpNRSolver::initialize(bool continuePrevious) {
    // This method is called once on entering run()
    // This is the right place to set up vectors

    // Clear OP NR solver error
    clearError();

    // Clear flags
    clearFlags();

    // If not in continue mode set current states to 0
    if (!continuePrevious) {
        // Zero states
        states.zero();
    }

    // If bypass is enabled, prepare space for previous device states
    // Need to do this here because the user might sweep nr_bypass, but
    // the minimum requirement for calling rebuild() is that mapping 
    // changes. But nr_bypass does not affect mapping. 
    if (circuit.simulatorOptions().core().nr_bypass) {
        deviceStates.resize(circuit.deviceStatesCount());
    }
    
    // Set vectors for building linear system
    bool computeMaxResidualContribution = settings.residualCheck;

    loadSetup_.resistiveResidual = delta.data();
    loadSetup_.linearizedResistiveRhsResidual = delta.data();
    loadSetup_.maxResistiveResidualContribution = computeMaxResidualContribution ? maxResidualContribution_.data() : nullptr;

    // Set up tolerance reference value for solution
    auto& options = circuit.simulatorOptions().core();
    if (options.relrefsol==SimulatorOptions::relrefPointLocal) {
        globalSolRef = false;
        historicSolRef = false;
    } else if (options.relrefsol==SimulatorOptions::relrefLocal) {
        globalSolRef = false;
        historicSolRef = true;
    } else if (options.relrefsol==SimulatorOptions::relrefPointGlobal) {
        globalSolRef = true;
        historicSolRef = false;
    } else if (options.relrefsol==SimulatorOptions::relrefGlobal) {
        globalSolRef = true;
        historicSolRef = true;
    } else if (options.relrefsol==SimulatorOptions::relrefRelref) {
        if (options.relref == SimulatorOptions::relrefAlllocal) {
            globalSolRef = false;
            historicSolRef = true;
        } else if (options.relref == SimulatorOptions::relrefSigglobal) {
            globalSolRef = true;
            historicSolRef = true;
        } else if (options.relref == SimulatorOptions::relrefAllglobal) {
            globalSolRef = true;
            historicSolRef = true;
        } else {
            lastError = Error::BadSolReference;
            return false;
        }
    } else {
        lastError = Error::BadSolReference;
        return false;
    }

    // Set up tolerance reference value for residual
    if (options.relrefres==SimulatorOptions::relrefPointLocal) {
        globalResRef = false;
        historicResRef = false;
    } else if (options.relrefres==SimulatorOptions::relrefLocal) {
        globalResRef = false;
        historicResRef = true;
    } else if (options.relrefres==SimulatorOptions::relrefPointGlobal) {
        globalResRef = true;
        historicResRef = false;
    } else if (options.relrefres==SimulatorOptions::relrefGlobal) {
        globalResRef = true;
        historicResRef = true;
    } else if (options.relrefres==SimulatorOptions::relrefRelref) {
        if (options.relref == SimulatorOptions::relrefAlllocal) {
            globalResRef = false;
            historicResRef = true;
        } else if (options.relref == SimulatorOptions::relrefSigglobal) {
            globalResRef = false;
            historicResRef = true;
        } else if (options.relref == SimulatorOptions::relrefAllglobal) {
            globalResRef = true;
            historicResRef = true;
        } else {
            lastError = Error::BadResReference;
            return false;
        }
    } else {
        lastError = Error::BadResReference;
        return false;
    }
    
    return true;
}

bool OpNRSolver::preIteration(bool continuePrevious) {
    // Clear maximal residual contribution
    zero(maxResidualContribution_);
    // Pass iteration number to Verilog-A models
    commons.iteration = iteration;
    return true;    
}

bool OpNRSolver::postSolve(bool continuePrevious) {
    auto& acct = circuit.tables().accounting();
    acct.acctNew.bpinst += evalSetup_.bypassableInstances;
    acct.acctNew.bpopport += evalSetup_.bypassOpportunuties;
    acct.acctNew.bpbypassed += evalSetup_.bypassedInstances;
    
    return true;
}

bool OpNRSolver::postConvergenceCheck(bool continuePrevious) {
    // If algorithm converged we are going to exit next and states must be 
    // rotated because the new state belongs to the current solution
    // We also rotate states if no convergence yet
    // beacuse we must prepare for next iteration. 
    states.rotate();

    // Print debug information on convergence
    if (settings.debug) {
        std::stringstream ss;
        ss << std::scientific << std::setprecision(2);
        Simulator::dbg() << "Iteration " << std::to_string(iteration) << (preventedConvergence ? ", convergence not allowed" : "");
        if (!preventedConvergence) {
            Simulator::dbg() << (iterationConverged ? ", converged" : "");
            if (settings.residualCheck) {
                ss.str(""); ss << maxResidual;
                Simulator::dbg() << ", worst residual=" << ss.str() << " @ " << (maxResidualNode ? maxResidualNode->name() : "(unknown)");
            }
            if (iteration>1) {
                ss.str(""); ss << maxDelta;
                Simulator::dbg() << ", worst delta=" << ss.str() << " @ " << (maxDeltaNode ? maxDeltaNode->name() : "(unknown)");
            }
        }
        Simulator::dbg() << "\n";
    }
    return true;
}

bool OpNRSolver::postIteration(bool continuePrevious) {
    return true;
}

void OpNRSolver::loadShunts(double gshunt, bool loadJacobian) {
    // Now load gshunt if it is greater than 0.0
    // Gshunt current (and its residual contribution) is
    //   gshunt * x
    double* xprev = solution.data();
    auto nUnknowns = circuit.unknownCount();
    for(decltype(nUnknowns) i=1; i<=nUnknowns; i++) {
        auto ptr = diagPtrs[i];
        if (!isFlow[i]) {
            // Jacobian
            if (loadJacobian) {
                *ptr += gshunt;
            }
            // Residual
            delta[i] += gshunt*xprev[i];
        }
    }
}

bool OpNRSolver::evalAndLoadWrapper(EvalSetup& evalSetup, LoadSetup& loadSetup) {
    lastError = Error::OK;
    evalSetup.requestHighPrecision = highPrecision;
    if (!circuit.evalAndLoad(commons, &evalSetup, &loadSetup, nullptr)) {
        // Load error
        lastError = Error::EvalAndLoad;
        if (settings.debug>2) {
            Simulator::dbg() << "Evaluation error.\n";
        }
        return false;
    }

    // Store Abort, Finish, and Stop flag
    if (evalSetup.requests.abort) {
        setFlags(Flags::Abort);
    }
    if (evalSetup.requests.finish) {
        setFlags(Flags::Finish);
    }
    if (evalSetup.requests.stop) {
        setFlags(Flags::Stop);
    }
    
    // Handle abort right now, finish and stop are handled outside NR loop
    if (checkFlags(Flags::Abort)) {
        if (settings.debug>2) {
            Simulator::dbg() << "Abort requested during evaluation.\n";
        }
        return false;
    }

    return true;
}

void OpNRSolver::setNodesetAndIcFlags(bool continuePrevious) {
    auto nsiter = circuit.simulatorOptions().core().op_nsiter;

    // Set nodesetEnabled flag in esSystem
    // Forces slot 0 (continuation nodesets) and 1 (user nodesets) 
    // are used in ordinary OP analysis
    // Nodesets are enabled in iterations 1..op_nsiter 
    // if continuePrevious is false. 
    // They are also enabled if slot 1 (user nodesets) is active. 
    evalSetup_.nodesetEnabled = (iteration<=nsiter) && (continuePrevious==false);
    
    // Set icEnabled flag in esSystem
    // Slot 2 holds permanent forces for computing initial conditions
    // when OP analysis is invoked from tran core. 
    // Whenever this slot is active transient forces are enables 
    // and we are applying initial conditions. 
    evalSetup_.icEnabled = forcesEnabled.size()>2 && forcesEnabled[2];
}

std::tuple<bool, bool> OpNRSolver::buildSystem(bool continuePrevious) {
    lastError = Error::OK;

    auto n = circuit.unknownCount();

    // Remove forces originating from nodesets after nsiter iterations
    auto nsiter = circuit.simulatorOptions().core().op_nsiter;
    // Do this only at nsiter+1 (first iteration has index 1)
    if (iteration==nsiter+1) {
        // Continuation nodesets
        enableForces(0, false);
        // User-specified nodesets
        enableForces(1, false);
    }

    // Set nodeset and IC flags
    setNodesetAndIcFlags(continuePrevious); 
    
    // Raw arrays
    double* xprev = solution.data();
    
    // Inject value for limiting debugging purposes
    // if (iteration==1) {
    //     xprev[1]=1.0;
    // }

    // Init limits if not in continue mode and iteration is 1
    evalSetup_.initializeLimiting = !continuePrevious && (iteration==1);
    
    // Inactive element bypass enabled, take it into account at evaluation time
    if (circuit.simulatorOptions().core().nr_bypass) {
        // In continue mode bypass is allowed in first iteration
        // Otherwise we allow it starting with the second one 
        if (continuePrevious || iteration>1) {
            evalSetup_.allowBypass = true;
        } else {
            evalSetup_.allowBypass = false;
        }
        // For bypass check
        evalSetup_.deviceStates = deviceStates.data();
    }

    // Force instance evaluation bypass if requested
    evalSetup_.forceBypass = commons.requestForcedBypass;
    
    // Clear bypass forcing request
    commons.requestForcedBypass = false;

    // Evaluate and load
    auto evalSt = evalAndLoadWrapper(evalSetup_, loadSetup_);
    if (!evalSt) {
        lastError = Error::EvalAndLoad;
        errorIteration = iteration;
        return std::make_tuple(false, evalSetup_.limitingApplied);
    }
    delta[0] = 0.0;

    // Now load gshunt if it is greater than 0.0
    auto gshunt = commons.gshunt;
    if (gshunt>0) {
        loadShunts(gshunt);
    }

    // Add forced values to the system
    if (haveForces() && !loadForces(true)) {
        if (settings.debug) {
            Simulator::dbg() << "Failed to load forced values at iteration " << iteration << "\n";
        }
        lastOpNRError = OpNRSolverError::LoadForces;
        errorIteration = iteration;
        std::make_tuple(false, evalSetup_.limitingApplied);
    }

    // Prevent convergence if limiting was applied
    return std::make_tuple(true, evalSetup_.limitingApplied); 
}

bool OpNRSolver::loadForces(bool loadJacobian) {
    // Are any forces enabled? 
    auto nf = forcesList.size();
    
    // Get row norms
    jac.rowMaxNorm(dataWithoutBucket(rowNorm));

    // Load forces
    auto n = jac.nRow();
    double* xprev = solution.data();
    for(decltype(nf) iForce=0; iForce<nf; iForce++) {
        // Skip disabled force lists
        if (!forcesEnabled[iForce]) {
            continue;
        }

        // First, handle forced unknowns
        auto& enabled = forcesList[iForce].unknownForced_;
        auto& force = forcesList[iForce].unknownValue_;
        auto nForceNodes = force.size();
        // Load only if the number of forced unknowns matches 
        // the number of unknowns in the circuit including ground
        if (nForceNodes==n+1) {
            for(decltype(nForceNodes) i=1; i<=n; i++) {
                if (enabled[i]) {
                    double factor = rowNorm[i]*settings.forceFactor;
                    if (factor==0.0) {
                        factor = 1.0;
                    }
                    // Jacobian entry: factor
                    // Residual: factor * x_i - factor * nodeset_i
                    auto ptr = diagPtrs[i];
                    if (ptr) {
                        // Jacobian
                        if (loadJacobian) {
                            *ptr += factor;
                        }
                        // Residual
                        delta[i] += factor * xprev[i] - factor * force[i];
                    }
                }
            }
        }

        // Second, handle forced deltas
        auto& extraDiagPtrs = extraDiags[iForce]; 
        auto& deltas = forcesList[iForce].deltaValue_;
        auto nDeltas = deltas.size();
        auto& uPairs = forcesList[iForce].deltaIndices_;
        // Load only if number of extradiagonal pointer pairs matches
        // the number of forced deltas
        if (extraDiagPtrs.size()==nDeltas) {
            for(decltype(nDeltas) i=0; i<nDeltas; i++) {
                auto [u1, u2] = uPairs[i];
                auto [extraDiagPtr1, extraDiagPtr2] = extraDiagPtrs[i];

                double factor1 = rowNorm[u1]*settings.forceFactor;
                double factor2 = rowNorm[u2]*settings.forceFactor;

                double contrib1 = factor1 * (xprev[u1] - xprev[u2]) - factor1 * deltas[i]; 
                double contrib2 = factor2 * (xprev[u2] - xprev[u1]) + factor2 * deltas[i]; 

                // Jacobian entry: 
                //         u1        u2
                //   u1    factor1  -factor1
                //   u2   -factor2   factor2
                // 
                // Residual at KCL u1: factor1 * (u1-u2) - factor1 * nodeset
                // Residual at KCL u2: factor2 * (u2-u1) + factor2 * nodeset
                *(diagPtrs[u1]) += factor1;
                *extraDiagPtr1 += -factor1;

                *(diagPtrs[u2]) += factor2;
                *extraDiagPtr2 += -factor2;
                
                delta[u1] += contrib1;
                delta[u2] += contrib2;
            }
        }
    }

    return true;
}

std::tuple<bool, bool> OpNRSolver::checkResidual() {
    // Options
    auto& options = circuit.simulatorOptions().core();

    // Compute norms only in debug mode
    bool computeNorms = settings.debug;

    // In residual we have the residual at previous solution
    // We are going to check that residual
    
    // Number of unknowns (vector length includes a bucket at index 0)
    auto n = circuit.unknownCount();

    // Results
    maxResidual = 0.0;
    maxNormResidual = 0.0;
    l2normResidual2 = 0.0;
    maxResidualNode = nullptr;
    
    // Assume residual is OK
    bool residualOk = true;
    
    // Get point maximum for each residual nature
    zero(pointMaxResidualContribution_); 
    for(decltype(n) i=1; i<=n; i++) {
        double c = std::fabs(maxResidualContribution_[i]);
        // Get residual nature index
        auto ndx = commons.residual_natureIndex[i];
        if (c>pointMaxResidualContribution_[ndx]) {
            pointMaxResidualContribution_[ndx] = c;
        }
    }
    
    // Go through all variables (except ground)
    for(decltype(n) i=1; i<=n; i++) {
        // Representative node, associated flow nature index
        auto rn = circuit.reprNode(i);
        // Skip internal device nodes
        if (rn->checkFlags(Node::Flags::InternalDeviceNode)) {
            continue;
        }
        // Get residual nature index
        auto ndx = commons.residual_natureIndex[i];
        
        // Compute tolerance reference
        // Point local reference by default
        // Compute tolerance reference, start with previous value of the i-th unknown
        double tolref = std::fabs(maxResidualContribution_[i]);

        // Account for global and historic references
        if (historicResRef) {
            if (globalResRef) {
                // Historic global reference, ndx is the nature index
                tolref = std::max(tolref, globalMaxResidualContribution_[ndx]);
            } else {
                // Historic local reference, i is the index of unknown
                tolref = std::max(tolref, historicMaxResidualContribution_[i]);
            }
        } else if (globalResRef) {
            // Point global reference, ndx is the nature index
            tolref = std::max(tolref, pointMaxResidualContribution_[ndx]);
        }

        // Residual tolerance (Designer's Guide to Spice and Spectre, chapter 2.2.2)
        auto tol = std::max(std::fabs(tolref*options.reltol), commons.residual_abstol[i]);

        // TODO: internal nodes when NR converges have a very low residual contribution
        //       because a single OSDI instance provides all the residual contribution to them 
        //       and the contribution is close to 0. 
        //       This means that absolute tolerances may be too low and prevent convergence forever. 
        //       This is somehow remedied if relrefres is set to pointglobal or global. 

        // Residual component
        double rescomp = fabs(delta[i]);
    
        // Normalized residual component
        double normResidual = rescomp/tol;

        if (computeNorms) {
            l2normResidual2 += normResidual*normResidual;
            // Update largest normalized component
            if (i==1 || normResidual>maxNormResidual) {
                maxResidual = rescomp;
                maxNormResidual = normResidual;
                maxResidualNode = rn;
            }
        }

        // See if residual component exceeds tolerance
        if (rescomp>tol) {
            residualOk = false;
            // Can exit if not computing norms
            if (!computeNorms) {
                return std::make_tuple(true, residualOk); 
            }
        }
    }
    
    return std::make_tuple(true, residualOk); 
}

std::tuple<bool, bool> OpNRSolver::checkDelta() {
    // Options
    auto& options = circuit.simulatorOptions().core();

    // Compute norms only in debug mode
    bool computeNorms = settings.debug;

    // In delta we have the solution change
    // Check it for convergence
    
    // Number of unknowns (vector length includes a bucket at index 0)
    auto n = circuit.unknownCount();

    // Raw arrays
    double* xprev = solution.data();
    
    maxDelta = 0.0;
    maxNormDelta = 0.0;
    maxDeltaNode = nullptr;
    
    // Check convergence (see if delta is small enough), 
    // but only if this is iteration 2 or later
    // In iteration 1 assume we did not converge
    
    // Assume we converged
    bool deltaOk = true;
    
    double* xdelta = delta.data();

    // Get point maximum for each solution nature
    zero(pointMaxSolution_);
    auto xold = solution.data();
    for(decltype(n) i=1; i<=n; i++) {
        double c = std::fabs(xold[i]);
        // Get unknown nature index
        auto ndx = commons.unknown_natureIndex[i];
        if (c>pointMaxSolution_[ndx]) {
            pointMaxSolution_[ndx] = c;
        }
    }

    // Use 1-based index (with bucket) because same indexing is used for variables
    for(decltype(n) i=1; i<=n; i++) {
        // Get unknown nature index
        auto ndx = commons.unknown_natureIndex[i];

        // Compute tolerance reference
        // Point local reference by default
        // Compute tolerance reference, start with previous value of the i-th unknown
        double tolref = std::fabs(xprev[i]);
        
        // Account for global and historic references
        if (historicSolRef) {
            if (globalSolRef) {
                // Historic global reference, ndx is the nature index
                tolref = std::max(tolref, globalMaxSolution_[ndx]);
            } else {
                // Historic local reference, i is the index of unknown
                tolref = std::max(tolref, historicMaxSolution_[i]);
            }
        } else if (globalSolRef) {
            // Point global reference, ndx is the nature index
            tolref = std::max(tolref, pointMaxSolution_[ndx]);
        }
        
        // Compute tolerance
        double tol = std::max(std::fabs(tolref*options.reltol), commons.unknown_abstol[i]);

        // Absolute solution change 
        double deltaAbs = fabs(delta[i]);
        
        if (computeNorms) {
            double normDelta = deltaAbs/tol;
            if (i==1 || normDelta>maxNormDelta) {
                maxDelta = deltaAbs;
                maxNormDelta = normDelta;
                maxDeltaNode = circuit.reprNode(i);
            }
        }

        // Check tolerance
        if (deltaAbs>tol) {
            // Did not converge
            deltaOk = false;
            
            // Can exit if not computing norms
            if (!computeNorms) {
                return std::make_tuple(true, deltaOk);
            }
        }
    }
    
    return std::make_tuple(true, deltaOk);
}

void OpNRSolver::resetMaxima() {
    zero(historicMaxSolution_);
    zero(globalMaxSolution_);
    zero(pointMaxSolution_); 
    zero(historicMaxResidualContribution_);
    zero(globalMaxResidualContribution_);
    zero(pointMaxResidualContribution_); 
}  

void OpNRSolver::initializeMaxima(OpNRSolver& other) {
    historicMaxSolution_ = other.historicMaxSolution_;
    globalMaxSolution_ = other.globalMaxSolution_;
    historicMaxResidualContribution_ = other.historicMaxResidualContribution_;
    globalMaxResidualContribution_ = other.globalMaxResidualContribution_;
}

void OpNRSolver::updateMaxima() {
    auto n = circuit.unknownCount();
    auto* x = solution.data();
    auto* mrc = maxResidualContribution_.data();
    for(decltype(n) i=1; i<=n; i++) {
        double c;
        UnknownIndex ndx;
        
        // Historic local and historic global maximal solution
        // Voltage nodes -> potential nature is voltage (0) 
        // Flow nodes -> potential nature is current (1) 
        // Unknown nature index
        ndx = commons.unknown_natureIndex[i];
        c = std::fabs(x[i]);
        if (c>historicMaxSolution_[i]) {
            historicMaxSolution_[i] = c;
        }
        if (c>globalMaxSolution_[ndx]) {
            globalMaxSolution_[ndx] = c;
        }
        
        // Historic local and historic global maximal residual contribution
        // Voltage nodes -> flow nature is current (1) 
        // Flow nodes -> flow nature is voltage (0) 
        // Residual nature index
        ndx = commons.residual_natureIndex[i];
        c = std::fabs(mrc[i]);
        if (c>historicMaxResidualContribution_[i]) {
            historicMaxResidualContribution_[i] = c;
        }
        if (c>globalMaxResidualContribution_[ndx]) {
            globalMaxResidualContribution_[ndx] = c;
        }
    }
}

bool OpNRSolver::formatError(Status& s, NameResolver* resolver) const {
    // Error in NRSolver
    if (lastError!=NRSolver::Error::OK) {
        NRSolver::formatError(s, resolver);
        return false;
    }

    switch (lastOpNRError) {
        case OpNRSolverError::ConflictNode:
            s.set(Status::Force, "Conflicting forces for node '"+std::string(errorNode1->name())+"'.");
            return false;
        case OpNRSolverError::ConflictDelta:
            s.set(Status::Force, "Forcing delta on node pair ('"
                        +std::string(errorNode1->name())+"', '"
                        +std::string(errorNode2->name())
                        +"') conflicts previous forces."
                    );
            return false;
        case OpNRSolverError::LoadForces:
            s.set(Status::Force, "Failed to load forces.");
            return false;
        default:
            return true;
    }
}

void OpNRSolver::dumpSolution(std::ostream& os, double* solution, const char* prefix) {
    auto n = circuit.unknownCount();
    std::ios original_state(nullptr);
    original_state.copyfmt(os);
    os << std::setprecision(15);
    for(decltype(n) i=1; i<=n; i++) {
        auto rn = circuit.reprNode(i);
        os << prefix << rn->name() << " : " << solution[i] << "\n";
    }
    os.copyfmt(original_state);
}

}
