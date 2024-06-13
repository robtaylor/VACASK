#include "nrsolver.h"
#include "simulator.h"
#include "common.h"
#include <iomanip>

namespace NAMESPACE {

// We are solving a nonlinear system of equations of the form
// 
//   f(x) = 0
// 
// where f(x) is a vector-valued function (R^n->R^n) of vector x. 
// f(x) is also referred to as the residual. 
// The matrix of derivatives of residual wrt. components of x 
// is teh Jacobian of f, also denoted by J. 
// 
// Suppose we have an approximate solution x_i. 
// x_{i+1} is found by solving
//
//   J(x_i) (x_{i+1} - x_i) = -f(x_i)
//
// We reorganize this into 
// 
//   J(x_i) (dnx_i) = f(x_i)
//
// where dnx_i = x_i - x_{i+1}. Now we can so,ve for dnx_i and 
// compute x_{i+1} as
//
//   x_{i+1} = x_i - dnx_i
//
// In this way we can directly use the Jacobian and residual loader 
// functions without changing the sign of the RHS. 

// NRSolver force slots that contain forced deltas (branch forces) 
// need to be populated before bebuild is called. 
// Slots containing branch forces affect the circuit topology. 
// This needs to be taken into account in Analysis::setSweepState() and 
// Analysis::setAnalysisOptions() before rebuild() is called. 
// Slots can be activated/deactivated,. 
NRSolver::NRSolver(
    Circuit& circuit, KluRealMatrix& jac, 
    VectorRepository<double>& states, VectorRepository<double>& solution, 
    NRSettings& settings
) : circuit(circuit), jac(jac), states(states), solution(solution), settings(settings), 
    iteration(0) {
}

bool NRSolver::rebuild() {
    // Allocate space in vectors
    auto n = circuit.unknownCount();
    diagPtrs.resize(n+1);
    isFlow.resize(n+1);
    delta.resize(n+1);
    rowNorm.resize(n+1);
    historicMaxSolution_.resize(n+1);
    historicMaxResidualContribution_.resize(n+1);
    resetMaxima();
    
    // Bind diagonal matrix elements
    // Needed for forcing unknown values and setting gshunts
    for(decltype(n) i=1; i<=n; i++) {
        auto* node = circuit.reprNode(i);
        diagPtrs[i] = jac.elementPtr(i, i, Component::RealPart);
        isFlow[i] = node->maskedFlags(Node::Flags::NodeTypeMask)==Node::Flags::PotentialNode;
    }

    // Bind extradiagonal matrix entries for forced deltas
    extraDiags.resize(forcesList.size());
    auto nForces = forcesList.size();
    for(decltype(nForces) iForce=0; iForce<nForces; iForce++) {
        auto& deltaIndices = forcesList[iForce].deltaIndices(); 
        auto nDelta = deltaIndices.size();
        auto& ptrs = extraDiags[iForce];
        ptrs.clear();
        for(decltype(nDelta) i=0; i<nDelta; i++) {
            auto [u1, u2] = deltaIndices[i];
            ptrs.push_back(
                std::make_tuple(
                    jac.elementPtr(u1, u2, Component::RealPart),
                    jac.elementPtr(u2, u1, Component::RealPart)
                )
            );
        }
    }

    return true;
}

bool NRSolver::loadForces(bool loadJacobian) {
    // Are any forces enabled? 
    auto nf = forcesList.size();
    
    // Get row norms
    jac.rowMaxNorm(dataWithoutBucket(rowNorm));

    // Load forces
    auto n = circuit.unknownCount();
    double* xprev = solution.data();
    for(decltype(nf) iForce=0; iForce<nf; iForce++) {
        // Skip disabled force lists
        if (!forcesEnabled[iForce]) {
            continue;
        }

        // First, handle forced unknowns
        auto& enabled = forcesList[iForce].unknownForced();
        auto& force = forcesList[iForce].unknownValue();
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
                    // Jacobian entry: -factor
                    // Residual: - factor * x_i + factor * nodeset_i
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
        auto& deltas = forcesList[iForce].deltaValue();
        auto nDeltas = deltas.size();
        auto& uPairs = forcesList[iForce].deltaIndices();
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

void stateDump(Vector<double>& states) {
    auto n=states.size();
    for(decltype(n) i=0; i<n; i++) {
        std::cout << states[i] << " ";
    }
}

void NRSolver::resizeForces(Int n) {
    forcesList.resize(n);
    forcesEnabled.resize(n, false);
}

Forces& NRSolver::forces(Int ndx) {
    return forcesList.at(ndx);
} 

void NRSolver::resetMaxima() {
    zero(historicMaxSolution_);
    zero(historicMaxResidualContribution_);
    globalMaxSolution_ = 0;
    globalMaxResidualContribution_ = 0;
    pointMaxSolution_ = 0;
    pointMaxResidualContribution_ = 0;
}  

void NRSolver::initializeMaxima(NRSolver& other) {
    historicMaxSolution_ = other.historicMaxSolution();
    globalMaxSolution_ = other. globalMaxSolution();
    historicMaxResidualContribution_ = other.historicMaxResidualContribution();
    globalMaxResidualContribution_ = other.globalMaxResidualContribution();
}

void NRSolver::updateMaxima() {
    auto n = circuit.unknownCount();
    auto* x = solution.data();
    auto* mrc = maxResidualContribution_.data();
    for(decltype(n) i=1; i<=n; i++) {
        double c;
        c = std::fabs(x[i]);
        if (c>historicMaxSolution_[i]) {
            historicMaxSolution_[i] = c;
        }
        if (c>globalMaxSolution_) {
            globalMaxSolution_ = c;
        }
        c = std::fabs(mrc[i]);
        if (c>historicMaxResidualContribution_[i]) {
            historicMaxResidualContribution_[i] = c;
        }
        if (c>globalMaxResidualContribution_) {
            globalMaxResidualContribution_ = c;
        }
    }
}

std::tuple<bool, double, double, Node*> NRSolver::checkDelta(bool* deltaOk, bool computeNorms) {
    // In delta we have the solution change
    // Check it for convergence
    
    // Number of unknowns (vector length includes a bucket at index 0)
    auto n = circuit.unknownCount();

    // Raw arrays
    double* xprev = solution.data();
    
    double maxDelta = 0.0;
    double maxNormDelta = 0.0;
    Node* maxDeltaNode = nullptr;
    
    // Check convergence (see if delta is small enough), 
    // but only if this is iteration 2 or later
    // In iteration 1 assume we did not converge
    
    // Assume we converged
    if (deltaOk) {
        *deltaOk = true;
    }

    double* xdelta = delta.data();

    // Get point maximum
    pointMaxSolution_ = 0;
    for(decltype(n) i=1; i<=n; i++) {
        double c = std::fabs(xprev[i]);
        if (c>pointMaxSolution_) {
            pointMaxSolution_ = c;
        }
    }

    // Use 1-based index (with bucket) because same indexing is used for variables
    for(decltype(n) i=1; i<=n; i++) {
        auto rn = circuit.reprNode(i);

        // Compute tolerance reference
        double tolref = std::fabs(xprev[i]);
        
        // Cannot account for new solution because damping has not been performed yet

        // Account for global and historic references
        if (settings.historicSolRef) {
            if (settings.globalSolRef) {
                tolref = std::max(tolref, globalMaxSolution_);
            } else {
                tolref = std::max(tolref, historicMaxSolution_[i]);
            }
        } else if (settings.globalSolRef) {
            tolref = std::max(tolref, pointMaxSolution_);
        }
        
        // Compute tolerance
        double tol = circuit.solutionTolerance(rn, tolref);

        // Absolute solution change 
        double deltaAbs = fabs(delta[i]);
        
        if (computeNorms) {
            double normDelta = deltaAbs/tol;
            if (i==1 || normDelta>maxNormDelta) {
                maxDelta = deltaAbs;
                maxNormDelta = normDelta;
                maxDeltaNode = rn;
            }
        }

        // Check tolerance
        if (deltaAbs>tol) {
            // Did not converge
            if (deltaOk) {
                *deltaOk = false;
            }
            // Can exit if not computing norms
            if (!computeNorms) {
                break;
            }
        }
    }
    
    return std::make_tuple(true, maxDelta, maxNormDelta, maxDeltaNode);
}

std::tuple<bool, double, double, double, Node*> NRSolver::checkResidual(bool* residualOk, bool computeNorms) {
    // In residual we have the residual at previous solution
    // We are going to check that residual
    
    // Number of unknowns (vector length includes a bucket at index 0)
    auto n = circuit.unknownCount();

    // Results
    double maxResidual = 0.0;
    double maxNormResidual = 0.0;
    double l2normResidual2 = 0.0;
    Node* maxResidualNode = nullptr;
    
    // Assume residual is OK
    if (residualOk) {
        *residualOk = true;
    }

    // Get point maximum
    pointMaxResidualContribution_ = 0;
    for(decltype(n) i=1; i<=n; i++) {
        double c = std::fabs(maxResidualContribution_[i]);
        if (c>pointMaxResidualContribution_) {
            pointMaxResidualContribution_ = c;
        }
    }
    
    // Go through all variables (except ground)
    for(decltype(n) i=1; i<=n; i++) {
        // Get representative node for i-th variable
        auto rn = circuit.reprNode(i);

        // Compute tolerance reference
        double tolref = std::fabs(maxResidualContribution_[i]);

        // Account for global and historic references
        if (settings.historicResRef) {
            if (settings.globalResRef) {
                tolref = std::max(tolref, globalMaxResidualContribution_);
            } else {
                tolref = std::max(tolref, historicMaxResidualContribution_[i]);
            }
        } else if (settings.globalResRef) {
            tolref = std::max(tolref, pointMaxResidualContribution_);
        }

        // Residual tolerance (Designer's Guide to Spice and Spectre, chapter 2.2.2)
        auto tol = circuit.residualTolerance(rn, tolref);

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
            if (residualOk) {
                *residualOk = false;
                // Can exit if not computing norms
                if (!computeNorms) {
                    break;
                }
            }
        }
    }
    
    return std::make_tuple(true, maxResidual, maxNormResidual, l2normResidual2, maxResidualNode); 
}

bool NRSolver::run(bool continuePrevious) {
    auto t0 = Accounting::wclk();
    circuit.tables().accounting().acctNew.nrcall++;

    jac.setAccounting(circuit.tables().accounting());

    // Clear error
    clearError();

    // Number of unknowns (vector length includes a bucket at index 0)
    auto n = solution.length()-1;
    
    // Iteration limit and damping
    int itlim;
    if (continuePrevious) {
        itlim = settings.itlimCont;
    } else {
        itlim = settings.itlim;
    }

    // If not in continue mode set current solution and state to 0
    if (!continuePrevious) {
        // Zero current solution and states
        solution.zero();
        states.zero();
    }

    if (settings.debug) {
        Simulator::dbg() << "Starting NR algorithm " << (continuePrevious ? "with given initial solution" : "with zero initial solution") << ".\n";
    }

    // Main loop
    bool deltaOk;
    bool residualOk;
    double maxResidual;
    double maxNormResidual;
    double l2normResidual2;
    Node* maxResidualNode; 
    double maxDelta;
    double maxNormDelta; 
    Node* maxDeltaNode;
    iteration = 0;
    Int convIter = 0;
    bool converged = false;

    clearError();

    // Initialize structures
    if (!initialize(continuePrevious)) {
        // Assume initialize() has set lastError
        circuit.tables().accounting().acctNew.tnr += Accounting::wclkDelta(t0);
        return false;
    }

    // Minimal damping factor, needed for adjusting the iteration limit
    double minimalDampingFactor = 1.0;

    bool exitNrLoop = false;
    do {
        // Simulator::dbg() << "NR it=" << iteration << " : at=" << solution.position() << "/" << solution.size()
        //     << " state at=" << states.position() << "/" << states.size()
        //     << " current data ptr=" << size_t(solution.data()) << "\n";

        // Assume no convergence
        bool iterationConverged = false;

        // Iteration counter (first iteration is 1)
        iteration++;

        // Pass iteration number to Verilog-A models
        circuit.simulatorInternals().iteration = iteration;
        
        // Zero matrix, new solution, new states, and residual/delta
        jac.zero();
        solution.zeroFuture();
        states.zeroFuture();
        zero(delta);
        
        double* xprev = solution.data();
        double* xnew = solution.futureData();
        double* xdelta = delta.data();

        if (settings.debug>=4) {
            std::cout << "Old solution at iteration " << iteration << "\n";
            circuit.dumpSolution(std::cout, solution.data(), "  ");
            std::cout << "\n";
        }
        // std::cout << "New solution at iteration " << iteration << "\n";
        // circuit.dumpSolution(std::cout, solution.futureArray(), "  ");
        // std::cout << "\n";
        
        // Clear maximal residual contribution
        zero(maxResidualContribution_);
        
        auto [buildOk, preventConvergence] = buildSystem(continuePrevious);
        if (!buildOk) {
            // Load error or abort
            break;
        }
        xdelta[0] = 0.0; // Set RHS bucket to 0

        // Add forced values to the system
        if (haveForces() && !loadForces(true)) {
            if (settings.debug) {
                Simulator::dbg() << "Failed to load forced values at iteration " << iteration << "\n";
            }
            break;
        }

        // Check if system is finite
        if (settings.matrixCheck && !jac.isFinite(true, true)) {
            lastError = Error::LinearSolver;
            errorIteration = iteration;
            if (settings.debug) {
                Simulator::dbg() << "A matrix entry is not finite. Solver failed.\n";
            }
            break;
        }

        if (settings.rhsCheck && !jac.isFinite(dataWithoutBucket(delta), true, true)) {
            lastError = Error::LinearSolver;
            errorIteration = iteration;
            if (settings.debug) {
                Simulator::dbg() << "An RHS entry is not finite. Solver failed.\n";
            }
            break;
        }

        // std::cout << "Old solution before residual check at iteration " << iteration << "\n";
        // circuit.dumpSolution(std::cout, solution.array(), "  ");
        // std::cout << "\n";
        // std::cout << "New solution before residual check at iteration " << iteration << "\n";
        // circuit.dumpSolution(std::cout, solution.futureArray(), "  ");
        // std::cout << "\n";

        // Residual norms are needed if we have dynamic damping or debug is set
        bool computeResidualNorms = settings.dampingSteps>0 || settings.debug;
        // Assume residual is OK
        bool residualOk = true;
        // Call residualCheck() if 
        // - check is required and convergence was not prevented or
        // - residual norms are needed
        if (settings.residualCheck || computeResidualNorms) {
            bool residualCheckOk;
            std::tie(residualCheckOk, maxResidual, maxNormResidual, l2normResidual2, maxResidualNode) = 
                checkResidual(settings.residualCheck ? &residualOk : nullptr, computeResidualNorms);
            if (!residualCheckOk) { 
                if (settings.debug) {
                    Simulator::dbg() << "Residual check error.\n";
                }
                break;
            }
        }

        // std::cout << "Old solution after residual check at iteration " << iteration << "\n";
        // circuit.dumpSolution(std::cout, solution.array(), "  ");
        // std::cout << "\n";
        // std::cout << "New solution after residual check at iteration " << iteration << "\n";
        // circuit.dumpSolution(std::cout, solution.futureArray(), "  ");
        // std::cout << "\n";
        
        // std::cout << "States 0: ";
        // stateDump(states.vector(0));
        // std::cout << "\n";
        // 
        // std::cout << "States 1: ";
        // stateDump(states.vector(1));
        // std::cout << "\n";
        // 
        // std::cout << "Future states: ";
        // stateDump(states.futureVector());
        // std::cout << "\n" << size_t(states.futureData()) << "\n";

        if (settings.debug>=2) {
            Simulator::dbg() << "Linear system at iteration " << iteration << "\n";
            jac.dump(Simulator::dbg(), dataWithoutBucket(delta));
            Simulator::dbg() << "\n";
        }
        
        // jac.dumpSparsityTables(std::cout);
        // std::cout << std::endl;
        // jac.dumpEntries(std::cout);
        // std::cout << std::endl;
        
        // Weird... KLU does not detect singular matrix for two dangling serially connected resistors
        //          Probably due to numerical errors (1e-16 is taken as a valid pivot)
        //          if there is a nonzero excitation present the iteration won't converge 
        //          because it fails the residual check

        // Factorization
        bool forceFullFactorization = false;        
        if (jac.isFactored()) {
            // Refactor (if possible)
            if (!jac.refactor()) {
                // Failed, try again by fully factoring
                forceFullFactorization = true;
            } 
        }
        if (forceFullFactorization || !jac.isFactored()) {
            // Full factorization
            if (!jac.factor()) {
                // Failed, give up
                lastError = Error::LinearSolver;
                errorIteration = iteration;
                if (settings.debug) {
                    Simulator::dbg() << "LU factorization failed.\n";
                }
                break;
            }
        }

        // Solve, use vector without ground component (+1)
        if (!jac.solve(dataWithoutBucket(delta))) {
            lastError = Error::LinearSolver;
            errorIteration = iteration;
            if (settings.debug) {
                Simulator::dbg() << "Failed to solve factored system.\n";
            }
            break;
        }

        circuit.tables().accounting().acctNew.nriter++;

        if (settings.solutionCheck && !jac.isFinite(dataWithoutBucket(delta), true, true)) {
            lastError = Error::SolutionError;
            errorIteration = iteration;
            if (settings.debug) {
                Simulator::dbg() << "A solution entry is not finite. Solver failed.\n";
            }
            break;
        }

        // std::cout << "Old solution after solve at iteration " << iteration << "\n";
        // circuit.dumpSolution(std::cout, solution.array(), "  ");
        // std::cout << "\n";
        // std::cout << "New solution after solve at iteration " << iteration << "\n";
        // circuit.dumpSolution(std::cout, solution.futureArray(), "  ");
        // std::cout << "\n";
        
        // Set solution delta and new solution buckets to 0. 
        xdelta[0] = 0.0;
        xnew[0] = 0.0; 

        // std::cout << "Negative solution delta at iteration " << iteration << "\n";
        // circuit.dumpSolution(std::cout, delta.data(), "  ");
        // std::cout << "\n";

        // Delta norms are computed in debug mode
        bool computeDeltaNorms = settings.debug;
        // Check convergence (see if delta is small enough)
        // Assume delta did not converge
        bool deltaOk = false;
        // Check is not performed in first iteration
        if (iteration>1) {
            bool deltaCheckOk;
            std::tie(deltaCheckOk, maxDelta, maxNormDelta, maxDeltaNode) =
                checkDelta(&deltaOk, computeDeltaNorms);
            if (!deltaCheckOk) {
                if (settings.debug) {
                    Simulator::dbg() << "Solution delta check error.\n";
                }
                break;
            }
        }

        // Check if this iteration converged
        iterationConverged = !preventConvergence && deltaOk && residualOk;
        
        // Print debug messages
        if (settings.debug) {
            bool iterationConverged = deltaOk && residualOk;
            std::stringstream ss;
            ss << std::scientific << std::setprecision(2);
            Simulator::dbg() << "Iteration " << std::to_string(iteration) << (preventConvergence ? ", convergence not allowed" : "");
            if (!preventConvergence) {
                Simulator::dbg() << (iterationConverged ? ", converged" : "");
                if (computeResidualNorms) {
                    ss.str(""); ss << maxResidual;
                    Simulator::dbg() << ", worst residual=" << ss.str() << " @ " << maxResidualNode->name();
                }
                if (iteration>1) {
                    ss.str(""); ss << maxDelta;
                    Simulator::dbg() << ", worst delta=" << ss.str() << " @ " << maxDeltaNode->name();
                }
            }
            Simulator::dbg() << "\n";
        }

        // Count consecutive convergent iterations
        if (iterationConverged) {
            convIter++;
        } else {
            convIter = 0;
        }

        // Set converged flag
        converged = convIter >= settings.convIter;
        
        // std::cout << "Iteration " << iteration << " residualOK=" << residualOk << " deltaOk=" << deltaOk << "\n";

        // Exit if converged
        if (converged) {
            // Rotate states because the new state belongs to the current solution
            states.rotate();

            // std::cout << "On NR exit:\n";
            // std::cout << "States 0: ";
            // stateDump(states.vector(0));
            // std::cout << "\n";
            // 
            // std::cout << "States 1: ";
            // stateDump(states.vector(1));
            // std::cout << "\n";

            break;
        }

        // std::cout << "Old solution before update at iteration " << iteration << "\n";
        // circuit.dumpSolution(std::cout, solution.array(), "  ");
        // std::cout << "\n";
        // std::cout << "New solution before update at iteration " << iteration << "\n";
        // circuit.dumpSolution(std::cout, solution.futureArray(), "  ");
        // std::cout << "\n";
        
        // Not converged yet, compute new solution 
        if (settings.dampingSteps<=0) {
            // Static damping
            // Update pointMaxSolution_
            pointMaxSolution_ = 0;
            for(decltype(n) i=1; i<=n; i++) {
                xnew[i] = xprev[i] - xdelta[i]*settings.dampingFactor;
                double c = std::fabs(xnew[i]);
                if (c>pointMaxSolution_) {
                    pointMaxSolution_ = c;
                }
            }
        } else {
            // Dynamic damping, note that maxResidualContribVec and l2normResidual2 are available
            // Start at op_nrdamping
            double dampedMaxNormResidual;
            double dampedL2normResidual2;
            Int cnt=0;
            bool dampingResidualOk = false;
            Int dampingSteps = 0;
            double dampingFactor = settings.dampingFactor;
            do {
                // Compute new solution
                // Update pointMaxSolution_
                pointMaxSolution_ = 0;
                for(decltype(n) i=1; i<=n; i++) {
                    xnew[i] = xprev[i] - xdelta[i]*dampingFactor;
                    double c = std::fabs(xnew[i]);
                    if (c>pointMaxSolution_) {
                        pointMaxSolution_ = c;
                    }
                }

                // Zero residual/delta
                zero(delta);

                // Compute residual
                auto [residualComputationOk, dummy_] = computeResidual(continuePrevious);
                if (!residualComputationOk) {
                    // Residual computation error or abort
                    exitNrLoop = true;
                    break;
                }

                // Add forced values
                if (haveForces() && !loadForces(false)) {
                    if (settings.debug) {
                        Simulator::dbg() << "Failed to load forced values in iteration " << iteration << ", damping step " << (dampingSteps+1) << "\n";
                    }
                    break;
                }

                dampingSteps++;
                
                auto [residualCheckOk, dampedMaxResidual, dampedMaxNormResidual, dampedL2normResidual2, dummyNode] = 
                    checkResidual(nullptr, true);
                if (!residualCheckOk) {
                    if (settings.debug) {
                        Simulator::dbg() << "Residual norm computation failed.\n";
                    }
                    exitNrLoop = true;
                    break;
                }
                
                // Damping is successful if 
                // - max norm residual is below 1 (all components below tol) or
                // - squared L2 norm residual is not greater than squared L2 norm residual at old solution
                if (settings.debug) {
                    std::stringstream ss;
                    ss << std::scientific << std::setprecision(2);
                    ss.str(""); ss << dampingFactor;
                    Simulator::dbg() << "NR damping step " << dampingSteps << ", damping=" << ss.str();
                    if (dampedMaxNormResidual<1) { 
                        Simulator::dbg() << ", residual within tolerance.\n"; 
                    } else {
                        ss.str(""); ss << std::sqrt(dampedL2normResidual2);
                        Simulator::dbg() << ", L2 normalized residual " << ss.str();
                        ss.str(""); ss << std::sqrt(l2normResidual2);
                        if (dampedL2normResidual2<=l2normResidual2) {
                            Simulator::dbg() << " <= " << ss.str() << " (ok).\n"; 
                        } else {
                            Simulator::dbg() << " > " << ss.str() << " (too big).\n"; 
                        }
                    }
                }
                // Damped point residual is below tolerance or L2 norm is smaller than at old point
                dampingResidualOk = dampedMaxNormResidual<1 || dampedL2normResidual2<=l2normResidual2;
                if (dampingResidualOk) {
                    // Damping successfull, exit loop
                    break;
                }
                dampingFactor *= settings.dampingStep;
            } while (dampingSteps<settings.dampingSteps); // Damping loop

            // Update minimal damping factor
            if (dampingFactor<minimalDampingFactor) {
                minimalDampingFactor = dampingFactor;
                if (settings.debug) {
                   Simulator::dbg() << "Minimal damping decreased, adjusting NR iteration limit to " << size_t(itlim/minimalDampingFactor) << ".\n";  
                }
            }

            // What to do if damping fails to reduce residual?
            // This happens when 
            // - damping loop is exited due to abort or residual evaluation error
            // - avalable damping iterations are exhausted without finding a sufficiently small residual
            // At the moment just keep the shortened step. 
            if (!dampingResidualOk) {
                if (settings.debug) {
                    Simulator::dbg() << "Damping failed, keeping shortened step.\n";
                }
            }
        } 

        // std::cout << "Old solution after update at iteration " << iteration << "\n";
        // circuit.dumpSolution(std::cout, solution.array(), "  ");
        // std::cout << "\n";

        if (settings.debug>=3) {
            Simulator::dbg() << "New solution in iteration " << iteration << "\n";
            circuit.dumpSolution(Simulator::dbg(), solution.futureData(), "  ");
            Simulator::dbg() << "\n";
        }

        // Rotate RHS vectors and state (swap current and future)
        solution.rotate();
        states.rotate();

        // std::cout << "Old solution after rotation " << iteration << "\n";
        // circuit.dumpSolution(std::cout, solution.array(), "  ");
        // std::cout << "\n";
        // std::cout << "New solution after rotation " << iteration << "\n";
        // circuit.dumpSolution(std::cout, solution.futureArray(), "  ");
        // std::cout << "\n";

    } while (iteration<(itlim/minimalDampingFactor) && !exitNrLoop); // NR loop

    if (settings.debug) {
        Simulator::dbg() << "NR algorithm " << (converged ? "converged in " : "failed to converge in ") << iteration << " iteration(s).\n";
    }

    if (!converged && lastError==Error::OK) {
        lastError = Error::Convergence;
        errorIteration = iteration;
    }

    circuit.tables().accounting().acctNew.tnr += Accounting::wclkDelta(t0);
    return converged;
}

bool NRSolver::formatError(Status& s, NameResolver* resolver) const {
    auto matrixError = jac.formatError(s, resolver);
    switch (lastError) {
        case Error::ForcesIndex:
            s.set(Status::Range, "Force index out of range.");
            return false;
        case Error::EvalAndLoad:
            s.set(Status::NonlinearSolver, "Evaluation/load error.");
            s.extend("Leaving core NR loop in iteration "+std::to_string(errorIteration)+"."); 
            return false;
        case Error::LinearSolver: 
            s.extend("Leaving core NR loop in iteration "+std::to_string(errorIteration)+"."); 
            return false;
        case Error::SolutionError:
            s.extend("Solution component is not finite."); 
            s.extend("Leaving core NR loop in iteration "+std::to_string(errorIteration)+"."); 
            return false;
        case Error::Convergence:
            s.set(Status::NonlinearSolver, "NR solver failed to converge.");
            s.extend("Leaving core NR loop in iteration "+std::to_string(errorIteration)+"."); 
            return false;
        case Error::BadSolReference: 
            s.set(Status::NonlinearSolver, "Unsupported relrefsol value.");
            break;
        case Error::BadResReference: 
            s.set(Status::NonlinearSolver, "Unsupported relrefres value.");
            break;
    }
    return matrixError;
}


}
