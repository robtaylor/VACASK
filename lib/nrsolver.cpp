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
// is the Jacobian of f, also denoted by J. 
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
    Accounting& acct, KluRealMatrix& jac, 
    VectorRepository<double>& solution, 
    NRSettings& settings
) : acct(acct), jac(jac), solution(solution), settings(settings), 
    iteration(0) {
}

bool NRSolver::rebuild() {
    // Allocate space in vectors
    // Jacobian is already built, get number of unknowns excluding ground
    auto n = jac.nRow();
    diagPtrs.resize(n+1);
    delta.resize(n+1);
    rowNorm.resize(n+1);
    
    // Bind diagonal matrix elements
    // Needed for forcing unknown values and setting gshunts
    for(decltype(n) i=0; i<n; i++) {
        diagPtrs[i+1] = jac.elementPtr(MatrixEntryPosition(i+1, i+1), Component::Real);
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
                    jac.elementPtr(MatrixEntryPosition(u1, u2), Component::Real),
                    jac.elementPtr(MatrixEntryPosition(u2, u1), Component::Real)
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
    auto n = jac.nRow();
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

void NRSolver::resizeForces(Int n) {
    forcesList.resize(n);
    forcesEnabled.resize(n, false);
}

Forces& NRSolver::forces(Int ndx) {
    return forcesList.at(ndx);
} 

bool NRSolver::run(bool continuePrevious) {
    auto t0 = Accounting::wclk();
    acct.acctNew.nrcall++;

    jac.setAccounting(acct);

    // Number of unknowns (vector length includes a bucket at index 0)
    auto n = solution.length()-1;
    
    // Iteration limit and damping
    int itlim;
    if (continuePrevious) {
        itlim = settings.itlimCont;
    } else {
        itlim = settings.itlim;
    }

    // Main loop setup
    bool deltaOk;
    bool residualOk;
    double maxResidual;
    double maxNormResidual;
    double l2normResidual2;
    Id maxResidualNodeId; 
    double maxDelta;
    double maxNormDelta; 
    Id maxDeltaNodeId;
    iteration = 0;
    Int convIter = 0;
    converged = false;

    // Allow lower precision
    highPrecision = false;
    
    // Clear error
    clearError();

    // Initialize structures
    if (!initialize(continuePrevious)) {
        // Assume initialize() has set lastError
        acct.acctNew.tnr += Accounting::wclkDelta(t0);
        return false;
    }

    // If not in continue mode set current solution to 0
    if (!continuePrevious) {
        // Zero current solution
        solution.zero();
    }

    if (settings.debug) {
        Simulator::dbg() << "Starting NR algorithm " << (continuePrevious ? "with given initial solution" : "with zero initial solution") << ".\n";
    }
    
    // Main loop
    do {
        // Assume no convergence
        iterationConverged = false;

        // Iteration counter (first iteration is 1)
        iteration++;

        // Zero matrix, new solution, and residual/delta
        jac.zero();
        solution.zeroFuture();
        zero(delta);
        
        double* xprev = solution.data();
        double* xnew = solution.futureData();
        double* xdelta = delta.data();

        if (settings.debug>=4) {
            std::cout << "Old solution at iteration " << iteration << "\n";
            dumpSolution(std::cout, solution.data(), "  ");
            std::cout << "\n";
        }
        // std::cout << "New solution at iteration " << iteration << "\n";
        // dumpSolution(std::cout, solution.futureArray(), "  ");
        // std::cout << "\n";
        
        // Pre-iteration tasks
        if (!preIteration(continuePrevious)) {
            if (settings.debug) {
                Simulator::dbg() << "Pre-iteration step failed.\n";
            }
            converged = false;
            break;
        }
        
        auto [buildOk, preventedConvergence] = buildSystem(continuePrevious);
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
        // dumpSolution(std::cout, solution.array(), "  ");
        // std::cout << "\n";
        // std::cout << "New solution before residual check at iteration " << iteration << "\n";
        // dumpSolution(std::cout, solution.futureArray(), "  ");
        // std::cout << "\n";

        // Assume residual is OK
        bool residualOk = true;
        // Call residualCheck() if check is required
        if (settings.residualCheck) {
            bool residualCheckOk;
            std::tie(residualCheckOk, residualOk) = checkResidual();
            if (!residualCheckOk) { 
                if (settings.debug) {
                    Simulator::dbg() << "Residual check error.\n";
                }
                break;
            }
        }

        // std::cout << "Old solution after residual check at iteration " << iteration << "\n";
        // dumpSolution(std::cout, solution.array(), "  ");
        // std::cout << "\n";
        // std::cout << "New solution after residual check at iteration " << iteration << "\n";
        // dumpSolution(std::cout, solution.futureArray(), "  ");
        // std::cout << "\n";
        
        if (settings.debug>=2) {
            Simulator::dbg() << "Linear system at iteration " << iteration << "\n";
            jac.dump(Simulator::dbg(), dataWithoutBucket(delta));
            Simulator::dbg() << "\n";
        }
        
        // jac.dumpSparsityTables(std::cout);
        // std::cout << "\n";
        // jac.dumpEntries(std::cout);
        // std::cout << "\n";
        
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

        // Post-solve tasks (instance convergence check)
        if (!postSolve(continuePrevious)) {
            if (settings.debug) {
                Simulator::dbg() << "Post-solve step failed.\n";
            }
            converged = false;
            break;
        }

        acct.acctNew.nriter++;

        if (settings.solutionCheck && !jac.isFinite(dataWithoutBucket(delta), true, true)) {
            lastError = Error::SolutionError;
            errorIteration = iteration;
            if (settings.debug) {
                Simulator::dbg() << "A solution entry is not finite. Solver failed.\n";
            }
            break;
        }

        // std::cout << "Old solution after solve at iteration " << iteration << "\n";
        // dumpSolution(std::cout, solution.array(), "  ");
        // std::cout << "\n";
        // std::cout << "New solution after solve at iteration " << iteration << "\n";
        // dumpSolution(std::cout, solution.futureArray(), "  ");
        // std::cout << "\n";
        
        // Set solution delta and new solution buckets to 0. 
        xdelta[0] = 0.0;
        xnew[0] = 0.0; 

        // std::cout << "Negative solution delta at iteration " << iteration << "\n";
        // dumpSolution(std::cout, delta.data(), "  ");
        // std::cout << "\n";

        // Check convergence (see if delta is small enough)
        // Assume delta did not converge
        deltaOk = false;
        // Check is not performed in the first iteration
        if (iteration>1) {
            bool deltaCheckOk;
            std::tie(deltaCheckOk, deltaOk) = checkDelta();
            if (!deltaCheckOk) {
                if (settings.debug) {
                    Simulator::dbg() << "Solution delta check error.\n";
                }
                break;
            }
        }

        // If delta is OK, but residual is not, request increased precision
        if (deltaOk && !residualOk) { 
            highPrecision = true;
        } else {
            highPrecision = false;
        }


        // Check if this iteration converged
        iterationConverged = !preventedConvergence && deltaOk && residualOk;
        
        // Count consecutive convergent iterations
        if (iterationConverged) {
            convIter++;
        } else {
            convIter = 0;
        }

        // Set converged flag based on how many consecutive iterations converged
        converged = convIter >= settings.convIter;
        
        // std::cout << "Iteration " << iteration << " residualOK=" << residualOk << " deltaOk=" << deltaOk << "\n";

        // Post convergence check tasks
        if (!postConvergenceCheck(continuePrevious)) {
            if (settings.debug) {
                Simulator::dbg() << "Post convergence check step failed.\n";
            }
            converged = false;
            break;
        }

        // Exit if converged
        if (converged) {
            break;
        }

        // std::cout << "Old solution before update at iteration " << iteration << "\n";
        // dumpSolution(std::cout, solution.array(), "  ");
        // std::cout << "\n";
        // std::cout << "New solution before update at iteration " << iteration << "\n";
        // dumpSolution(std::cout, solution.futureArray(), "  ");
        // std::cout << "\n";
        
        // Not converged yet, compute new solution 
        
        // Compute new solution, use static damping
        for(decltype(n) i=1; i<=n; i++) {
            xnew[i] = xprev[i] - xdelta[i]*settings.dampingFactor;
        }

        // std::cout << "Old solution after update at iteration " << iteration << "\n";
        // dumpSolution(std::cout, solution.array(), "  ");
        // std::cout << "\n";

        if (settings.debug>=3) {
            Simulator::dbg() << "New solution in iteration " << iteration << "\n";
            dumpSolution(Simulator::dbg(), solution.futureData(), "  ");
            Simulator::dbg() << "\n";
        }

        // Rotate RHS vectors (swap current and future)
        solution.rotate();

        // Post-iteration tasks
        if (!postIteration(continuePrevious)) {
            if (settings.debug) {
                Simulator::dbg() << "Post-iteration step failed.\n";
            }
            converged = false;
            break;
        }
    
        // std::cout << "Old solution after rotation " << iteration << "\n";
        // dumpSolution(std::cout, solution.array(), "  ");
        // std::cout << "\n";
        // std::cout << "New solution after rotation " << iteration << "\n";
        // dumpSolution(std::cout, solution.futureArray(), "  ");
        // std::cout << "\n";

    } while (iteration<itlim); // NR loop

    if (settings.debug) {
        Simulator::dbg() << "NR algorithm " << (converged ? "converged in " : "failed to converge in ") << iteration << " iteration(s).\n";
    }

    if (!converged && lastError==Error::OK) {
        lastError = Error::Convergence;
        errorIteration = iteration;
    }

    acct.acctNew.tnr += Accounting::wclkDelta(t0);
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
        case Error::ConvergenceCheck:
            s.set(Status::NonlinearSolver, "Instance convergence check error.");
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
