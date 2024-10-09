#include "corehbnr.h"
#include "simulator.h"
#include "densematrix.h"
#include "common.h"

namespace NAMESPACE {

HBNRSolver::HBNRSolver(
        Circuit& circuit, 
        KluBlockSparseRealMatrix& jacColoc, 
        KluBlockSparseRealMatrix& bsjac, 
        VectorRepository<double>& solution, 
        Vector<Complex>& solutionFD, 
        Vector<Real>& frequencies, 
        Vector<Real>& timepoints, 
        DenseMatrix<Real>& DDT, 
        DenseMatrix<Real>& DDTcolMajor, 
        DenseMatrix<double>& APFT, 
        NRSettings& settings
) : circuit(circuit), jacColoc(jacColoc), bsjac(bsjac), solutionFD(solutionFD), 
    frequencies(frequencies), timepoints(timepoints), DDT(DDT), DDTcolMajor(DDTcolMajor), APFT(APFT), 
    NRSolver(circuit.tables().accounting(), bsjac, solution, settings) {
    // Slot 0 is for sweep continuation and homotopy
    // Do not need slot 1 as we do not support explicit nodesets
    resizeForces(1);

    // For constructing the linearized system in NR loop
    evalSetup_ = EvalSetup {
        // Inputs
        .solution = &oldSolutionAtTk, 
        .dummyStates = &dummyStates, 

        // Signal this is not a static DC analysis
        // Evaluation is in time domain so effectively 
        // we are doing the same thing as transient analysis. 
        .staticAnalysis = false, 
        .dcAnalysis = false, 
        .tranAnalysis = true, 

        // Evaluation
        // - no limiting
        // - resistive and reactive
        // - no opvars for now - maybe later we can collect their time-domain points
        //   and transform them to frequency domain before dumping them
        .enableLimiting = false, 
        .evaluateResistiveJacobian = true, 
        .evaluateReactiveJacobian = true, 
        .evaluateResistiveResidual = true, 
        .evaluateReactiveResidual = true, 
        .evaluateOpvars = true, 
    };

    loadSetup_ = LoadSetup {
        .states = nullptr, 
        .loadResistiveJacobian = true, 
        .loadReactiveJacobian = true, 
        // Used for loading with offset 0
        .reactiveJacobianFactor = 1.0, 
    };
}

bool HBNRSolver::setForces(Int ndx, const AnnotatedSolution& solution, bool abortOnError) {
    // Get forces
    auto& f = forces(ndx);

    // Clear forced values
    f.clear();
    
    // Number of unknowns
    auto n = circuit.unknownCount();

    // Number of components per unknown
    auto nt = timepoints.size();

    // Make space for variable forces (also include bucket)
    f.resizeUnknownForces(n*nt+1);

    // Number of frequencies in solution and solver
    auto nfSolution = solution.auxData().size();
    auto nfSolver = frequencies.size();

    // Prepare frequency translator between solution and solver
    // Translator stores the solver frequency index for each solutiuon frequency index
    std::vector<int> xlat;
    // Assume no frequency can be trasnlated (negative index)
    xlat.resize(nfSolution, -1);

    // Translate DC (it is always present)
    if (nfSolver>0 && nfSolution>0) {
        xlat[0] = 0;
    }

    // Translate the rest
    decltype(nfSolver) ndxSolver = 1;
    decltype(nfSolution) ndxSolution = 1;
    for(; ndxSolver<nfSolver && ndxSolution<nfSolution;) {
        auto fSolver = frequencies[ndxSolver];
        auto fSolution = solution.auxData()[ndxSolution];
        if (std::abs(fSolver-fSolution)<=std::max(std::abs(fSolver), std::abs(fSolution))*1e-14) {
            // Frequencies are almost the same, store translator
            xlat[ndxSolution] = ndxSolver;
            // Advance both indices
            ndxSolver++;
            ndxSolution++;
        } else if (fSolver<fSolution) {
            // Solver frequency is lower, advance solver index
            ndxSolver++;
        } else {
            // Solution frequency is lower, advance solution index
            ndxSolution++;
        }
    }

    bool error = false;

    // Go hrough annotated solution. fill APFT spectrum 
    // Use resistive residual vector for APFT spectrum
    auto& forcesFD = resistiveResidual;
    // Set APFT spectrum to 0
    forcesFD.resize(n*nt, 0);

    // Force values at colocation points are stored in reactive residual vector
    auto& forcesTD = reactiveResidual;
    forcesTD.resize(n*nt);

    // Solution spectrum
    auto& solSpec = solution.values();

    // Copy matching frequency components
    if (nfSolution>0) {
        // Copy DC from APFT spectrum of stored solution
        forcesFD[0] = solSpec[0];
        // Copy the rest
        for(decltype(nfSolution) i=1; i<nfSolution; i++) {
            auto xlfi = xlat[i];
            if (xlfi>=0) {
                // Translation exists
                auto ndx = (i-1)*2+1;
                auto destNdx = (xlfi-1)*2+1;
                forcesFD[destNdx] = solSpec[ndx];
                forcesFD[destNdx+1] = solSpec[ndx+1];
            }
        }
    }

    // Perform inverse APFT
    
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

        if (!f.setForceOnUnknown(node, value)) {
            error = true;
            if (abortOnError) {
                return false;
            }
        }
    }

    return !error;
}

bool HBNRSolver::rebuild() {
    // Call parent's rebuild
    if (!NRSolver::rebuild()) {
        // Assume parent has set the error flag
        return false;
    }

    // Old states at one timepoint (dummy vector of zeros because we do no limiting)
    dummyStates.resize(circuit.statesCount());
    zero(dummyStates);
    
    // Jacobian is sized in core or analysis.
    // solution is sized in core. 
    // delta is resized by NRSolver::rebuild() based on Jacobian size. 

    // Because vector lengths and Jacobian size may change 
    // due to different number of frequency components 
    // we size them in initialize() before first iteration. 
    // Analysis asks cores if they request a rebuild. 
    // HB core replies that it does if the set of frequencies changes. 
    
    return true;
}

bool HBNRSolver::initialize(bool continuePrevious) {
    // Number fo frequency components and timepoints
    auto nt = timepoints.size();
    
    // DDT, APFT, and IAPFT are already set up

    // Number of nodes
    auto n = circuit.unknownCount();

    // Number of block rows
    auto nb = bsjac.nBlockElementRows();

    // Old solution and derivative wrt time at all timepoints
    resistiveResidual.resize(n*nt);
    reactiveResidual.resize(n*nt);

    // Old solution and resistive residual at one timepoint
    // Includes ground node because it is used by evalAndLoad()
    oldSolutionAtTk.upsize(1, n+1);
    resistiveResidualAtTk.resize(n+1);
    reactiveResidualAtTk.resize(n+1);

    // Maximum residual contribution at single timepoint
    maxResidualContributionAtTk_.resize(n+1);

    // Maximum residual contribution for each equation at each timepoint
    maxResidualContribution_.resize(n*nt);

    // Maximum across all equations at given timepoint for each nature
    // Computed in checkResidual()
    pointMaxResidualContribution_.resize(2, nb);

    // Maximum across all timepoints for each nature
    // Computed in checkDelta()
    pointMaxSolution_.resize(2, nb);

    // Set up loading
    // Resistive residual
    loadSetup_.resistiveResidual = resistiveResidualAtTk.data();
    // Reactive residual
    loadSetup_.reactiveResidual = reactiveResidualAtTk.data();

    // Maximal residual contribution computed by evalAndLoad()
    loadSetup_.maxResistiveResidualContribution = maxResidualContributionAtTk_.data();

    // Zero states
    zero(dummyStates);

    // Set up tolerance reference value for solution
    auto& options = circuit.simulatorOptions().core();
    if (options.relrefsol==SimulatorOptions::relrefPointLocal) {
        globalSolRef = false;
        // historicSolRef = false;
    } else if (options.relrefsol==SimulatorOptions::relrefLocal) {
        globalSolRef = false;
        // historicSolRef = true;
    } else if (options.relrefsol==SimulatorOptions::relrefPointGlobal) {
        globalSolRef = true;
        // historicSolRef = false;
    } else if (options.relrefsol==SimulatorOptions::relrefGlobal) {
        globalSolRef = true;
        // historicSolRef = true;
    } else if (options.relrefsol==SimulatorOptions::relrefRelref) {
        if (options.relref == SimulatorOptions::relrefAlllocal) {
            globalSolRef = false;
            // historicSolRef = true;
        } else if (options.relref == SimulatorOptions::relrefSigglobal) {
            globalSolRef = true;
            // historicSolRef = true;
        } else if (options.relref == SimulatorOptions::relrefAllglobal) {
            globalSolRef = true;
            // historicSolRef = true;
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
    } else if (options.relrefres==SimulatorOptions::relrefLocal) {
        globalResRef = false;
    } else if (options.relrefres==SimulatorOptions::relrefPointGlobal) {
        globalResRef = true;
    } else if (options.relrefres==SimulatorOptions::relrefGlobal) {
        globalResRef = true;
    } else if (options.relrefres==SimulatorOptions::relrefRelref) {
        if (options.relref == SimulatorOptions::relrefAlllocal) {
            globalResRef = false;
        } else if (options.relref == SimulatorOptions::relrefSigglobal) {
            globalResRef = false;
        } else if (options.relref == SimulatorOptions::relrefAllglobal) {
            globalResRef = true;
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

bool HBNRSolver::preIteration(bool continuePrevious) {
    return true;
}

bool HBNRSolver::postSolve(bool continuePrevious) {
    // Nothing to do - we have no bypassing
    return true;
}

bool HBNRSolver::postConvergenceCheck(bool continuePrevious) {
    // Print debug information on convergence
    if (settings.debug) {
        std::stringstream ss;
        ss << std::scientific << std::setprecision(2);
        Simulator::dbg() << "Iteration " << std::to_string(iteration) << (preventedConvergence ? ", convergence not allowed" : "");
        if (!preventedConvergence) {
            Simulator::dbg() << (iterationConverged ? ", converged" : "");
            if (settings.residualCheck) {
                ss.str(""); 
                ss << maxResidual;
                Simulator::dbg() << ", worst residual=" << ss.str() << " @ " << (maxResidualNode ? maxResidualNode->name() : "(unknown)")
                                 << "~t" << maxResidualTimepointIndex << "=" << timepoints[maxResidualTimepointIndex];
            }
            if (iteration>1) {
                ss.str(""); ss << maxDelta;
                Simulator::dbg() << ", worst delta=" << ss.str() << " @ " << (maxDeltaNode ? maxDeltaNode->name() : "(unknown)")
                                 << "~t" << maxDeltaTimepointIndex << "=" << timepoints[maxDeltaTimepointIndex];
            }
        }
        Simulator::dbg() << "\n";
    }
    return true;
}

bool HBNRSolver::postIteration(bool continuePrevious) {
    return true;
}

bool HBNRSolver::postRun(bool continuePrevious) {
    if (converged) {
        // If converged, convert solution from TD to FD
        auto n = circuit.unknownCount();
        auto nf = frequencies.size();
        auto nt = timepoints.size();
        solutionFD.resize(n*nf); // no bucket

        for(decltype(n) i=0; i<n; i++) {
            auto cxSpecPtr = solutionFD.data()+nf*i;
            // APFT computes spectrum as complex values, with the exception of DC which is stored as real. 
            // We write APFT output starting at the imaginary part of the DC complex magnitude. 
            auto inPtr = solution.data()+1+i*nt;
            auto outPtr = reinterpret_cast<double*>(cxSpecPtr)+1;
            auto outVec = VectorView<Real>(outPtr, nt, 1);
            APFT.multiply(VectorView<Real>(inPtr, nt, 1), outVec);
            // Move DC from imaginary to real part of DC complex magnitude. 
            *cxSpecPtr = cxSpecPtr->imag();
        }
    }
    return true;
}

bool HBNRSolver::evalAndLoadWrapper(EvalSetup& evalSetup, LoadSetup& loadSetup) {
    lastError = Error::OK;
    evalSetup.requestHighPrecision = highPrecision;
    if (!circuit.evalAndLoad(&evalSetup, &loadSetup, nullptr)) {
        // Load error
        lastError = Error::EvalAndLoad;
        if (settings.debug>2) {
            Simulator::dbg() << "Evaluation error.\n";
        }
        return false;
    }

    // Store Abort, Finish, and Stop flag
    if (evalSetup_.requests.abort) {
        setFlags(Flags::Abort);
    }
    if (evalSetup_.requests.finish) {
        setFlags(Flags::Finish);
    }
    if (evalSetup_.requests.stop) {
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

std::tuple<bool, bool> HBNRSolver::buildSystem(bool continuePrevious) {
    // Jacobian values at colocation points are stored in jacColoc with dense 
    // blocks of size nt x 2, where nt is the number of colocation points. 
    // Resistive Jacobian is bound to 0-based subentry (0, 0) of each dense block. 
    // Reactive Jacobian is bound to 0-based subentry (0, 1) of each dense block. 
    // As Jacobian load offset goes from 0..nb-1 (nb=ntimepoints)
    // the first two columns of each block are loaded with Jacobian values, 
    // i.e. (k, 0) with resistive and (k, 1) with reactive Jacobian values 
    // at times coresponding to timepoints tk, k=0..nb-1 because KLU matrices are 
    // stored in column major order. 
    // 
    // Let Jr_ijk and Jc_ijk denote the resistive and reactive Jacobian value 
    // from block with 1-based position (i+1, j+1) at timepoint with index k. 
    // i, j and k are all 0-based. 
    // After nt evalAndLoad() calls the two columns of each block in jacColoc 
    // are filled with Jr_ijk and Jc_ijk. 
    // 
    // The unknowns and the equations are in time domain, i.e. we formulate 
    // HB in time domain. Let Gij denote the diagonal nb x nb matrix holding 
    // the resistive Jacobian values for the block at position (i+1, j+1). 
    // Similarly Cij is the diagonal nb x nb matrix holding the reactive 
    // Jacobian values for the block. The HB Jacobian block can then be 
    // expressed as
    //   Gij + Gamma^-1 Omega Gamma Cij
    // where Gamma and Gamma^-1 are the APFT and inverse APFT matrices, 
    // while Omega is the time-derivative operator in frequency domain. 
    // Matrix DDT holds Gamma^-1 Omega Gamma. 
    // The HB Jacobian block at (i+1, j+1) is constructed by scaling columns 
    // of DDT with reactive Jacobian values Jc_ijk and then adding the 
    // resistive Jacoban values Jr_ijk to the diagonal. 
    // 
    // Because blocks are stored in column major order the innermost loop 
    // iterates over values of k. In this way good cache locality is achieved. 
    
    // Get sizes
    auto n = bsjac.nBlockRows();
    auto nb = bsjac.nBlockElementRows();

    // Clear maximal residual contribution
    zero(maxResidualContribution_);

    // Old solution is in time domain. Get it. 
    auto solTD = solution.data();

    // Clear Jacobian at colocation points
    jacColoc.zero();
    
    // Loop through timepoints 0..nb-1
    for(decltype(nb) k=0; k<nb; k++) {
        // We read old solution starting at index 1+k
        // (skip bucket, k-th unknown, first timepoint)
        // Vector length n, stride nb
        // We write to the vector of old solutions at timepoint t_k, 
        // start at index 1 (skip bucket), length n, stride 1
        VectorView(oldSolutionAtTk.vector(), 1, n, 1) = VectorView(solution.vector(), 1+k, n, nb);

        // Zero residual vectors where evalAndLoad() will load the residuals at t_k
        zero(resistiveResidualAtTk);
        zero(reactiveResidualAtTk);

        // Zero maximal residual contribution at timepoint
        zero(maxResidualContributionAtTk_);
        
        // Set time and offset
        evalSetup_.time = timepoints[k];
        loadSetup_.jacobianLoadOffset = k;

        // This is a kludge, implement real source stepping
        decltype(settings.itlim) ramp;
        if (continuePrevious) {
            ramp = settings.itlimCont*0.1;
        } else {
            ramp = settings.itlim*0.1;
        }
        if (ramp<1) {
            ramp = 1;
        }
        // NR starts with iteration=1
        circuit.simulatorInternals().sourcescalefactor = (iteration-1)<ramp ? (iteration-1)*1.0/ramp : 1.0;
    
        // For k-th timepoint (t_k) load 
        // - resistive and reactive Jacobian at t_k with offset i, 
        // - resistive residuals for all equations at t_k
        // - reactive residuals for all equations at t_k
        // Values are stored in jacColoc. 
        auto ok = evalAndLoadWrapper(evalSetup_, loadSetup_);

        // Put resistive residuals at t_k in the residuals vector
        VectorView(resistiveResidual, k, n, nb) = VectorView(resistiveResidualAtTk, 1, n, 1);
        VectorView(reactiveResidual, k, n, nb) = VectorView(reactiveResidualAtTk, 1, n, 1);

        // Put maximal resistive residual contribution at t_k into maxResidualContribution_
        VectorView(maxResidualContribution_, k, n, nb) = 
            VectorView(maxResidualContributionAtTk_, 1, n, 1);
    }

    // For each block (ordered in column major order)
    for(auto& pos : circuit.sparsityMap().positions()) {
        // Get block position (for debugging), make position 0-based
        auto [i, j] = pos;
        i--;
        j--;

        // Get dense block with Jacobian values at colocation points
        auto [colocBlock, found1] = jacColoc.block(pos);

        // Get HB Jacobian dense block
        auto [block, found2] = bsjac.block(pos);

        // Get g_ijk and c_ijk columns from block (column elements are indexed by k)
        auto gCol = colocBlock.column(0);
        auto cCol = colocBlock.column(1);

        // Scan columns in block in from 2 to nb-1 
        for(decltype(nb) l=0; l<nb; l++) {   
            // Get target column
            auto targetColumn = block.column(l);

            // Get DDT column
            auto ddtColumn = DDTcolMajor.column(l);

            // Write scaled Jc_ijk
            targetColumn.writeScaled(ddtColumn, cCol[l]);

            // Add diagonal Jr_ijk
            targetColumn[l] += gCol[l];
        }

        // Now handle residuals
        // delta is zeroed at the beginning of each iteration by NRSolver

        // Block-transform reactive residual with DDT, store it in delta
        auto resPtr = resistiveResidual.data();
        auto reacPtr = reactiveResidual.data();
        auto maxResPtr = maxResidualContribution_.data();
        // Skip bucket
        auto deltaPtr = delta.data()+1;
        for(decltype(n) i=0; i<n; i++) {
            // Perform DDT on reactive residual block
            VectorView dest(deltaPtr, nb, 1);
            DDT.multiply(
                // length nb, stride 1
                VectorView(reacPtr, nb, 1), 
                dest
            );
            
            // Take values from dest and update maximal residual contribution with them
            VectorView maxRes(maxResPtr, nb, 1);
            for(decltype(nb) k=0; k<nb; k++) {
                auto c = std::abs(dest[k]);
                if (c>maxRes[k]) {
                    maxRes[k] = c;
                }
            }

            // Add resistive residual to block
            dest.add(VectorView(resPtr, nb, 1));

            // Move on to next block
            resPtr += nb;
            reacPtr += nb;
            deltaPtr += nb;
            maxResPtr += nb;
        }

        // Bucket
        delta[0] = 0.0;
    }

    // OK, do not prevent convergence
    return std::make_tuple(true, false); 
}

std::tuple<bool, bool> HBNRSolver::checkResidual() {
    // Compute norms only in debug mode
    bool computeNorms = settings.debug;

    // In residual we have the residual at previous solution
    // We are going to check that residual
    
    // Number of unknowns excluding ground
    auto n = circuit.unknownCount();

    // Number of timepoints
    auto nt = timepoints.size();

    // Results
    maxResidual = 0.0;
    maxNormResidual = 0.0;
    l2normResidual2 = 0.0;
    maxResidualNode = nullptr;
    maxResidualTimepointIndex = 0;
    
    // Assume residual is OK
    bool residualOk = true;
    
    // Get point maximum for each residual nature
    pointMaxResidualContribution_.zero(); 
    // Loop through all nodes
    auto compPtr = maxResidualContribution_.data();
    for(decltype(n) i=0; i<n; i++) {
        // Get representative node (1-based index) and nature index
        auto rn = circuit.reprNode(i+1);
        bool isPotential = ((rn->flags() & Node::Flags::PotentialNode) == Node::Flags::PotentialNode); 
        size_t ndx = isPotential ? 1 : 0;
        // Loop through all timepoints
        for(decltype(nt) k=0; k<nt; k++) {
            double c = std::fabs(*compPtr);
            if (c>pointMaxResidualContribution_.at(ndx, k)) {
                pointMaxResidualContribution_.at(ndx, k) = c;
            }
            compPtr++;
        }
    }
    
    // Go through all variables (except ground)
    for(decltype(n) i=1; i<=n; i++) {
        // Representative node (1-based index), associated flow nature index
        auto rn = circuit.reprNode(i);
        bool isPotential = ((rn->flags() & Node::Flags::PotentialNode) == Node::Flags::PotentialNode); 
        size_t ndx = isPotential ? 1 : 0;
        
        // Go through all timepoints
        for(decltype(nt) k=0; k<nt; k++) {
            // Compute tolerance reference
            // Point local reference by default
            // Compute tolerance reference, start with previous value of the i-th unknown
            double tolref = std::fabs(maxResidualContribution_[(i-1)*nt+k]);
            
            // Account for global references
            if (globalResRef) {
                // Point global reference, ndx is the nature index
                tolref = std::max(tolref, pointMaxResidualContribution_.at(ndx, k));
            }

            // Residual tolerance (Designer's Guide to Spice and Spectre, chapter 2.2.2)
            auto tol = circuit.residualTolerance(rn, tolref);

            // Residual component
            double rescomp = fabs(delta[1+(i-1)*nt+k]);

            // Normalized residual component
            double normResidual = rescomp/tol;

            if (computeNorms) {
                l2normResidual2 += normResidual*normResidual;
                // Update largest normalized component
                if (i==0 || normResidual>maxNormResidual) {
                    maxResidual = rescomp;
                    maxNormResidual = normResidual;
                    maxResidualNode = rn;
                    maxResidualTimepointIndex = k;
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
    }
    
    return std::make_tuple(true, residualOk); 
}

std::tuple<bool, bool> HBNRSolver::checkDelta() {
    // Compute norms only in debug mode
    bool computeNorms = settings.debug;

    // In delta we have the solution change
    // Check it for convergence
    
    // Number of unknowns (vector length includes a bucket at index 0)
    auto n = circuit.unknownCount();
    auto nt = timepoints.size();

    maxDelta = 0.0;
    maxNormDelta = 0.0;
    maxDeltaNode = nullptr;
    maxDeltaTimepointIndex = 0;
    
    // Check convergence (see if delta is small enough), 
    // but only if this is iteration 2 or later
    // In iteration 1 assume we did not converge
    
    // Assume we converged
    bool deltaOk = true;
    
    // Get point maximum for each solution nature
    auto xold = solution.data();
    // Skip bucket
    for(decltype(n) i=1; i<=n; i++) {
        // Representative node, associated potential nature index
        auto rn = circuit.reprNode(i);
        bool isPotential = ((rn->flags() & Node::Flags::PotentialNode) == Node::Flags::PotentialNode); 
        size_t ndx = isPotential ? 0 : 1;
        for(decltype(nt) k=0; k<nt; k++) {
            double c = std::fabs(xold[i]);
            // Rows are natures, columns are frequency components (DC, f1, f2, ...)
            if (c>pointMaxSolution_.at(ndx, k)) {
                pointMaxSolution_.at(ndx, k) = c;
            }
        }
    }

    // Use 1-based index (with bucket) because same indexing is used for variables
    auto xdelta = delta.data();
    for(decltype(n) i=1; i<=n; i++) {
        // Representative node, associated potential nature index
        auto rn = circuit.reprNode(i);
        bool isPotential = ((rn->flags() & Node::Flags::PotentialNode) == Node::Flags::PotentialNode); 
        size_t ndx = isPotential ? 0 : 1;
        for(decltype(nt) j=0; j<nt; j++) {
            // Compute tolerance reference
            // Point local reference by default
            // Compute tolerance reference, start with previous value of the i-th unknown at j-th frequency
            double tolref = xold[i];
            
            // Account for global references, no historic reference because we are in frequency domain
            if (globalSolRef) {
                // Point global reference, ndx is the nature index
                tolref = std::max(tolref, pointMaxSolution_.at(ndx, j));
            }
            
            // Compute tolerance
            double tol = circuit.solutionTolerance(rn, tolref);

            // Absolute solution change 
            double deltaAbs = fabs(xdelta[i]);;

            if (computeNorms) {
                double normDelta = deltaAbs/tol;
                if (i==1 || normDelta>maxNormDelta) {
                    maxDelta = deltaAbs;
                    maxNormDelta = normDelta;
                    maxDeltaNode = rn;
                    maxDeltaTimepointIndex = j;
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
    }
    
    return std::make_tuple(true, deltaOk);
}

void HBNRSolver::dumpSolution(std::ostream& os, double* solution, const char* prefix) {
    auto n = circuit.unknownCount();
    auto nt = timepoints.size();
    for(decltype(n) i=1; i<=n; i++) {
        auto rn = circuit.reprNode(i);
        for(decltype(nt) k=0; k<nt; k++) {
            os << prefix << rn->name() << "@t" << k << " : " << solution[1+(i-1)*nt+k] << "\n";
        }
    }
}

}
