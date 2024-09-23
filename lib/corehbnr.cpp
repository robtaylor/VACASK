#include "corehbnr.h"
#include "simulator.h"
#include "densematrix.h"
#include "common.h"

namespace NAMESPACE {

HBNRSolver::HBNRSolver(
    Circuit& circuit, KluBlockSparseRealMatrix& bsjac, 
        VectorRepository<double>& solution, 
        Vector<Real>& spectrum, 
        Vector<Real>& timepoints, 
        DenseMatrix<Real>& XF, 
        DenseMatrix<Real>& XFdot, 
        NRSettings& settings
) : circuit(circuit), bsjac(bsjac) , 
    spectrum(spectrum), timepoints(timepoints), XF(XF), XFdot(XFdot), 
    NRSolver(circuit.tables().accounting(), bsjac, solution, settings) {
    resizeForces(0);

    // For constructing the linearized system in NR loop
    evalSetup_ = EvalSetup {
        // Inputs
        .solution = &oldSolutionTD, 
        .states = &dummyStatesRepo, 

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
        .states = &dummyStatesRepo, 
        .loadResistiveJacobian = true, 
        .loadReactiveJacobian = true, 
    };
}

bool HBNRSolver::rebuild() {
    // Call parent's rebuild
    if (!NRSolver::rebuild()) {
        // Assume parent has set the error flag
        return false;
    }

    // Old states at one timepoint (dummy vector of zeros because we do no limiting)
    dummyStatesRepo.upsize(1, circuit.statesCount());
    dummyStatesRepo.zero();
    
    // Jacobian is sized in core or analysis.
    // solution is sized in core. 
    // delta is resized by NRSolver::rebuild() based on Jacobian size. 

    // Because vector lengths and Jacobian size may change 
    // due to different number of frequency components 
    // we size them in initialize() before first iteration. 
    // Analysis asks cores if they request a rebuild. 
    // HB core replies that it does if the spectrum changes. 
    
    return true;
}

bool HBNRSolver::initialize(bool continuePrevious) {
    // Number fo frequency components and timepoints
    auto nf = spectrum.size();
    auto nt = timepoints.size();
    // nt = 2*nf-1

    // XF and XFdot are already set up

    // Number of nodes
    auto n = circuit.unknownCount();

    // Number of block rows
    auto nb = bsjac.nBlockElementRows();

    // Resize vectors and matrices
    // XF and XFdot row maximum
    XFrowMax.resize(nt);
    XFdotRowMax.resize(nt);
    
    // Old solution and derivative wrt time at all timepoints
    // Bucket for ground node is needed because these vectors 
    // are used by NRSolver which assumes a bucket of length 1. 
    oldSolutionTD.upsize(1, n*nt+1);
    oldSolutionTDDot.resize(n*nt+1);

    // Old solution and resistive residual at one timepoint
    // Includes ground node because it is used by evalAndLoad()
    oldSolutionTDtk.upsize(1, n+1);
    resistiveResidualAtTk.resize(n+1);

    // Maximum residual contribution at single timepoint
    maxResidualContributionAtTk_.resize(n+1);

    // Maximum residual contribution for each equation at each timepoint
    maxResidualContribution_.resize(n*nt+1);

    // Maximum across all timepoints for each equation
    historicMaxResidualContribution_.resize(n+1);

    // Maximum across all timepoints and equations for each nature
    globalMaxResidualContribution_.resize(2);

    // Maximum across all equations at given timepoint for each nature
    // Computed in checkResidual()
    pointMaxResidualContribution_.resize(2);

    // Maximum across all (complex) unknowns for each nature and frequency
    // Computed in checkDelta()
    pointMaxSolution_.resize(2, nf);

    // Set up loading
    // Resistive residual
    loadSetup_.resistiveResidual = resistiveResidualAtTk.data();
    // Maximal residual contribution computed by evalAndLoad()
    loadSetup_.maxResistiveResidualContribution = maxResidualContributionAtTk_.data();

    // Zero states
    dummyStatesRepo.zero();
    
    // Compute maximal row value for XF and XFdot
    XFrowMax.resize(nt);
    XFdotRowMax.resize(nt);
    zero(XFrowMax);
    zero(XFdotRowMax);
    for(decltype(nt) i=0; i<nt; i++) {
        auto XFrow = XF.row(i);
        auto XFdotRow = XFdot.row(i);
        for(decltype(nt) j=0; j<nt; j++) {
            auto c = std::fabs(XFrow[j]);
            if (c>XFrowMax[i]) {
                XFrowMax[i] = c;
            }
            auto cdot = std::fabs(XFdotRow[j]);
            if (cdot>XFdotRowMax[i]) {
                XFdotRowMax[i] = cdot;
            }
        }
    }

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

bool HBNRSolver::preIteration(bool continuePrevious) {
    // Clear maximal residual contribution
    zero(maxResidualContribution_);

    // Zero historic values (we compute them at each iteration)
    zero(historicMaxResidualContribution_);
    zero(globalMaxResidualContribution_);

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
                ss.str(""); ss << maxResidual;
                Simulator::dbg() << ", worst residual=" << ss.str() << " @ " << (maxResidualNode ? maxResidualNode->name() : "(unknown)")
                                 << ", t" << maxResidualTimepointIndex << "=" << timepoints[maxResidualTimepointIndex];
            }
            if (iteration>1) {
                ss.str(""); ss << maxDelta;
                Simulator::dbg() << ", worst delta=" << ss.str() << " @ " << (maxDeltaNode ? maxDeltaNode->name() : "(unknown)")
                                 << ", f" << maxDeltaFreqIndex << "=" << spectrum[maxDeltaFreqIndex];
            }
        }
        Simulator::dbg() << "\n";
    }
    return true;
}

bool HBNRSolver::postIteration(bool continuePrevious) {
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
    // Resistive Jacobian is bound to 1-based subentry (1, 1) of each dense block. 
    // Reactive Jacobian is bound to 1-based subentry (1, 2) of each dense block. 
    // As Jacobian load offset goes from 0..nb-1 (nb=ntimepoints)
    // the first two columns are loaded with Jacobian values coresponding 
    // to timepoints because KLU matrices are stored in column major order. 
    // Let Jr_ijk and Jc_ijk denote the resistive and reactive Jacobian value 
    // from block with 1-based position (i+1, j+1) at timepoint with index k. 
    // i, j and k are all 0-based. 
    // After evalAndLoad() the first column of each block is filled Jr_ijk 
    // and the second column of each block is filled with Jc_ijk. k is the 
    // timepoint index and corresponds to the 0-based row index within the block. 
    // 
    // At a particular timepoint with index k the dense block row k is 
    // obtained as 
    //   Jr_ijk * XF_kl + Jc * XFdot_kl
    // where l and k are the 0-based row and column index within the block. 
    //
    // Let x_(i*nb+k) and xdot_(i*nb+k) denote the time-domain old solution 
    // and its derivative wrt time at timepoint t_k. 
    // 
    // The resistive residual at t_k for unknown i is added to delta_(i*nb+k)
    // immediately after evalAndLoad() for t_k. 
    // 
    // Because the first column of XF is 1 and the first column of XFdot is 0 
    // the first column already contains Jr_ijk which is what should be there 
    // anyway. We fill each dense block in the following way:
    //   loop l = nb-1..1 // loop over columns in reverse
    //     loop k = 0..nb-1 // loop over rows
    //       g = Jr_ijk = block_k0
    //       c = Jr_ijk = block_k1
    //       block_kl = g * XF_kl + c * XFdot_kl
    //       if l==1 
    //         // add reactive residual to delta (total residual)
    //         delta_(i*nb+k) += c * xdot_(j*nb+k)
    //       
    // The inner loop iterates j from 0 to nb-1. This way continuous vectors
    // (i.e. columns of dense blocks are processed in sequence) and good 
    // cache locality is obtained because dense blocks are stored as column
    // major order with spacings between consecutive columns due to KLU matrix 
    // structure. 
    // 
    // XF and XFdot are iterated column-wise in the inner loop. That is not 
    // optimal because they are stored as row major matrices. Fortunately 
    // they are fairly small so we expect them to be fully cached. 

    // TODO: When computing residual instead of 
    //       dq(x(t))/dt = dq/dx dx/dt use
    //       q(x(t)) -> freq domain -> compute derivative -> time domain = dq(x(t))/dt
    //       This will correctly handle elements that depend on time. 

    // Get sizes
    auto n = bsjac.nBlockRows();
    auto nb = bsjac.nBlockElementRows();

    // Transform old solution from frequency domain (x) to time domain (oldSolutionTd). 
    // Compute derivative of old solution wrt time to obtain oldSolutionTDdot. 
    // Work in blocks of length nb. 
    auto solFD = solution.data();
    auto solTD = oldSolutionTD.data();
    auto solTDDot = oldSolutionTDDot.data();
    auto solFDPtr = solFD;
    auto solTDPtr = solTD;
    auto solTDDotPtr = solTDDot;

    // Fill bucket with 0, skip bucket
    *solTDPtr = 0;
    *solTDDotPtr = 0;
    solTDPtr++;
    solTDDotPtr++;
    
    // Compute all unknowns, do not compute ground because it is not in the vector
    for(decltype(n) i=0; i<n; i++) {
        // Transform i-th unknown from frequency domain to time domain
        auto vvFD = VectorView<double>(solFDPtr, nb, 1);
        auto vvTD = VectorView<double>(solTDPtr, nb, 1);
        auto vvTDdot = VectorView<double>(solTDDotPtr, nb, 1);
        XF.vecMul(vvFD, vvTD);
        XFdot.vecMul(vvFD, vvTDdot);
        solFDPtr += nb;
        solTDPtr += nb;
        solTDDotPtr += nb;
    }

    // delta is zeroed at the beginning of each iteration by NRSolver
    
    // Loop through timepoints
    for(decltype(nb) k=0; k<nb; k++) {
        // Collect old solution vector at i-th timepoint
        // This vector is used by evalAndLoad()
        auto oldSolutionTDtkPtr = oldSolutionTDtk.data();
        
        // Ground node (bucket), set and skip
        *oldSolutionTDtkPtr = 0;
        oldSolutionTDtkPtr++;

        // Now, collect data for all other unknowns at t_k
        solTDPtr = solTD+1;
        for(decltype(n) l=0; l<n; l++) {
            *oldSolutionTDtkPtr = *solTDPtr;
            oldSolutionTDtkPtr++;
            solTDPtr += nb;
        }

        // Zero residual vector where evalAndLoad() will load the 
        // resistive residual at t_k
        zero(resistiveResidualAtTk);

        // Zero maximal residual contribution at timepoint
        zero(maxResidualContributionAtTk_);
        
        // Set time and offset
        evalSetup_.time = timepoints[k];
        loadSetup_.jacobianLoadOffset = k;
    
        // For i-th timepoint load 
        // - resistive and reactive Jacobian at t_k with offset i, 
        // - resistive residuals for all equations at t_k
        auto ok = evalAndLoadWrapper(evalSetup_, loadSetup_);

        // Put resistive residuals at t_k in residual vector (delta)
        for(decltype(n) i=0; i<n; i++) {
            // 1-based            1-based
            delta[1+i*nb+k] = resistiveResidualAtTk[1+i];
        }

        // Handle maximal resistive residual contribution at this timepoint
        // Loop through nodes
        for(decltype(n) i=0; i<n; i++) {
            auto c = maxResidualContributionAtTk_[1+i];
            // 1-based                              // 1-based
            maxResidualContribution_[1+i*nb+k] = maxResidualContributionAtTk_[1+i];
            
            // Update historic max for node
            if (c>historicMaxResidualContribution_[1+i]) {
                historicMaxResidualContribution_[1+i] = c;
            }

            // Update global historic max for nature (representative node index is 1-based)
            auto rn = circuit.reprNode(1+i);
            bool isPotential = ((rn->flags() & Node::Flags::PotentialNode) == Node::Flags::PotentialNode); 
            size_t ndx = isPotential ? 1 : 0;
            if (c>globalMaxResidualContribution_[ndx]) {
                globalMaxResidualContribution_[ndx] = c;
            }
        }
    }

    // For each block (ordered in column major order)
    auto oldSolutionTDtkPtr = oldSolutionTDtk.data();
    for(auto& pos : circuit.sparsityMap().positions()) {
        // Get dense block
        auto [block, found] = bsjac.block(pos);

        // Get block position, make position 0-based
        auto [i, j] = pos;
        i--;
        j--;

        // Get g and c columns from block
        auto gCol = block.column(0);
        auto cCol = block.column(1);
        
        // Scan columns in block in reverse from nb-1 to 1 
        // l=0 is for DC which should hold Jr_ijk*1 + Jc_ijk*0 
        // and is already filled with correct values (Jr_ijk). 
        for(decltype(nb) l=nb-1; l>0; l--) {   
            // Get target column
            auto targetColumn = block.column(l);

            // Get XF and XFdot columns
            auto XFCol = XF.column(l);
            auto XFdotCol = XFdot.column(l);

            // Get solution derivative for unknown j in time-domain 
            auto solDot = VectorView<double>(solTDDot+1+j*nb, nb, 1);

            // Get delta (residual) subvector for equation i
            auto deltaSub = VectorView<double>(delta.data()+1+i*nb, nb, 1);

            // When l==1, update max residual contribution with reactive cotribution
            if (l==1) {
                // Scan rows in block
                for(decltype(nb) k=0; k<nb; k++) {
                    auto c = cCol[k]*solDot[k];
                    if (c>maxResidualContribution_[1+i*nb+k]) {
                        maxResidualContribution_[1+i*nb+k] = c;
                    }
                }
            }

            // Scan rows in block
            for(decltype(nb) k=0; k<nb; k++) {
                // Get resistive and reactive Jacobian of block (i, j) at t_k
                // before gCol is filled with contributions when l=1. 
                auto g = gCol[k]; 
                auto c = cCol[k]; 
                
                // block_kl = Jc_ijk*XFdot_kl
                targetColumn[k] = c * XFdotCol[k];

                // block_kl += Jr_ijk*XF_kl
                targetColumn[k] += g * XFCol[k];

                // Only once in the end (when l=1) add reactive residual contribution to delta 
                if (l==1) {
                    // Jc_ijk*xdot_(j*nb+k)
                    // 1-based           1-based
                    deltaSub[k] += c*solDot[k];
                }
            }
        }
    }

    // OK, do not prevent convergence
    return std::make_tuple(true, false); 
}

std::tuple<bool, bool> HBNRSolver::checkResidual() {
    // Compute norms only in debug mode
    bool computeNorms = settings.debug;

    // In residual we have the residual at previous solution
    // We are going to check that residual
    
    // Number of unknowns (vector length includes a bucket at index 0)
    auto n = circuit.unknownCount();

    // Number of timepoints
    auto nt = timepoints.size();

    // Results
    double maxResidual = 0.0;
    double maxNormResidual = 0.0;
    double l2normResidual2 = 0.0;
    Node* maxResidualNode = nullptr;
    
    // Assume residual is OK
    bool residualOk = true;
    
    // Get point maximum for each residual nature
    zero(pointMaxResidualContribution_); 
    // Loop through all nodes
    for(decltype(n) i=0; i<n; i++) {
        // Get representative node (1-based index) and nature index
        auto rn = circuit.reprNode(i+1);
        bool isPotential = ((rn->flags() & Node::Flags::PotentialNode) == Node::Flags::PotentialNode); 
        size_t ndx = isPotential ? 1 : 0;
        // Loop through all timepoints
        for(decltype(nt) k=0; k<nt; k++) {
            double c = std::fabs(maxResidualContribution_[1+i*nt+k]);
            if (c>pointMaxResidualContribution_[ndx]) {
                pointMaxResidualContribution_[ndx] = c;
            }
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
            double tolref = std::fabs(maxResidualContribution_[1+(i-1)*nt+k]);
            
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
    auto nf = spectrum.size();

    maxDelta = 0.0;
    maxNormDelta = 0.0;
    Node* maxDeltaNode = nullptr;
    maxDeltaFreqIndex = 0;
    
    // Check convergence (see if delta is small enough), 
    // but only if this is iteration 2 or later
    // In iteration 1 assume we did not converge
    
    // Assume we converged
    bool deltaOk = true;
    
    // Get point maximum for each solution nature
    auto compPtr = solution.data();
    // Skip bucket
    compPtr++;
    for(decltype(n) i=1; i<=n; i++) {
        // Representative node, associated potential nature index
        auto rn = circuit.reprNode(i);
        bool isPotential = ((rn->flags() & Node::Flags::PotentialNode) == Node::Flags::PotentialNode); 
        size_t ndx = isPotential ? 0 : 1;
        for(decltype(nf) j=0; j<nf; j++) {
            double c;
            if (j==0) {
                // DC
                c = std::fabs(*compPtr);
                // Advance by 1
                compPtr++;
            } else {
                // Not DC, it is complex
                c = std::abs(*reinterpret_cast<Complex*>(compPtr));
                // Advance by 2 (real, imaginary)
                compPtr += 2;
            }
            // Rows are natures, columns are frequency components (DC, f1, f2, ...)
            if (c>pointMaxSolution_.at(ndx, j)) {
                pointMaxSolution_.at(ndx, j) = c;
            }
        }
    }

    // Use 1-based index (with bucket) because same indexing is used for variables
    compPtr = solution.data();
    auto deltaPtr = delta.data();
    // Skip bucket
    compPtr++;
    deltaPtr++;
    for(decltype(n) i=1; i<=n; i++) {
        // Representative node, associated potential nature index
        auto rn = circuit.reprNode(i);
        bool isPotential = ((rn->flags() & Node::Flags::PotentialNode) == Node::Flags::PotentialNode); 
        size_t ndx = isPotential ? 0 : 1;
        for(decltype(nf) j=0; j<nf; j++) {
            // Compute tolerance reference
            // Point local reference by default
            // Compute tolerance reference, start with previous value of the i-th unknown at j-th frequency
            double tolref;
            if (j==0) {
                // DC
                tolref = std::fabs(*compPtr);
                compPtr++;
            } else {
                // Complex
                tolref = std::abs(*reinterpret_cast<Complex*>(compPtr));
                compPtr += 2;
            }
            
            // Account for global references, no historic reference because we are in frequency domain
            if (globalSolRef) {
                // Point global reference, ndx is the nature index
                tolref = std::max(tolref, pointMaxSolution_.at(ndx, j));
            }
            
            // Compute tolerance
            double tol = circuit.solutionTolerance(rn, tolref);

            // Absolute solution change 
            double deltaAbs;
            if (j==0) {
                // DC
                deltaAbs = fabs(*deltaPtr);
                deltaPtr++;
            } else {
                // Complex
                tolref = std::abs(*reinterpret_cast<Complex*>(deltaPtr));
                deltaPtr += 2;
            }

            if (computeNorms) {
                double normDelta = deltaAbs/tol;
                if (i==1 || normDelta>maxNormDelta) {
                    maxDelta = deltaAbs;
                    maxNormDelta = normDelta;
                    maxDeltaNode = rn;
                    maxDeltaFreqIndex = j;
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


}

