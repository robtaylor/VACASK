#include "corehbnr.h"
#include "simulator.h"
#include "densematrix.h"
#include "common.h"

namespace NAMESPACE {

HBNRSolver::HBNRSolver(
        Circuit& circuit, 
        CommonData& commons, 
        KluBlockSparseRealMatrix& jacColoc, 
        KluBlockSparseRealMatrix& bsjac, 
        VectorRepository<double>& solution, 
        Vector<Complex>& solutionFD, 
        Vector<Real>& frequencies, 
        Vector<Real>& timepoints, 
        DenseMatrix<Real>& DDT, 
        DenseMatrix<Real>& DDTcolMajor, 
        DenseMatrix<double>& APFT, 
        DenseMatrix<double>& IAPFT, 
        NRSettings& settings
) : circuit(circuit), commons(commons), jacColoc(jacColoc), bsjac(bsjac), solutionFD(solutionFD), 
    frequencies(frequencies), timepoints(timepoints), DDT(DDT), DDTcolMajor(DDTcolMajor), 
    APFT(APFT), IAPFT(IAPFT), 
    NRSolver(circuit.tables().accounting(), bsjac, solution, settings) {
    // Slot 0 is for sweep continuation and homotopy (set via CoreStateStorage object)
    // Slot 1 is for nodesets that are read from stored results. 
    resizeForces(2);

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

bool HBNRSolver::setForces(Int ndx, const AnnotatedSolution& storedSolution, bool abortOnError) {
    // Get forces
    auto& f = forces(ndx);

    // Clear forced values, set number of forces to 0
    f.clear();
    
    // Number of unknowns
    auto n = circuit.unknownCount();

    // Number of components per unknown
    auto nf = frequencies.size(); // number of frequencies per unknown
    auto blockSize = 2*nf-1; // number of timepoints per unknown

    // Make space for variable forces (also include bucket)
    f.unknownValue_.resize(n*blockSize+1);
    // By default turn off all forces
    f.unknownForced_.resize(n*blockSize+1, false);
    
    // Number of frequencies in solution and solver
    auto nfSolution = storedSolution.auxData().size();
    auto nfSolver = frequencies.size();

    // Prepare frequency translator between solution and solver
    // Translator stores the solver frequency index for each solutiuon frequency index
    // Assume no frequency can be translated into solution frequency (negative index)
    std::vector<int> xlat(nfSolver, -1);
    
    // DC can be translated as new 0 -> old 0 (it is always present)
    if (nfSolver>0 && nfSolution>0) {
        xlat[0] = 0;
    }

    // Translate the rest
    decltype(nfSolver) ndxSolver = 1;
    decltype(nfSolution) ndxSolution = 1;
    for(; ndxSolver<nfSolver && ndxSolution<nfSolution;) {
        auto fSolver = frequencies[ndxSolver];
        auto fSolution = storedSolution.auxData()[ndxSolution];
        if (std::abs(fSolver-fSolution)<=std::max(std::abs(fSolver), std::abs(fSolution))*1e-14) {
            // Frequencies are almost the same, store translator
            xlat[ndxSolver] = ndxSolution;
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

    // Go through annotated solution. fill APFT spectrum 
    // Use resistive residual vector for APFT spectrum
    auto& forcesFD = resistiveResidual;
    // Resistive residual can hold n APFT spectra, we need only the first one
    // Resize it, just in case
    forcesFD.resize(n*blockSize);

    // Solution spectrum
    auto& solSpec = storedSolution.cxValues();
    auto& solNames = storedSolution.names();

    // Check if we have solution name annotations
    // with matching length. 
    bool checkNames;
    if (solNames.size()==n+1) {
        // Yes, check names
        checkNames = true;
    } else if (solNames.size()==0 && solSpec.size()==nfSolver*n) {
        // No annotations, solutions vector has correct length
        checkNames = false;
    } else {
        // Cannot apply stored solution, no names nor matching length vector
        lastHBNRError = HBNRSolverError::ForcesError;
        return false;
    }

    // Go through all unknowns. 
    for(decltype(n) i=1; i<=n; i++) {
        Node* node;
        if (checkNames) {
            // Stored solution has name annotations, get node by the name from the solution
            node = circuit.findNode(solNames[i]);
        } else {
            // Stored solution is coherent, simply take the representative node of i-th unknown
            node = circuit.reprNode(i);
        }
        if (!node) {
            // Node not found. No forces will be applied to this unknown. 
            continue;
        }
        // Copy spectrum for one node
        auto ui = node->unknownIndex();
        // Origin index in complex spectrum vector (no bucket)
        auto srcOrigin = (ui-1)*nf;
        // Origin index in destination vector of TD values (no bucket)
        auto destOrigin = (ui-1)*blockSize;
        
        // Copy DC (one real value)
        forcesFD[0] = solSpec[srcOrigin].real();
        // Scan all nonzero frequencies of solver's spectrum
        for(decltype(nf) k=1; k<nf; k++) {
            // Translate solver frequency into solution frequency
            auto xlf = xlat[k];
            // Index of real component (DC is stored as a single real number)
            auto ndx = 1+2*(k-1);
            if (xlf>=0) {
                // Translation exists, copy solution component
                forcesFD[ndx] = solSpec[srcOrigin+k].real();
                forcesFD[ndx+1] = solSpec[srcOrigin+k].imag();
            } else {
                // No translation, fill with zeros
                forcesFD[ndx] = 0;
                forcesFD[ndx+1] = 0;
            }
        }
        // Inverse APFT, store in forces vector
        auto fd = VectorView<double>(forcesFD.data(), blockSize, 1);
        auto td = VectorView<double>(f.unknownValue_.data()+1+destOrigin, blockSize, 1);
        IAPFT.multiply(fd, td);
        // After IAPFT the resulting timepoints are all valid forces, even if not all spectral components were copied
        // Mark all forces for this unknown as set. 
        for(decltype(nf) k=0; k<blockSize; k++) {
            f.unknownForced_[1+destOrigin+k] = true;
        }
    }

    // std::cout << "Set forces:\n";
    // f.dump(circuit, std::cout);
    
    // Ignore errors (conflicting forces are overwritten by newer value)
    // Error checking makes sense in case of manual forces (nodeset, ic). 
    // Therefore we ignore abortOnError. 
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

    // Get diagonal pointers for forces
    auto n = circuit.unknownCount();
    auto nt = timepoints.size();
    diagPtrs.resize(n*nt+1);
    
    // Bind diagonal matrix elements
    // Needed for forcing unknown values
    for(decltype(n) i=0; i<n; i++) {
        for(decltype(nt) j=0; j<nt; j++) {
            // We know the matrix type so we can use the elementPtr() non-virtual function
            diagPtrs[1+i*nt+j] = bsjac.elementPtr(MatrixEntryPosition(i+1, i+1), Component::Real, MatrixEntryPosition(j, j));
        }
    }

    return true;
}

bool HBNRSolver::initialize(bool continuePrevious) {
    // Clear HB NR solver error
    clearError();

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
        // If converged, convert solution from TD to FD, store as complex spectrum
        auto n = circuit.unknownCount();
        auto nf = frequencies.size();
        auto nt = timepoints.size();
        solutionFD.resize(n*nf); // no bucket

        for(decltype(n) i=0; i<n; i++) {
            auto cxSpecPtr = solutionFD.data()+nf*i;
            // APFT computes spectrum as complex values, with the exception of DC which is stored as a real value. 
            // We write APFT output starting at the imaginary part of the DC complex magnitude. 
            // This way all complex values will be in the right place, except for the DC value which 
            // will be placed in the DC solution's imaginary part. 
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
    if (!circuit.evalAndLoad(commons, &evalSetup, &loadSetup, nullptr)) {
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

    // Remove forces originating from nodesets after nsiter iterations
    auto nsiter = circuit.simulatorOptions().core().op_nsiter;
    // Do this only at nsiter+1 (first iteration has index 1)
    if (iteration==nsiter+1) {
        // Continuation nodesets
        enableForces(0, false);
        // User-specified nodesets
        enableForces(1, false);
    }

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

    // Add forced values to the system
    if (haveForces() && !loadForces(true)) {
        if (settings.debug) {
            Simulator::dbg() << "Failed to load forced values at iteration " << iteration << "\n";
        }
        lastHBNRError = HBNRSolverError::LoadForces;
        errorIteration = iteration;
        std::make_tuple(false, evalSetup_.limitingApplied);
    }


    // OK, do not prevent convergence
    return std::make_tuple(true, false); 
}

bool HBNRSolver::loadForces(bool loadJacobian) {
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
    }

    return true;
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
        // Skip internal device nodes
        if (rn->checkFlags(Node::Flags::InternalDeviceNode)) {
            continue;
        }
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
            double tolref = xold[1+(i-1)*nt+j];
            
            // Account for global references, no historic reference because we are in frequency domain
            if (globalSolRef) {
                // Point global reference, ndx is the nature index
                tolref = std::max(tolref, pointMaxSolution_.at(ndx, j));
            }
            
            // Compute tolerance
            double tol = circuit.solutionTolerance(rn, tolref);

            // Absolute solution change 
            double deltaAbs = std::fabs(xdelta[1+(i-1)*nt+j]);;

            if (computeNorms) {
                double normDelta = deltaAbs/tol;
                if ((i==1 && j==0) || normDelta>maxNormDelta) {
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

bool HBNRSolver::formatError(Status& s, NameResolver* resolver) const {
    // Error in NRSolver
    if (lastError!=NRSolver::Error::OK) {
        NRSolver::formatError(s, resolver);
        return false;
    }

    switch (lastHBNRError) {
        case HBNRSolverError::ForcesError:
            s.set(Status::Force, "Failed to apply forces.");
            return false;            
        case HBNRSolverError::LoadForces:
            s.set(Status::Force, "Failed to load forces.");
            return false;
        default:
            return true;
    }
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
