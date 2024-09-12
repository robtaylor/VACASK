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

    // Allocate space in vectors
    auto n = circuit.unknownCount();
    auto nf = spectrum.size();
    auto nt = timepoints.size();
    auto ncomp = nf*2-1;

    // Old states at one timepoint (dummy vector of zeros because we do no limiting)
    dummyStatesRepo.upsize(1, circuit.statesCount());
    dummyStatesRepo.zero();
    
    // Old solution and derivative wrt time at all timepoints
    // Bucket for ground node is needed because these vectors 
    // are used by NRSolver which assumes a bucket of length 1. 
    oldSolutionTD.upsize(1, n*nt+1);
    oldSolutionTDDot.upsize(1, n*nt+1);

    // Old solution and resistive residual at one timepoint
    // Includes ground node because it is used by evalAndLoad()
    oldSolutionTDtk.upsize(1, n+1);
    resistiveResidualTk.resize(n+1);

    // Jacobian is sized is core or analysis.
    // solution is sized in core. 
    // delta is resized by NRSolver::rebuild() based on Jacobian size. 

    return true;
}

bool HBNRSolver::initialize(bool continuePrevious) {
    return true;
}

bool HBNRSolver::preIteration(bool continuePrevious) {
    return true;
}

bool HBNRSolver::postSolve(bool continuePrevious) {
    return true;
}

bool HBNRSolver::postConvergenceCheck(bool continuePrevious) {
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

    // Get sizes
    auto n = bsjac.blocksInColumn();
    auto nb = bsjac.blockRows();

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

        // Clear residual vector where evalAndLoad() will load the 
        // resistive residual at t_k
        resistiveResidualTk.clear();

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
            delta[1+i*nb+k] = resistiveResidualTk[1+i];
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
        
        // Scan columns in reverse from nb-1 to 1 
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
            auto deltaSub =VectorView<double>(delta.data()+1+i*nb, nb, 1);

            // Scan rows
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

}

