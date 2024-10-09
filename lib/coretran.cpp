#include "coretran.h"
#include "an.h"
#include "simulator.h"
#include "common.h"
#include <filesystem>
#include <algorithm>

namespace NAMESPACE {

// Application of general implicit integration formulae
// 
// q(t_k)       ... terminal's reactive residual at t_k
// qdot(t_k)    ... time derivative of terminal's reactive residual at t_k (current at t_k)
// x_k          ... circuit solution at t_k
// h_k = t_{k+1} - t_k .. timestep to new timepoint
// 
// Value of reactive residual's time derivative at t_{k+1}
//   qdot(t_{k+1}) = 
//       1 / (h_k b_{-1}) * q(t_{k+1})                 // nonlinear function of new solution
//     - sum_{i=0}^{numX-1} a_i/(h_k b_{-1}) * q_{k-i} // past values of node's reactive residual
//     - sum_{i=0}^{numXdot-1} b_i/b_{-1} * qdot_{k-i} // past values of node's reactive residual's time derivative
//
// Superscripts
//   j   .. solution of previous NR iteration (known)
//   j+1 .. solution of current NR iteration (to be solved for)
// 
// Contribution of reactive residual of a terminal at t_{k+1} in iteration j+1 of NR algorithm
//   i_{k+1}^{j+1} = 
//       1 / (h_k b_{-1}) * q(x_{k+1}^{j+1})           // nonlinear function of new solution
//     - sum_{i=0}^{numX-1} a_i/(h_k b_{-1}) * q_{k-i} // past values of node's reactive residual
//     - sum_{i=0}^{numXdot-1} b_i/b_{-1} * qdot_{k-i} // past values of node's reactive residual's time derivative
//
// di |                   1      dq |
// -- |             = ---------- -- |               reactive Jacobian contribution to matrix of coefficients
// dx | t_{k+1}, j    h_k b_{-1} dx | t_{k+1}, j
//                    \________/ \_____________/
//                         |            |
//                         |      reactive Jacobian
//                         |     
//                     alpha in load_jacobian_react() and load_jacobian_tran()


/*
// Solution of I-L circuit, sinusoidal excitation current, V0=0.5, f=50
double globalDbg[][2] = { // time v
{ 0.000000e+00,	 0.000000e+00 },	
{ 2.500000e-05,	 1.570780e-01 },	
{ 2.560286e-05,	 1.570747e-01 },	
{ 2.615840e-05,	 1.570744e-01 },	
{ 2.726947e-05,	 1.570738e-01 },	
{ 2.949161e-05,	 1.570730e-01 },	
{ 3.393589e-05,	 1.570706e-01 },	
{ 4.282445e-05,	 1.570657e-01 },	
{ 6.060157e-05,	 1.570517e-01 },	
{ 9.615581e-05,	 1.570107e-01 },	
{ 1.672643e-04,	 1.568731e-01 },	
{ 3.094813e-04,	 1.563796e-01 },	
{ 5.094813e-04,	 1.551326e-01 },	
{ 7.094813e-04,	 1.532346e-01 },	
{ 9.094813e-04,	 1.507706e-01 },	
{ 1.109481e-03,	 1.476730e-01 },	
{ 1.309481e-03,	 1.440312e-01 },	
{ 1.509481e-03,	 1.397823e-01 },	
{ 1.709481e-03,	 1.350204e-01 },	
{ 1.909481e-03,	 1.296869e-01 },	
{ 2.109481e-03,	 1.238804e-01 },	
{ 2.309481e-03,	 1.175462e-01 },	
{ 2.509481e-03,	 1.107869e-01 },	
{ 2.709481e-03,	 1.035516e-01 },	
{ 2.909481e-03,	 9.594633e-02 },	
{ 3.109481e-03,	 8.792374e-02 },	
{ 3.309481e-03,	 7.959282e-02 },	
{ 3.509481e-03,	 7.090912e-02 },	
{ 3.709481e-03,	 6.198423e-02 },	
{ 3.909481e-03,	 5.277607e-02 },	
{ 4.109481e-03,	 4.339827e-02 },	
{ 4.309481e-03,	 3.381055e-02 },	
{ 4.509481e-03,	 2.412805e-02 },	
{ 4.709481e-03,	 1.431167e-02 },	
{ 4.909481e-03,	 4.477461e-03 },	
{ 5.109481e-03,	-5.413069e-03 },	
{ 5.309481e-03,	-1.524358e-02 },	
{ 5.509481e-03,	-2.505259e-02 },	
{ 5.709481e-03,	-3.472408e-02 },	
{ 5.909481e-03,	-4.429717e-02 },	
{ 6.109481e-03,	-5.365680e-02 },	
{ 6.309481e-03,	-6.284331e-02 },	
{ 6.509481e-03,	-7.174316e-02 },	
{ 6.709481e-03,	-8.039853e-02 },	
{ 6.909481e-03,	-8.869795e-02 },	
{ 7.109481e-03,	-9.668597e-02 },	
{ 7.309481e-03,	-1.042538e-01 },	
{ 7.509481e-03,	-1.114488e-01 },	
{ 7.709481e-03,	-1.181653e-01 },	
{ 7.909481e-03,	-1.244541e-01 },	
{ 8.109481e-03,	-1.302131e-01 },	
{ 8.309481e-03,	-1.354969e-01 },	
{ 8.509481e-03,	-1.402072e-01 },	
{ 8.709481e-03,	-1.444029e-01 },	
{ 8.909481e-03,	-1.479901e-01 },	
{ 9.109481e-03,	-1.510318e-01 },	
{ 9.309481e-03,	-1.534388e-01 },	
{ 9.509481e-03,	-1.552790e-01 },	
{ 9.709481e-03,	-1.564677e-01 },	
{ 9.909481e-03,	-1.570775e-01 },	
{ 1.010948e-02,	-1.570287e-01 },	
{ 1.030948e-02,	-1.563989e-01 },	
{ 1.050948e-02,	-1.551132e-01 },	
{ 1.070948e-02,	-1.532540e-01 },	
{ 1.090948e-02,	-1.507513e-01 },	
{ 1.110948e-02,	-1.476923e-01 },	
{ 1.130948e-02,	-1.440118e-01 },	
{ 1.150948e-02,	-1.398016e-01 },	
{ 1.170948e-02,	-1.350010e-01 },	
{ 1.190948e-02,	-1.297063e-01 },	
{ 1.210948e-02,	-1.238610e-01 },	
{ 1.230948e-02,	-1.175656e-01 },	
{ 1.250948e-02,	-1.107675e-01 },	
{ 1.270948e-02,	-1.035709e-01 },	
{ 1.290948e-02,	-9.592698e-02 },	
{ 1.310948e-02,	-8.794309e-02 },	
{ 1.330948e-02,	-7.957347e-02 },	
{ 1.350948e-02,	-7.092847e-02 },	
{ 1.370948e-02,	-6.196489e-02 },	
{ 1.390948e-02,	-5.279541e-02 },	
{ 1.410948e-02,	-4.337893e-02 },	
{ 1.430948e-02,	-3.382990e-02 },	
{ 1.450948e-02,	-2.410870e-02 },	
{ 1.470948e-02,	-1.433101e-02 },	
{ 1.490948e-02,	-4.458116e-03 },	
{ 1.510948e-02,	 5.393724e-03 },	
{ 1.530948e-02,	 1.526293e-02 },	
{ 1.550948e-02,	 2.503325e-02 },	
{ 1.570948e-02,	 3.474342e-02 },	
{ 1.590948e-02,	 4.427783e-02 },	
{ 1.610948e-02,	 5.367614e-02 },	
{ 1.630948e-02,	 6.282397e-02 },	
{ 1.650948e-02,	 7.176251e-02 },	
{ 1.670948e-02,	 8.037918e-02 },	
{ 1.690948e-02,	 8.871729e-02 },	
{ 1.710948e-02,	 9.666662e-02 },	
{ 1.730948e-02,	 1.042731e-01 },	
{ 1.750948e-02,	 1.114294e-01 },	
{ 1.770948e-02,	 1.181846e-01 },	
{ 1.790948e-02,	 1.244347e-01 },	
{ 1.810948e-02,	 1.302324e-01 },	
{ 1.830948e-02,	 1.354775e-01 },	
{ 1.850948e-02,	 1.402266e-01 },	
{ 1.870948e-02,	 1.443836e-01 },	
{ 1.890948e-02,	 1.480094e-01 },	
{ 1.910948e-02,	 1.510125e-01 },	
{ 1.930948e-02,	 1.534582e-01 },	
{ 1.950948e-02,	 1.552596e-01 },	
{ 1.970948e-02,	 1.564870e-01 },	
{ 1.990948e-02,	 1.570581e-01 },	
{ 2.010948e-02,	 1.570481e-01 },	
{ 2.030948e-02,	 1.563796e-01 },	
{ 2.050948e-02,	 1.551326e-01 },	
{ 2.070948e-02,	 1.532346e-01 },	
{ 2.090948e-02,	 1.507706e-01 },	
{ 2.110948e-02,	 1.476730e-01 },	
{ 2.130948e-02,	 1.440312e-01 },	
{ 2.150948e-02,	 1.397823e-01 },	
{ 2.170948e-02,	 1.350204e-01 },	
{ 2.190948e-02,	 1.296869e-01 },	
{ 2.210948e-02,	 1.238804e-01 },	
{ 2.230948e-02,	 1.175462e-01 },	
{ 2.250948e-02,	 1.107869e-01 },	
{ 2.270948e-02,	 1.035516e-01 },	
{ 2.290948e-02,	 9.594633e-02 },	
{ 2.310948e-02,	 8.792374e-02 },	
{ 2.330948e-02,	 7.959282e-02 },	
{ 2.350948e-02,	 7.090912e-02 },	
{ 2.370948e-02,	 6.198423e-02 },	
{ 2.390948e-02,	 5.277607e-02 },	
{ 2.410948e-02,	 4.339827e-02 },	
{ 2.430948e-02,	 3.381055e-02 },	
{ 2.450948e-02,	 2.412805e-02 },	
{ 2.470948e-02,	 1.431167e-02 },	
{ 2.490948e-02,	 4.477461e-03 },	
{ 2.510948e-02,	-5.413069e-03 },	
{ 2.530948e-02,	-1.524358e-02 },	
{ 2.550948e-02,	-2.505259e-02 },	
{ 2.570948e-02,	-3.472408e-02 },	
{ 2.590948e-02,	-4.429717e-02 },	
{ 2.610948e-02,	-5.365680e-02 },	
{ 2.630948e-02,	-6.284331e-02 },	
{ 2.650948e-02,	-7.174316e-02 },	
{ 2.670948e-02,	-8.039853e-02 },	
{ 2.690948e-02,	-8.869795e-02 },	
{ 2.710948e-02,	-9.668597e-02 },	
{ 2.730948e-02,	-1.042538e-01 },	
{ 2.750948e-02,	-1.114488e-01 },	
{ 2.770948e-02,	-1.181653e-01 },	
{ 2.790948e-02,	-1.244541e-01 },	
{ 2.810948e-02,	-1.302131e-01 },	
{ 2.830948e-02,	-1.354969e-01 },	
{ 2.850948e-02,	-1.402072e-01 },	
{ 2.870948e-02,	-1.444029e-01 },	
{ 2.890948e-02,	-1.479901e-01 },	
{ 2.910948e-02,	-1.510318e-01 },	
{ 2.930948e-02,	-1.534388e-01 },	
{ 2.950948e-02,	-1.552790e-01 },	
{ 2.970948e-02,	-1.564677e-01 },	
{ 2.990948e-02,	-1.570775e-01 },	
{ 3.000000e-02,	-1.570395e-01 },	
};
int nGlobalDbg = sizeof(globalDbg)/(sizeof(double)*2);
*/

Id TranCore::icmodeOp = Id::createStatic("op");
Id TranCore::icmodeUic = Id::createStatic("uic"); 

TranParameters::TranParameters() {
    // Turn off output for op analysis
    opParams.writeOutput = 0; 
    icmode = TranCore::icmodeOp;
}

template<> int Introspection<TranParameters>::setup() {
    registerMember(step);
    registerMember(stop);
    registerMember(start);
    registerMember(maxstep);
    registerMember(icmode);
    registerNamedMember(opParams.nodeset, "nodeset");
    registerMember(ic);
    registerMember(store);
    
    return 0;
}
instantiateIntrospection(TranParameters);


TranCore::TranCore(
    OutputDescriptorResolver& parentResolver, TranParameters& params, OperatingPointCore& opCore, 
    Circuit& circuit, 
    KluRealMatrix& jacobian, VectorRepository<double>& solution, VectorRepository<double>& states
) : AnalysisCore(parentResolver, circuit), params(params), outfile(nullptr), opCore_(opCore), 
    jacobian(jacobian), solution(solution), states(states), 
    nrSolver(circuit, jacobian, states, solution, nrSettings, integCoeffs) { 
    // Slots 0 (current) and -1 (future) are used for the NR solver
    // Slots 1, 2, ... correspond to past values (at t_{k}, t_{k-1}, ...)
    // Therefore historyOffset needs to be set to 1 when calling 
    // preparePredictorHistory() and prepareDifferentiatorHistory(). 
    
    // Make another forces slot in opCore_'s solver
    opCore_.solver().resizeForces(3);

    // Set analysis type for the initial operating point analysis
    auto& esSystem = opCore_.solver().evalSetup();
    esSystem.staticAnalysis = true;
    esSystem.dcAnalysis = false;
    esSystem.tranAnalysis = true;
}

TranCore::~TranCore() {
    delete outfile;
}

bool TranCore::addDefaultOutputDescriptors() {
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    if (savesCount==0) {
        return addAllUnknowns(PTSave(Loc::bad, "default", Id(), Id()));
    }
    return true;
}

bool TranCore::resolveOutputDescriptors(bool strict) {
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
        case OutdTime:
            outputSources.emplace_back(&(nrSolver.evalSetup().time));
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

bool TranCore::addCoreOutputDescriptors() {
    clearError();
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    
    if (!addOutputDescriptor(OutputDescriptor(OutdTime, "time"))) {
        lastError = Error::Descriptor;
        errorId = "time";
        return false;
    }
    return true;
}

std::tuple<bool, bool> TranCore::preMapping(Status& s) {
    // Go through nodesets. Decode, set, and check delta part. 
    // No need to check unknowns part because diagonal entries are
    // always allocated by Circuit. 
    auto& icParam = params.ic;
    if (icParam.type()==Value::Type::String) {
        // It is a string. 
        // Will retrieve DC solution and set it as nodeset. 
        // Nothing to do here because DC solutions do not introduce delta foces. 
        return std::make_tuple(true, false);
    } else if (icParam.type()==Value::Type::ValueVec) {
        // It is a list
        // Retrieve it and check if we need extra matrix entries
        auto& valVec = icParam.val<ValueVector>();
        auto [ok, needsMapping] = preprocessedIc.set(circuit, valVec, s);
        if (!ok) {
            s.extend("Failed to preprocess initial conditions.");
        }
        return std::make_tuple(ok, needsMapping);
    } else {
        // Unsupported type, signal an error.
        s.set(Status::BadArguments, "Initial conditions must be a list or a string.");
        return std::make_tuple(false, false);
    }
}

bool TranCore::populateStructures(Status& s) {
    // Go through node pairs and add entries to sparsity map
    for(auto& pair : preprocessedIc.nodePairs) {
        auto [node1, node2] = pair;
        if (!node1 | !node2) {
            // One of the two nodes was not found, ignore pair
            continue;
        }
        
        // IC forces are all resistive
        if (auto [_, ok] = circuit.createJacobianEntry(node1, node2, EntryFlags::Resistive, s); !ok) {
            return false;
        }
        
        if (auto [_, ok] = circuit.createJacobianEntry(node2, node1, EntryFlags::Resistive, s); ok) {
            return false;
        }
    }

    return true;
}

bool TranCore::rebuild(Status& s) {
    // We are using the same Jacobian as operating point analysis
    
    // Bind Jacobian entries
    // OperatingPointCore has bound the resistive part of the Jacobian
    // Let's bind the reactive part 
    if (!circuit.bind(nullptr, Component::Real, std::nullopt, &jacobian, Component::Real, std::nullopt, s)) {
        return false;
    }
    
    // Prepare NR solver settings
    auto& options = circuit.simulatorOptions().core();
    nrSettings = NRSettings {
        .debug = options.nr_debug, 
        .itlim = options.tran_itl, 
        .itlimCont = options.tran_itl, 
        .convIter = options.nr_conviter, 
        .residualCheck = bool(options.nr_residualcheck),  
        .dampingFactor = options.nr_damping, 
        .matrixCheck = bool(options.matrixcheck), 
        .rhsCheck = bool(options.rhscheck), 
        .solutionCheck = bool(options.solutioncheck), 
        .forceFactor = options.nr_force, 
    };

    // Create Forces from preprocessed initial conditions in opCore's 
    // NR solver force slot 2. 
    // This has to be done now for forces that include delta forces
    // (e.g. user nodesets and ics) because nrSolver.rebuild() 
    // (which is called before nrSolver.run()) requires 
    // nrSolver.forcesList[].deltaIndices() to be set up. 
    // Setting up is done by nrSolver.forces().set(). 
    // We also handle ordinary initial conditions from solution repository. 
    // The enclosing Tran class must first call TranCore::rebuild()
    // and after that OperatingPointCore::rebuild() because the 
    // latter one requires slot 2 to be populated in order to gather 
    // pointers to matrix entries. 
    auto strictforce = circuit.simulatorOptions().core().strictforce; 
    auto& icParam = params.ic;
    if (icParam.type()==Value::Type::String) {
        String& solutionName = icParam.val<String>();
        if (solutionName.length()>0) {
            // Get solution from repository
            auto solPtr = circuit.storedSolution("dc", solutionName);    
            if (!solPtr) {
                // No initial conditions 
                opCore_.solver().forces(2).clear();
                if (params.icmode==icmodeOp) {
                    Simulator::wrn() << "Warning, solution '"+solutionName+"' not found. No initial conditions applied.\n";
                } else {
                    Simulator::wrn() << "Warning, solution '"+solutionName+"' not found. Zero initial conditions applied.\n";
                }
            } else {
                // Initial conditions from solution repository
                // if (!opCore_.solver().forces(2).set(circuit, *solPtr, strictforce)) {
                if (!opCore_.solver().setForces(2, *solPtr, strictforce)) {
                    // Abort on error if strictforce is set
                    if (strictforce) {
                        opCore_.solver().forces(2).formatError(s);
                        return false;
                    }
                }
            }
        } else {
            // No initial conditions, clear slot
            opCore_.solver().forces(2).clear();
        }
    } else if (icParam.type()==Value::Type::ValueVec) {
        // A list with possibly delta forces
        // Set slot 2, use op mode for setting up initial conditions
        bool uicMode = false;
        if (!opCore_.solver().forces(2).set(circuit, preprocessedIc, uicMode, strictforce)) {
            // Abort on error if strictforce is set
            if (strictforce) {
                opCore_.solver().forces(2).formatError(s);
                return false;
            }
        }
    } else {
        // Error
        s.set(Status::BadArguments, "Initial conditions must be a list or a string.");
        return false;
    }

    // Rebuild transient NR solver structures
    if (!nrSolver.rebuild()) {
        s.set(Status::NonlinearSolver, "Failed to rebuild internal structures of nonlinear solver.");
        return false;
    }

    // opCore_'s solver is rebuilt later (done by enclosing class). 

    return true;
}

bool TranCore::initializeOutputs(Id name, Status& s) {
    // Create output file if not created yet
    if (!outfile) {
        outfile = new OutputRawfile(
            name, outputDescriptors, outputSources,
            (circuit.simulatorOptions().core().rawfile==SimulatorOptions::rawfileBinary ? OutputRawfile::Flags::Binary : OutputRawfile::Flags::None) |
                OutputRawfile::Flags::Padded);
        outfile->setTitle(circuit.title());
        outfile->setPlotname("Transient Analysis");
    }
    outfile->prologue();

    return true;
}

bool TranCore::finalizeOutputs(Status& s) {
    outfile->epilogue();
    delete outfile;
    outfile = nullptr;
    
    // Write DC solution to repository if analysis is OK
    if (finished && params.store.length()>0) {
        auto sol = circuit.newStoredSolution("dc", params.store);
        sol->setNames(circuit);
        sol->values() = solution.vector();
    }

    return true;
}

bool TranCore::deleteOutputs(Id name, Status& s) {
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

bool TranCore::evalAndLoadWrapper(EvalSetup& evalSetup, LoadSetup& loadSetup) {
    clearError();
    if (!circuit.evalAndLoad(&evalSetup, &loadSetup, nullptr)) {
        // Load error
        setError(TranError::EvalAndLoad);
        if (circuit.simulatorOptions().core().tran_debug>1) {
            Simulator::dbg() << "  Evaluation error.\n";
        }
        return false;
    }
    
    // Handle abort right now, finish and stop are handled outside NR loop
    if (evalSetup.requests.abort) {
        if (circuit.simulatorOptions().core().tran_debug>1) {
            Simulator::dbg() << "  Abort requested during evaluation.\n";
        }
        return false;
    }

    // By default ignore Finish and Stop. Will be handled in timestep loop. 
    // Outside timestep loop they are ignored. 

    return true;
}

Id TranCore::methodAM = Id::createStatic("am");
Id TranCore::methodBDF = Id::createStatic("bdf");
Id TranCore::methodGear = Id::createStatic("gear");
Id TranCore::methodEuler = Id::createStatic("euler");
Id TranCore::methodTrapezoidal = Id::createStatic("trap");
Id TranCore::methodBDF2 = Id::createStatic("bdf2");
Id TranCore::methodGear2 = Id::createStatic("gear2");

CoreCoroutine TranCore::coroutine(bool continuePrevious) {
    clearError();
    auto& options = circuit.simulatorOptions().core();
    auto& internals = circuit.simulatorInternals(); 
    auto debug = options.tran_debug;

    auto n = circuit.unknownCount();

    Int maxOrder;
    Int order;

    finished = false;

    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);

    // Check parameters
    if (params.step<=0) {
        setError(TranError::Tstep);
        co_yield CoreState::Aborted;
    }

    if (params.stop<=0) {
        setError(TranError::Tstop);
        co_yield CoreState::Aborted;
    }

    if (params.start>=params.stop) {
        setError(TranError::Tstart);
        co_yield CoreState::Aborted;
    }

    // Set up integration method
    if (options.tran_method == methodAM) {
        maxOrder = options.tran_maxord;
        integCoeffs.setMethod(IntegratorCoeffs::Method::AdamsMoulton, maxOrder, options.tran_xmu);
    } else if (options.tran_method == methodBDF || options.tran_method == methodGear) {
        maxOrder = options.tran_maxord;
        integCoeffs.setMethod(IntegratorCoeffs::Method::BDF, maxOrder, options.tran_xmu);
    } else if (options.tran_method == methodEuler) {
        maxOrder = 1;
        integCoeffs.setMethod(IntegratorCoeffs::Method::AdamsMoulton, maxOrder, options.tran_xmu);
    } else if (options.tran_method == methodTrapezoidal) {
        maxOrder = 2;
        integCoeffs.setMethod(IntegratorCoeffs::Method::AdamsMoulton, maxOrder, options.tran_xmu);
    } else if (options.tran_method == methodBDF2 || options.tran_method == methodGear2) {
        maxOrder = 2;
        integCoeffs.setMethod(IntegratorCoeffs::Method::BDF, maxOrder, options.tran_xmu);
    } else {
        setError(TranError::Method);
        errorId = options.tran_method;
        co_yield CoreState::Aborted;
    }

    if (progressReporter) {
        progressReporter->setValueFormat(ProgressReporter::ValueFormat::Scientific, 6);
        progressReporter->setValueDecoration("", "s");    
    }
    initProgress(params.stop, 0);
    
    // Set up predictor
    predictorCoeffs.setMethod(IntegratorCoeffs::Method::PolynomialExtrapolation, maxOrder, options.tran_xmu);

    // Make space in history, need two extra slots for NR algorithm when computing future solution
    // Predictor works with solutions
    solution.upsize(predictorCoeffs.pastStatesNeeded()+2, n+1);
    // Integrator works with states
    states.upsize(integCoeffs.pastStatesNeeded()+2, circuit.statesCount());
    // Filtered past solutions (for use with trap ringing filter)
    // Need max(3, past points required by predictor) + last accepted solution
    if (options.tran_trapltefilter) {
        auto filteredSolutionLength = std::max(3, predictorCoeffs.pastStatesNeeded()) + 1;
        filteredSolution.upsize(filteredSolutionLength, n+1);
    }
    
    // Make space in predicted solution and scaled LTE
    predictedSolution.resize(n+1);
    scaledLte.resize(n+1);
    pastTimesteps.upsize(n+1);
    
    // breakPoints.at(1) <= tk < tSolve <= breakPoints.at(0)
    //   Entry 1 ... break points before and including last accepted solution
    //   Entry 0 ... next break point reported at last computed solution (tSolve)
    breakPoints.upsize(3);
    
    // Compute initial solution at t=0
    // Regardless of uic, op analysis is used with ic forcing
    // Solve OP analysis, in continuation mode old solution is used for starting NR iteration
    // via old solution when topologies match or nodesets (temporary forces) when topologies 
    // do not match. 
    tk = 0.0; 
    nrSolver.evalSetup().time = tk; 
    nPoints = 0; 
    if (params.icmode==icmodeOp) {
        if (debug>0) {
            Simulator::dbg() << "Solving initial conditions with OP analysis.\n";
        }
        // Use the OP analysis logic to set up force slots 0 (continuation nodesets)
        // and 1 (user-specified nodesets). Slot 2 contains permanent forces (i.e. 
        // initial condition). 
        // If ic parameter is given, activate slot 2. 
        // Forces in slot 2 were set during last call to rebuild(). 
        opCore_.solver().enableForces(2, true);
        
        // Run op analysis
        if (!opCore_.run(continuePrevious)) {
            setError(TranError::OperatingPointError);
            co_yield CoreState::Aborted;
        }
    } else if (params.icmode==icmodeUic) {
        if (debug>0) {
            Simulator::dbg() << "Setting initial conditions.\n";
        }
        // Prepare RHS with values based on ic parameter
        auto strictforce = circuit.simulatorOptions().core().strictforce; 
        // First, build the RHS values using force mechanism
        bool uicMode = true;
        if (!uicForces.set(circuit, preprocessedIc, uicMode, strictforce)) {
            // Abort on error if strictforce is set
            if (strictforce) {
                setError(TranError::UicForces);
                co_yield CoreState::Aborted;
            }
        }
        // Copy values to RHS
        solution.vector() = uicForces.unknownValue();
    } else {
        setError(TranError::IcMode);
        co_yield CoreState::Aborted;
    }
    // opCore_.dump(Simulator::dbg()); Simulator::dbg() << "\n";
    nPoints++; 

    // Write results at t=0, but only if tstart=0
    if (params.start<=0) {
        outfile->addPoint();
    }

    setProgress(0, 0);

    // Initialize reactive residual state of instances
    // Compute next breakpoint
    // Compute maximal frequency
    EvalSetup esInit = {
        .solution = &solution, 
        .states = &states, 

        .staticAnalysis = false, 
        .dcAnalysis = false, 
        .tranAnalysis = true, 
        .nodesetEnabled = false, 
        .icEnabled = true, 

        .evaluateReactiveResidual = true, 
        .evaluateLinearizedReactiveRhsResidual = true, 

        .storeReactiveState = true, 
        
        .computeBoundStep = true, 
        .computeNextBreakpoint = true, 
        .computeMaxFreq = true, 
    };
    LoadSetup lsInit = {
        .states = &states, 
    };

    // States will be loaded into the future slot
    // Make a copy of states from OP analysis 
    // in the future slot (-1)
    states.futureVector() = states.vector();
    // Evaluate and load residual and residual derivative in future slot
    // allowBypass is false by default so no bypass takes place
    if (!evalAndLoadWrapper(esInit, lsInit)) {
        // Error was already set in the wrapper
        co_yield CoreState::Aborted;
    }

    // Check for Stop/Finish at t=0 
    // In OP IC mode check opCore_ and esInit
    // In UIC IC mode check only esInit
    bool stopFlag = false;
    bool finishFlag = false;
    if (params.icmode==icmodeOp) {
        finishFlag |= opCore_.solver().evalSetup().requests.finish;
        stopFlag |= opCore_.solver().evalSetup().requests.stop;
    }
    finishFlag |= esInit.requests.finish;
    stopFlag |= esInit.requests.stop;
    if (finishFlag) {
        co_yield CoreState::Finished;
    } else if (stopFlag) {
        co_yield CoreState::Stopped;
    }

    // Clear Converged, Bypassed, and HasDeviceHistory flags
    // to make sure the device history will be initilialized 
    // after the first timepoint is computed. 
    circuit.applyInstanceFlags(
        Instance::Flags::HasDeviceHistory |
        Instance::Flags::Converged |
        Instance::Flags::Bypassed, 
        Instance::NoFlags
    );
    
    // Rotate states -1 and 0 so that future (-1) becomes current (0)
    states.rotate();

    // Set up break points
    // We just accepted solution at tsolve=0. 
    // Therefore t=0 is the first past breakpoint. 
    breakPoints.add(0.0);

    // Retrieve next breakpoint, take into account stop time and start time
    // Ignore breakpoints <=0
    double nextBreakPoint = esInit.nextBreakPoint;
    updateBreakPoint(nextBreakPoint, params.stop, 0.0);
    updateBreakPoint(nextBreakPoint, params.start, 0.0);
    // Due to stop time we definitely have a finite break point after 0.0
    breakPoints.add(nextBreakPoint);

    // Need to store last accepted point's boundStep value
    // so that we can apply it when timepoint is rejected
    acceptedBoundStep = esInit.boundStep;

    // Initial timestep
    auto h0 = params.step;
    // Limited by maxstep
    if (params.maxstep>0) {
        h0 = std::min(h0, params.maxstep);
    }
    // Limit by maxFreq (tran_ffmax*period/2)
    // Maybe get rid of this
    if (options.tran_ffmax>0 && esInit.maxFreq>0) { 
        h0 = std::min(h0, options.tran_ffmax/(2*esInit.maxFreq));
    }

    // Need to store last accepted point's hmax value
    // so that we can apply it when timepoint is rejected
    acceptedHmax = h0;
    
    // Limit by distance between last and next break point times tran_fbr
    auto breakDelta = breakPoints.at(0) - breakPoints.at(1);
    if (breakDelta>0) {
        h0 = std::min(h0, options.tran_fbr*breakDelta);
    }
    // Limit by initial boundStep
    if (esInit.boundStep>0) {
        h0 = std::min(h0, esInit.boundStep);
    }

    // Scale by tran_fs
    h0 = options.tran_fs*h0;

    // Set timestep
    auto hk = h0;

    // Compute tSolve as tk + hk, except when tk+hk=breakpoint. 
    // At breakpoints tSolve=breakpoint to avoid rounding errors. 
       
    // Set next timepoint
    auto tSolve = tk+hk;
    
    // Set current algorithm order
    order = 1;
    
    // Predictor and corrector expect t_k (in our case IC) in slot 1
    // Currently IC is in slot 0
    // Advance solution and states by 1 slot so that current (0) becomes t_k slot (1)
    solution.advance(1);
    states.advance(1); 

    // Slots summary: 
    // Slots 1, 2, 3, ... are used for past timepoints
    // Slot 0 is used for the previous solution in the NR loop
    // Slot -1 is used for the next solution in the NR loop
    // Currently IC solution and state are in slot 0
    
    // Timestep loop
    // Number of past points since last discontinuity
    size_t pointsSinceLastDiscontinuity;
    if (params.icmode==icmodeOp) {
        // Initially we have one almost valid past point
        // Operating point incorrectly computes 
        // - the current through a capacitor when it is driven by a voltage source
        // - the voltage across an inductor when it is driven by a current source 
        pointsSinceLastDiscontinuity = 0;
    } else {
        // In UIC mode the first point at t=0 is set from initial conditions 
        // specified by the user and may not be consistent. 
        pointsSinceLastDiscontinuity = 0;
    }

    // Number of consecutive points computed with trapezoidal integration
    size_t trapHistory = 0;
    
    // Initialize maximal past solution and residual contribution
    if (params.icmode==icmodeOp) {
        nrSolver.initializeMaxima(opCore_.solver());
    } else {
        nrSolver.resetMaxima();
    }
    
    while (true) {
        // NR will be applied at tSolve
        nrSolver.evalSetup().time = tSolve;

        if (debug>1) {
            ss.str(""); ss << tSolve;
            Simulator::dbg() << "  Solving at t="+ss.str();
            ss.str(""); ss << hk;
            Simulator::dbg() << " with hk="+ss.str()+"\n";
        }

        // Simulator::dbg() << "tran loop start at=" << solution.position() << "/" << solution.size() 
        //     << " state at=" << states.position() << "/" << states.size() 
        //     << " current data ptr=" << size_t(solution.data()) << "\n";

        // Compute and scale coefficients
        // Coefficients need to be recomputed when 
        // - timestep history changes
        // - current timestep changes
        // Basically this means we have to recompute them at each timestep
        predictorCoeffs.setOrder(order);
        integCoeffs.setOrder(order);
        bool havePredictor = pointsSinceLastDiscontinuity>=predictorCoeffs.minimalPredictorHistory(); // ???
        if (havePredictor && !(predictorCoeffs.compute(pastTimesteps, hk) && predictorCoeffs.scalePredictor(hk))) {
            setError(TranError::Predictor);
            co_yield CoreState::Aborted;
        }
        // Integrator coeffs must be scaled
        if (!(integCoeffs.compute(pastTimesteps, hk) && integCoeffs.scaleDifferentiator(hk))) {
            setError(TranError::Corrector);
            co_yield CoreState::Aborted;
        }

        // Compute predictor, use it to start NR solver
        // Put predicted value in slot 0
        if (!havePredictor) {
            // Not enough points to compute predictor
            // Copy solution at t_k
            predictedSolution = solution.vector(1);
        } else {
            // Prepare history for predictor, first past state is at offset 1 (we are using slots 1, 2, ...)
            if (
                options.tran_trapltefilter && 
                integCoeffs.method()==IntegratorCoeffs::Method::AdamsMoulton && 
                integCoeffs.order()==2
            ) {
                // Use filtered history if we are using trapezoidal algorithm
                predictorCoeffs.preparePredictorHistory(filteredSolution, 1);
            } else {
                predictorCoeffs.preparePredictorHistory(solution, 1);
            }
            // Predict solution
            predictorCoeffs.predict(predictedSolution);
            // for(int i=1; i<=n; i++) {
            //     auto manual = solution.at(1)[i]+(solution.at(1)[i]-solution.at(2)[i])/pastTimesteps.at(0)*hk;
            //     std::cout << i << " " << predictedSolution[i] << " manual " << manual << "\n";
            // }
        }
        
        // Turn off predictor for debugging purposes
        // zero(predictedSolution);
        
        // Initial NR solution and states at t_k go to slot 0
        if (options.tran_predictor) {
            // Use predicted solution
            solution.vector() = predictedSolution;
        } else {
            // Use previous solution
            solution.vector() = solution.vector(1);
        }
        states.vector() = states.vector(1);

        // Prepare differentiator history, first past state is at offset 1 (we are using slots 1, 2, ...)
        integCoeffs.prepareDifferentiatorHistory(states, 1);

        // Solve
        auto solutionOk = nrSolver.run(true);
        // Simulator::out() << "  Solver iterations: " << nrSolver.iterations() << "\n";
        if (!solutionOk) {
            // Solver failed. Did we have an Abort request? 
            if (nrSolver.checkFlags(TranNRSolver::Flags::Abort)) {
                co_yield CoreState::Aborted;
            }
            // No Abort, solver failed
        }
        
        if (solutionOk && options.tran_trapltefilter) {
            // Solver success, apply trap ringing filter
            // Do this only if method we were using was trapezoidal
            if (integCoeffs.method()==IntegratorCoeffs::Method::AdamsMoulton && integCoeffs.order()==2) {
                // Need at least 3 past points
                if (trapHistory>=3) {
                    filteredSolution.vector()[0] = 0;
                    for(UnknownIndex i=1; i<=n; i++) {
                        // Latest point
                        auto x0 = solution.vector()[i];
                        // Past points
                        auto x1 = solution.vector(1)[i];
                        auto x2 = solution.vector(2)[i];
                        auto x3 = solution.vector(3)[i];
                        // Past timesteps
                        auto hk1 = pastTimesteps.at(0);
                        auto hk2 = pastTimesteps.at(1);
                        // Slope of envelope defined by x1 and x3
                        auto k13 = (x1-x3)/(hk1+hk2);
                        // Slope of envelope defined by x2 and x0
                        auto k02 = (x0-x2)/(hk1+hk);
                        // Find crossing, decide if we do correction
                        auto kdelta = k13-k02;
                        bool correct = false;
                        if (kdelta==0) {
                            correct = true;
                        } else {
                            auto hcross = (x2+k02*hk1-x1)/kdelta;
                            // Crossing after x0 or before x3
                            if (hcross>hk || hcross<-(hk1+hk2)) {
                                correct = true;
                            }
                        }
                        // Always correct
                        correct = true;
                        if (correct) {
                            // Do correction
                            auto deltaEnvelopeAt2 = x2 - (x3+k13*hk2);
                            // Corrected x0 with half of envelope width
                            filteredSolution.vector()[i] = x0 - deltaEnvelopeAt2/2;
                        } else {
                            // No correction
                            filteredSolution.vector()[i] = x0;
                        }
                    }
                } else {
                    // Not enough points yet, copy point
                    filteredSolution.vector() = solution.vector();    
                }
                // Increase timepoint counter
                if (trapHistory<10) {
                    trapHistory++;
                }
            } else {
                // Not trapezoidal algorithm, copy point
                filteredSolution.vector() = solution.vector();
                // reset timepoint counter
                trapHistory = 0;
            }
        }

        // Handle Finish and Stop
        if (nrSolver.checkFlags(TranNRSolver::Flags::Finish)) {
            if (debug>0) {
                Simulator::dbg() << "Finish requested during transient analysis.\n";
            }
            // Mark analysis as finished
            finished = true;
            break;
        }
        if (nrSolver.checkFlags(TranNRSolver::Flags::Stop)) {
            if (debug>0) {
                Simulator::dbg() << "Stop requested during transient analysis.\n";
            }
            // Analysis is not finished
            break;
        }

        // Next break point assuming tSolve will be accepted,
        nextBreakPoint = nrSolver.evalSetup().nextBreakPoint;
        auto boundStep = nrSolver.evalSetup().boundStep;
        auto discontinuity = nrSolver.evalSetup().discontinuity;
        // Update breakpoint with stop and start time, ignore breakpoints <= tSolve
        // tsolve=tk+hk>0 because tk>=0 and hk>0. 
        updateBreakPoint(nextBreakPoint, params.stop, tSolve);
        updateBreakPoint(nextBreakPoint, params.start, tSolve);
        
        // Maximal timestep 
        double hmax;
        hmax = params.stop - params.start;
        // Limit to tran_rmax*step if tran_rmax>=1
        if (options.tran_rmax>=1) {
            hmax = std::min(hmax, options.tran_rmax*params.step);
        } 
        // Limit to (stop-start)/tran_minpts if tran_minpts>0
        if (options.tran_minpts>0) {
            // distance to stop
            hmax = std::min(hmax, (params.stop - params.start)/options.tran_minpts); // - tSolve;
        }
        // Limit by maxstep if given
        if (params.maxstep>0) {
            hmax = std::min(hmax, params.maxstep);
        }
        // Limit by maxFreq (tran_ffmax*period/2)
        // Maybe get rid of this
        // Use esInit's computed maxFreq value
        if (options.tran_ffmax>0 && esInit.maxFreq>0) { 
            hmax = std::min(hmax, options.tran_ffmax/(2*esInit.maxFreq));
        }

        bool accept;
        double hkNew;
        Int newOrder = order;
        if (!solutionOk) {
            // Solver failed
            accept = false;
            // New timestep will be shorter
            hkNew = hk*options.tran_ft;
            // Reduce integration order to 1
            newOrder = 1;
            if (debug>0) {
                if (debug>0) {
                    Status tmps;
                    auto nr = UnknownNameResolver(circuit);
                    nrSolver.formatError(tmps, &nr);
                    Simulator::out() << tmps.message() << "\n";
                }
                Simulator::dbg() << "  Solver failed, will reject point.\n";
            }
        } else if (!havePredictor) {
            // Cannot use LTE timestep control, not even with order=1, 
            // because there are not enough points in the history for predictor. 
            // Accept timestep, use the same timestep again. 
            accept = true;
            hkNew = hk; 
            if (debug>1) {
                Simulator::dbg() << "  Cannot estimate LTE. Will accept the point and keep the step unchanged.\n";
            }
        } else {
            // Solution OK, solver iterations count small enough, assume the timestep will be accepted
            // We have a predicted value because havePredicotr==false was handled in previous else if. 
            accept = true;
            // Compute (corrector - predictor) and scale to obtain LTE
            // factor is always<1 because errorCoeff() of polynomial predictor is negative while 
            // the errorCoeff() of integrator is positive
            double factor = integCoeffs.errorCoeff() / (integCoeffs.errorCoeff()-predictorCoeffs.errorCoeff());
            // Assume 2*hk as new timestep
            hkNew = 2*hk;

            // Spectre reference value for reltol
            // (from Designer's guide to Spice and Spectre)
            // global      .. max abs over time and nodes 
            // local       .. max abs over time, individual for each node 
            // pointlocal  .. abs, individual for each node 
            //
            // pointglobal .. max abs over nodes, individual for each time (an extra not availbale in Spectre)
            //
            // relref       solution delta      lte criterion   residual criterion      errpreset (Spectre)
            //              (2.9)               (4.57)          (2.10)
            // allglobal    global              global          global                  liberal
            // sigglobal    global              global          local                   moderate
            // alllocal     local               local           local                   conservative
            // pointlocal   pointlocal          pointlocal      pointlocal

            // Compute global reference values across past
            // Compute maximum across all unknowns at this point 
            double pointMax = 0;
            for(decltype(n) i=1; i<=n; i++) {
                double c = std::fabs(solution.vector()[i]);
                if (c>pointMax) {
                    pointMax = c;
                }
            }

            // Go through all unknowns, except for the ground
            bool haveRatio = false;
            double maxRatio = 0.0;
            for(decltype(n) i=1; i<=n; i++) {
                // Representative node, associated potential nature index
                auto rn = circuit.reprNode(i);
                bool isPotential = ((rn->flags() & Node::Flags::PotentialNode) == Node::Flags::PotentialNode); 
                size_t ndx = isPotential ? 0 : 1;

                double lte;
                if (
                    options.tran_trapltefilter && 
                    integCoeffs.method()==IntegratorCoeffs::Method::AdamsMoulton && 
                    integCoeffs.order()==2
                ) {
                    // Use filtered history if the algorithm we are using is trapezoidal
                    lte = factor*(filteredSolution.vector()[i] - predictedSolution[i]);
                } else {
                    lte = factor*(solution.vector()[i] - predictedSolution[i]);
                }
                
                // Compute tolerance
                // Scale tolerance by lteratio>1
                // Greater lteratio results in greater LTE tolerance and greater timestep

                double tol;
                if (options.relreflte == SimulatorOptions::relrefGlobal) {
                    // Maximum over all unknowns, maximum over past timepoints
                    tol = circuit.solutionTolerance(rn, nrSolver.globalMaxSolution()[ndx]);
                } else if (options.relreflte == SimulatorOptions::relrefPointGlobal) {
                    // Maximum over all unknowns, for each timepoint
                    tol = circuit.solutionTolerance(rn, nrSolver.pointMaxSolution()[ndx]);
                } else if (options.relreflte == SimulatorOptions::relrefLocal) {
                    // For each unknown, maximum over past timepoints
                    tol = circuit.solutionTolerance(rn, nrSolver.historicMaxSolution()[i]);
                } else if (options.relreflte == SimulatorOptions::relrefPointLocal) {
                    // For each unknown, for each timepoint
                    tol = circuit.solutionTolerance(rn, solution.vector()[i]);
                } else if (options.relreflte == SimulatorOptions::relrefRelref) {
                    if (options.relref == SimulatorOptions::relrefAlllocal) {
                        // For each unknown, maximum over past timepoints
                        tol = circuit.solutionTolerance(rn, nrSolver.historicMaxSolution()[i]);
                    } else if (options.relref == SimulatorOptions::relrefSigglobal) {
                        // Maximum over all unknowns, maximum over past timepoints
                        tol = circuit.solutionTolerance(rn, nrSolver.globalMaxSolution()[ndx]);
                    } else if (options.relref == SimulatorOptions::relrefAllglobal) {
                        // Maximum over all unknowns, maximum over past timepoints
                        tol = circuit.solutionTolerance(rn, nrSolver.globalMaxSolution()[ndx]);
                    } else {
                        setError(TranError::BadLteReference);
                        co_yield CoreState::Aborted;
                    }
                } else {
                    setError(TranError::BadLteReference);
                    co_yield CoreState::Aborted;
                }

                // std::cout << i << " predicted=" << predictedSolution[i] 
                //     << " obtained=" << solution.vector()[i]
                //     << " delta=" << (solution.vector()[i] - predictedSolution[i])
                //     << " lte=" << lte << " tol=" << tol << "\n";// ratio>1 means we need to decrease timestep
                
                double ratio;
                if (options.tran_spicelte) {
                    // SPICE forgets to divide tol by (order+1)!
                    ratio = std::abs(lte)/(tol*options.tran_lteratio*IntegratorCoeffs::ffactorial(order+1));
                } else {
                    ratio = std::abs(lte)/(tol*options.tran_lteratio);
                }
                // Looking for largest ratio (worst LTE)
                if (ratio>maxRatio) {
                    maxRatio = ratio;
                }
            }
            // Update timestep only if maxRatio>0 (0 means no LTE, so step can become infinite)
            // LTE grows with (order+1)-th power of the timestep
            if (maxRatio>0) {
                double hkFactor;
                if (options.tran_spicelte) {
                    hkFactor = std::pow(maxRatio, -1.0/order);
                } else {
                    hkFactor = std::pow(maxRatio, -1.0/(order+1));
                }
                // auto hkFactor = std::pow(maxRatio, -1.0/(order));
                hkNew = std::min(hkNew, hk * hkFactor);
            }
            if (debug>1) {
                ss.str(""); ss << maxRatio;
                Simulator::dbg() << "  Maximal LTE/tol="+ss.str();
                ss.str(""); ss << hkNew;
                Simulator::dbg() << " suggests dt="+ss.str()+".\n";
            }
            
            // hkNew is now the smallest of these two
            // - 2*hk
            // - maximal timestep from t_k to t_{k+1} that keeps LTE small enough
            auto hkRatio = hk/hkNew;
            if (options.tran_redofactor>0 && hkRatio > options.tran_redofactor) {
                // The timestep we used for reaching t_{k+1} is greater than tran_redofactor*hkNew. 
                // LTE is too large, need to reject timepoint
                accept = false;
                if (debug>0) {
                    ss.str(""); ss << hkNew;
                    Simulator::dbg() << "  Timestep too large, hk/hknew=" << hkRatio << ". Will reject point.\n";
                }
            } else {
                // LTE at t_{k+1} is small enough, accept timestep, 
                // use hkNew as next timestep
                accept = true;
                // Increase order
                if (newOrder<maxOrder) {
                    newOrder = newOrder + 1;
                    if (debug>1) {
                        Simulator::dbg() << "  Increasing order to "+std::to_string(newOrder)+".\n";
                    }
                }
            }
        }

        // 
        // Beyond this point we are no longer allowed to change accept
        // 

        // Timestep computation and breakpoint handling
        
        // Origin from which timestep cutting due to break point will take place
        double cutOrigin;
        if (accept) {
            // Accepting the timepoint

            // Update maxima
            nrSolver.updateMaxima();

            circuit.tables().accounting().acctNew.accepted++;
            
            // Note that this is no way to compare floats, 
            // but we set tSolve to breakPoints.at(0) when we 
            // reached breakPoints.at(0) so it should be OK. 
            if (tSolve == breakPoints.at(0)) {
                // At last stored breakpoint
                // Store next break point
                breakPoints.add(nextBreakPoint);

                // We are at a discontinuity
                discontinuity = 1;
            } else {
                // Last stored break point not reached yet
                // Is next break point after tSolve
                if (nextBreakPoint>tSolve) {
                    // Replace last stored break point
                    breakPoints.at(0) = nextBreakPoint;
                }
            }
            // Update maximal timestep due to break points
            hmax = std::min(hmax, options.tran_fbr*(breakPoints.at(0)-breakPoints.at(1))); 

            // Store last accepted boundStep and hmax value
            acceptedBoundStep = boundStep;
            acceptedHmax = hmax;

            // Timestep cutting origin is tSolve
            cutOrigin = tSolve;
        } else {
            // Rejecting the timepoint
            circuit.tables().accounting().acctNew.rejected++;
            // Fallback point tk should be between breakPoints.at(1) and breakPoints.at(0)
            // If not, panic
            if (tk<breakPoints.at(1) || tk>breakPoints.at(0)) {
                DBGCHECK(true, "Internal breakpoint handling error at rejection, t="+std::to_string(tk)+".");
                setError(TranError::BreakPointPanic);
                co_yield CoreState::Aborted;
            }
            // Is next break point after tk
            if (nextBreakPoint>tk) {
                // Replace last stored break point
                breakPoints.at(0) = nextBreakPoint;
            }
            // Update maximal timestep due to break points
            hmax = std::min(acceptedHmax, options.tran_fbr*(breakPoints.at(0)-breakPoints.at(1))); 

            // Use boundStep of last accepted point
            boundStep = acceptedBoundStep;

            // Timestep cutting origin is tk
            cutOrigin = tk;
        }

        // Limit timestep due to requests from instances
        if (boundStep>0) {
            if (debug>1 && boundStep<hkNew) {
                ss.str(""); ss << boundStep;
                Simulator::dbg() << "  Instance(s) limit timestep to dt="+ss.str()+".\n";
            }
            hkNew = std::min(hkNew, boundStep);
        }

        // Limit timestep due to hmax
        if (hmax>0) {
            if (hmax<hkNew && debug>1) {
                ss.str(""); ss << hmax;
                Simulator::dbg() << "  Timestep limited by hmax to dt="+ss.str()+".\n";
            }
            hkNew = std::min(hkNew, hmax);
        }

        // Timepoint (break point) to which hkNew should be cut
        double cutPoint = breakPoints.at(0);
        if (debug>1) {
            ss.str(""); ss << cutPoint;
            Simulator::dbg() << "  Next break point is at t="+ss.str()+".\n";
        }

        // Timestep cutting (next point crosses breakpoint) 
        // and shortening (next point close to, but before breakpoint)
        auto tNext = cutOrigin + hkNew;
        double tSolveNew;
        if (tNext>cutPoint) {
            // Next point crosses the cut point, cut the next step
            // Make sure next solve point is exactly at breakpoint
            hkNew = cutPoint - cutOrigin;
            tSolveNew = cutPoint;
            if (debug>1) {
                ss.str(""); ss << hkNew;
                Simulator::dbg() << "  Timestep cut due to break point, dt="+ss.str()+".\n";
            }
        } else if (cutPoint-tNext<0.1*hkNew) {
            // Cut point is not crossed, but next point is close to the cut point
            // Shorten timestep, i.e. divide hkNew by two so that the timestep after hkNew will not be too small
            hkNew /= 2;
            tSolveNew = cutOrigin + hkNew;
            if (debug>1) {
                ss.str(""); ss << hkNew;
                Simulator::dbg() << "  Next timepoint is close to a break point. Setting dt="+ss.str()+".\n";
            }
        } else {
            // No cutting or shortening
            tSolveNew = tNext;
        }
        
        // // Force timepoints
        // if (nPoints+1<nGlobalDbg) {
        //     hkNew=globalDbg[nPoints+1][0]-tSolve;
        //     accept=true;
        // }
        // 
        // // Check point 
        // auto sol = solution.vector()[1];
        // auto refSol = globalDbg[nPoints][1];
        // if (std::abs(sol-refSol)>1e-4) {
        //     int a=1;
        // }
        
        // Accept or reject
        if (accept) {
            nPoints++; 
            
            // If we are accepting a point at a discontinuty
            if (discontinuity>=0) {
                if (debug>1) {
                    Simulator::dbg() << "  Discontinuity reached, setting order to 1.\n";
                }
                
                // Set order to 1
                newOrder = 1;

                // Reset points since last discontinuity counter
                pointsSinceLastDiscontinuity = 1;
            } else {
                // Not at a discontinuity, increase counter
                pointsSinceLastDiscontinuity++;
            }
            
            if (debug>0) {
                ss.str(""); ss << tSolve;
                Simulator::dbg() << "Point #"+std::to_string(nPoints)+" accepted at t="+ss.str() << ", dt=";
                ss.str(""); ss << hkNew;
                Simulator::dbg() << ss.str();
                Simulator::dbg() << ", order=" << newOrder << ".\n";
            }
            
            // Write results, starting at tSolve=params.start
            if (tSolve>=params.start-timeRelativeTolerance*tk) {
                outfile->addPoint();
            }

            setProgress(tSolve, tSolve);

            // Check if we are done
            if (tSolve>=params.stop) {
                // Reached stop time, mark analysis as finished
                finished = true;
                break;
            }

            // Store timestep
            pastTimesteps.add(hk);
            
            // Advance time
            tk = tSolve; 

            // Check Finish and Stop
            // Verilog-AMS LRM states that Finish and Stop should be taken into account
            // at converged iterations (we assume that this means accepted timepoints in transient analysis). 
            if (nrSolver.evalSetup().requests.finish) {
                co_yield CoreState::Finished;
            } else if (nrSolver.evalSetup().requests.stop) {
                co_yield CoreState::Stopped;
            }
            
            // Advance history so that t_{k+1} (slot 0) becomes t_k (slot 1)
            solution.advance();
            states.advance();
            if (options.tran_trapltefilter) {
                filteredSolution.advance();
            }
        } else {
            if (debug>0) {
                ss.str(""); ss << tSolve;
                Simulator::dbg() << "Point rejected at t="+ss.str()+", dt=";
                ss.str(""); ss << hkNew;
                Simulator::dbg() << ss.str();
                Simulator::dbg() << ", order=" << newOrder << ".\n";
            }
            
            // Nothing to do, t_k slot (1) remains at the same place
        }

        // Set new hk, tSolve, and order
        hk = hkNew;
        tSolve = tSolveNew;
        order = newOrder;
        
        // Check for timestep too small
        if (hk<tSolve*timeRelativeTolerance) {
            setError(TranError::TimestepTooSmall);
            co_yield CoreState::Aborted;
        }

        // If we just accepted a timepoint, we still have points to compute, and tran_contbypass=1
        // request full bypass in next NR iteration because 
        // all instance cores are already evaluated at the inputs corresponding to the last NR iteration. 
        if (accept && options.nr_contbypass) {
            internals.requestForcedBypass = true; 
        }
    } // while

    if (debug) {
        Simulator::dbg() << "Transient analysis completed.\n";
    }
    
    co_yield CoreState::Finished;
}

bool TranCore::run(bool continuePrevious) {
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

bool TranCore::formatError(Status& s) const {
    auto nr = UnknownNameResolver(circuit);
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4);
    
    // First, handle AnalysisCore errors
    if (lastError!=Error::OK) {
        AnalysisCore::formatError(s);
        return false;
    }
    
    // Then handle TranCore errors
    switch (lastTranError) {
        case TranError::EvalAndLoad:
            s.set(Status::Analysis, "Initial state evaluation failed.");
            break;
        case TranError::MatrixError:
            jacobian.formatError(s, &nr);
            break;
        case TranError::OperatingPointError:
            opCore_.formatError(s);
            break;
        case TranError::UicForces:
            uicForces.formatError(s);
            break;
        case TranError::Tstep:
            s.set(Status::Analysis, "Transient time step must be greater than zero.");
            break;
        case TranError::Tstop: 
            s.set(Status::Analysis, "Transient stop time must be greater than zero.");
            break;
        case TranError::Tstart: 
            s.set(Status::Analysis, "Transient recording start time is after stop time.");
            break;
        case TranError::Method: 
            s.set(Status::Analysis, "Unknown integration method '"+std::string(errorId)+"'.");
            break;
        case TranError::IcMode: 
            s.set(Status::Analysis, "Unknown initial conditions mode.");
            break;
        case TranError::Predictor:
            s.set(Status::Analysis, "Failed to compute predictor coefficients.");
            break;
        case TranError::Corrector: 
            s.set(Status::Analysis, "Failed to compute integrator coefficients.");
            break;
        case TranError::TimestepTooSmall: 
            s.set(Status::Analysis, "Timestep too small. Transient analysis aborted.");
            break;
        case TranError::BadLteReference: 
            s.set(Status::Analysis, "Unsupported relreflte value.");
            break;
        case TranError::BreakPointPanic: 
            s.set(Status::Analysis, "Panic in breakpoint handling.");
            break;
        default:
            return true;
    }
    return false;
}


void TranCore::dump(std::ostream& os) const {
    AnalysisCore::dump(os);
    os << "  Points accepted: "+std::to_string(nPoints)+"\n";
    os << "  Last accepted point at: "+std::to_string(tk)+"\n"; 
}

}
