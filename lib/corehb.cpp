#include <vector>
#include <algorithm>
#include "core.h"
#include "corehb.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

Id HBCore::truncateRaw = Id::createStatic("raw");
Id HBCore::truncateBox = Id::createStatic("box");
Id HBCore::truncateDiamond = Id::createStatic("diamond");

Id HBCore::sampleUniform = Id::createStatic("uniform");
Id HBCore::sampleRandom = Id::createStatic("random");

HBParameters::HBParameters() {
    truncate = HBCore::truncateDiamond;
    sample = HBCore::sampleRandom;
}

template<> int Introspection<HBParameters>::setup() {
    registerMember(freq);
    registerMember(nharm);
    registerMember(immax);
    registerMember(truncate);
    registerMember(samplefac);
    registerMember(nper);
    registerMember(sample);
    
    return 0;
}
instantiateIntrospection(HBParameters);


HBCore::HBCore(
    OutputDescriptorResolver& parentResolver, HBParameters& params, Circuit& circuit, 
    KluBlockSparseRealMatrix& jacColoc, KluBlockSparseRealMatrix& jacobian, VectorRepository<double>& solution
) : AnalysisCore(parentResolver, circuit), params(params), outfile(nullptr), jacColoc(jacColoc), 
    nrSolver(circuit, jacColoc, jacobian, solution, solutionFD, frequencies, timepoints, DDT, DDTcolMajor, APFT, nrSettings), 
    bsjac(jacobian), solution(solution), firstBuild(true) {
};

HBCore::~HBCore() {
    delete outfile;
}

bool HBCore::addCoreOutputDescriptors() {
    clearError();
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    if (!addOutputDescriptor(OutputDescriptor(OutdFrequency, "frequency"))) {
        lastError = Error::Descriptor;
        errorId = "frequency";
        return false;
    }
    return true;
}

bool HBCore::addDefaultOutputDescriptors() {
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    if (savesCount==0) {
        return addAllUnknowns(PTSave(Loc::bad, "default", Id(), Id()));
    }
    return true;
}

bool HBCore::resolveOutputDescriptors(bool strict, Status& s) {
    // Clear output sources
    outputSources.clear();
    // Resolve output descriptors
    bool ok = true;
    for (auto it = outputDescriptors.cbegin(); it != outputDescriptors.cend(); ++it) {
        Node *node;
        Instance *inst;
        // TODO: handle opvars someday
        switch (it->type) {
        case OutdSolComponent:
            ok = addComplexVarOutputSource(strict, it->id, outputPhasors); 
            break;
        case OutdFrequency:
            outputSources.emplace_back(&outputFreq);
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

bool HBCore::initializeOutputs(Id name, Status& s) {
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    // Create output file if not created yet
    if (!outfile) {
        outfile = new OutputRawfile(
            name, outputDescriptors, outputSources,
            (circuit.simulatorOptions().core().rawfile==SimulatorOptions::rawfileBinary ? OutputRawfile::Flags::Binary : OutputRawfile::Flags::None) |
                OutputRawfile::Flags::Padded | OutputRawfile::Flags::Complex);
        outfile->setTitle(circuit.title());
        outfile->setPlotname("Harmonic Balance Analysis");
    }
    outfile->prologue();

    return true;
}

bool HBCore::finalizeOutputs(Status& s) {
    if (outfile) {
        outfile->epilogue();
        delete outfile;
        outfile = nullptr;
    }
    return true;
}

bool HBCore::deleteOutputs(Id name, Status& s) {
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

size_t HBCore::stateStorageSize() const { 
    return analysisStateRepository.size(); 
}

size_t HBCore::allocateStateStorage(size_t n) { 
    auto nOld = analysisStateRepository.size();
    analysisStateRepository.resize(nOld+n); 
    
    // Initially no state is coherent nor valid
    for(decltype(nOld) i=nOld; i<nOld+n; i++) {
        analysisStateRepository[i].coherent = false;
        analysisStateRepository[i].valid = false;
    }
    
    return nOld;
}

void HBCore::deallocateStateStorage(size_t n) {
    auto nOld = analysisStateRepository.size();
    if (n==0 || n>nOld) {
        n = nOld;
    }
    analysisStateRepository.resize(nOld-n);
}

bool HBCore::storeState(size_t ndx, bool storeDetails) {
    auto& repo = analysisStateRepository.at(ndx);
    // Store current solution as annotated solution
    if (storeDetails) {
        repo.solution.setNames(circuit);
    }
    
    // solutionFD is a complex spectrum, we need to convert it to an APFT spectrum
    auto nt = timepoints.size();
    auto& sol = repo.solution.values();
    sol.resize(nt);
    // DC
    sol[0] = solutionFD[0].real();
    // f>0
    auto nf = solutionFD.size();
    for(decltype(nf) i=1; i<nf; i++) {
        sol[1+2*(i-1)] = solutionFD[i].real();
        sol[1+2*(i-1)+1] = solutionFD[i].imag();
    }
    
    // Store frequencies
    repo.solution.auxData() = frequencies;
    
    // Stored state is coherent and valid
    repo.coherent = true;
    repo.valid = true;
    return true;
}

bool HBCore::restoreState(size_t ndx) {
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

void HBCore::makeStateIncoherent(size_t ndx) {
    analysisStateRepository.at(ndx).coherent = false;
}

// Analysis asks cores if they request a rebuild. 
// HB core replies that it does if the set of frequencies changes. 
// Along with changed set of frequencies this function recomputes
// - colocation points
// - transform matrices
std::tuple<bool, bool> HBCore::requestsRebuild(Status& s) {
    // First build, nothing to compare to
    if (firstBuild) {
        return std::make_tuple(true, true);
    }

    // Did parameters that affect the set of frequencies and the colocation timepoints change
    bool needsRebuild = 
        oldParams.freq != params.freq ||
        oldParams.nharm != params.nharm ||
        oldParams.immax != params.immax ||
        oldParams.truncate != params.truncate || 
        oldParams.samplefac != params.samplefac ||
        oldParams.nper != params.nper || 
        oldParams.sample != params.sample;
    oldParams = params;
    return std::make_tuple(true, needsRebuild);
}

bool HBCore::rebuild(Status& s) {
    clearError();

    auto& options = circuit.simulatorOptions().core();
    nrSettings = NRSettings {
        .debug = options.nr_debug, 
        .itlim = options.hb_itl, 
        .itlimCont = options.hb_itlcont, 
        .convIter = options.nr_conviter, 
        .residualCheck = bool(options.nr_residualcheck),  
        .dampingFactor = options.nr_damping, 
        .matrixCheck = bool(options.matrixcheck), 
        .rhsCheck = bool(options.rhscheck), 
        .solutionCheck = bool(options.solutioncheck), 
        .forceFactor = options.nr_force, 
    };

    // Compute set of frequencies
    if (!buildGrid(s)) {
        return false;
    }

    // Compute colocation
    if (!buildColocation(s)) {
        return false;
    }

    // Recompute transforms
    if (!buildAPFT(s)) {
        return false;
    }

    // Number of colocation points
    auto nt = timepoints.size();

    // Jacobian entries at colocation points
    jacColoc.rebuild(circuit.sparsityMap(), circuit.unknownCount(), nt, 2);

    // HB Jacobian
    if (!bsjac.rebuild(circuit.sparsityMap(), circuit.unknownCount(), nt, nt)) {
        setError(HBError::MatrixError);
        return false;
    }

    // Bind resistive residuals to 0-based subelement (0,0) 
    // Bind reactive residuals to 0-based subelement (0,1) 
    
    // Resistive Jacobian entries remain bound to OP Jacobian, 
    // reactive parts will be bound to imaginary entries of acMatrix
    if (!circuit.bind(
        &jacColoc, Component::Real, MatrixEntryPosition(0, 0), 
        &jacColoc, Component::Real, MatrixEntryPosition(0, 1), 
        s
    )) {
        return false;
    }

    // Rebuild NR solver structures
    if (!nrSolver.rebuild()) {
        s.set(Status::NonlinearSolver, "Failed to rebuild internal structures of nonlinear solver.");
        return false;
    }
    
    firstBuild = false;
    return true;
}

std::tuple<bool, bool> HBCore::runSolver(bool continuePrevious) {
    auto& options = circuit.simulatorOptions().core();
    bool hasInitialStates;
    // Handle continuation
    if (continuePrevious) {
        // Continue mode
        if (continueState &&
            continueState->valid && continueState->coherent &&
            continueState->solution.values().size()==circuit.unknownCount()*timepoints.size() && 
            continueState->solution.auxData().size()==freq.size()
        ) {
            // Continue a state
            // State is valid, coherent, and its lengths match those of the solver vectors
            // Restore current state
            solution.vector() = continueState->solution.values();
            hasInitialStates = true;
            // No nodesets applied
            nrSolver.enableForces(0, false);
            nrSolver.enableForces(1, false);
            if (options.op_debug>1) {
                Simulator::dbg() << "OP using ordinary continue mode with stored state.\n";
            }
            // Use forced bypass if allowed
            circuit.simulatorInternals().requestForcedBypass = circuit.simulatorInternals().allowContinueStateBypass;
        } else if (continueState && continueState->valid) {
            // Stored analysis state is valid, but not coherent with current circuit, 
            // its lengths may not match those of the solver vectors. 
            // Use nodesets to continue, but set no initial states vector nor initial solution. 
            // Ignore nodeset conflicts arising from stored solution. 
            // There should be no such conflicts as we are applying nodesets to nodes only, not node deltas. 
            // TODO: nrSolver.forces(0).set(circuit, continueState->solution, false);
            nrSolver.enableForces(0, true);
            // Disable user-specified nodesets
            nrSolver.enableForces(1, false);
            if (options.op_debug>1) {
                Simulator::dbg() << "OP using nodeset continue mode with stored state.\n";
            }
            // Forced bypass is not allowed
            circuit.simulatorInternals().requestForcedBypass = false;
        } else {
            // No valid state, continue with whatever is in solution and states vector
            hasInitialStates = true;
            // No nodesets applied
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

    auto converged = nrSolver.run(hasInitialStates);
    auto abort = nrSolver.checkFlags(HBNRSolver::Flags::Abort);
    return std::make_tuple(converged, abort);
}

Int HBCore::iterations() const {
    return nrSolver.iterations();
}

Int HBCore::iterationLimit(bool continuePrevious) const {
    return continuePrevious ? nrSettings.itlimCont : nrSettings.itlim;
}

CoreCoroutine HBCore::coroutine(bool continuePrevious) {
    initProgress(1, 0);

    clearError();
    
    auto& options = circuit.simulatorOptions().core();
    auto& internals = circuit.simulatorInternals();
    converged_ = false;
    auto debug = options.hb_debug;
    auto n = circuit.unknownCount(); 
    auto nb = timepoints.size();
    auto nf = freq.size();

    // Make sure structures are large enough
    solution.upsize(2, n*nb+1);
    
    if (debug>0) {
        Simulator::dbg() << "Starting HB analysis.\n";
    }

    // Run solver (time domain formulation)
    converged_ = nrSolver.run(continuePrevious);
    if (!converged_) {
        setError(HBError::SolverError);
    }

    if (converged_ && params.writeOutput) {
        // Collect results
        outputPhasors.upsize(1, n+1);
        auto outvec = outputPhasors.data();
        for(decltype(nf) k=0; k<nf; k++) {
            outputFreq = freq[k].f;
            for(decltype(n) i=0; i<n; i++) {
                outvec[i+1] = solutionFD[i*nf+k];
            }
            
            // Dump to output
            outfile->addPoint();
        }
    }
    
    setProgress(1);

    // HB analysis can only Abort or Finish
    if (converged_) {
        co_yield CoreState::Finished;
    } else {
        co_yield CoreState::Aborted;
    }
}

bool HBCore::run(bool continuePrevious) {
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


class HbUnknownNameResolver : public NameResolver {
public:
    HbUnknownNameResolver(Circuit& circuit, size_t nb) : circuit(circuit), nb(nb) {};

    virtual Id operator()(MatrixEntryIndex u) {
        return circuit.reprNode(u/nb+1)->name();
    };

private:
    Circuit& circuit;
    size_t nb;
};


bool HBCore::formatError(Status& s) const {
    auto nb = timepoints.size();
    auto nr = HbUnknownNameResolver(circuit, nb);
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4);

    // First, handle AnalysisCore errors
    if (lastError!=Error::OK) {
        AnalysisCore::formatError(s);
        return false;
    }
    
    // Then handle OperatingPointCore errors
    switch (lastHbError) {
        case HBError::MatrixError:
            bsjac.formatError(s, &nr);
            return false;
        case HBError::SolverError:
            nrSolver.formatError(s, &nr);
            return false;
    }
    return true;
}

void HBCore::dump(std::ostream& os) const {
    AnalysisCore::dump(os);
    os << "  Results\n";
    auto n = circuit.unknownCount();
    auto nf = frequencies.size();
    for(decltype(n) i=1; i<=n; i++) {
        auto rn = circuit.reprNode(i);
        for(decltype(nf) k=0; k<nf; k++) {
            auto c = solutionFD[i];
            os << "    " << rn->name() << "@" << frequencies[k] << "Hz : " << c.real();
            if (c.imag()>=0) {
                os << "+";
            }
            os << c.imag();
            os << "i\n";
        }
    }
}

bool HBCore::test() {
    Status s;
    bool ok = true;

    HBParameters p;
    p.freq = {1000, 100000};
    p.nharm = 4;
    p.truncate = "diamond";
    p.sample = "random";
    p.samplefac = 4;

    // Dummy strutures
    OutputDescriptorResolver dummyResolver;
    KluBlockSparseRealMatrix jacColoc;
    KluBlockSparseRealMatrix bsjac;
    VectorRepository<double> sol;
    ParserTables tab;
    Circuit dummyCircuit(tab);

    HBCore hb(dummyResolver, p, dummyCircuit, jacColoc, bsjac, sol);

    if (ok && !hb.buildGrid(s)) {
        ok = false;
        std::cout << "Failed to build grid: " << s.message() << "\n";
    }
    
    if (ok && !hb.buildColocation(s)) {
        ok = false;
        std::cout << "Failed to select colocation points: " << s.message() << "\n";
    } 
    
    if (!hb.buildAPFT(s)) {
        ok = false;
        std::cout << "Failed to build APFT: " << s.message() << "\n";
    }

    if (ok) {
        double delta;
        auto n = hb.timepoints.size();
        DenseMatrix<double> result(n, n);

        // APFT*IAPFT
        hb.APFT.multiply(hb.IAPFT, result);
        DenseMatrix<double> I(n, n);
        I.identity();
        result.subtract(I, result);
        delta = result.maxAbs();
        std::cout << "APFT * IAPFT - I :: delta = " << delta << "\n";
        std::cout << "\n";
        if (delta>1e-12) {
            ok = false;
            std::cout << "IAPFT inverse failed\n";
        }

        // APFT of first non zero frequency cosine
        std::vector<double> v(n, 0.0);
        std::vector<double> vres(n, 0.0);
        auto f = hb.freq[1].f;
        auto mag = 10;
        for(size_t i=0; i<n; i++) {
            auto t = hb.timepoints[i];
            v[i] = mag*std::cos(2*std::numbers::pi*f*t);
        }

        auto vv = VectorView<double>(v);
        auto vvres = VectorView<double>(vres);
        hb.APFT.multiply(vv, vvres);
        auto norm = vvres.maxAbs();
        std::cout << "APFT of cosine at f1\n";
        vvres.dump(std::cout);
        std::cout << "\n";
        for(size_t i=0; i<n; i++) {
            if (
                i==1 && std::abs(vvres[i]-mag)/norm>1e-12 ||
                i!=1 && std::abs(vvres[i])/norm>1e-12
            ) {
                ok = false;
                std::cout << "APFT failed\n";
                break;
            }
        }

        // Derivative of first nonzero frequency cosine wrt time
        hb.DDT.multiply(vv, vvres);
        norm = vvres.maxAbs();
        std::cout << "APFT of DDT of cosine at f1\n";
        auto spec = std::vector<double>(n, 0);
        auto vspec = VectorView<double>(spec);
        hb.APFT.multiply(vvres, vspec);
        vspec.dump(std::cout);
        std::cout << "\n";
        for(size_t i=0; i<n; i++) {
            auto t = hb.timepoints[i];
            auto exact = -mag*2*std::numbers::pi*f*std::sin(2*std::numbers::pi*f*t);
            if (std::abs(vvres[i] - exact)/norm>1e-12) {
                ok = false;
                std::cout << "DDT failed\n";
                break;
            }
        }

    }

    if (!ok) {
        std::cout << "HB core test failed.\n";
    }

    return ok;
}

}
