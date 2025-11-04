#include <vector>
#include <algorithm>
#include "core.h"
#include "corehb.h"
#include "hmtpsrc.h"
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
    registerMember(write);
    registerMember(nodeset);
    registerMember(store);
    
    return 0;
}
instantiateIntrospection(HBParameters);


HBCore::HBCore(
    OutputDescriptorResolver& parentResolver, HBParameters& params, Circuit& circuit, CommonData& commons, 
    KluBlockSparseRealMatrix& jacColoc, KluBlockSparseRealMatrix& jacobian, VectorRepository<double>& solution
) : AnalysisCore(parentResolver, circuit, commons), params(params), outfile(nullptr), jacColoc(jacColoc), 
    nrSolver(circuit, commons, jacColoc, jacobian, solution, solutionFD, frequencies, timepoints, DDT, DDTcolMajor, APFT, IAPFT, nrSettings), 
    bsjac(jacobian), solution(solution), firstBuild(true), continueState(nullptr) {
};

HBCore::~HBCore() {
    delete outfile;
}

bool HBCore::addCoreOutputDescriptors() {
    clearError();
    // If output is suppressed, skip all this work
    if (!params.write || Simulator::noOutput()) {
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
    if (!params.write || Simulator::noOutput()) {
        return true;
    }
    if (savesCount==0) {
        return addAllUnknowns(PTSave("default", Id(), Id()));
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
    if (!params.write || Simulator::noOutput()) {
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

    // Write DC solution to repository if analysis is OK
    if (converged_ && params.store.length()>0) {
        auto sol = circuit.newStoredSolution("hb", params.store);
        sol->setNames(circuit);
        sol->setCxValues(solutionFD);
        sol->setAuxData(frequencies);
    }
    return true;
}

bool HBCore::deleteOutputs(Id name, Status& s) {
    if (!params.write || Simulator::noOutput()) {
        return true;
    }

    // Cannot assume outfile is available
    auto fname = std::string(name)+".raw";
    if (std::filesystem::exists(fname)) {
        std::filesystem::remove(fname);
    }
    return true;
}

bool HBCore::storeState(size_t ndx, bool storeDetails) {
    auto& repo = coreStates.at(ndx);
    // Store current solution as annotated solution
    if (storeDetails) {
        repo.solution.setNames(circuit);
    } else {
        repo.solution.clearNames();
    }
    
    // Store solution in frequency domain (complex spectrum)
    repo.solution.setCxValues(solutionFD);
    
    // Store frequencies
    repo.solution.setAuxData(frequencies);
    
    // Stored state is coherent and valid
    repo.coherent = true;
    repo.valid = true;
    return true;
}

bool HBCore::restoreState(size_t ndx) {
    auto& state = coreStates.at(ndx);
    if (state.valid) {
        // State is valid
        continueState = &state;
        return true;
    } else { 
        // Nothing to restore, do not use continuation mode
        return false;
    }
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

    // Prepare nodesets
    auto strictforce = circuit.simulatorOptions().core().strictforce; 
    String& solutionName = params.nodeset;
    if (solutionName.length()>0) {
        // Get solution from repository
        auto solPtr = circuit.storedSolution("hb", solutionName);
        if (!solPtr) {
            // No nodesets
            nrSolver.forces(1).clear();
            Simulator::wrn() << "Warning, solution '"+solutionName+"' not found. No user nodesets applied.\n";
        } else {
            // Nodesets from solution repository
            if (!nrSolver.setForces(1, *solPtr, strictforce)) {
                // Abort if strictforce is set
                if (strictforce) {
                    nrSolver.formatError(s);
                    return false;
                }
            }
        }
    } else {
        // No nodesets, clear slot
        nrSolver.forces(1).clear();
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
    // Assume no initial state given, start with standard initial point. 
    // Coherence information is set by an.cpp and homotopy. 
    // They are responsible for detecting topology changes/rebuilds. 
    // Homotopy sweeps parameters that do not cause topology changes/rebuilds. 
    bool runInContinueMode = false;
    // Handle continuation
    if (continuePrevious) {
        // Continue mode
        if (continueState &&
            continueState->valid && continueState->coherent &&
            continueState->solution.cxValues().size()==circuit.unknownCount()*timepoints.size() && 
            continueState->solution.auxData().size()==freq.size()
        ) {
            // Continue a state
            // State is valid, coherent, and its lengths match those of the solver vectors
            // Restore current state
            solution.vector() = continueState->solution.values();
            runInContinueMode = true;
            // No forces applied
            nrSolver.enableForces(0, false);
            nrSolver.enableForces(1, false);
            if (options.hb_debug>1) {
                Simulator::dbg() << "HB using ordinary continue mode with stored analysis state.\n";
            }
            // Forced bypass is not allowed
            commons.requestForcedBypass = false;
        } else if (continueState && continueState->valid) {
            // Continue a state
            // Stored analysis state is valid, but not coherent with current circuit, 
            // its lengths may not match those of the solver vectors. 
            // Use forces to continue, but set no initial states vector nor initial solution. 
            // Ignore forces conflicts arising from stored solution. 
            // There should be no such conflicts as we are applying forces to nodes only, not node deltas. 
            nrSolver.setForces(0, continueState->solution, false);
            nrSolver.enableForces(0, true);
            // Disable user-specified forces (nodesets)
            nrSolver.enableForces(1, false);
            if (options.hb_debug>1) {
                Simulator::dbg() << "HB using forced continue mode with stored analysis state.\n";
            }
            // Forced bypass is not allowed
            commons.requestForcedBypass = false;
        } else {
            // Do not continue a state (either not provided or not valid)
            // Continue with whatever is in solution and states vector
            runInContinueMode = true;
            // No forces applied
            nrSolver.enableForces(0, false);
            nrSolver.enableForces(1, false);
            if (options.hb_debug>1) {
                Simulator::dbg() << "HB using ordinary continue mode with previous solution.\n";
            }
            // Forced bypass is not allowed
            commons.requestForcedBypass = false;
        }
        // Continue state is spent after first use
        continueState = nullptr;
    } else {
        // Continue mode not requested
        // Disable continuation forces in slot 0
        nrSolver.enableForces(0, false);

        // Apply forces specified by user in slot 1 (nodeset parameter)
        nrSolver.enableForces(1, true); 

        // Forced bypass is not allowed
        commons.requestForcedBypass = false;
        
        if (options.hb_debug>1) {
            Simulator::dbg() << "HB using standard initial solution with forced nodesets.\n";
        }
    }

    auto converged = nrSolver.run(runInContinueMode);
    auto abort = nrSolver.checkFlags(HBNRSolver::Flags::Abort);
    if (!converged_ || abort) {
        setError(HBError::SolverError);
    }

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
    converged_ = false;
    bool leave = false;
    bool tried = false;
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
    // Initial plain HB
    auto skipinitial = options.hb_skipinitial;
    if (!skipinitial) {
        tried = true;
        std::tie(converged_, leave) = runSolver(continuePrevious);
        if (!converged_) {
            setError(HBError::InitialHB);
        }
        if (debug>0) {
            if (converged_) {
                Simulator::dbg() << "HB core algorithm converged in " << std::to_string(nrSolver.iterations()) << " NR iteration(s).\n";
            } else {
                Simulator::dbg() << "HB core algorithm failed to converge in " << std::to_string(nrSolver.iterations()) << " NR iteration(s).\n";
            }
        }
    }

    // Try homotopy
    if (!converged_ && !leave && options.hb_homotopy.size()>0) {
        Homotopy* homotopy;
        for(auto it : options.hb_homotopy) {
            if (it==Homotopy::src) {
                if (debug>0) {
                    Simulator::dbg() << "Trying source stepping.\n";
                }
                homotopy = new SourceStepping(circuit, *this);
            } else {
                if (debug>0) {
                    Simulator::dbg() << "Unknown homotopy '"+std::string(it)+"'.\n";
                }
                homotopy = nullptr;
            }
            if (!homotopy) {
                continue;
            }
            // Run
            tried = true;
            std::tie(converged_, leave) = homotopy->run();
            if (debug>0) {
                if (converged_) {
                    Simulator::dbg() << "Homotopy converged in " << std::to_string(homotopy->stepCount()) << " step(s).\n";
                } else {
                    Simulator::dbg() << "Homotopy failed to converge in " << std::to_string(homotopy->stepCount()) << " step(s).\n";
                }
            }
            delete homotopy;
            if (leave || converged_) {
                break;
            }
        }
        if (!converged_) {
            setError(HBError::Homotopy);
        }
    }

    if (!leave) {
        // Did not leave early
        if (!tried) {
            // No algorithm tried
            setError(HBError::NoAlgorithm);
        } else if (converged_) {
            // Tried and converged, write results
            if (outfile && params.write) {
                // Collect results for one frequency
                outputPhasors.upsize(1, n+1);
                auto outvec = outputPhasors.data();
                for(decltype(nf) k=0; k<nf; k++) {
                    outputFreq = freq[k].f;
                    for(decltype(n) i=0; i<n; i++) {
                        outvec[i+1] = solutionFD[i*nf+k];
                    }                    
                    // Dump values at current frequency to output
                    outfile->addPoint();
                }
            }
        }
    } else {
        // Leaving early, did not converge
        // Add a status message one level higher
        converged_ = false;
    }
    
    setProgress(1);

    // std::cout << "spectrum\n";
    // int i=0;
    // for(auto it : solutionFD) {
    //     std::cout << "  " << i << " " << it << "\n";
    //     i++;
    // }

    // nrSolver.dumpSolution(std::cout, solution.data(), "  ");

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
    
    // Then handle HBCore errors
    switch (lastHbError) {
        case HBError::InitialHB:
            s.extend("Initial HB analysis failed.");
            return false;
        case HBError::Homotopy:
            s.set(Status::Analysis, "Homotopy failed.");
            return false;
        case HBError::NoAlgorithm:
            s.set(Status::Analysis, "No HB algorithm tried."); 
            return false;
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
    CommonData dummyCommons;

    HBCore hb(dummyResolver, p, dummyCircuit, dummyCommons, jacColoc, bsjac, sol);

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
