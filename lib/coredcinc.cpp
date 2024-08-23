#include <iomanip>
#include <cmath>
#include <filesystem>
#include "coredcinc.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

// Default parameters
DcIncrParameters::DcIncrParameters() {
    opParams.writeOutput = 0; 
}

template<> int Introspection<DcIncrParameters>::setup() {
    registerMember(dumpop);
    registerNamedMember(opParams.nodeset, "nodeset");
    registerNamedMember(opParams.store, "store");

    return 0;
}
instantiateIntrospection(DcIncrParameters);


DcIncrementalCore::DcIncrementalCore(
    Analysis& analysis, DcIncrParameters& params, OperatingPointCore& opCore, Circuit& circuit, 
    KluRealMatrix& jacobian, Vector<double>& incrementalSolution
) : AnalysisCore(analysis, circuit), params(params), outfile(nullptr), opCore_(opCore), 
    jacobian(jacobian), incrementalSolution(incrementalSolution) {

    // Set analysis type for the initial operating point analysis
    auto& elsSystem = opCore_.solver().evalSetupSystem();
    elsSystem.staticAnalysis = true;
    elsSystem.dcAnalysis = false;
    elsSystem.acAnalysis = true;
}

DcIncrementalCore::~DcIncrementalCore() {
    delete outfile;
}

// Implement this in every derived class so that calls to 
// resolveOutputDescriptor() will be inlined. 
bool DcIncrementalCore::resolveOutputDescriptors(bool strict) {
    // Clear output sources
    outputSources.clear();
    // Resolve output descriptors
    bool ok = true; 
    for (auto it = outputDescriptors.cbegin(); it != outputDescriptors.cend(); ++it) {
        Node *node;
        Instance *inst;
        switch (it->type) {
        case OutdSolComponent:
            ok = addRealVarOutputSource(strict, it->id, incrementalSolution);
            break;
        default:
            // Delegate to parent
            ok = analysis.resolveOutputDescriptor(*it, outputSources, strict);
            break;
        }
        if (!ok) {
            break;
        }
    }
    return ok;
}

bool DcIncrementalCore::addDefaultOutputDescriptors() {
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    if (savesCount==0) {
        return addAllUnknowns(PTSave(Loc::bad, "default", Id(), Id()));
    }
    return true;
}

bool DcIncrementalCore::initializeOutputs(Id name) {
    // Create output file if not created yet
    if (!outfile) {
        outfile = new OutputRawfile(
            name, outputDescriptors, outputSources,
            (circuit.simulatorOptions().core().rawfile==SimulatorOptions::rawfileBinary ? OutputRawfile::Flags::Binary : OutputRawfile::Flags::None) |
                OutputRawfile::Flags::Padded);
        outfile->setTitle(circuit.title());
        outfile->setPlotname("DC Incremental Response Analysis");
    }
    outfile->prologue();

    return true;
}

bool DcIncrementalCore::finalizeOutputs() {
    outfile->epilogue();
    delete outfile;
    outfile = nullptr;
    return true;
}

bool DcIncrementalCore::deleteOutputs(Id name) {
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
    
bool DcIncrementalCore::rebuild(Status& s) {
    return true;
}

// System of equations is 
//   G(x) dx = dJ
CoreCoroutine DcIncrementalCore::coroutine(bool continuePrevious) {
    initProgress(1, 0);

    jacobian.setAccounting(circuit.tables().accounting());

    clearError();
    auto n = circuit.unknownCount();
    // Make sure structures are large enough
    incrementalSolution.resize(n+1);
    
    auto opOk = opCore_.run(continuePrevious);
    if (!opOk) {
        setError(DcIncrError::OpError);
        co_yield CoreState::Aborted;
    }

    LoadSetup lsRhs { 
        // Outputs
        .dcIncrementResidual = incrementalSolution.data()
    };

    auto& options = circuit.simulatorOptions().core();
    Int debug = options.smsig_debug;

    if (debug>0) {
        Simulator::dbg() << "Starting DC incremental analysis.\n";
    }

    // Jacobian is already factored (done by op core)
    
    // Prepare RHS (add excitations given by delta parameter)
    zero(incrementalSolution);
    auto filter = [](Device* device) { return device->checkFlags(Device::Flags::GeneratesDCIncremental); };
    if (!circuit.evalAndLoad(nullptr, &lsRhs, filter)) {
        // Load error
        setError(DcIncrError::EvalAndLoad);
        if (debug>0) {
            Simulator::dbg() << "Error in DC incremental excitation load.\n";
        }
        co_yield CoreState::Aborted;
    }

    // Change sign of residual because it is on the RHS 
    // and we need the small signal response with the correct sign
    for(decltype(n) i=0; i<=n; i++) {
        incrementalSolution[i] = -incrementalSolution[i];
    }

    if (debug>=100) {
        Simulator::dbg() << "Linear system\n";
        jacobian.dump(Simulator::dbg(), dataWithoutBucket(incrementalSolution)); 
        Simulator::dbg() << "\n";
    }

    // We don't need max residual contribution because we do not check residual

    // Solve 
    if (!jacobian.solve(dataWithoutBucket(incrementalSolution))) {
        setError(DcIncrError::MatrixError);
        co_yield CoreState::Aborted;
    }

    // Set solution bucket to 0
    incrementalSolution[0] = 0.0;

    if (options.solutioncheck && !jacobian.isFinite(dataWithoutBucket(incrementalSolution), true, true)) {
        setError(DcIncrError::SolutionError);
        if (options.smsig_debug) {
            Simulator::dbg() << "A solution entry is not finite. Solver failed.\n";
        }
        co_yield CoreState::Aborted;
    }
    
    if (debug>0) {
        Simulator::dbg() << "DC incremental analysis finished.\n";
    }

    // Dump solution
    if (params.writeOutput && outfile) {
        outfile->addPoint();
    }
    
    setProgress(1);

    co_yield CoreState::Finished;
}

bool DcIncrementalCore::run(bool continuePrevious) {
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

bool DcIncrementalCore::formatError(Status& s) const {
    auto nr = UnknownNameResolver(circuit);
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4);
    
    // First, handle AnalysisCore errors
    if (lastError!=Error::OK) {
        AnalysisCore::formatError(s);
        return false;
    }
    
    // Then handle AcCore errors
    switch (lastDcIncrError) {
        case DcIncrError::EvalAndLoad:
            s.set(Status::Analysis, "Jacobian evaluation failed.");
            break;
        case DcIncrError::MatrixError:
            jacobian.formatError(s, &nr);
            break;
        case DcIncrError::SolutionError:
            jacobian.formatError(s, &nr);
            s.extend("Solution component is not finite.");
            break;
        case DcIncrError::OpError:
            opCore_.formatError(s);
            break;
        default:
            return true;
    }
    s.extend("Leaving DC incremental analysis.");
    return false;
}

void DcIncrementalCore::dump(std::ostream& os) const {
    AnalysisCore::dump(os);
    os << "  Results" << std::endl;
    circuit.dumpSolution(os, incrementalSolution.data(), "    ");
}


}
