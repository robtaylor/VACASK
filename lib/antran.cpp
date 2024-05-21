#include "antran.h"
#include "common.h"


namespace NAMESPACE {

Tran::Tran(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : Analysis(name, circuit, ptAnalysis), 
      opCore(*this, params.core().opParams, circuit, jac, solution, states), 
      tranCore(*this, params.core(), opCore, circuit, jac, solution, states) { 
}

Tran::~Tran() {
}

Analysis* Tran::create(PTAnalysis& ptAnalysis, Circuit& circuit, Status& s) {
    auto* an = new Tran(ptAnalysis.name(), circuit, ptAnalysis);
    return an;
}

void Tran::clearOutputDescriptors() {
    // This function is called once before analysis/sweep is started
    // and before any output descriptors are added. 
    // Disable op analysis output. 
    params.core().opParams.writeOutput = 0;
    
    opCore.clearOutputDescriptors();
    tranCore.clearOutputDescriptors();
}

bool Tran::addCommonOutputDescriptor(const OutputDescriptor& desc) {
    // Any error causes immediate exit
    bool s1 = opCore.addOutputDescriptor(desc);
    bool s2 = tranCore.addOutputDescriptor(desc);
    return s1 && s2;
}

bool Tran::addCoreOutputDescriptors(Status& s) {
    // False is returned if something goes wrong
    if (!opCore.addCoreOutputDescriptors()) {
        opCore.formatError(s);
        return false;
    }
    if (!tranCore.addCoreOutputDescriptors()) {
        tranCore.formatError(s);
        return false;
    }
    return true;
}

bool Tran::resolveOutputDescriptors(bool strict, Status& s) {
    // Any error causes immediate exit if strict is true
    // Before exit an error message is formatted and status is set
    if (!opCore.resolveOutputDescriptors(strict)) {
        if (strict) {
            opCore.formatError(s);
            return false;
        }
    }
    if (!tranCore.resolveOutputDescriptors(strict)) {
        if (strict) {
            tranCore.formatError(s);
            return false;
        }
    }
    return true;
}

bool Tran::resolveSave(const PTSave& save, bool verify, Status& s) {
    // Tran saves
    static const auto idOpDefault = Id("default");
    static const auto idOpFull = Id("full");
    static const auto idV = Id("v");
    static const auto idI = Id("i");
    static const auto idP = Id("p");

    bool st = true;
    if (save.typeName() == idOpDefault) {
        st = tranCore.addAllUnknowns(save);
    } else if (save.typeName() == idOpFull) {
        st = tranCore.addAllNodes(save);
    } else if (save.typeName() == idV) {
        st = tranCore.addNode(save);
    } else if (save.typeName() == idI) {
        st = tranCore.addFlow(save);
    } else if (save.typeName() == idP) {
        st = tranCore.addInstanceOpvar(save);
    } else {
        // Report error only if verification is required
        if (verify) {
            s.set(Status::Save, std::string("Analysis does not support save directive."));
            s.extend(save.location());
            return false;
        } else {
            // No verification required, OK
            return true;
        }
    }
    if (verify && !st) {
        // Format error
        tranCore.formatError(s);
        s.extend(save.location());
        return false;
    }

    // No error
    return true;
}

bool Tran::addDefaultOutputDescriptors() {
    // Must be invoked on all cores regardless of return value
    auto s1 = opCore.addDefaultOutputDescriptors();
    auto s2 = tranCore.addDefaultOutputDescriptors();
    return s1 && s2;
}

bool Tran::initializeOutputs(Status& s) {
    // Any error exits immediately
    if (!opCore.initializeOutputs(std::string(name_)+".op")) {
        opCore.formatError(s);
        return false;
    }
    if (!tranCore.initializeOutputs(name_)) {
        tranCore.formatError(s);
        return false;
    }
    return true;
}

bool Tran::finalizeOutputs(Status& s) {
    // Finalization has to be performed on all cores, regardless of errors
    auto ok1 = opCore.finalizeOutputs();
    auto ok2 = tranCore.finalizeOutputs();
    if (!ok1) {
        opCore.formatError(s);
    }
    if (!ok2) {
        // Error in tranCore will mask the error in op core
        tranCore.formatError(s);
    }
    return ok1 && ok2;
}

bool Tran::deleteOutputs(Status& s) {
    // Output needs to be deleted for all cores
    auto ok1 = opCore.deleteOutputs(std::string(name_)+".op");
    auto ok2 = tranCore.deleteOutputs(name_);
    if (!ok1) {
        opCore.formatError(s);
    }
    if (!ok2) {
        // Error in tranCore will mask the error in op core
        tranCore.formatError(s);
    }
    return ok1 && ok2;
}

bool Tran::rebuildCores(Status& s) {
    // Create Jacobian - it is common to both cores, so we need to rebuild it here
    if (!jac.rebuild(circuit.sparsityMap(), circuit.unknownCount())) {
        jac.formatError(s);
        return false;
    }

    // std::cout << "Sparsity pattern" << std::endl;
    // jac->dumpSparsityTables(std::cout);
    // std::cout << std::endl;
    // jac->dumpSparsity(std::cout);
    // std::cout << std::endl;

    // First rebuild the tranCore because its rebuild function stores ICs 
    // in slot 2 of opCore's nrSolver's forces. 
    if (!tranCore.rebuild(s)) {
        return false;
    }
    if (!opCore.rebuild(s)) {
        return false;
    }
    return true;
}

size_t Tran::analysisStateStorageSize() const { 
    // Only op core has storage
    return opCore.stateStorageSize();
}

void Tran::resizeAnalysisStateStorage(size_t n) { 
    // Only op core has storage
    opCore.resizeStateStorage(n);
}

bool Tran::storeState(size_t ndx) {
    // Only op core has storage
    return opCore.storeState(ndx);
}

bool Tran::restoreState(size_t ndx) {
    // Only op core has storage
    return opCore.restoreState(ndx);
}

void Tran::makeStateIncoherent(size_t ndx) {
    opCore.makeStateIncoherent(ndx);
}

std::tuple<bool, bool> Tran::preMapping(Status& s) {
    auto [ok, needsMapping] = opCore.preMapping(s);
    if (!ok) {
        return std::make_tuple(false, needsMapping);
    }
    auto [ok1, map1] = tranCore.preMapping(s);
    return std::make_tuple(ok&&ok1, needsMapping||map1);
    
}

bool Tran::populateStructures(Status& s) {
    auto ok = opCore.populateStructures(s);
    if (!ok) {
        return false;
    }
    return tranCore.populateStructures(s);
}

bool Tran::runCores(bool continuePrevious, Status& s) {
    // Tran core will run op core
    if (!tranCore.run(continuePrevious)) {
        tranCore.formatError(s);
        return false;
    }
    return true;
}

void Tran::dump(std::ostream& os) const {
    Analysis::dump(os);
    os << "Analysis type: transient analysis"<< std::endl;
    os << "OP analysis core:" << std::endl;
    opCore.dump(os);
    os << std::endl;
    os << "Transient analysis core:" << std::endl;
    tranCore.dump(os);
}

}
