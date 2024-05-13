#include "anac.h"
#include "common.h"


namespace NAMESPACE {

Ac::Ac(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : Analysis(name, circuit, ptAnalysis), 
      opCore(*this, params.core().opParams, circuit, jac, solution, states), 
      acCore(*this, params.core(), opCore, circuit, jac, solution, states, acMatrix, acSolution) {
}

Ac::~Ac() {
}

Analysis* Ac::create(PTAnalysis& ptAnalysis, Circuit& circuit, Status& s) {
    auto* an = new Ac(ptAnalysis.name(), circuit, ptAnalysis);
    return an;
}

void Ac::clearOutputDescriptors() {
    // This function is called once before analysis/sweep is started
    // and before any output descriptors are added. 
    // This is the place to enable op core output if dumpop=1. 
    // This is done only once. If dumpop changes due to changed global 
    // parameters, the change is ignored. 
    params.core().opParams.writeOutput = params.core().dumpop;
    
    opCore.clearOutputDescriptors();
    acCore.clearOutputDescriptors();
}

bool Ac::addCommonOutputDescriptor(const OutputDescriptor& desc) {
    // Any error causes immediate exit
    return opCore.addOutputDescriptor(desc) &&
        acCore.addOutputDescriptor(desc);
}

bool Ac::addCoreOutputDescriptors(Status& s) {
    // Any error causes immediate exit
    return opCore.addCoreOutputDescriptors(s) &&
        acCore.addCoreOutputDescriptors(s);
}

bool Ac::resolveOutputDescriptors(bool strict, Status& s) {
    // Any error causes immediate exit
    return opCore.resolveOutputDescriptors(strict, s) &&
        acCore.resolveOutputDescriptors(strict, s);
}

bool Ac::resolveSave(const PTSave& save, bool verify, Status& s) {
    // DC incremental saves
    static const auto idDefault = Id("default");
    static const auto idFull = Id("full");
    static const auto idDv = Id("dv");
    static const auto idDi = Id("di");

    // OP saves
    static const auto idOpDefault = Id("opdefault");
    static const auto idOpFull = Id("opfull");
    static const auto idV = Id("v");
    static const auto idI = Id("i");
    static const auto idP = Id("p");

    if (save.typeName() == idDefault) {
        return acCore.addAllUnknowns(save, verify, s);
    } else if (save.typeName() == idFull) {
        return acCore.addAllNodes(save, verify, s);
    } else if (save.typeName() == idDv) {
        return acCore.addNode(save, verify, s);
    } else if (save.typeName() == idDi) {
        return acCore.addFlow(save, verify, s);
    } else if (save.typeName() == idOpDefault) {
        return opCore.addAllUnknowns(save, verify, s);
    } else if (save.typeName() == idOpFull) {
        return opCore.addAllNodes(save, verify, s);
    } else if (save.typeName() == idV) {
        return opCore.addNode(save, verify, s);
    } else if (save.typeName() == idI) {
        return opCore.addFlow(save, verify, s); 
    } else if (save.typeName() == idP) {
        return opCore.addInstanceOpvar(save, verify, s);
    } else {
        // Report error only if verification is required
        if (verify) {
            s.set(Status::Save, std::string("Analysis does not support save directive."));
            s.extend(save.location());
            return false;
        }
    }
    return true;
}

bool Ac::addDefaultOutputDescriptors(Status& s) {
    // Any error causes immediate exit
    return opCore.addDefaultOutputDescriptors(s) && 
        acCore.addDefaultOutputDescriptors(s);
}

bool Ac::initializeOutputs(Status& s) {
    return opCore.initializeOutputs(std::string(name_)+".op", s) && 
        acCore.initializeOutputs(name_, s);
}

bool Ac::finalizeOutputs(Status& s) {
    auto ok1 = opCore.finalizeOutputs(s);
    auto ok2 = acCore.finalizeOutputs(s);
    return ok1 && ok2;
}

bool Ac::deleteOutputs(Status& s) {
    auto ok1 = opCore.deleteOutputs(std::string(name_)+".op", s);
    auto ok2 = acCore.deleteOutputs(name_, s);
    return ok1 && ok2;
}

bool Ac::rebuildCores(Status& s) {
    // Create Jacobian - it is common to both cores, so we need to rebuild it here
    if (!jac.rebuild(circuit.sparsityMap(), circuit.unknownCount())) {
        jac.formatError(s);
        return false;
    }

    return opCore.rebuild(s) && 
        acCore.rebuild(s);
}

size_t Ac::analysisStateStorageSize() const { 
    // Only op core has storage
    return opCore.stateStorageSize();
}

void Ac::resizeAnalysisStateStorage(size_t n) { 
    // Only op core has storage
    opCore.resizeStateStorage(n);
}

bool Ac::storeState(size_t ndx) {
    // Only op core has storage
    return opCore.storeState(ndx);
}

bool Ac::restoreState(size_t ndx) {
    // Only op core has storage
    return opCore.restoreState(ndx);
}

void Ac::makeStateIncoherent(size_t ndx) {
    opCore.makeStateIncoherent(ndx);
}

std::tuple<bool, bool> Ac::preMapping(Status& s) {
    return opCore.preMapping(s);
}

bool Ac::populateStructures(Status& s) {
    return opCore.populateStructures(s);
}

bool Ac::runCores(bool continuePrevious, Status& s) {
    // Ac core will run op core
    return acCore.run(continuePrevious, s);
}

void Ac::dump(std::ostream& os) const {
    Analysis::dump(os);
    os << "Analysis type: AC small-signal"<< std::endl;
    os << "OP analysis core:" << std::endl;
    opCore.dump(os);
    os << std::endl;
    os << "AC small-signal analysis core:" << std::endl;
    acCore.dump(os);
}

}