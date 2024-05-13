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
    return opCore.addOutputDescriptor(desc) &&
         tranCore.addOutputDescriptor(desc);
}

bool Tran::addCoreOutputDescriptors(Status& s) {
    // Any error causes immediate exit
    return opCore.addCoreOutputDescriptors(s) &&
         tranCore.addCoreOutputDescriptors(s);
}

bool Tran::resolveOutputDescriptors(bool strict, Status& s) {
    // Any error causes immediate exit
    return opCore.resolveOutputDescriptors(strict, s) &&
         tranCore.resolveOutputDescriptors(strict, s);
}

bool Tran::resolveSave(const PTSave& save, bool verify, Status& s) {
    // Tran saves
    static const auto idOpDefault = Id("default");
    static const auto idOpFull = Id("full");
    static const auto idV = Id("v");
    static const auto idI = Id("i");
    static const auto idP = Id("p");

    if (save.typeName() == idOpDefault) {
        return tranCore.addAllUnknowns(save, verify, s);
    } else if (save.typeName() == idOpFull) {
        return tranCore.addAllNodes(save, verify, s);
    } else if (save.typeName() == idV) {
        return tranCore.addNode(save, verify, s);
    } else if (save.typeName() == idI) {
        return tranCore.addFlow(save, verify, s); 
    } else if (save.typeName() == idP) {
        return tranCore.addInstanceOpvar(save, verify, s);
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

bool Tran::addDefaultOutputDescriptors(Status& s) {
    // Any error causes immediate exit
    return opCore.addDefaultOutputDescriptors(s) && 
         tranCore.addDefaultOutputDescriptors(s);
}

bool Tran::initializeOutputs(Status& s) {
    return opCore.initializeOutputs(std::string(name_)+".op", s) && 
         tranCore.initializeOutputs(name_, s);
}

bool Tran::finalizeOutputs(Status& s) {
    auto ok1 = opCore.finalizeOutputs(s);
    auto ok2 = tranCore.finalizeOutputs(s);
    return ok1 && ok2;
}

bool Tran::deleteOutputs(Status& s) {
    auto ok1 = opCore.deleteOutputs(name_, s);
    auto ok2 = tranCore.deleteOutputs(name_, s);
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
    return tranCore.rebuild(s) && 
             opCore.rebuild(s);
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
    return tranCore.preMapping(s);
}

bool Tran::populateStructures(Status& s) {
    return tranCore.populateStructures(s);
}

bool Tran::runCores(bool continuePrevious, Status& s) {
    // Tran core will run op core
    return tranCore.run(continuePrevious, s);
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
