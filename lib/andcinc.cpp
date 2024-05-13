#include "andcinc.h"
#include "common.h"


namespace NAMESPACE {

DcIncremental::DcIncremental(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : Analysis(name, circuit, ptAnalysis), 
      opCore(*this, params.core().opParams, circuit, jac, solution, states), 
      dcIncCore(*this, params.core(), opCore, circuit, jac, incrementalSolution) {
}

DcIncremental::~DcIncremental() {
}

Analysis* DcIncremental::create(PTAnalysis& ptAnalysis, Circuit& circuit, Status& s) {
    auto* an = new DcIncremental(ptAnalysis.name(), circuit, ptAnalysis);
    return an;
}

void DcIncremental::clearOutputDescriptors() {
    // This function is called once before analysis/sweep is started
    // and before any output descriptors are added. 
    // This is the place to enable op core output if dumpop=1. 
    // This is done only once. If dumpop changes due to changed global 
    // parameters, the change is ignored. 
    params.core().opParams.writeOutput = params.core().dumpop;
    
    opCore.clearOutputDescriptors();
    dcIncCore.clearOutputDescriptors();
}

bool DcIncremental::addCommonOutputDescriptor(const OutputDescriptor& desc) {
    // Any error causes immediate exit
    return opCore.addOutputDescriptor(desc) &&
        dcIncCore.addOutputDescriptor(desc);
}

bool DcIncremental::addCoreOutputDescriptors(Status& s) {
    // Any error causes immediate exit
    return opCore.addCoreOutputDescriptors(s) &&
        dcIncCore.addCoreOutputDescriptors(s);
}

bool DcIncremental::resolveOutputDescriptors(bool strict, Status& s) {
    // Any error causes immediate exit
    return opCore.resolveOutputDescriptors(strict, s) &&
        dcIncCore.resolveOutputDescriptors(strict, s);
}

bool DcIncremental::resolveSave(const PTSave& save, bool verify, Status& s) {
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
        return dcIncCore.addAllUnknowns(save, verify, s);
    } else if (save.typeName() == idFull) {
        return dcIncCore.addAllNodes(save, verify, s);
    } else if (save.typeName() == idDv) {
        return dcIncCore.addNode(save, verify, s);
    } else if (save.typeName() == idDi) {
        return dcIncCore.addFlow(save, verify, s);
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

bool DcIncremental::addDefaultOutputDescriptors(Status& s) {
    // Any error causes immediate exit
    return opCore.addDefaultOutputDescriptors(s) && 
        dcIncCore.addDefaultOutputDescriptors(s);
}

bool DcIncremental::initializeOutputs(Status& s) {
    return opCore.initializeOutputs(std::string(name_)+".op", s) && 
        dcIncCore.initializeOutputs(name_, s);
}

bool DcIncremental::finalizeOutputs(Status& s) {
    auto ok1 = opCore.finalizeOutputs(s); 
    auto ok2 = dcIncCore.finalizeOutputs(s); 
    return ok1 && ok2;
}

bool DcIncremental::deleteOutputs(Status& s) {
    auto ok1 = opCore.deleteOutputs(std::string(name_)+".op", s); 
    auto ok2 = dcIncCore.deleteOutputs(name_, s);
    return ok1 && ok2;
}

bool DcIncremental::rebuildCores(Status& s) {
    // Create Jacobian - it is common to both cores, so we need to rebuild it here
    if (!jac.rebuild(circuit.sparsityMap(), circuit.unknownCount())) {
        jac.formatError(s);
        return false;
    }

    return opCore.rebuild(s) &&
        dcIncCore.rebuild(s);
}

size_t DcIncremental::analysisStateStorageSize() const { 
    // Only op core has storage
    return opCore.stateStorageSize();
}

void DcIncremental::resizeAnalysisStateStorage(size_t n) { 
    // Only op core has storage
    opCore.resizeStateStorage(n);
}

bool DcIncremental::storeState(size_t ndx) {
    // Only op core has storage
    return opCore.storeState(ndx);
}

bool DcIncremental::restoreState(size_t ndx) {
    // Only op core has storage
    return opCore.restoreState(ndx);
}

void DcIncremental::makeStateIncoherent(size_t ndx) {
    opCore.makeStateIncoherent(ndx);
}

std::tuple<bool, bool> DcIncremental::preMapping(Status& s) {
    return opCore.preMapping(s);
}

bool DcIncremental::populateStructures(Status& s) {
    return opCore.populateStructures(s);
}

bool DcIncremental::runCores(bool continuePrevious, Status& s) {
    // Dc incremental core will run op core
    return dcIncCore.run(continuePrevious, s);
}

void DcIncremental::dump(std::ostream& os) const {
    Analysis::dump(os);
    os << "Analysis type: DC incremental"<< std::endl;
    os << "OP analysis core:" << std::endl;
    opCore.dump(os);
    os << std::endl;
    os << "DC incremental analysis core:" << std::endl;
    dcIncCore.dump(os);
}

}
