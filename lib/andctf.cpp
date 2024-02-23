#include "andctf.h"
#include "common.h"


namespace NAMESPACE {

DcTf::DcTf(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : Analysis(name, circuit, ptAnalysis), 
      opCore(*this, params.core().opParams, circuit, jac, solution, states), 
      dcTfCore(*this, params.core(), opCore, sourceIndex, circuit, jac, incrementalSolution, sources, tf, yin, zin) {
    jac.setResolver(&resolver); 
}

DcTf::~DcTf() {
}

Analysis* DcTf::create(PTAnalysis& ptAnalysis, Circuit& circuit, Status& s) {
    auto* an = new DcTf(ptAnalysis.name(), circuit, ptAnalysis);
    return an;
}

void DcTf::clearOutputDescriptors() {
    // This function is called once before analysis/sweep is started
    // and before any output descriptors are added. 
    // This is the place to enable op core output if dumpop=1. 
    // This is done only once. If dumpop changes due to changed global 
    // parameters, the change is ignored. 
    params.core().opParams.writeOutput = params.core().dumpop;
    
    opCore.clearOutputDescriptors();
    dcTfCore.clearOutputDescriptors();
}

bool DcTf::addCommonOutputDescriptor(const OutputDescriptor& desc) {
    // Any error causes immediate exit
    return opCore.addOutputDescriptor(desc) &&
        dcTfCore.addOutputDescriptor(desc);
}

bool DcTf::addCoreOutputDescriptors(Status& s) {
    // Any error causes immediate exit
    return opCore.addCoreOutputDescriptors(s) &&
        dcTfCore.addCoreOutputDescriptors(s);
}

bool DcTf::resolveOutputDescriptors(bool strict, Status& s) {
    // Any error causes immediate exit
    return opCore.resolveOutputDescriptors(strict, s) &&
        dcTfCore.resolveOutputDescriptors(strict, s);
}

bool DcTf::resolveSave(const PTSave& save, bool verify, Status& s) {
    // DC TF saves
    static const auto idDefault = Id("default");
    static const auto idTf  = Id("tf");
    static const auto idZin = Id("zin");
    static const auto idYin = Id("yin");

    // OP saves
    static const auto idOpDefault = Id("opdefault");
    static const auto idOpFull = Id("opfull");
    static const auto idV = Id("v");
    static const auto idI = Id("i");
    static const auto idP = Id("p");

    if (save.typeName() == idDefault) {
        return dcTfCore.addAllTfZin(save, verify, sourceIndex, s);
    } else if (save.typeName() == idTf) {
        return dcTfCore.addTf(save, verify, sourceIndex, s);
    } else if (save.typeName() == idZin) {
        return dcTfCore.addZin(save, verify,sourceIndex,  s);
    } else if (save.typeName() == idYin) {
        return dcTfCore.addYin(save, verify, sourceIndex, s);
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

bool DcTf::addDefaultOutputDescriptors(Status& s) {
    // Any error causes immediate exit
    return opCore.addDefaultOutputDescriptors(s) && 
        dcTfCore.addDefaultOutputDescriptors(s);
}

bool DcTf::initializeOutputs(Status& s) {
    return opCore.initializeOutputs(std::string(name_)+".op", s) && 
        dcTfCore.initializeOutputs(name_, s);
}

bool DcTf::finalizeOutputs(Status& s) {
    auto ok1 = opCore.finalizeOutputs(s);
    auto ok2 = dcTfCore.finalizeOutputs(s);
    return ok1 && ok2;
}

bool DcTf::deleteOutputs(Status& s) {
    auto ok1 = opCore.deleteOutputs(std::string(name_)+".op", s);
    auto ok2 = dcTfCore.deleteOutputs(name_, s);
    return ok1 && ok2;
}

bool DcTf::rebuildCores(Status& s) {
    // Create Jacobian - it is common to both cores, so we need to rebuild it here
    if (!jac.rebuild(circuit.sparsityMap(), circuit.unknownCount(), s)) {
        return false;
    }

    return opCore.rebuild(s) && 
        dcTfCore.rebuild(s);
}

size_t DcTf::analysisStateStorageSize() const { 
    // Only op core has storage
    return opCore.stateStorageSize();
}

void DcTf::resizeAnalysisStateStorage(size_t n) { 
    // Only op core has storage
    opCore.resizeStateStorage(n);
}

bool DcTf::storeState(size_t ndx) {
    // Only op core has storage
    return opCore.storeState(ndx);
}

bool DcTf::restoreState(size_t ndx) {
    // Only op core has storage
    return opCore.restoreState(ndx);
}

void DcTf::makeStateIncoherent(size_t ndx) {
    opCore.makeStateIncoherent(ndx);
}

std::tuple<bool, bool> DcTf::preMapping(Status& s) {
    return opCore.preMapping(s);
}

bool DcTf::populateStructures(Status& s) {
    return opCore.populateStructures(s);
}

bool DcTf::runCores(bool continuePrevious, Status& s) {
    // DcTf core will run op core
    return dcTfCore.run(continuePrevious, s);
}

void DcTf::dump(std::ostream& os) const {
    Analysis::dump(os);
    os << "Analysis type: DC TF"<< std::endl;
    os << "OP analysis core:" << std::endl;
    opCore.dump(os);
    os << std::endl;
    os << "DC TF analysis core:" << std::endl;
    dcTfCore.dump(os);
}

}