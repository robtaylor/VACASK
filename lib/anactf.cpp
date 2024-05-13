#include "anactf.h"
#include "common.h"


namespace NAMESPACE {

AcTf::AcTf(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : Analysis(name, circuit, ptAnalysis), 
      opCore(*this, params.core().opParams, circuit, jac, solution, states), 
      acTfCore(*this, params.core(), opCore, sourceIndex, circuit, jac, solution, states, acMatrix, acSolution, sources, tf, yin, zin) {
}

AcTf::~AcTf() {
}

Analysis* AcTf::create(PTAnalysis& ptAnalysis, Circuit& circuit, Status& s) {
    auto* an = new AcTf(ptAnalysis.name(), circuit, ptAnalysis);
    return an;
}

void AcTf::clearOutputDescriptors() {
    // This function is called once before analysis/sweep is started
    // and before any output descriptors are added. 
    // This is the place to enable op core output if dumpop=1. 
    // This is done only once. If dumpop changes due to changed global 
    // parameters, the change is ignored. 
    params.core().opParams.writeOutput = params.core().dumpop;
    
    opCore.clearOutputDescriptors();
    acTfCore.clearOutputDescriptors();
}

bool AcTf::addCommonOutputDescriptor(const OutputDescriptor& desc) {
    // Any error causes immediate exit
    return opCore.addOutputDescriptor(desc) &&
        acTfCore.addOutputDescriptor(desc);
}

bool AcTf::addCoreOutputDescriptors(Status& s) {
    // Any error causes immediate exit
    return opCore.addCoreOutputDescriptors(s) &&
        acTfCore.addCoreOutputDescriptors(s);
}

bool AcTf::resolveOutputDescriptors(bool strict, Status& s) {
    // Any error causes immediate exit
    return opCore.resolveOutputDescriptors(strict, s) &&
        acTfCore.resolveOutputDescriptors(strict, s);
}

bool AcTf::resolveSave(const PTSave& save, bool verify, Status& s) {
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
        return acTfCore.addAllTfZin(save, verify, sourceIndex, s);
    } else if (save.typeName() == idTf) {
        return acTfCore.addTf(save, verify, sourceIndex, s);
    } else if (save.typeName() == idZin) {
        return acTfCore.addZin(save, verify,sourceIndex,  s);
    } else if (save.typeName() == idYin) {
        return acTfCore.addYin(save, verify, sourceIndex, s);
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

bool AcTf::addDefaultOutputDescriptors(Status& s) {
    // Any error causes immediate exit
    return opCore.addDefaultOutputDescriptors(s) && 
        acTfCore.addDefaultOutputDescriptors(s);
}

bool AcTf::initializeOutputs(Status& s) {
    return opCore.initializeOutputs(std::string(name_)+".op", s) && 
        acTfCore.initializeOutputs(name_, s);
}

bool AcTf::finalizeOutputs(Status& s) {
    bool ok1 = opCore.finalizeOutputs(s);
    bool ok2 = acTfCore.finalizeOutputs(s);
    return ok1 && ok2;
}

bool AcTf::deleteOutputs(Status& s) {
    auto ok1 = opCore.deleteOutputs(std::string(name_)+".op", s);
    auto ok2 = acTfCore.deleteOutputs(name_, s);
    return ok1 && ok2;
}

bool AcTf::rebuildCores(Status& s) {
    // Create Jacobian - it is common to both cores, so we need to rebuild it here
    if (!jac.rebuild(circuit.sparsityMap(), circuit.unknownCount())) {
        jac.formatError(s);
        return false;
    }

    return opCore.rebuild(s) && 
        acTfCore.rebuild(s);
}

size_t AcTf::analysisStateStorageSize() const { 
    // Only op core has storage
    return opCore.stateStorageSize();
}

void AcTf::resizeAnalysisStateStorage(size_t n) { 
    // Only op core has storage
    opCore.resizeStateStorage(n);
}

bool AcTf::storeState(size_t ndx) {
    // Only op core has storage
    return opCore.storeState(ndx);
}

bool AcTf::restoreState(size_t ndx) {
    // Only op core has storage
    return opCore.restoreState(ndx);
}

void AcTf::makeStateIncoherent(size_t ndx) {
    opCore.makeStateIncoherent(ndx);
}

std::tuple<bool, bool> AcTf::preMapping(Status& s) {
    return opCore.preMapping(s);
}

bool AcTf::populateStructures(Status& s) {
    return opCore.populateStructures(s);
}

bool AcTf::runCores(bool continuePrevious, Status& s) {
    // AcTf core will run op core
    return acTfCore.run(continuePrevious, s);
}

void AcTf::dump(std::ostream& os) const {
    Analysis::dump(os);
    os << "Analysis type: AC TF"<< std::endl;
    os << "OP analysis core:" << std::endl;
    opCore.dump(os);
    os << std::endl;
    os << "AC TF analysis core:" << std::endl;
    acTfCore.dump(os);
}

}
