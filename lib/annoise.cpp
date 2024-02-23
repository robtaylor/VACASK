#include "annoise.h"
#include "common.h"


namespace NAMESPACE {

Noise::Noise(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : Analysis(name, circuit, ptAnalysis), 
      opCore(*this, params.core().opParams, circuit, jac, solution, states), 
      noiseCore(*this, params.core(), opCore, contributionOffset, circuit, jac, solution, states, acMatrix, acSolution, results, powerGain, outputNoise) {
    jac.setResolver(&resolver); 
    acMatrix.setResolver(&resolver);
}

Noise::~Noise() {
}

Analysis* Noise::create(PTAnalysis& ptAnalysis, Circuit& circuit, Status& s) {
    auto* an = new Noise(ptAnalysis.name(), circuit, ptAnalysis);
    return an;
}

void Noise::clearOutputDescriptors() {
    // This function is called once before analysis/sweep is started
    // and before any output descriptors are added. 
    // This is the place to enable op core output if dumpop=1. 
    // This is done only once. If dumpop changes due to changed global 
    // parameters, the change is ignored. 
    params.core().opParams.writeOutput = params.core().dumpop;
    
    opCore.clearOutputDescriptors();
    noiseCore.clearOutputDescriptors();
}

bool Noise::addCommonOutputDescriptor(const OutputDescriptor& desc) {
    // Any error causes immediate exit
    return opCore.addOutputDescriptor(desc) &&
        noiseCore.addOutputDescriptor(desc);
}

bool Noise::addCoreOutputDescriptors(Status& s) {
    // Any error causes immediate exit
    return opCore.addCoreOutputDescriptors(s) &&
        noiseCore.addCoreOutputDescriptors(s);
}

bool Noise::resolveOutputDescriptors(bool strict, Status& s) {
    // Any error causes immediate exit
    return opCore.resolveOutputDescriptors(strict, s) &&
        noiseCore.resolveOutputDescriptors(strict, s);
}

bool Noise::resolveSave(const PTSave& save, bool verify, Status& s) {
    // DC TF saves
    static const auto idDefault = Id("default");
    static const auto idFull = Id("full");
    static const auto idN  = Id("n");
    static const auto idNc = Id("nc");

    // n(instance) ... total noise contributed by instance
    // n(instance, contrib) ... noise contribution of a noise source within an isntance
    // nc(instance) ... total noise contributed by instance and all partial noise source contributions
    
    // OP saves
    static const auto idOpDefault = Id("opdefault");
    static const auto idOpFull = Id("opfull");
    static const auto idV = Id("v");
    static const auto idI = Id("i");
    static const auto idP = Id("p");

    if (save.typeName() == idDefault) {
        return noiseCore.addAllNoiseContribInst(save, verify, false, s);
    } else if (save.typeName() == idFull) {
        return noiseCore.addAllNoiseContribInst(save, verify, true, s);
    } else if (save.typeName() == idN) {
        return noiseCore.addNoiseContribInst(save, verify, false, s);
    } else if (save.typeName() == idNc) {
        return noiseCore.addNoiseContribInst(save, verify, true,  s);
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

bool Noise::addDefaultOutputDescriptors(Status& s) {
    // Any error causes immediate exit
    return opCore.addDefaultOutputDescriptors(s) && 
        noiseCore.addDefaultOutputDescriptors(s);
}

bool Noise::initializeOutputs(Status& s) {
    return opCore.initializeOutputs(std::string(name_)+".op", s) && 
        noiseCore.initializeOutputs(name_, s);
}

bool Noise::finalizeOutputs(Status& s) {
    bool ok1 = opCore.finalizeOutputs(s);
    bool ok2 = noiseCore.finalizeOutputs(s);
    return ok1 && ok2;
}

bool Noise::deleteOutputs(Status& s) {
    auto ok1 = opCore.deleteOutputs(std::string(name_)+".op", s);
    auto ok2 = noiseCore.deleteOutputs(name_, s);
    return ok1 && ok2;
}

bool Noise::rebuildCores(Status& s) {
    // Create Jacobian - it is common to both cores, so we need to rebuild it here
    if (!jac.rebuild(circuit.sparsityMap(), circuit.unknownCount(), s)) {
        return false;
    }

    return opCore.rebuild(s) && 
        noiseCore.rebuild(s);
}

size_t Noise::analysisStateStorageSize() const { 
    // Only op core has storage
    return opCore.stateStorageSize();
}

void Noise::resizeAnalysisStateStorage(size_t n) { 
    // Only op core has storage
    opCore.resizeStateStorage(n);
}

bool Noise::storeState(size_t ndx) {
    // Only op core has storage
    return opCore.storeState(ndx);
}

bool Noise::restoreState(size_t ndx) {
    // Only op core has storage
    return opCore.restoreState(ndx);
}

void Noise::makeStateIncoherent(size_t ndx) {
    opCore.makeStateIncoherent(ndx);
}

std::tuple<bool, bool> Noise::preMapping(Status& s) {
    return opCore.preMapping(s);
}

bool Noise::populateStructures(Status& s) {
    return opCore.populateStructures(s);
}

bool Noise::runCores(bool continuePrevious, Status& s) {
    // Noise core will run op core
    return noiseCore.run(continuePrevious, s);
}

void Noise::dump(std::ostream& os) const {
    Analysis::dump(os);
    os << "Analysis type: noise analysis"<< std::endl;
    os << "OP analysis core:" << std::endl;
    opCore.dump(os);
    os << std::endl;
    os << "Noise analysis core:" << std::endl;
    noiseCore.dump(os);
}

}
