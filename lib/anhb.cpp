#include <cmath>
#include "anhb.h"
#include "introspection.h"
#include "devbase.h"
#include <iomanip>
#include <cmath>
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

HB::HB(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : Analysis(name, circuit, ptAnalysis), 
      jacColoc(true), // need a large bucket for loading with offsets
      jac(false), // a small bucket will do 
      core(*this, params.core(), circuit, jacColoc, jac, solution) {
};

HB::~HB() {
}

Analysis* HB::create(PTAnalysis& ptAnalysis, Circuit& circuit, Status& s) {
    auto* an = new HB(ptAnalysis.name(), circuit, ptAnalysis);
    return an;
}

void HB::clearOutputDescriptors() {
    core.clearOutputDescriptors();
}

bool HB::addCommonOutputDescriptor(const OutputDescriptor& desc) {
    return core.addOutputDescriptor(desc);
}

bool HB::addCoreOutputDescriptors(Status& s) {
    if (!core.addCoreOutputDescriptors()) {
        core.formatError(s);
        return false;
    }
    return true;
}

bool HB::resolveOutputDescriptors(bool strict, Status& s) {
    // Trigger resolving in core analyses
    return core.resolveOutputDescriptors(strict, s);
}

bool HB::resolveSave(const PTSave& save, bool verify, Status& s) {
    static const auto idDefault = Id("default");
    static const auto idFull = Id("full");
    static const auto idV = Id("v");
    static const auto idI = Id("i");
    static const auto idP = Id("p");

    bool st = true;
    // TODO: handle opvars someday
    if (save.typeName() == idDefault) {
        st = core.addAllUnknowns(save);
    } else if (save.typeName() == idFull) {
        st = core.addAllNodes(save);
    } else if (save.typeName() == idV) {
        st = core.addNode(save);
    } else if (save.typeName() == idI) {
        st = core.addFlow(save);
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
        core.formatError(s);
        s.extend(save.location());
        return false;
    }

    // No error
    return true;
}

bool HB::addDefaultOutputDescriptors() {
    return core.addDefaultOutputDescriptors();
}

bool HB::initializeOutputs(Status& s) {
    return core.initializeOutputs(name_, s);
}

bool HB::finalizeOutputs(Status& s) {
    return core.finalizeOutputs(s);
}

bool HB::deleteOutputs(Status& s) {
    return core.deleteOutputs(name_, s);
}

bool HB::rebuildCores(Status& s) {
    // Jacobian will be built by the core
    return core.rebuild(s);
}

size_t HB::analysisStateStorageSize() const { 
    return core.stateStorageSize();
}

void HB::resizeAnalysisStateStorage(size_t n) { 
    core.resizeStateStorage(n);
}

bool HB::storeState(size_t ndx) {
    return core.storeState(ndx);
}

bool HB::restoreState(size_t ndx) {
    return core.restoreState(ndx);
}

void HB::makeStateIncoherent(size_t ndx) {
    core.makeStateIncoherent(ndx);
}

std::tuple<bool, bool> HB::preMapping(Status& s) {
    return core.preMapping(s);
}

bool HB::populateStructures(Status& s) {
    return core.populateStructures(s);
}

void HB::dump(std::ostream& os) const {
    Analysis::dump(os);
    os << "Analysis type: harmonic balance"<< std::endl;
    os << "HB analysis core:" << std::endl;
    core.dump(os);
}

}
