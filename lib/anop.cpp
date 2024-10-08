#include <cmath>
#include "anop.h"
#include "introspection.h"
#include "devbase.h"
#include <iomanip>
#include <cmath>
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

OperatingPoint::OperatingPoint(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : Analysis(name, circuit, ptAnalysis), 
      core(*this, params.core(), circuit, jac, solution, states) {
};

OperatingPoint::~OperatingPoint() {
}

Analysis* OperatingPoint::create(PTAnalysis& ptAnalysis, Circuit& circuit, Status& s) {
    auto* an = new OperatingPoint(ptAnalysis.name(), circuit, ptAnalysis);
    return an;
}

void OperatingPoint::clearOutputDescriptors() {
    core.clearOutputDescriptors();
}

bool OperatingPoint::addCommonOutputDescriptor(const OutputDescriptor& desc) {
    return core.addOutputDescriptor(desc);
}

bool OperatingPoint::addCoreOutputDescriptors(Status& s) {
    if (!core.addCoreOutputDescriptors()) {
        core.formatError(s);
        return false;
    }
    return true;
}

bool OperatingPoint::resolveOutputDescriptors(bool strict, Status& s) {
    // Trigger resolving in core analyses
    return core.resolveOutputDescriptors(strict, s);
}

bool OperatingPoint::resolveSave(const PTSave& save, bool verify, Status& s) {
    static const auto idDefault = Id("default");
    static const auto idFull = Id("full");
    static const auto idV = Id("v");
    static const auto idI = Id("i");
    static const auto idP = Id("p");

    bool st = true;
    if (save.typeName() == idDefault) {
        st = core.addAllUnknowns(save);
    } else if (save.typeName() == idFull) {
        st = core.addAllNodes(save);
    } else if (save.typeName() == idV) {
        st = core.addNode(save);
    } else if (save.typeName() == idI) {
        st = core.addFlow(save);
    } else if (save.typeName() == idP) {
        st = core.addInstanceOpvar(save);
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

bool OperatingPoint::addDefaultOutputDescriptors() {
    return core.addDefaultOutputDescriptors();
}

bool OperatingPoint::initializeOutputs(Status& s) {
    return core.initializeOutputs(name_, s);
}

bool OperatingPoint::finalizeOutputs(Status& s) {
    return core.finalizeOutputs(s);
}

bool OperatingPoint::deleteOutputs(Status& s) {
    return core.deleteOutputs(name_, s);
}

bool OperatingPoint::rebuildCores(Status& s) {
    // Create Jacobian - it is common to both cores in small-signal analyses
    // Therefore we build it outside the core before the core is rebuilt. 
    if (!jac.rebuild(circuit.sparsityMap(), circuit.unknownCount())) {
        jac.formatError(s);
        return false;
    }

    return core.rebuild(s);
}

size_t OperatingPoint::analysisStateStorageSize() const { 
    return core.stateStorageSize();
}

size_t OperatingPoint::allocateAnalysisStateStorage(size_t n) { 
    return core.allocateStateStorage(n);
}

void OperatingPoint::deallocateAnalysisStateStorage(size_t n) { 
    core.deallocateStateStorage(n);
}

bool OperatingPoint::storeState(size_t ndx, bool storeDetails) {
    return core.storeState(ndx, storeDetails);
}

bool OperatingPoint::restoreState(size_t ndx) {
    return core.restoreState(ndx);
}

void OperatingPoint::makeStateIncoherent(size_t ndx) {
    core.makeStateIncoherent(ndx);
}

std::tuple<bool, bool> OperatingPoint::preMapping(Status& s) {
    return core.preMapping(s);
}

bool OperatingPoint::populateStructures(Status& s) {
    return core.populateStructures(s);
}

void OperatingPoint::dump(std::ostream& os) const {
    Analysis::dump(os);
    os << "Analysis type: operating point"<< std::endl;
    os << "OP analysis core:" << std::endl;
    core.dump(os);
}

}
