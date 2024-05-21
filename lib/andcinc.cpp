#include "andcinc.h"
#include "common.h"


namespace NAMESPACE {

template<> SmallSignal<DcIncrementalCore, DcIncrData>::SmallSignal(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : Analysis(name, circuit, ptAnalysis), 
      opCore(*this, params.core().opParams, circuit, jac, solution, states), 
      smsigCore(*this, params.core(), opCore, circuit, jac, incrementalSolution) {
}

template<> bool SmallSignal<DcIncrementalCore, DcIncrData>::resolveSave(const PTSave& save, bool verify, Status& s) {
    // DC incremental saves
    static const auto idDefault = Id("default");
    static const auto idFull = Id("full");
    static const auto idDv = Id("dv");
    static const auto idDi = Id("di");

    bool st = true;
    bool handled = true;
    if (save.typeName() == idDefault) {
        st = smsigCore.addAllUnknowns(save);
    } else if (save.typeName() == idFull) {
        st = smsigCore.addAllNodes(save);
    } else if (save.typeName() == idDv) {
        st = smsigCore.addNode(save);
    } else if (save.typeName() == idDi) {
        st = smsigCore.addFlow(save);
    } else {
        // Handle OP saves
        std::tie(st, handled) = resolveOpSave(save, verify, s); 
        // Not handled error was formatted by resolveOpSave()
        // Also all op errors were formatted
        if (verify) {
            // Verification required, return status
            return st;
        } else {
            // No verification required, OK
            return true;
        }
    }

    // Handled save via smsigCore, check error if verification required
    if (verify && !st) {
        // Format error
        smsigCore.formatError(s);
        s.extend(save.location());
        return false;
    } 
    
    // No error
    return true;
}

template<> void SmallSignal<DcIncrementalCore, DcIncrData>::dump(std::ostream& os) const {
    Analysis::dump(os);
    os << "Analysis type: DC incremental"<< std::endl;
    os << "OP analysis core:" << std::endl;
    opCore.dump(os);
    os << std::endl;
    os << "DC incremental analysis core:" << std::endl;
    smsigCore.dump(os);
}

}
