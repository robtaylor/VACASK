#include "anac.h"
#include "common.h"


namespace NAMESPACE {

template<> SmallSignal<AcCore, AcData>::SmallSignal(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : Analysis(name, circuit, ptAnalysis), 
      opCore(*this, params.core().opParams, circuit, commons, jac, solution, states), 
      smsigCore(*this, params.core(), opCore, circuit, commons, jac, solution, states, acMatrix, acSolution) {
}

template<> bool SmallSignal<AcCore, AcData>::resolveSave(const PTSave& save, bool verify, Status& s) {
    // AC saves
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

template<> void SmallSignal<AcCore, AcData>::dump(std::ostream& os) const {
    Analysis::dump(os);
    os << "Analysis type: AC small-signal\n";
    os << "OP analysis core:\n";
    opCore.dump(os);
    os << "\n";
    os << "AC small-signal analysis core:\n";
    smsigCore.dump(os);
}

}