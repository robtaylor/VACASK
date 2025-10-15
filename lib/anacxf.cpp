#include "anacxf.h"
#include "common.h"


namespace NAMESPACE {

template<> SmallSignal<ACXFCore, ACXFData>::SmallSignal(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : Analysis(name, circuit, ptAnalysis), 
      opCore(*this, params.core().opParams, circuit, commons, jac, solution, states), 
      smsigCore(*this, params.core(), opCore, sourceIndex, circuit, commons, jac, solution, states, acMatrix, acSolution, sources, tf, yin, zin) {
}

template<> bool SmallSignal<ACXFCore, ACXFData>::resolveSave(const PTSave& save, bool verify, Status& s) {
    // ACXF saves
    static const auto idDefault = Id("default");
    static const auto idTf  = Id("tf");
    static const auto idZin = Id("zin");
    static const auto idYin = Id("yin");

    bool st = true;
    bool handled = true;
    if (save.typeName() == idDefault) {
        st = smsigCore.addAllTfZin(save, sourceIndex);
    } else if (save.typeName() == idTf) {
        st = smsigCore.addTf(save, sourceIndex);
    } else if (save.typeName() == idZin) {
        st = smsigCore.addZin(save, sourceIndex);
    } else if (save.typeName() == idYin) {
        st = smsigCore.addYin(save, sourceIndex);
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

template<> void SmallSignal<ACXFCore, ACXFData>::dump(std::ostream& os) const {
    Analysis::dump(os);
    os << "Analysis type: AC TF"<< std::endl;
    os << "OP analysis core:" << std::endl;
    opCore.dump(os);
    os << std::endl;
    os << "AC TF analysis core:" << std::endl;
    smsigCore.dump(os);
}

}
