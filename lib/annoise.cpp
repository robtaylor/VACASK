#include "annoise.h"
#include "common.h"


namespace NAMESPACE {

template<> SmallSignal<NoiseCore, NoiseData>::SmallSignal(Id name, Circuit& circuit, PTAnalysis& ptAnalysis) 
    : Analysis(name, circuit, ptAnalysis), 
      opCore(*this, params.core().opParams, circuit, jac, solution, states), 
      smsigCore(*this, params.core(), opCore, contributionOffset, circuit, jac, 
      solution, states, acMatrix, acSolution, results, powerGain, outputNoise) {
}

template<> bool SmallSignal<NoiseCore, NoiseData>::resolveSave(const PTSave& save, bool verify, Status& s) {
    // Noise saves
    static const auto idDefault = Id("default");
    static const auto idFull = Id("full");
    static const auto idN  = Id("n");
    static const auto idNc = Id("nc");

    bool st = true;
    bool handled = true;
    if (save.typeName() == idDefault) {
        st = smsigCore.addAllNoiseContribInst(save, false);
    } else if (save.typeName() == idFull) {
        st = smsigCore.addAllNoiseContribInst(save, true);
    } else if (save.typeName() == idN) {
        st = smsigCore.addNoiseContribInst(save, false);
    } else if (save.typeName() == idNc) {
        st = smsigCore.addNoiseContribInst(save, true);
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

template<> void SmallSignal<NoiseCore, NoiseData>::dump(std::ostream& os) const {
    Analysis::dump(os);
    os << "Analysis type: noise analysis"<< std::endl;
    os << "OP analysis core:" << std::endl;
    opCore.dump(os);
    os << std::endl;
    os << "Noise analysis core:" << std::endl;
    smsigCore.dump(os);
}

}
