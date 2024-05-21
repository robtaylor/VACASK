#ifndef __ANNOISE_DEFINED
#define __ANNOISE_DEFINED

#include "ansmsig.h"
#include "coreop.h"
#include "corenoise.h"
#include "parameterized.h"
#include "common.h"


namespace NAMESPACE {

// Noise analysis data
class NoiseData {
protected:
    KluComplexMatrix acMatrix; 
    Vector<Complex> acSolution;

    std::unordered_map<std::pair<Id, Id>, size_t> contributionOffset; 
    Vector<double> results;  
    double powerGain;
    double outputNoise;
};

// Constructor specialization
template<> SmallSignal<NoiseCore, NoiseData>::SmallSignal(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);

// Resolve save specialization
template<> bool SmallSignal<NoiseCore, NoiseData>::resolveSave(const PTSave& save, bool verify, Status& s);

// Dump specialization
template<> void SmallSignal<NoiseCore, NoiseData>::dump(std::ostream& os) const;

// Typedef Ac
typedef SmallSignal<NoiseCore, NoiseData> Noise;

}

#endif
