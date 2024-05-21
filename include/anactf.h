#ifndef __ANACTF_DEFINED
#define __ANACTF_DEFINED

#include "ansmsig.h"
#include "coreop.h"
#include "coreactf.h"
#include "parameterized.h"
#include "common.h"


namespace NAMESPACE {

// AcTf analysis data
class AcTfData {
protected:
    KluComplexMatrix acMatrix; 
    Vector<Complex> acSolution;

    std::unordered_map<Id,size_t> sourceIndex;
    std::vector<Instance*> sources;
    Vector<Complex> tf;
    Vector<Complex> yin;
    Vector<Complex> zin;
};

// Constructor specialization
template<> SmallSignal<AcTfCore, AcTfData>::SmallSignal(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);

// Resolve save specialization
template<> bool SmallSignal<AcTfCore, AcTfData>::resolveSave(const PTSave& save, bool verify, Status& s);

// Dump specialization
template<> void SmallSignal<AcTfCore, AcTfData>::dump(std::ostream& os) const;

// Typedef AcTf
typedef SmallSignal<AcTfCore, AcTfData> AcTf;

}

#endif
