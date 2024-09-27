#ifndef __ANACXF_DEFINED
#define __ANACXF_DEFINED

#include "ansmsig.h"
#include "coreop.h"
#include "coreacxf.h"
#include "parameterized.h"
#include "common.h"


namespace NAMESPACE {

// ACXF analysis data
class ACXFData {
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
template<> SmallSignal<ACXFCore, ACXFData>::SmallSignal(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);

// Resolve save specialization
template<> bool SmallSignal<ACXFCore, ACXFData>::resolveSave(const PTSave& save, bool verify, Status& s);

// Dump specialization
template<> void SmallSignal<ACXFCore, ACXFData>::dump(std::ostream& os) const;

// Typedef ACXF
typedef SmallSignal<ACXFCore, ACXFData> ACXF;

}

#endif
