#ifndef __ANDCTF_DEFINED
#define __ANDCTF_DEFINED

#include "ansmsig.h"
#include "coreop.h"
#include "coredctf.h"
#include "parameterized.h"
#include "common.h"


namespace NAMESPACE {

// Ac analysis data
class DcTfData {
protected:
    std::unordered_map<Id,size_t> sourceIndex;
    std::vector<Instance*> sources; 
    Vector<double> incrementalSolution; // Incremental solution
    Vector<double> tf;
    Vector<double> yin;
    Vector<double> zin;
};

// Constructor specialization
template<> SmallSignal<DcTfCore, DcTfData>::SmallSignal(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);

// Resolve save specialization
template<> bool SmallSignal<DcTfCore, DcTfData>::resolveSave(const PTSave& save, bool verify, Status& s);

// Dump specialization
template<> void SmallSignal<DcTfCore, DcTfData>::dump(std::ostream& os) const;

// Typedef Ac
typedef SmallSignal<DcTfCore, DcTfData> DcTf;

}

#endif
