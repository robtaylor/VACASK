#ifndef __ANDCINC_DEFINED
#define __ANDCINC_DEFINED

#include "ansmsig.h"
#include "coreop.h"
#include "coredcinc.h"
#include "parameterized.h"
#include "common.h"


namespace NAMESPACE {

// DcIncr analysis data
class DcIncrData {
protected:
    Vector<double> incrementalSolution; // Incremental solution
};

// Constructor specialization
template<> SmallSignal<DcIncrementalCore, DcIncrData>::SmallSignal(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);

// Resolve save specialization
template<> bool SmallSignal<DcIncrementalCore, DcIncrData>::resolveSave(const PTSave& save, bool verify, Status& s);

// Dump specialization
template<> void SmallSignal<DcIncrementalCore, DcIncrData>::dump(std::ostream& os) const;

// Typedef DcIncremental
typedef SmallSignal<DcIncrementalCore, DcIncrData> DcIncremental;

}

#endif
