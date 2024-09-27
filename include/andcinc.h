#ifndef __ANDCINC_DEFINED
#define __ANDCINC_DEFINED

#include "ansmsig.h"
#include "coreop.h"
#include "coredcinc.h"
#include "parameterized.h"
#include "common.h"


namespace NAMESPACE {

// DcIncr analysis data
class DCIncrementalData {
protected:
    Vector<double> incrementalSolution; // Incremental solution
};

// Constructor specialization
template<> SmallSignal<DCIncrementalCore, DCIncrementalData>::SmallSignal(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);

// Resolve save specialization
template<> bool SmallSignal<DCIncrementalCore, DCIncrementalData>::resolveSave(const PTSave& save, bool verify, Status& s);

// Dump specialization
template<> void SmallSignal<DCIncrementalCore, DCIncrementalData>::dump(std::ostream& os) const;

// Typedef DCIncremental
typedef SmallSignal<DCIncrementalCore, DCIncrementalData> DCIncremental;

}

#endif
