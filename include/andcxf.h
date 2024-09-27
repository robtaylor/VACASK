#ifndef __ANDCXF_DEFINED
#define __ANDCXF_DEFINED

#include "ansmsig.h"
#include "coreop.h"
#include "coredcxf.h"
#include "parameterized.h"
#include "common.h"


namespace NAMESPACE {

// AC analysis data
class DCXFData {
protected:
    std::unordered_map<Id,size_t> sourceIndex;
    std::vector<Instance*> sources; 
    Vector<double> incrementalSolution; // Incremental solution
    Vector<double> tf;
    Vector<double> yin;
    Vector<double> zin;
};

// Constructor specialization
template<> SmallSignal<DCXFCore, DCXFData>::SmallSignal(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);

// Resolve save specialization
template<> bool SmallSignal<DCXFCore, DCXFData>::resolveSave(const PTSave& save, bool verify, Status& s);

// Dump specialization
template<> void SmallSignal<DCXFCore, DCXFData>::dump(std::ostream& os) const;

// Typedef AC
typedef SmallSignal<DCXFCore, DCXFData> DCXF;

}

#endif
