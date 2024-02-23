#ifndef __ANCORE_DEFINED
#define __ANCORE_DEFINED

#include <unordered_map>
#include "circuit.h"
#include "output.h"
#include "common.h"

namespace std {
template <>
    class hash <std::pair<NAMESPACE::Id,NAMESPACE::Id>> {
    public :
        size_t operator()(const pair<NAMESPACE::Id,NAMESPACE::Id> &x) const {
            size_t h = std::hash<int>()(x.first) ^ std::hash<int>()(x.second);
            return  h;
        };
    };
}

namespace NAMESPACE {

class Analysis;

// Analysis core, one analysis can have multiple analysis cores (i.e. op, tran, ...)
class AnalysisCore {
public: 
    AnalysisCore(Analysis& analysis, Circuit& circuit);
    
    AnalysisCore           (const AnalysisCore&)  = delete;
    AnalysisCore           (      AnalysisCore&&) = delete;
    AnalysisCore& operator=(const AnalysisCore&)  = delete;
    AnalysisCore& operator=(      AnalysisCore&&) = delete;

    // Clear output descriptors
    void clearOutputDescriptors();

    // Add an output descriptor to descriptors list of the core analysis
    bool addOutputDescriptor(const OutputDescriptor& descr);
    bool addOutputDescriptor(OutputDescriptor&& descr);

    // Add output descriptors that are not based on saves but are specific 
    // to analysis core (e.g. frequency, time). By default add nothing. 
    bool addCoreOutputDescriptors(Status& s=Status::ignore);

    // Add default output descriptors if no save has been provided
    bool addDefaultOutputDescriptors(Status& s=Status::ignore) { return true; };
    
    // Resolve all output descriptors into output sources
    // Delegate resolving of unknown decriptors to analysis
    bool resolveOutputDescriptors(bool strict, Status &s) { return true; };

    // Core functionality

    // The following two are called before output descriptors are added or resolved
    
    // Check if we need to add sparsity map or states vector entries 
    // Return value: ok, need mapping
    std::tuple<bool, bool> preMapping(Status& s=Status::ignore) { return std::make_tuple(true, false); }; 

    // Add sparsity map and states vector entries, set up forces on NR solver that require extradiagonal elements
    bool populateStructures(Status& s=Status::ignore) { return true; };

    // Called before core is run
    // - calls rebuild() for Jacobians
    // - binds instances to Jacobian entries 
    // - calls NRSolver::rebuild() (if a NR solver is used)
    bool rebuild(Status& s=Status::ignore) { return true; }; 
    
    // Called before core is run (and once per sweep) to initalize output files
    bool initializeOutputs(Id name, Status& s=Status::ignore) { return true; };

    // Runs the core
    bool run(bool continuePrevious, Status& s=Status::ignore) { return true; };

    // Called after core is run (and once per sweep) to finalieze and close output files
    bool finalizeOutputs(Status& s=Status::ignore) { return true; };

    // Called if analysis fails to remove output files
    bool deleteOutputs(Id name, Status& s=Status::ignore) { return true; };

    // Core state storage (used by continuation in sweeps), do nothing by default
    size_t stateStorageSize() const { return 0; };
    void resizeStateStorage(size_t n) { return; };
    bool storeState(size_t ndx) { return true; };
    bool restoreState(size_t ndx) { return true; };
    void makeStateIncoherent(size_t ndx) { return; };

    // Dump internals
    void dump(std::ostream& os) const;

    // Common handlers for save directive -> output descriptor(s) 
    bool addAllUnknowns(const PTSave& save, bool verify, Status& s=Status::ignore);
    bool addAllNodes(const PTSave& save, bool verify, Status& s=Status::ignore);
    bool addNode(const PTSave& save, bool verify, Status& s=Status::ignore);
    bool addFlow(const PTSave& save, bool verify, Status& s=Status::ignore);
    bool addInstanceOpvar(const PTSave& save, bool verify, Status& s=Status::ignore);
    bool addAllTfZin(const PTSave& save, bool verify, std::unordered_map<Id,size_t>& nameMap, Status& s=Status::ignore);
    bool addTf(const PTSave& save, bool verify, std::unordered_map<Id,size_t>& nameMap, Status& s=Status::ignore);
    bool addZin(const PTSave& save, bool verify, std::unordered_map<Id,size_t>& nameMap, Status& s=Status::ignore);
    bool addYin(const PTSave& save, bool verify, std::unordered_map<Id,size_t>& nameMap, Status& s=Status::ignore);
    bool addAllNoiseContribInst(const PTSave& save, bool verify, bool details, Status& s=Status::ignore);
    bool addNoiseContribInst(const PTSave& save, bool verify, bool details, Status& s=Status::ignore);
    
    // Common handlers for output descriptor -> output source
    bool addRealVarOutputSource(bool strict, Id name, const Vector<double>& solution, Status& s=Status::ignore);
    bool addRealVarOutputSource(bool strict, Id name, const VectorRepository<double>& solution, Status& s=Status::ignore);
    bool addComplexVarOutputSource(bool strict, Id name, const Vector<Complex>& solution, Status& s=Status::ignore);
    bool addComplexVarOutputSource(bool strict, Id name, const VectorRepository<Complex>& solution, Status& s=Status::ignore);
    bool addOpvarOutputSource(bool strict, Id instance, Id opvar, Status& s=Status::ignore);

protected:
    std::tuple<bool, UnknownIndex, UnknownIndex> getOutput(Value& v, Status& s=Status::ignore);
    std::tuple<bool, Instance*> getInput(Id name, Status& s=Status::ignore);
    Analysis& analysis;
    Circuit& circuit;
    Output::DescriptorList outputDescriptors;
    Output::SourcesList outputSources;
    std::unordered_map<Id,size_t> outputDescriptorIndices;
    // Number of save directives that produced at least one output descriptor
    size_t savesCount; 
};


}

#endif
