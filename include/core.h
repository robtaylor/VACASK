#ifndef __ANCORE_DEFINED
#define __ANCORE_DEFINED

#include <unordered_map>
#include "circuit.h"
#include "output.h"
#include "hash.h"
#include "common.h"

namespace NAMESPACE {

typedef struct IdPairHash {
    size_t operator()(const std::pair<Id,Id> &x) const {
        return hash_val(x.first.id(), x.second.id());
    }
} IdPairHash;

class Analysis;

// Analysis core, one analysis can have multiple analysis cores (i.e. op, tran, ...)
class AnalysisCore {
public: 
    enum class Error {
        OK, 
        Arguments, 
        NodeNotFound, 
        OpvarNotFound, 
        InstanceNotFound, 
        OutputSpec, 
        OutputType, 
        InstanceNotSource, 
        Descriptor, 
    };

    AnalysisCore(Analysis& analysis, Circuit& circuit);
    
    AnalysisCore           (const AnalysisCore&)  = delete;
    AnalysisCore           (      AnalysisCore&&) = delete;
    AnalysisCore& operator=(const AnalysisCore&)  = delete;
    AnalysisCore& operator=(      AnalysisCore&&) = delete;

    // Clear error
    void clearError() { lastError = Error::OK; }; 

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore) const; 

    // Clear output descriptors
    void clearOutputDescriptors();

    // Add an output descriptor to descriptors list of the core analysis
    bool addOutputDescriptor(const OutputDescriptor& descr);
    bool addOutputDescriptor(OutputDescriptor&& descr);

    // Add output descriptors that are not based on saves but are specific 
    // to analysis core (e.g. frequency, time). By default add nothing. 
    bool addCoreOutputDescriptors() { return true; };

    // Add default output descriptors if no save has been provided
    bool addDefaultOutputDescriptors() { return true; };
    
    // Resolve all output descriptors into output sources
    // Delegate resolving of unknown decriptors to analysis
    bool resolveOutputDescriptors(bool strict) { return true; };

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
    bool addAllUnknowns(const PTSave& save);
    bool addAllNodes(const PTSave& save);
    bool addNode(const PTSave& save);
    bool addFlow(const PTSave& save);
    bool addInstanceOpvar(const PTSave& save);
    bool addAllTfZin(const PTSave& save, std::unordered_map<Id,size_t>& nameMap);
    bool addTf(const PTSave& save, std::unordered_map<Id,size_t>& nameMap);
    bool addZin(const PTSave& save, std::unordered_map<Id,size_t>& nameMap);
    bool addYin(const PTSave& save, std::unordered_map<Id,size_t>& nameMap);
    bool addAllNoiseContribInst(const PTSave& save, bool details);
    bool addNoiseContribInst(const PTSave& save, bool details);
    
    // Common handlers for output descriptor -> output source
    // TODO: handle error formatting in callers
    bool addRealVarOutputSource(bool strict, Id name, const Vector<double>& solution);
    bool addRealVarOutputSource(bool strict, Id name, const VectorRepository<double>& solution);
    bool addComplexVarOutputSource(bool strict, Id name, const Vector<Complex>& solution);
    bool addComplexVarOutputSource(bool strict, Id name, const VectorRepository<Complex>& solution);
    bool addOpvarOutputSource(bool strict, Id instance, Id opvar);

protected:
    enum Error lastError;
    Int errorExpectedArgCount;
    Id errorId;
    Id errorId2;

    std::tuple<bool, UnknownIndex, UnknownIndex> getOutput(Value& v);
    std::tuple<bool, Instance*> getInput(Id name);
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
