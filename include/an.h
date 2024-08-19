#ifndef __ANALYSIS_DEFINED
#define __ANALYSIS_DEFINED

#include <unordered_map>
#include "circuit.h"
#include "parameterized.h"
#include "parseroutput.h"
#include "output.h"
#include "answeep.h"
#include "generator.h"
#include "common.h"

namespace NAMESPACE {

// Generic analysis
class Analysis {
public:
    typedef Analysis* (*AnalysisFactory)(PTAnalysis& ptAnalysis, Circuit& circuit, Status& s);

    // Aborted .. run() -> do nothing, return Aborted
    // Finished .. run() -> do nothing, return Finished
    // Stopped .. run() -> continue, return status
    enum class State { Uninitilized, Ready, Aborted, Stopped, Finished };
    
    Analysis(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);
    virtual ~Analysis();

    Analysis           (const Analysis&)  = delete;
    Analysis           (      Analysis&&) = delete;
    Analysis& operator=(const Analysis&)  = delete;
    Analysis& operator=(      Analysis&&) = delete;

    State state() const { return state_; };

    Id name() const { return name_; };
    IStruct<SimulatorOptions>& simulatorOptions() { return simOptions; };

    // Converts an output descriptor to output source and stores it in the given output sources list
    // This handles output descriptors that are not specific for an analysis core, 
    // i.e. it is called by a core when resolving a descriptor is delegated to the analysis
    bool resolveOutputDescriptor(const OutputDescriptor& descr, Output::SourcesList& srcs, bool strict);

    // Sweep API
    size_t sweepCount() const { return ptAnalysis.sweeps().data().size(); };
    bool updateSweeper(Status& s=Status::ignore);
    
    // Setup api
    void setSaves(PTSavesVector* commonSaves);
    void setParametrization(const PTParameterMap* optionsMap);

    // Interface method, does sweeping if needed
    Generator<Analysis::State> run(Status& s=Status::ignore);

    static bool registerFactory(Id name, AnalysisFactory factory);

    static Analysis* create(
        PTAnalysis& ptAnalysis, 
        PTSavesVector* commonSaves, 
        const PTParameterMap* optionsMap, 
        Circuit& circuit, Status& s=Status::ignore
    );
    
    virtual Parameterized& parameters() = 0; 
    virtual const Parameterized& parameters() const = 0; 

    // Updates analysis parameters that are defined as expressions
    // Return value: ok, changed
    std::tuple<bool, bool> updateParameterExpressions(Status& s=Status::ignore); 

    // Mechanism for adding entries to sparsity map and states vector
    // (part of setSweepState() and setAnalysisOptions())
    // used by analyses (and their cores)
    // 
    // 1) check if entries are needed 
    //    sparsity map is pre-filled by circuit elaborate() and/or previous analysis run
    //    Return value: ok, need mapping
    virtual std::tuple<bool, bool> preMapping(Status& s=Status::ignore) { return std::make_tuple(true, false); };

    // 2) fill sparsity map and states vector with entries specific to analysis
    virtual bool populateStructures(Status& s=Status::ignore) { return true; };
    
    virtual void dump(std::ostream& os) const;

protected:
    Id name_;
    Circuit& circuit;
    ParameterSweeper sweeper;
    IStruct<SimulatorOptions> originalSimOptions;
    IStruct<SimulatorOptions> simOptions;
    PTAnalysis& ptAnalysis;
    
    // Analysis::run() steps:

    // Call Sweeper::bind() 
    
    // Call Sweeper::storeState() before sweep starts to store circuit state

    // Call setSweepState() once per sweep point to set parameters
    // If not sweeping, call setAnalysisOptions()

    // Adds and optionally checks output descriptors (if strictsave simulator option is set)
    // Called by toplevel analysis function (an.cpp) when setting data sources 
    // before analysis is run.  
    // - clears all output descriptors by calling clearOutputDescriptors()
    // - adds common output descriptors (save variables) by calling addOutputDescriptors()
    // - adds core output descriptors (for time, frequency) by calling addCoreOutputDescriptors()
    // - interprets save directives by calling resolveSave()
    // - adds default output descriptors to cores if saves produce no output descriptors
    //   by calling addDefaultOutputDescriptors()
    bool addOutputDescriptors(Status& s=Status::ignore);

    // Clears all output descriptors
    // Called by addOutputDescriptors() before output descriptors are added. 
    // Calls clearOutputDescriptors() method of all analysis cores. 
    virtual void clearOutputDescriptors() = 0;
    
    // Adds an output descriptor that is common to all analysis cores (e.g. sweep variables)
    // Called by toplevel analysis function (an.cpp) when setting up sweep variable output. 
    // Calls addOutputDescriptor() method of all analysis cores. 
    virtual bool addCommonOutputDescriptor(const OutputDescriptor& desc) = 0; 

    // Adds core output descriptors of analysis cores (e.g. time, frequency). 
    // Called by addOutputDescriptors() before save directives are interpreted. 
    // Calls addCoreOutputDescriptors() method of all analysis cores. 
    virtual bool addCoreOutputDescriptors(Status& s=Status::ignore) = 0;

    // Resolves a save directive into one or more output descriptors and
    // stores them in output descriptor lists of analysis cores 
    // Called by addOutputDescriptors() when interpreting a save directive. 
    // Reimplemented by each analysis so that the dispatch of saves to cores 
    // can be customized. 
    virtual bool resolveSave(const PTSave& save, bool verify, Status& s=Status::ignore) = 0;

    // Adds default output descriptors to cores for which save directives 
    // produced no output descriptors. 
    // Called by addOutputDescriptors() after all save directives were interpreted. 
    // Calls addDefaultOutputDescriptors() method of all analysis cores. 
    virtual bool addDefaultOutputDescriptors() = 0;
    
    // Resolve output descriptors of all analysis cores
    // Called by toplevel analysis function (an.cpp) when setting data sources 
    // before analysis is run.  
    // Calls resolveOutputDescriptors() method of all analysis cores. 
    virtual bool resolveOutputDescriptors(bool strict, Status& s=Status::ignore) = 0; 

    // Rebuild internal structures
    virtual bool rebuildCores(Status& s=Status::ignore) = 0;

    // Initialize core output (once per sweep run)
    virtual bool initializeOutputs(Status& s=Status::ignore) = 0;

    // Call Sweeper::restoreState() after initializeOutputs() if a sweep wraparound happened
    
    // Core analysis method
    virtual bool runCores(bool continuePrevious, Status& s=Status::ignore) = 0;

    // Call Sweeper::advance() after runCores()

    // Call Sweeper::storeState()

    // Finalize core output (once per sweep run)
    virtual bool finalizeOutputs(Status& s=Status::ignore) = 0;

    // Delete output files 
    // Called after analysis fails if strictoutput options is set
    virtual bool deleteOutputs(Status& s=Status::ignore) = 0;

    // Call setSweepState() after finalizeOtuputs() and deleteOutputs() to restore circuit state
    // If not sweeping, call setAnalysisOptions() 

    // Analysis state repository size/resize
    // Used for homotopy in sweeps
    virtual size_t analysisStateStorageSize() const = 0;
    virtual void resizeAnalysisStateStorage(size_t n) = 0;
    
    // Store analysis state in internal repository 
    // Used for homotopy in sweeps
    virtual bool storeState(size_t ndx) = 0;

    // Restore analysis state from internal repository 
    // Used for homotopy in sweeps
    virtual bool restoreState(size_t ndx) = 0;

    // Make analysis state incoherent with current topology
    // Used for homotopy in sweeps
    // In OP analysis sweeps this forces the use of nodesets 
    // instead of previous solution. 
    virtual void makeStateIncoherent(size_t ndx) = 0;

    PTParameterMap parameterizedOptions;

    // Generator that runs all involved cores, replaces run
    // Analysis::runCoresEngine(Status& s=Status::ignore) 
    
    // Generator that runs a core
    // Core::run(Status& s=Status::ignore)

    
private:
    State state_;
    static std::unordered_map<Id,AnalysisFactory>& getRegistry() {
        static std::unordered_map<Id,AnalysisFactory> factoryMap;
        return factoryMap;
    };
    const PTSavesVector* commonSaves;
    Int advancedSweepIndex;
};

}

#endif
