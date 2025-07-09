#ifndef __DEVBASE_DEFINED
#define __DEVBASE_DEFINED

#include "identifier.h"
#include "flags.h"
#include "status.h"
#include "parseroutput.h"
#include "node.h"
#include "ansupport.h"
#include "parameterized.h"
#include "output.h"
#include "klumatrix.h"
#include "parseroutput.h"
#include "options.h"
#include "flags.h"
#include "coretrancoef.h"
#include "elsetup.h"
#include "value.h"
#include "common.h"
#include <cinttypes>


namespace NAMESPACE {
// Parameters are enumerated with 0..parameterCount-1
// Multiple Ids can map to the same parameter index
// Instance and model parameters have two independent zero-based enumerations

class Circuit;

class Model;
class Instance;

enum class DeviceFlags : uint16_t {
    None = 0, 
    // Device is valid
    IsValid = 2, 
    // Device generates DC incremental excitation
    GeneratesDCIncremental = 4, 
    // Device generates AC small-signal excitation
    GeneratesAC = 8, 
    // Allows bypass
    Bypassable = 16, 
};
DEFINE_FLAG_OPERATORS(DeviceFlags);

// The API tries to avoid many virtual method calls during mapping/setup/eval/load

class Device : public FlagBase<DeviceFlags> {
public:
    friend class Model;

    Device(Id name, const Loc& location=Loc::bad); 
    virtual ~Device();

    Device           (const Device&)  = delete;
    Device           (      Device&&) = default;
    Device& operator=(const Device&)  = delete;
    Device& operator=(      Device&&) = default;

    // Checks if this device is the same as some other device
    virtual bool operator==(const Device& other) const = 0;

    // Device location in input files
    const Loc location() const { return loc; };

    // Device name
    Id name() const { return name_; };
    
    // Does this device have subinstances/submodels (internal nodes do not count)? 
    // By default a device is not hierarchical. 
    virtual bool isHierarchical() const { return false; }; 

    // Clears the list of all models of this device. 
    void clearModelList() { models_.clear(); instanceCount_ = 0; };

    // Adds a model to model list, makes sure only models of this device are added
    bool addModel(Model* model); 

    // Expose models storage
    std::vector<Model*>& models() { return models_; };
    
    // Number of models
    size_t modelCount() const { return models_.size(); };

    // Returns a model
    Model* model(size_t ndx) { return models_[ndx]; };

    // Number of instances
    size_t instanceCount() const { return instanceCount_; };

    // Sets parameter defaults (model, instance), computes node collapsing (instance)
    // Return value: ok, unknowns changed, sparsity changed
    virtual std::tuple<bool, bool, bool> setup(Circuit& cir, bool force, DeviceRequests* devReq, Status& s=Status::ignore) { return std::make_tuple(true, false, false); };
    
    // Collapses nodes
    // Return value: ok, changed
    virtual bool collapseNodes(Circuit& cir, Status& s=Status::ignore) { return true; };
    
    // Populates Jacobian with entries, allocates slots in states vector
    virtual bool populateStructures(Circuit& cir, Status& s=Status::ignore) { return true; };

    // Does computations immediately before analysis is run
    // The computations are based on other instances and models that need to be in set up state
    virtual bool preAnalysis(Circuit& circuit, Status& s=Status::ignore) { return true; };
    
    // Stores pointers to variables where the Jacobian values will be written at load
    // Uses virtual methods to correctly handle real and complex matrices. 
    virtual bool bind(
        Circuit& cir, 
        KluMatrixAccess* matResist, Component compResist, const std::optional<MatrixEntryPosition>& mepResist, 
        KluMatrixAccess* matReact, Component compReact, const std::optional<MatrixEntryPosition>& mepReact, 
        Status& s=Status::ignore
    ) { return true; };

    // Evaluate instances and load results
    virtual bool eval(Circuit& circuit, EvalSetup& evalSetup) { return true; };
    
    // Load instance evaluation results into linear system
    virtual bool load(Circuit& circuit, LoadSetup& loadSetup) { return true; };

    // Evaluate instances and load results
    virtual bool evalAndLoad(Circuit& circuit, EvalSetup* evalSetup, LoadSetup* loadSetup) { return true; };
    
    // A model created with this method is owned by the circuit. 
    // No need to delete it manually, it will get deleted when the circuit is deleted. 
    // The evaluator's latest context must be the parent instance's context. 
    virtual Model* createModel(Circuit& cir, Instance* parentInstance, RpnEvaluator& evaluator, const PTModel& parsedModel, Status& s=Status::ignore) = 0;

    // Source API
    // Is this device an independent source? 
    virtual bool isSource() const { return false; };

    // Is it a voltage source? 
    virtual bool isVoltageSource() const { return false; };
    
    // Dumps the device
    virtual void dump(int indent, std::ostream& os) const {};

    // For local timing
    double tovh;
    size_t novh;

protected:
    Id name_;
    const Loc& loc;
    std::vector<Model*> models_;
    size_t instanceCount_;
};


enum class ModelFlags : uint8_t {
    // Model needs to be set up
    NeedsSetup = 1, 
    // Model is valid
    IsValid = 2, 
};
DEFINE_FLAG_OPERATORS(ModelFlags);

// Auxiliary data used in hierarchical instantiation.
// For now it is used only for detecting hierachical recursion.
typedef struct InstantiationData {
public:
    InstantiationData();
    InstantiationData(Instance*);

    InstantiationData           (const InstantiationData&)  = delete;
    InstantiationData           (      InstantiationData&&) = default;
    InstantiationData& operator=(const InstantiationData&)  = delete;
    InstantiationData& operator=(      InstantiationData&&) = default;

    bool addAncestor(Instance* inst);
    bool removeAncestor(Instance* inst);
    bool isAncestor(Model* model);
    Id translateDefinition(Id definitionName);

private:
    std::unordered_set<Model*> ancestorModels_;
    std::vector<std::string> ancestorPathStack_;
} InstantiationData;


class Model : public Parameterized, public FlagBase<ModelFlags> {
public:
    Model(Device* device, Id name, Instance* parent, const PTModel& parsedModel);
    virtual ~Model();

    Model           (const Model&)  = delete;
    Model           (      Model&&) = default;
    Model& operator=(const Model&)  = delete;
    Model& operator=(      Model&&) = default;

    // Location of this model in input files
    const Loc location() const { return parsedModel.location(); };
    
    // Name of this model
    Id name() const { return name_; };
    
    // Clears the list of all instances of this model. 
    void clearInstanceList() { instances_.clear(); };

    // Adds an instance to instances list, makes sure only instances of this model are added
    bool addInstance(Instance* instance); 

    // Expose instances storage
    std::vector<Instance*>& instances() { return instances_; };
    
    // Number of instances
    size_t instanceCount() const { return instances_.size(); };

    // Returns an instance
    Instance* instance(size_t ndx) { return instances_[ndx]; };

    // Returns the generic device to which this model belongs. 
    Device* device() { return device_; };
    const Device* device() const { return device_; };

    // Returns the parent instance of this model, i.e. 
    // where in the hierarchy is this model positioned. 
    Instance* parent() { return parent_; };

    // Abstract members inherited from Parameterized
    //   parameterCount()
    //   parameterIndex()
    //   parameterName()
    //   parameterType() .. from parameter index
    //   getParameter()  .. by parameter index
    //   setParameter()  .. by parameter index
    //   parameterGiven(); .. by parameter index
    
    // Checks if a model parameter is not specified as an expression. 
    bool parameterIsFree(Id name); 
    
    // Sets parameter defaults (model, instance), computes node collapsing (instance)
    // Return value: ok, unknowns changed, sparsity changed
    virtual std::tuple<bool, bool, bool> setup(Circuit& cir, bool force, DeviceRequests* devReq, Status& s=Status::ignore) { return std::make_tuple(true, false, false); };

    // An instance created with this method is owned by the circuit. 
    // No need to delete it manually, it will get deleted when the circuit is deleted. 
    // The evaluator's latest context must be the parent instance's context. 
    // If externalContext is given it is inserted into the context stack and 
    // becomes the context in which the expressions of dependent parameters, subinstances,
    // and submodels are evaluated. External context is not destroyed on context exit and
    // is never reallocated or moved (the pointer remains valid). 
    // Otherwise a new disposable context is created. 
    virtual Instance* createInstance(Circuit& cir, Instance* parentInstance, RpnEvaluator& evaluator, Context* externalContext, const PTInstance& parsedInstance, InstantiationData& idata, Status& s=Status::ignore) = 0;
    
    // Dumps the model
    virtual void dump(int indent, std::ostream& os) const {};

protected:
    Id name_;
    Device* device_;
    Instance* parent_;
    Flags flags_;
    std::vector<Instance*> instances_;
    const PTModel& parsedModel;
};


enum class InstanceFlags : uint8_t {
    // Parameters changed, need to propagate them
    ParamsChanged = 1, 
    // Instance needs to be set up
    NeedsSetup = 2, 
    // Instance was setup up before, we can compare node collapsing to previous state
    HasSetupHistory = 4, 
    // Instance is valid
    IsValid = 8, 
    // Limiting applied
    LimitingApplied = 16, 
    // Has device history, we can check for convergence
    HasDeviceHistory = 32, 
    // Device output is converged
    OutputConverged = 64, 
    // Device is bypassed
    Bypassed = 128, 
};
DEFINE_FLAG_OPERATORS(InstanceFlags);

class Instance : public Parameterized, public FlagBase<InstanceFlags> {
public:
    // This would be much more elegant with coroutines. 
    // Operator ++
    // - for hierarchical instance go to first child
    // - for ordinary instance go to next peer
    // - at end of peers go up until next peer is found
    // Start iterator = [[Instance*], 0]
    // End iterator = stack = []
    struct HierarchicalIterator {
        using position_type     = std::pair<std::vector<Instance*>*,size_t>;
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = size_t;
        using value_type        = Instance;
        using pointer           = Instance*;
        using reference         = Instance&;

        // Construct end iterator
        HierarchicalIterator();
        // Construct begin iterator
        HierarchicalIterator(pointer instance);
        // *it
        reference operator*() const;
        // it->
        pointer operator->();
        // ++it
        HierarchicalIterator& operator++();
        // it++
        HierarchicalIterator operator++(int);

        friend bool operator==(const HierarchicalIterator& a, const HierarchicalIterator& b);
        friend bool operator!=(const HierarchicalIterator& a, const HierarchicalIterator& b);

        void stopDescent();

        // Hierarchy traversal stack
        std::vector<position_type> stack;
        // Dummy entry for the parent of starting instance
        std::vector<Instance*> dummy;
        // Force move to peer at next increment
        bool forcePeer;
    };

    Instance(Model* model, Id name, Instance* parent, const PTInstance& parsedInstance);
    virtual ~Instance();

    Instance           (const Instance&)  = delete;
    Instance           (      Instance&&) = default;
    Instance& operator=(const Instance&)  = delete;
    Instance& operator=(      Instance&&) = default;

    // Location of this instance in input files
    const Loc location() const { return parsedInstance.location(); };
    
    // Name of this instance
    Id name() const { return name_; };
    
    // Returns the generic model to which this instance belongs. 
    Model* model() { return model_; };
    const Model* model() const { return model_; };
    
    // Returns the parent instance of this instance, i.e. 
    // where in the hierarchy is this instance positioned. 
    Instance* parent() { return parent_; };
    
    // Abstract members inherited from Parameterized
    //   parameterCount()
    //   parameterIndex()
    //   parameterName()
    //   parameterType() .. from parameter index
    //   getParameter()  .. by parameter index
    //   setParameter()  .. by parameter index
    //   parameterGiven(); .. by parameter index
    
    // Principal parameter is the parameter which is used when no parameter is specified in sweep
    // By default there is none. 
    // Return value: principal parameter index, principal parameter exists
    virtual std::tuple<ParameterIndex, bool> principalParameterIndex() const { return std::make_tuple(0, false); };
    
    // Checks if an instance parameter is not specified as an expression. 
    bool parameterIsFree(Id name); 
    
    // Number of nodes of this instance >= number of terminals
    virtual TerminalIndex staticNodeCount() const = 0;

    // Number of terminals of this instance
    virtual TerminalIndex terminalCount() const = 0;

    // Node index
    virtual std::tuple<TerminalIndex, bool> nodeIndex(Id name) const = 0;

    // Node name
    virtual Id nodeName(TerminalIndex ndx) const = 0;
    
    // Bind a node to a terminal
    virtual bool bindTerminal(TerminalIndex n, Node* node, Status& s=Status::ignore) = 0;

    // Return a node bound to a terminal
    virtual Node* terminal(TerminalIndex n, Status& s=Status::ignore) const = 0;

    // Unbind all terminals from nodes
    virtual bool unbindTerminals(Circuit& cir, Status& s=Status::ignore) = 0;

    // Hierarchy API 
    // Returns a vector of generic child instances
    // Returns nullptr if this instance is not hierarchical
    virtual std::vector<Instance*>* childInstances() { return nullptr; };

    // Returns a vector of generic child models
    // Returns nullptr if this instance is not hierarchical
    virtual std::vector<Model*>* childModels() { return nullptr; };

    // Returns an iterator to this instance's hierarchy beginning
    // The first instance returned by the iterator is this instance. 
    HierarchicalIterator hierarchyBegin() { return HierarchicalIterator(this); };

    // Returns an iterator to this instance's hierarchy end
    HierarchicalIterator hierarchyEnd() { return HierarchicalIterator(); };

    // Translates subinstance/submodel name by prepending it with instance's name followed by ':'. 
    // Toplevel instance (parent is nullptr) performs no translation. 
    Id translate(Id child);

    // Translates a child node. By default a node is prepended with instance's name followed by ':'. 
    // Global nodes are not translated. 
    // Subcircuits must take into account terminal names which must translate to 
    // the nodes bound to them. Therefore this must be virtual. 
    virtual Id translateNode(Circuit& cir, Id nodeName);

    // Translates a peer object 
    // Prepends it with parent instance's name or nothing if parent is toplevel instance
    Id translatePeer(Id peer);

    // Adds a child instance/model to the list of child instances/model of a hierarchical instance. 
    // By default does nothing. 
    // Maybe we should raise an exception because we should never end up here. 
    virtual void addChild(Instance* child) { return; };
    virtual void addChild(Model* child) { return; };

    // Enter evaluation context of the instance. 
    // If externalContext is specified it is inserted into the context stack. 
    // otherwise a new disposable context is created. 
    // By default does nothing. 
    // Hierarchical devices must call this function before any parameter propagation or 
    // hierarchy change check, or hierarchy change is performed. 
    // The returned contextMarker can be used to revert the context stack back to the state 
    // before enterContext() is called. 
    // On failure the context stack must be left in the same state as before call. 
    // Return value: ok, contextMarker
    virtual std::tuple<bool, size_t> enterContext(Circuit& circuit, Context* externalContext, bool addToPath, bool rebuild, Status& s=Status::ignore);

    // Leave evaluation context of the instance. 
    // Can return to an arbitrary context given by stackMarker. 
    // If the context is external it is simply unlinked, but not destroyed. 
    // By default does nothing. 
    // Hierarchical devices must call this after all parameter propagation and 
    // topology change are finished. 
    static bool revertContext(Circuit& circuit, size_t contextMarker);
    
    // Check if current parameter values require the subhierarchy of the instance to be rebuilt. 
    // The context of the instance must be established with enterContext(). 
    // By default reports subhierarchy has not changed. 
    // Return value: statusOk, subhierarchyChanged
    virtual std::tuple<bool, bool> subhierarchyChanged(Circuit& circuit, Status& s=Status::ignore) { return std::make_tuple(true, false); }; 

    // Propagation by default does nothing and clears the parameters changed flag
    // It is supposed to propagate parameter values to immediate descendents (instances and models)
    // The evaluator must have an established global context. 
    // Once the instance's context is established the function checks if the subhierarchy 
    // has changed due to parameters change. If yes, subhierarchy is deleted and rebuilt. 
    // If subhierarchy did not change, parameters are propagated to immediate descendants 
    // (instances and models). 
    // If parameters of descendant instances change their ParamsChanged flag is set 
    // indicating that they require propagation. 
    // In the end the ParamsChanged flag of this instance is cleared. 
    // By default we only clear the ParamsChanged flag as this is not a hierarchical instance. 
    virtual bool propagateParameters(Circuit& circuit, RpnEvaluator& evaluator, Status& s=Status::ignore) { clearFlags(Flags::ParamsChanged); return true; };
    
    // Deletes the subhierarchy of an instance
    // Return value: ok
    virtual bool deleteHierarchy(Circuit& cir, Status& s=Status::ignore) { return true; }; 
    
    // Builds the instance's subhierarchy
    // Evaluator's latest context must be the parent instance's context or the global context
    // Return value: ok
    virtual bool buildHierarchy(Circuit& cir, RpnEvaluator& evaluator, InstantiationData& idata, Status& s=Status::ignore) { return false; };  
    
    // Source API
    // Returns the positive and negative equation for excitation produced by this instance
    // if this instance is an independent source. By default returns (0,0) which 
    // means that the excitation will go to the bucket. 
    virtual std::tuple<EquationIndex,EquationIndex> sourceExcitation(Circuit& cir) const  { return std::make_tuple(0, 0); };

    // Returns the positive and negative unknown for getting the response at this instance
    // if this instance is an independent source. By default returns (0,0) which 
    // means that the response will come from the ground and will be 0.0.  
    virtual std::tuple<UnknownIndex,UnknownIndex> sourceResponse(Circuit& cir) const { return std::make_tuple(0, 0); };

    // When excitaton is set to this value the response will correspond to
    // the case when the magnitude parameter of the source is set to 1. 
    virtual double scaledUnityExcitation() const { return 1.0; };

    // Factor by which the response measured at the source must be scaled to 
    // reflect the total response of all parallel instances combined. 
    virtual double responseScalingFactor() const { return 1.0; };
    
    // Opvar API
    // Returns the number of opvars
    virtual ParameterIndex opvarCount() const { return 0; };

    // Returns the opvar index corresponding to opvar name
    virtual std::tuple<ParameterIndex, bool> opvarIndex(Id name) const { return std::make_tuple(0, false); };

    // Returns the opvar name for given opvar index
    virtual Id opvarName(ParameterIndex ndx) const { return Id(); };

    // Returns the opvar type for given opvar index
    // Return value: type, opvar exists
    virtual std::tuple<Value::Type,bool> opvarType(ParameterIndex ndx, Status& s=Status::ignore) const { return std::make_tuple(Value::Type::Int, false); };
    
    // Returns the opvar type for given opvar name
    // Return value: type, opvar exists
    virtual std::tuple<Value::Type,bool> opvarType(Id name, Status& s=Status::ignore) const;
    
    // Returns the value of an opvar given by index
    virtual bool getOpvar(ParameterIndex ndx, Value& v, Status& s=Status::ignore) const { return false; };

    // Returns the value of an opvar given by name
    virtual bool getOpvar(Id name, Value& v, Status& s=Status::ignore) const;

    // Returns an output source corresponding to an opvar with given index
    // Return value: ok, output source
    virtual std::tuple<bool, OutputSource> opvarOutputSource(ParameterIndex ndx) const { return std::make_tuple(false, OutputSource()); };
    
    // Noise API
    virtual ParameterIndex noiseSourceCount() const { return 0; };
    virtual ParameterIndex uniqueNoiseSourceCount() const { return 0; };
    virtual Id noiseSourceName(ParameterIndex ndx) const { return Id(); };
    virtual std::tuple<ParameterIndex, bool> uniqueNoiseSourceIndex(ParameterIndex ndx) const { return std::make_tuple(0, false); };
    virtual std::tuple<ParameterIndex, bool> uniqueNoiseSourceIndex(Id name) const { return std::make_tuple(0, false); };

    virtual std::tuple<EquationIndex, EquationIndex> noiseExcitation(Circuit& cir, ParameterIndex ndx) const { return std::make_tuple(0, 0); };
    virtual bool loadNoise(Circuit& circuit, double freq, double* noiseDensity) { return true; };

    // Sets parameter defaults, computes node collapsing
    // Return value: ok, unknowns changed, sparsity changed
    virtual std::tuple<bool, bool, bool> setup(Circuit& cir, const SimulatorOptions& opt, SimulatorInternals& internals, bool force, Status& s=Status::ignore) { return std::make_tuple(true, false, false); };
    
    // Dumps the instance
    virtual void dump(int indent, const Circuit& cir, std::ostream& os) const {};

protected:
    // Create/get internal node, increase its reference count
    Node* getInternalNode(Circuit& circuit, const std::string& name, Node::Flags flags, Status& s=Status::ignore);
    
    Id name_; // Hierarchical name of instance
    Model* model_;
    Instance* parent_;
    Flags flags_;
    const PTInstance& parsedInstance;
};

}

#endif
