#ifndef __CIRCUIT_DEFINED
#define __CIRCUIT_DEFINED

#include <vector>
#include <unordered_map>
#include <memory>
#include "flags.h"
#include "parseroutput.h"
#include "rpneval.h"
#include "ansupport.h"
#include "ansolution.h"
#include "value.h"
#include "options.h"
#include "parameterized.h"
#include "klumatrix.h"
#include "hierdevice.h"
#include "answeep.h"
#include "node.h"
#include "pool.h"
#include "common.h"



namespace NAMESPACE {

// Device/instance/model pointer storage
// - The circuit has a vector of device pointers, it owns these devices
// - Circuit has a map from Id (name) to device index
// - Circuit has a map from Id (name) to model pointer, it owns the models. 
// - Circuit has a map from Id (name) to instance pointer, it owns the instances. 
// - Each device has a vector of model pointers, created by Circuit::buildEntityLists()
// - Each model has a vector of instance pointers, created by Circuit::buildEntityLists()
// - Each hierarchical instance has a vector of child model pointers 
//   listed in the same order as in the parsed subcircuit
// - Each hierarchical instance has a vector of child instance pointers 
//   listed in the same order as in the parsed subcircuit
//
// Each model pointer appears in 3 lists: 
// - under the corresponding device 
// - in the name-pointer map of the circuit
// - in the vector of child models of the parent hierarchical instance, 
//   an exception to this is the toplevel subcircuit definition, which has no parent instance. 
//   A pointer to it is stored by the circuit. 
//
// Each instance pointer appears in 3 lists
// - under the corresponding model
// - in the name-pointer map of the circuit
// - in the vector of child instances of the parent hierarchical instance, 
//   an exception to this is the toplevel subcircuit instance, which has no parent instance. 
//   A pointer to it is stored by the circuit. 

// Nodes are owned by nodePool
// A map from Id to node pointer is in nodeMap. 

// Data dependencies
// variables
//    |------------------|-----------------  
//    v                  v                v
// options      analysis parameters    all instances in the hierarchy
// 
// parameters of the default toplevel subcircuit instance
//   |
//   v
// all instances in the hierarchy
//
// parameters of additional toplevel subcircuit instances
//   |
//   v
// all instances under that particular toplevel subcircuit instance
//
// Contexts at parameter value computation (in order in which they are checked 
// when looking up a parameter). 
// 1) context of the hierarchical instance that contains the 
//    dependent parameter or particular instance/model whose parameter is being computed
// 2) context of the toplevel instance under which the particular instance is located
// 3) context of the default toplevel instance if it is not the same as 2)
// 4) circuit variables
// 5) constants

// At circuit construction
//   1) Remember parser tables
//   2) Set title
//   3) Create builtin devices
//   3) Load OSDI devices
//   4) Add default toplevel subcircuit definition
//   5) Add subcircuit definitions defined within toplevel subcircuit
//   6) Set reference to common saves 
//      (either default in no saves in tables or the saves in tables)
//   7) Set up global parameter name->index map 
//      and set global parameters (global params given as constants)
//   8) Update global context (import global params given as constants and params given as expressions)
// 
// At circuit elaboration
//   1) clear() old circuit
//   2) Set options
//   3) Create default subcircuit definition, create all additional toplevel definitions
//      also creates all nested definitions within them
//   4) Set up ground nodes and global nodes
//   5) Instantiate toplevel subcircuit instances
//      default subcircuit definition first, additional toplevel subcircuits 
//      in the same order as they are given
//   6) Circuit::buildEntityLists()
//        creates a list of models under each device
//        creates a list of instances under each model
//   7) Circuit::nodeOrdering()
//        sort nodes: primary criterion: ground nodes first, 
//                    secondary criterion: alphabetically by name
//        assign initial (temprary) unknown indices to nodes
//   8) Circuit::setup() (force full setup)
//      call Model::setup() and Instance::setup()
//      detect whether node collapsing changed (i.e. list of unknowns changed)
//            e.g. when a series resistance is 0 and causes two nodes to merge
//      clear Model::Flags::SetupNeeded and Instance::Flags::SetupNeeded
//      detect if sparsity changed (i.e. same uknowns but different matrix nonzero entries)
//            e.g. when controlling instance of a current-controlled source changes
//   9) Circuit::mapUnknowns()
//     ---- Enter here when unknowns change
//          At this point all instances, model, and nodes are created, 
//          instances and models are set up, 
//          but indices of unknowns corresponding to nodes are not yet known. 
//     9.1) Device::collapseNodes() ... calls Instance::collapseNodes() for all instances ...
//          Calls Circuit::collapseNodes() with node pairs to do actual collapsing. 
//          Collapsing creates node groups that share the same unknon. 
//     9.2) Assign actual unknown indices to nodes.
//     9.3) Assign representative nodes to unknowns (never an internal device node)
//  10) Circuit::buildSparsityAndStates()
//      ---- Enter here when sparsity pattern changes
//          At this point indices of unknowns corresponding to nodes are known, 
//          but the sparsity map and the states vector are not populated yet. 
//    10.1) Clear SparsityMap 
//          Create Jacobian entries and state vector entries by calling 
//          Device::populateStructures() ... calls Instance::populateStructures() ...
//          calls Circuit::createJacobianEntry() and Circuit::allocateStates()
//    10.2) Create diagonal entries by calling Circuit::createJacobianEntry()
//  11) Circuit::enumerateSystem()
//    11.1) Calls SparsityMap::enumerate() to assign indices in sparse Jacobian vector
//          to Jacobian entries

class Device;
class Model;
class Instance;

class AnnotatedSolution;
class Analysis;

// Compiles supported files, returns canonical path of compilation result
// Return value: ok, compiled file
class SourceCompiler {
public:
    virtual std::tuple<bool, bool> compile(const std::string& timeRefCanonicalPath, const std::string& fileName, const std::string& canonicalPath, std::string& outputCanonicalPath, Status& s=Status::ignore) = 0;
};

enum class CircuitFlags : uint8_t {
    VariablesChanged = 1, 
    HierarchyAffectingOptionsChanged = 2, 
    MappingAffectingOptionsChanged = 4, 
    
    // This one should be set manually whenever an instance/model paramater is changed
    // so that elaborateChanges() will do its job correctly and propagate changes. 
    HierarchyParametersChanged = 8,  

    Elaborated = 16,  

    Abort = 32,  // Exit analysis immediately, even in the middle of computing a point
    Finish = 64, // Wait until current point is computed to the end, then exit simulation
                 // i.e. for multipoint analyses (sweep, frequency sweep, time sweep) 
                 // wait until current point is computed, then exit
                 // Do not exit sweep. 
    Stop = 128,   // Stop analysis to possibly continue it later
                 // Exit sweep. 
    EvalFlags = Abort|Finish|Stop, // All eval flags
};
DEFINE_FLAG_OPERATORS(CircuitFlags);

// TODO: handle states count change after setup() is called
class Circuit : public FlagBase<CircuitFlags> {
public:
    Circuit(ParserTables& tab, SourceCompiler& compiler, Status& s=Status::ignore);
    ~Circuit();

    Circuit           (const Circuit&)  = delete;
    Circuit           (      Circuit&&) = delete;
    Circuit& operator=(const Circuit&)  = delete;
    Circuit& operator=(      Circuit&&) = delete;
    
    bool isValid() const { return valid; };

    // Clear circuit (return to the state that existed immediately after the constructor was called)
    void clear();

    // Devices
    size_t deviceCount() const { return devices.size(); };
    Device* device(size_t ndx) { return devices[ndx].get(); };

    // Variables API
    const Value* getVariable(Id name, Status& s=Status::ignore) const;
    bool setVariable(Id name, const Value& v, Status& s=Status::ignore);
    bool clearVariables(Status& s=Status::ignore);
    
    // Update global context
    bool updateGlobalContext(Status& s=Status::ignore);

    // Sets option, sets HierarchyAffectingOptionsChanged and MappingAffectingOptionsChanged if needed
    // Return value: ok, option changed
    std::tuple<bool, bool> setOption(Id name, const Value& v, Status& s=Status::ignore);

    // Sets options from another options structure 
    // Sets HierarchyAffectingOptionsChanged and MappingAffectingOptionsChanged if needed
    // Return value: options changed
    bool setOptions(IStruct<SimulatorOptions>& opt);

    // Sets options from parsed options
    // Sets HierarchyAffectingOptionsChanged and MappingAffectingOptionsChanged if needed
    // Return value: ok, options changed
    std::tuple<bool, bool> setOptions(const PTParameters& params, Status& s=Status::ignore);

    // Collapse nodes, assign unknowns to nodes, assign representative nodes to unknowns
    bool mapUnknowns(Status& s=Status::ignore);

    // Build sparsity pattern, allocate states vector entries
    bool buildSparsityAndStates(Status& s=Status::ignore);

    // Enumerate Jacobian entries
    bool enumerateSystem(Status& s=Status::ignore);

    // Build circuit from parsed description
    // Prefix is used for prefixing the definition name to obtain the toplevel instance name
    bool elaborate(
        const std::vector<Id>& toplevelDefinitions={}, 
        const std::string& topDefName="__topdef__", const std::string& topInstName="__topinst__", 
        SimulatorOptions* opt=nullptr, 
        Status& s=Status::ignore
    );

    // Title API
    const std::string& title() const { return title_; };
    void setTitle(const std::string& t) { title_ = t; };

    // Node API    
    // Node count
    NodeIndex nodeCount() const;
    // Lookup node by its index (can be used after nodeOrdering() is called)
    Node* node(NodeIndex i) { return nodeOrder[i]; };
    // Lookup node, create it if it does not exist, prohibit creation of a node with same name but different type
    // Increase node reference count
    Node* getNode(Id name, Node::Flags type=Node::Flags::PotentialNode, Status& s=Status::ignore);
    // Decrease node reference count, free it if ref counter reaches 0, return true if node was freed
    bool releaseNode(Node* node, Status& s=Status::ignore);
    // Check if a node is global (use globalNodes set)
    bool isGlobalNode(Id node) const;
    // Representative node corresponding to unknown
    Node* reprNode(UnknownIndex u) const { return unknownToReprNode[u]; };
    // Node value from solution vector
    double nodeValue(Node* node, const double* solution) const { return solution[node->unknownIndex()]; };
    Complex nodeValue(Node* node, const Complex* solution) const { return solution[node->unknownIndex()]; };
    // Add global node (does not create a node, just marks the node name as global)
    bool addGlobal(Id name, Status& s=Status::ignore);
    // Create ground node, there can be multiple ground nodes, all corresponding to the same variable (0)
    // Ground nodes are global
    bool addGround(Id name, Status& s=Status::ignore);

    // Adding entries to fast device, instance, and model access maps
    bool add(Device* dev, Status& s=Status::ignore);
    bool add(Model* mod, Status& s=Status::ignore);
    bool add(Instance* mod, Status& s=Status::ignore);

    // Remove istance from map and delete it
    bool remove(Instance* instance, Status& s=Status::ignore);
    // Remove model from map and delete it
    bool remove(Model* model, Status& s=Status::ignore);
    
    // Fast node, device, model, and instace by name
    Node* findNode(Id name);
    Device* findDevice(Id name, int* index=nullptr);
    Model* findModel(Id name);
    Instance* findInstance(Id name);
    
    std::vector<HierarchicalModel*>& toplevelModels() { return toplevelModels_; };
    std::vector<HierarchicalInstance*>& toplevelInstances() { return toplevelInstances_; };

    // Build lists of instances and models under models and devices
    bool buildEntityLists(Status& s=Status::ignore);

    // Retrieve evaluator, needed by 
    // - analyses to get the global context in which analysis parameters are evaluated
    RpnEvaluator& paramEvaluator() { return paramEvaluator_; };

    // Return variable evaluator
    RpnEvaluator& variableEvaluator() { return variableEvaluator_; };

    // Sparsity map
    SparsityMap& sparsityMap() { return sparsityMap_; };
    const SparsityMap& sparsityMap() const { return sparsityMap_; };
    EquationIndex unknownCount() const { return unknownCountExcludingGround; }; 

    // Helpers 
    // Collapse two nodes, if second node is nullptr, the first node is collapsed to ground
    bool collapseNodes(Node* n1, Node* n2, Status& s=Status::ignore);
    // Allocate Jacobian entry at row and column corresponding to nodes ne and nu
    // Returns pointer to integer where matrix entry index will be found and whether function succeeded
    std::tuple<MatrixEntryIndex*, bool> createJacobianEntry(Node* ne, Node* nu, Status& s=Status::ignore);
    // Allocate n entries in state vector, return global state index of first allocated entry
    StateIndex allocateStates(LocalStateIndex n);

    StateIndex statesCount() const { return statesCount_; };
    
    // Drivers
    // Return value: ok, unknowns changed, sparsity changed
    std::tuple<bool, bool, bool> setup(bool forceFull, Status& s=Status::ignore);
    bool preAnalysis(Status& s=Status::ignore);
    bool nodeOrdering(Status& s=Status::ignore);
    bool bind(
        KluRealMatrix* matResistReal, KluComplexMatrix* matResistCx, Component compResist, 
        KluRealMatrix* matReactReal, KluComplexMatrix* matReactCx, Component compReact, 
        Status& s=Status::ignore
    );
    // Return value: ok, hierarchy changed
    std::tuple<bool, bool> propagateDownHierarchy(Status& s=Status::ignore);
    
    bool evalAndLoad(EvalAndLoadSetup& els, bool (*deviceSelector)(Device*), Status& s=Status::ignore);
    void updateEvalFlags(EvalAndLoadSetup& els, Flags mask=Flags::EvalFlags);
    
    // Simulator options (Parameterized class with simulator options core)
    // IStruct<SimulatorOptions>& simulatorOptions() { return simOptions; };
    const IStruct<SimulatorOptions>& simulatorOptions() const { return simOptions; };

    // Simulator internals
    SimulatorInternals& simulatorInternals() { return simInternals; };
    
    // Dumpers for debugging
    void dumpDevices(int indent, std::ostream& os) const;
    void dumpModels(int indent, std::ostream& os) const;
    void dumpVariables(int indent, std::ostream& os) const;
    void dumpOptions(int indent, std::ostream& os) const;
    void dumpHierarchy(int indent, std::ostream& os) const;
    void dumpNodes(int indent, std::ostream& os) const;
    void dumpUnknowns(int indent, std::ostream& os) const;
    void dumpSparsity(int indent, std::ostream& os) const;
    void dumpSolution(std::ostream& os, const double* solution, const char* prefix="") const;
    void dumpSolution(std::ostream& os, const Complex* solution, const char* prefix="") const;

    // Tolerance computation 
    double solutionTolerance(Node* node, double oldSolution);
    double residualTolerance(Node* node, double maxRes);

    // Storing/restoring a solution
    bool storeDcSolution(Id name, Vector<double>& solution, Status& s=Status::ignore);
    const AnnotatedSolution* retrieveDcSolution(Id name, Status& s=Status::ignore) const;

    // Sets swept values, 
    // Applies common options expressions, analysis parameter expressions, and analysis options expressions
    // Propagates changes down hierarchy, applies changes to the hierarchy (if needed)
    // Rebuilds system (if needed)
    // Adds sparsity map entries and state vector slots requested by analysis (if needed)
    // Return value: ok, hierarchy changed, analysis binding needed
    std::tuple<bool, bool, bool> elaborateChanges(
        ParameterSweeper* sweeper, ParameterSweeper::WriteValues what, 
        Analysis* an, IStruct<SimulatorOptions>* opt, 
        PTParameterMap* optionsMap, 
        Status& s=Status::ignore
    );

private:
    // Remove node from map and delete it regardless of its ref count
    bool remove(Node* node, Status& s=Status::ignore);

    // Check if new options values require us to check if node collapsing changed
    bool mappingAffectingOptionsChanged(SimulatorOptions& opt);

    // Check if new options values require us to check if hierarchy changed
    bool hierarchyAffectingOptionsChanged(SimulatorOptions& opt);

    // Create a subcircuit definitions and all its subdefinitions
    HierarchicalModel* processSubcircuitDefinition(
        const PTSubcircuitDefinition& def, 
        const std::unordered_set<Id>* toplevelDefIds, 
        const std::string& topDefName, 
        const std::string& pathPrefix, 
        int depth, 
        Status& s=Status::ignore
    ); 

    // Build a toplevel subcircuit instance
    bool buildTopInstance(HierarchicalModel* model, Id name, Context& context, Status& s);
    
    bool valid; 

    ParserTables& tables_;
    
    std::string title_;
    std::vector<std::unique_ptr<Device>> devices;
    std::unordered_map<Id,size_t> deviceIndex;

    RecyclingTypedPoolAllocator<Node> nodePool;
    
    std::unordered_map<Id,std::unique_ptr<Instance>> instanceMap;
    std::unordered_map<Id,std::unique_ptr<Model>> modelMap;
    std::unordered_map<Id,Node*> nodeMap;

    // Global parameters
    Value* parameters;
    std::unordered_map<Id, ParameterIndex> parameterIndex_;

    // Hierarchical device (builtin)
    HierarchicalDevice* hdev;
    // Default toplevel model
    HierarchicalModel* defaultToplevelModel_;
    // Toplevel instance fake parser table entries
    std::vector<PTInstance> toplevelInstancesPT_;
    // Toplevel model on which the toplevel instance is based
    std::vector<HierarchicalModel*> toplevelModels_;
    // Toplevel instance
    std::vector<HierarchicalInstance*> toplevelInstances_;
    // Global hierarchy context of toplevel instances
    std::vector<Context> toplevelContext_;
    
    std::unordered_set<Id> globalNodes;
    
    // Node ordering
    std::vector<Node*> nodeOrder;

    // Mapping from unknown index to Node index. 
    // Each unknown can correspond tu multiple nodes. 
    // Unknowns are ordered by their index. 
    // Unknown=0 corresponds to ground node and is used as the RHS bucket. 
    std::unordered_multimap<UnknownIndex,NodeIndex> unknownToNodes;
    
    // Mapping from unknown to representative node
    std::vector<Node*> unknownToReprNode;

    // Ordered map of Jacobian entries
    // Ordering is by column (unknown) first
    SparsityMap sparsityMap_;
    EquationIndex unknownCountExcludingGround;

    // States count
    StateIndex statesCount_;

    // Evaluator for hierarchy
    RpnEvaluator paramEvaluator_;

    // Evaluator for variables
    RpnEvaluator variableEvaluator_;

    // Context holding the variables
    Context variables;
    
    // Simulator options corresponding to the current circuit state
    IStruct<SimulatorOptions> simOptions;

    // Simulator internals
    SimulatorInternals simInternals;
    
    // Annotated solutions (for nodesets and ics)
    std::unordered_map<Id, AnnotatedSolution> dcSolutionRepository;
};


class UnknownNameResolver : public NameResolver {
public:
    UnknownNameResolver(Circuit& circuit) : circuit(circuit) {};

    virtual Id operator()(MatrixEntryIndex u);

private:
    Circuit& circuit;
};


}

#endif
