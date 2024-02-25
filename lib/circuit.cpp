#include <filesystem>
#include <limits>
#include <algorithm>
#include "circuit.h"
#include "osdifile.h"
#include "osdidevice.h"
#include "osdimodel.h"
#include "filesystem.h"
#include "devbase.h"
#include "hierdevice.h"
#include "introspection.h"
#include "simulator.h"
#include "devbuiltins.h"
#include "common.h"


namespace NAMESPACE {

/*
2024/01/09

Suppose we have

    subckt sub1
        subckt subsub1
        ends

        x1 () subsub1
    ends

    subckt sub2
        subckt subsub2
            subckt subsubsub2
            ends

            x3 () subsubsub2
        ends

        x2 () subsub2
        x4 () sub1
    ends

    x5 () sub1

1) If circuit is built from the default toplevel definition we have the following hierarchical models
    __topdef__
    sub1
    sub1:subsub1
    sub2
    sub2:subsub2
    sub2:subsub2:subsubsub2

   There is only one toplevel instance based on __topdef__. 

   Definition lookup at instantiation of x1
     sub1:subsub1 - found
     
   Definition lookup at instantiation of x2
     sub2:subsub2 - found

   Definition lookup at instantiation of x3
     sub2:subsub2:subsubsub2 - found
     
   Definition lookup at instantiation of x4
     sub2:sub1 - not found found
     sub1 - found

   Definition lookup at instantiation of x5
     sub1 - found
   
   The following contexts are available
   1) constants context
      -- persistent, stored in the ContextStack as static member
   2) circuit variables context
      -- persistent, stored in the Circuit
   3) default toplevel instance parameters context (__topdef__)
      (accessible in all instances)
      -- persistent, stored in the Circuit
   4) hierarchical instance parameters context 
      (accessible only in that instance)
      -- not stored, built as required

2) If the circuit is built from sub2 we have the following hierarchical models
    __topdef__         - toplevel definition, does not add prefix
    sub1               - __topdef__ prefix omitted
    sub1:subsub1       - __topdef__ prefix omitted
    sub2               - toplevel definition, does not add prefix
    subsub2            - sub2 prefix omitted
    subsub2:subsubsub2 - sub2 prefix omitted
    
    There are as many toplevel instances as there were subcircuits involved in 
    building the circuit plus the toplevel instance based on __topdef__. 

    Subcircuit definitions defined inside the definition that was used for building 
    the toplevel circuit behave as if they were defined in the toplevel circuit. 

   Definition lookup at instantiation of x1 - takes place in sub1
     sub1:subsub1

   Definition lookup at instantiation of x2 - takes place in sub2
     subsub2 - found

   Definition lookup at instantiation of x3 - takes place in subsub2
     subsub2:subsubsub2 - found

   Definition lookup at instantiation of x4 - takes place in sub2
     sub1 - found
   
   Definition lookup at instantiation of x5 - takes place in __topdef__
     sub1 - found

   The following contexts are available
   1) constants context
      -- persistent, stored in the ContextStack as static member
   2) circuit variables context
      -- persistent, stored in the Circuit
   2) default toplevel instance parameters context (__topdef__)
      (accessible in all subcircuit instances)
      -- persistent, stored in the Circuit
   3) additional toplevel instance parameters context (sub2)
      (accessible only in the corresponding toplevel instance and its children)
      -- persistent, stored in the Circuit
   4) hierarchical instance parameters context 
      (accessible only in that instance)
      -- not stored, built as required

Toplevel instances have an order in which they are created. 
Suppose the order is
  1) __topdef__ (always first)
  2) sub1
  3) sub2

Changes in variables require propagation in all toplevel instances. 
Beside that toplevel instances must be checked for hierarchy change. 

Changes in __topdef__ params require propagation in all toplevel instances. 

Changes in sub2 require propagation only in sub2. 

Hierarchy change in __topdef__ requires rebuild of
of all subsequent toplevel circuits (sub1, sub2) includiong __topdef__
because toplevel models may have changed. 

Hierarchy change in a subinstance of __topdef__ does not require
a rebuild of other toplevel circuits. 

Hierarchy change in sub1 requires rebuild of sub1 and sub2.  
Hierarchy change in a subinstance of sub1 does not require
a rebuild of other toplevel circuits. 

Hierarchy change in sub2 requires rebuild of sub2. 

For all instances, but toplevel ones
1) propagate parameters on a particular hierarchical instance. 
2) check hierarchy change
3) delete hierarchy
4) rebuild hierarchy

For toplevel instances, process them starting with the default one
and proceeding in the order as additional toplevel instances were specified. 
1) Propagate circuit variables to toplevel instance's context
2) Check topology change. 
3) If topology changed, stop variable propagation, mark this instance x. 
   
4) Delete hierarchy of x and all subsequent toplevel instances in reverse order. 

5) Propagate parameters and optionally rebuild parts of subhierarchy for all toplevel instances 
   whose hierarchy has not been deleted starting with the default one and ending with the one before x. 

6) Rebuild hierarchy of toplevel instances whose hierarchy has been deleted, starting with x
   and ending with the last one specified at build time. 

API:
  Instance::enterContext()        - enters evaluation context of the instance
                                    must be called before any topology change or parameter propagation is started
  Instance::revertContext()       - leaves evaluation context of the instance
                                    must be called after all topology change or parameter propagation is finished
  Instance::topologyChanged()     - based on current state of parameters and stored 
                                    topology state report if topology changed
  Instance::propagateParameters() - propagate parameters to immediate descendants
  Instance::deleteHierarchy()     - delete subhierarchy
  Instance::buildHierarchy()      - build subhierarchy, implicitly also propagates parameters
*/

Circuit::Circuit(ParserTables& tab, SourceCompiler& compiler, Status& s) 
    : valid(false), tables_(tab), title_("Untitled"), 
      unknownCountExcludingGround(0), hdev(nullptr), 
      parameters(nullptr) {
    // Construct stuff that is common

    // Prepare variable evaluator
    variableEvaluator_.contextStack().enter(&variables);

    // Set title
    title_ = tables_.title();

    // Create hierarchical device (index 0)
    {
        hdev = new HierarchicalDevice("__hierarchical__", s);
        if (!hdev->checkFlags(Device::Flags::IsValid)) {
            // Destructor will delete hdev
            return;
        }
        // Add to devices 
        auto ndx = devices.size();
        devices.push_back(std::unique_ptr<Device>(hdev));
        deviceIndex[hdev->name()] = ndx;
    }

    // Create builtin devices
    std::vector<Device*> builtinDevices;
    createBuiltins(builtinDevices);
    for(auto& dev : builtinDevices) {
        if (!add(dev, s)) {
            s.extend("Failed to add builtin device.");
            for(auto d : builtinDevices) {
                delete d;
            }
            return;
        }
        dev = nullptr;
    }
    
    // Load OSDI devices
    // If device file is specified with absolute path
    //   look only at the specified absolute path
    // else
    //   if the source of the load directive is a file, not a string
    //     look in the directory where the file with the load directive is located
    //   look in the current working directory
    //   look in the module search path
    // if file is found and is not a compiled file
    //   set compiled file path to compiler output
    //   if file needs compiling
    //     compile it
    // else
    //   set compiled file path to found file
    // load compiled file
    for(auto it=tables_.loads().cbegin(); it!=tables_.loads().cend(); ++it) {
        auto fileName = it->file();
        auto module = it->module();
        auto asModule = it->asModule();
        auto extension = std::filesystem::path(fileName).extension();

        // Get the canonical path of the file where the load directive is located
        auto [fs, pos, line, offset] = it->location().data();
        auto loadDirectiveCanonicalPath = fs->canonicalName(pos);
            
        // Canonical path to osdi file
        std::string canonicalPath;
        bool found;
        
        if (std::filesystem::path(fileName).is_absolute()) {
            // Absolute path given
            found = findFile(fileName, canonicalPath);
        } else {
            // Relative path given
            
            // Try the directory of the file where the load directive was found
            if (loadDirectiveCanonicalPath.size()>0) {
                // Have canonical name of the file (input came from a file)
                auto directory = std::filesystem::path(loadDirectiveCanonicalPath).parent_path().string();
                found = findFile(fileName, canonicalPath, directory);
            }

            // Try current directory
            if (!found) {
                found = findFile(fileName, canonicalPath, "");
            }

            // Try osdi path if no success
            if (!found) {
                found = findFile(fileName, canonicalPath, Simulator::modulePath());
            }
        }
                
        // Not found
        if (!found) {
            s.set(Status::NotFound, std::string("File '")+fileName+"' not found.");
            s.extend(it->location());
            return;
        }

        // Found file, handle compilation of supported files
        // Compiler decides where to put the compiled file and returns its path
        std::string outputPath;
        auto [ok, compiled] = compiler.compile(loadDirectiveCanonicalPath, fileName, canonicalPath, outputPath, s);
        if (!ok) {
            s.extend(it->location());
            return;
        }
        if (compiled) {
            canonicalPath = outputPath;
        }

        // Open possibly compiled file
        auto* osf = OsdiFile::open(canonicalPath, it->location(), s);
        if (!osf) {
            s.set(Status::NotFound, std::string("Failed to open OSDI file '")+fileName+"'.");
            s.extend(it->location());
            return;
        }

        if (!module) {
            // Get all devices
            auto devCount = osf->deviceCount();
            for(decltype(devCount) i=0; i<devCount; i++) {
                auto* dev = osf->createDevice(i, Id::none, it->location(), s);
                if (!add(dev, s)) {
                    s.extend("Failed to add device.");
                    s.extend(it->location());
                    delete dev;
                    return;
                }
            }
        } else {
            // Get specific device
            auto* dev = osf->createDevice(module, asModule, it->location(), s);
            if (!dev) {
                s.set(Status::NotFound, std::string("Module '")+std::string(module)+"' not found in OSDI file '"+std::string(fileName)+"'.");
                s.extend(it->location());
                return;
            }
            if (!add(dev, s)) {
                s.extend("Failed to add device.");
                s.extend(it->location());
                delete dev;
                return;
            }
        }
    }

    // Create add variables context to paramEvaluator_ context stack
    if (!updateGlobalContext(s)) {
        return;
    }

    valid = true;
}

Circuit::~Circuit() {
    // Be careful about the order
    // Instances need models, which in turn need devices
    // Just clear structures that own objects
    instanceMap.clear();
    modelMap.clear();
    devices.clear();
    nodePool.clear();

    delete [] parameters;
}

void Circuit::clear() {
    // Be careful about the order
    // Instances need models, which in turn need devices
    instanceMap.clear();

    // Remove all model map entries
    modelMap.clear();

    // Remove all model map entries, except for the subcircuit definition entries
    // for(auto it=modelMap.begin(); it!=modelMap.end();) {
    //     auto next = it;
    //     ++next;
    //     if (it->second.get()->device() != hdev) {
    //         modelMap.erase(it);
    //     }
    //     it = next;
    // }

    // Clear lists of models under devices
    for(auto& it : devices) {
        it.get()->clearModelList();
    }

    // No need to clear entity lists under models
    // as they have been deleted with them. 
    // An exception to this are subcircuit models because they were not deleted. 
    // for(auto& it : modelMap) {
    //     it.second.get()->clearInstanceList();
    // }

    // Delete nodes
    nodePool.clear();
    nodeMap.clear();

    // Global nodes
    globalNodes.clear();

    // Node order
    nodeOrder.clear();

    // Unknown to nodes mapping
    unknownToNodes.clear();

    // Unknown to representative node
    unknownToReprNode.clear();

    // Sparsity map
    sparsityMap_.clear();

    // Misc
    toplevelInstancesPT_.clear();
    toplevelModels_.clear();
    toplevelInstances_.clear(); 
    toplevelContext_.clear();
    unknownCountExcludingGround = 0;
}

NodeIndex Circuit::nodeCount() const {
    return nodeMap.size();
}

Node* Circuit::getNode(Id name, Node::Flags type, Status& s) {
    // See if it is already there
    auto existingNode = findNode(name);
    if (existingNode) {
        if ((existingNode->maskedFlags(Node::Flags::NodeTypeMask)) == type) {
            existingNode->incRef();
            return existingNode;
        } else {
            s.set(Status::BadNode, "A node named '"+std::string(name)+"' exists but is of wrong type.");
            return nullptr;
        }
    }

    // Not there, create it
    auto node = nodePool.allocate(name, type);
    node->incRef();
    nodeMap.insert({name, node});
    return node;
}

bool Circuit::releaseNode(Node* node, Status& s) {
    if (!findNode(node->name())) {
        s.set(Status::Internal, "Node '"+std::string(node->name())+"' is not owned by this circuit.");
        return false;
    }
    if (node->decRef()==0) {
        // Release node
        return remove(node, s);
    }
    return true;
}

bool Circuit::isGlobalNode(Id node) const {
    return globalNodes.find(node)!=globalNodes.end();
}

bool Circuit::addGlobal(Id name, Status& s) {
    globalNodes.insert(name);
    return true;
}

bool Circuit::addGround(Id name, Status& s) {
    auto node = getNode(name, Node::Flags::PotentialNode, s);
    if (node) {
        // Unknown index is 0 for ground nodes
        node->setUnknownIndex(0);
        node->setFlags(Node::Flags::Ground);

        // Ground nodes are global
        if (!addGlobal(name, s)) {
            return false;
        }
    }
    return node != nullptr;
}

bool Circuit::add(Device* dev, Status& s) {
    auto [it, inserted] = deviceIndex.insert({dev->name(), devices.size()});
    
    // Insertion failed because there is a device present with the same name
    if (!inserted) {
        // Is it also the same device with different name
        auto other = devices[it->second].get();
        if (*dev==*other) {
            // Same type, same device, ignore (because we already have it)
            return true;
        }
        s.set(Status::Redefinition, std::string("A device with name '")+std::string(dev->name())+"' already exists.");
        if (other->location()) {
            s.extend("The existing device was first defined here");
            s.extend(other->location());
        }
        return false;
    }
    devices.push_back(std::unique_ptr<Device>(dev));
    deviceIndex[dev->name()] = devices.size()-1;
    return true;
}

bool Circuit::add(Model* mod, Status& s) {
    // Failure to insert will delete model
    auto [it, inserted] = modelMap.insert({mod->name(), nullptr});
    if (!inserted) {
        s.set(Status::Redefinition, "A model with name '"+std::string(mod->name())+"' already exists.");
        s.extend(mod->location());
        if (it->second->location()) {
            s.extend("The existing model was first defined here");
            s.extend(it->second->location());
        }
        return false;
    }
    it->second.reset(mod);
    return true;
}

bool Circuit::add(Instance* inst, Status& s) {
    // Failure to insert will delete instance
    auto [it, inserted] = instanceMap.insert({inst->name(), nullptr});
    if (!inserted) {
        s.set(Status::Redefinition, "An instance with name '"+std::string(inst->name())+"' already exists.");
        s.extend(inst->location());
        if (it->second->location()) {
            s.extend("The existing instance was first defined here");
            s.extend(it->second->location());
        }
        return false;
    }
    it->second.reset(inst);
    return true;
}

bool Circuit::remove(Node* node, Status& s) {
    auto name = node->name();
    auto checkNode = findNode(name);
    if (!checkNode || checkNode!=node) {
        s.set(Status::Internal, "Node '"+std::string(name)+"' is not owned by this circuit.");
        return false;
    }
    nodeMap.erase(name);
    nodePool.free(node);
    return true;
}

bool Circuit::remove(Instance* instance, Status& s) {
    auto name = instance->name();
    auto checkInstance = findInstance(name);
    if (!checkInstance || checkInstance!=instance) {
        s.set(Status::Internal, "Instance '"+std::string(name)+"' is not owned by this circuit.");
        return false;
    }
    instanceMap.erase(name);
    return true;
}

bool Circuit::remove(Model* model, Status& s) {
    auto name = model->name();
    auto checkModel = findModel(name);
    if (!checkModel || checkModel!=model) {
        s.set(Status::Internal, "Model '"+std::string(name)+"' is not owned by this circuit.");
        return false;
    }
    modelMap.erase(name);
    return true;
}

Node* Circuit::findNode(Id name) {
    auto it=nodeMap.find(name);
    if (it==nodeMap.end()) {
        return nullptr;
    }
    return it->second;
}

Device* Circuit::findDevice(Id name, int* index) {
    auto it = deviceIndex.find(name);
    if (it!=deviceIndex.end()) {
        if (index) {
            *index = it->second;
        }
        return devices[it->second].get();
    }
    return nullptr;
}

Model* Circuit::findModel(Id name) {
    auto it = modelMap.find(name);
    if (it!=modelMap.end()) {
        return it->second.get();
    }
    return nullptr;
}

Instance* Circuit::findInstance(Id name) {
    auto it = instanceMap.find(name);
    if (it!=instanceMap.end()) {
        return it->second.get();
    }
    return nullptr;
}

bool Circuit::buildEntityLists(Status& s) {
    for(auto& it : devices) {
        it.get()->clearModelList();
    }
    for(auto& it : modelMap) {
        auto modelPtr = it.second.get();
        modelPtr->clearInstanceList();
        modelPtr->device()->addModel(modelPtr);
    }
    for(auto& it : instanceMap) {
        auto instPtr = it.second.get();
        instPtr->model()->addInstance(instPtr);
    }
    return true;
}

HierarchicalModel* Circuit::processSubcircuitDefinition(
    const PTSubcircuitDefinition& def, 
    const std::unordered_set<Id>* toplevelDefIds, 
    const std::string& topDefName, 
    const std::string& pathPrefix, 
    int depth, 
    Status& s
) {
    auto name = def.name();
    
    // Do we treat this definition as a toplevel definition? 
    bool asToplevel;
    std::string childrenPathPrefix;
    if (depth==0) {
        // At depth 0 (default toplevel definition) we treat it as toplevel definition
        asToplevel = true;
        name = topDefName;
        // Definition names of children are not prefixed
        childrenPathPrefix = "";
    } else if (depth==1 && toplevelDefIds->contains(name)) {
        // At depth 1 we check if its name is in toplevelDefNames
        asToplevel = true;
        name = topDefName+"("+std::string(name)+")"; 
        // Definition names of children are not prefixed
        childrenPathPrefix = "";
    } else {
        // All others are treated as ordinary definitions and are prefixed by path
        asToplevel = false;
        // Add prefix to definition names of children
        childrenPathPrefix = std::string(name);
    }
    
    // Hierarchical definition name - prefix it with specified path
    std::string hierDefName;
    if (pathPrefix.size()>0) {
        // Have prefix
        hierDefName = pathPrefix + ":" + std::string(name);
    } else {
        // No prefix
        hierDefName = std::string(name);
    }

    // Create definition
    auto hierarchicalModel = new HierarchicalModel(hdev, hierDefName, nullptr, def, s);
    // std::cout << hierDefName << "\n";
    if (!hierarchicalModel->checkFlags(Model::Flags::IsValid)) {
        // Need to delete it manually because it is not stored in modelMap
        delete hierarchicalModel;
        return nullptr;
    }

    // Store it in modelMap
    if (!add(hierarchicalModel, s)) {
        return nullptr;
    }

    // Now process all definitions that were defined inside this definition
    for(auto& it : def.subDefs()) {
        if (!processSubcircuitDefinition(*(it.get()), toplevelDefIds, topDefName, childrenPathPrefix, depth+1, s)) {
            return nullptr;
        }
    }

    return hierarchicalModel;
}

bool Circuit::buildTopInstance(HierarchicalModel* model, Id name, Context& context, Status& s) {
    // Fake parser table entry
    toplevelInstancesPT_.push_back(
        std::move(
            PTInstance(
                Loc::bad, name, model->name(), 
                PTIdentifierList({}), PTParameters()
            )
        )
    );
    InstantiationData idata;
    auto inst = model->createInstance(*this, nullptr, paramEvaluator_, &context, toplevelInstancesPT_.back(), idata, s);
    if (!inst) {
        return false;
    }
    toplevelInstances_.push_back(static_cast<HierarchicalInstance*>(inst));
    return true;
}

bool Circuit::elaborate(
    const std::vector<Id>& toplevelDefinitions, 
    const std::string& topDefName, const std::string& topInstName, 
    SimulatorOptions* opt, 
    Status& s
) { 
    clear();
    clearFlags(Flags::Elaborated);

    // Set options
    if (opt) {
        // Options given
        simOptions.core() = *opt;
    } else {
        // No option given, reset to defaults
        simOptions.core() = SimulatorOptions();
    }

    // Convert vector to set
    std::unordered_set<Id> defIdSet(toplevelDefinitions.begin(), toplevelDefinitions.end());

    // This also creates definitions defined inside toplevel definition. 
    // Definitions that were specified in toplevelDefinitions are processed differently,  
    // that is definitions defined within them are not prefixed 
    // (i.e. behave as if they were defined in the toplevel definition). 
    defaultToplevelModel_ = processSubcircuitDefinition(tables_.defaultSubDef(), &defIdSet, topDefName, "", 0, s);
    if (!defaultToplevelModel_) {
        return false;
    }
    
    // // Propagate global parameters to options, skip propagation to analysis because we have none
    // IStruct<SimulatorOptions> opt;
    // if (auto [ok, changed] = opt.setParameters(tables_.options(), paramEvaluator_, s); !ok) {
    //     s.extend("Initial options computation failed.");
    //     return false;
    // }
    // 
    // // Set options
    // auto changed = setOptions(opt);

    // Now we are ready to build the circuit
    
    // Set up ground nodes
    for(auto it=tables_.groundNodes().cbegin(); it!=tables_.groundNodes().cend(); ++it) {
        if (!addGround(it->name(), s)) {
            s.set(Status::CreationFailed, "Failed to add ground node '"+std::string(it->name())+"'.");
            return false;
        }
    }

    // Set up global nodes
    for(auto& it : tables_.globalNodes()) {
        if (!addGlobal(it.name(), s)) { 
            s.set(Status::CreationFailed, "Failed to add global node '"+std::string(it.name())+"'.");
            return false;
        }
    }

    // Default toplevel circuit is always built
    toplevelModels_.push_back(defaultToplevelModel_);
    // First look up all definitions
    for(auto& defName : toplevelDefinitions) {
        auto actualName = topDefName + "(" + std::string(defName) +")";
        auto model = findModel(actualName);
        if (!model) {
            s.set(Status::BadArguments, "Subcircuit definition '"+std::string(defName)+"' not found.\nOnly subcircuit definitions from the default toplevel circuit can be used.");
            return false;
        }
        if (!model->device()->isHierarchical()) {
            s.set(Status::BadArguments, "'"+std::string(defName)+"' is not a subcircuit definition.");
            return false;
        }
        HierarchicalModel* hMod = static_cast<HierarchicalModel*>(model);
        // Add to toplevel models
        toplevelModels_.push_back(hMod);
    }
    // Toplevel context
    auto n = toplevelModels_.size();
    toplevelContext_.resize(n);
    for(auto it : toplevelContext_) {
        it.clear();
    }
    // Get stack marker to which we return in case of failure
    size_t contextMarker = paramEvaluator_.contextMarker();
    // Next, build all toplevel instances
    for(decltype(n) i=0; i<n; i++) {
        auto model = toplevelModels_[i];
        Id instanceName;
        if (i==0) {
            // Default subcircuit definition has instance name equal to topInstName
            instanceName = topInstName;
        } else {
            // Other subcircuit definitions have instance name topInstName(name)
            instanceName = topInstName+"("+std::string(toplevelDefinitions[i-1])+")";
        }
        HierarchicalModel* hMod = static_cast<HierarchicalModel*>(model);
        if (!buildTopInstance(hMod, instanceName, toplevelContext_[i], s)) {
            // Revert to context stack state before building
            Instance::revertContext(*this, contextMarker);
            return false;
        }
        // For second and all further toplevel instances the context
        // of default toplevel circuit is global. So insert it into 
        // the context stack and add it to the search path, but do not rebuild it. 
        // Do this only once for all toplevel instances when the first 
        // non-default toplevel instance is built. 
        if (i==0) {
            toplevelInstances_[0]->enterContext(*this, &(toplevelContext_[0]), true, false, s);
        }   
    }
    // Revert to context stack state before building
    Instance::revertContext(*this, contextMarker);
    
    // At this point all models and instances are created and all parameters are set. 

    // Build entity lists
    if (!buildEntityLists(s)) {
        s.extend("Initial entity list building failed.");
        return false;
    }

    // Order nodes
    if (!nodeOrdering(s)) {
        s.extend("Initial node ordering failed.");
        return false;
    }

    // We do a full setup
    auto [ok, unknownsChanged, sparsityChanged] = setup(true, s);
    if (!ok) {
        s.extend("Initial circuit setup failed.");
        return false;
    }

    // Assign unknowns to nodes
    if (!mapUnknowns(s)) {
        s.extend("Initial unknown mapping failed.");
        return false;
    }

    // Build sparsity pattern and states vector entries
    if (!buildSparsityAndStates(s)) {
        s.extend("Initial sparsity pattern creation and states allocation failed.");
        return false;
    }


    // We enumerate Jacobian entries
    if (!enumerateSystem(s)) {
        s.extend("System enumeration failed.");
        return false;
    }

    // Sanity checks

    // Check if number of nodes is nonzero
    if (nodeCount()<=0) {
        s.set(Status::Empty, "Circuit has no nodes.");
        return false;
    }

    // Check if number of unknowns is nonzero
    if (unknownCount()<=0) {
        s.set(Status::Empty, "Circuit has no unknowns.");
        return false;
    }

    // Check if there is at least one nonzero entry in sparsity map
    if (sparsityMap_.size()<=0) {
        s.set(Status::Empty, "Sparsity pattern has no nonzero entries.");
        return false;
    }
    
    // Circuit is in consistent state
    // Mark that by clearing flags
    clearFlags(Flags::VariablesChanged);
    clearFlags(Flags::HierarchyAffectingOptionsChanged);
    clearFlags(Flags::MappingAffectingOptionsChanged);
    clearFlags(Circuit::Flags::HierarchyParametersChanged);

    // Mark circuit as elaborated
    setFlags(Flags::Elaborated);

    return true;
}

std::tuple<bool, bool, bool> Circuit::setup(bool forceFull, Status& s) {
    // Do setup
    bool unknownsChanged = false;
    bool sparsityChanged = false;
    for(auto& dev : devices) {
        auto [ok, tmpUnknowns, tmpSparsity] = dev->setup(*this, forceFull, s);
        unknownsChanged |= tmpUnknowns;
        sparsityChanged |= tmpSparsity;
        if (!ok) {
            return std::make_tuple(false, unknownsChanged, sparsityChanged);
        }
    }
    return std::tuple(true, unknownsChanged, sparsityChanged);
}

bool Circuit::preAnalysis(Status& s) {
    // Do pre-analysis computations
    for(auto& dev : devices) {
        if (!dev->preAnalysis(*this, s)) {
            return false;
        }
    }
    return true;
}

// Driver function for collapsing two nodes, called by instances
bool Circuit::collapseNodes(Node* n1, Node* n2, Status& s) {
    // Collapsing node to itself, nothing to do
    if (n1==n2) {
        return true;
    }
    
    // Get unknowns
    auto u1 = n1->unknownIndex();
    // Assume node2 is ground
    decltype(u1) u2 = 0; 
    if (n2) {
        u2 = n2->unknownIndex();
    }

    // If the two unknowns are the same we have nothing to do
    if (u1==u2) {
        return true;
    }

    // Reorder so that u1<u2
    if (u2<u1) {
        auto tmp = u1;
        u1 = u2;
        u2 = tmp;
    }

    // std::cout << n1 << "<-" << n2 << " : " << u1 << "<-" << u2 << std::endl;

    // Relocate nodes corresponding to u2 to nodes corresponding to u1
    auto range = unknownToNodes.equal_range(u2);
    auto it = range.first;
    do {
        auto at = it;
        ++it;
        // Change node unknown to u1
        nodeOrder[at->second]->setUnknownIndex(u1);
        // Extract multimap node
        auto h = unknownToNodes.extract(at);
        // Change key (unknown) to u1
        h.key() = u1;
        // Reinsert
        auto itIns = unknownToNodes.insert(std::move(h));
    } while (it!=range.second);

    /*
    int lastu = -1;
    bool first = true;
    for(auto it=unknownToNodes.begin(); it!=unknownToNodes.end(); ++it) {
        if (first || it->first!=lastu) {
            lastu = it->first;
            std::cout << std::endl << "  u" << it->first << " : ";
            first = false;
        }
        std::cout << it->second << " ";
    }
    std::cout << std::endl;
    */
    return true;
}

std::tuple<MatrixEntryIndex*, bool> Circuit::createJacobianEntry(Node* ne, Node* nu, Status& s) {
    // Map nodes to equation, unknown pair
    auto e = ne->unknownIndex();
    auto u = nu->unknownIndex();

    // if e or u correspond to ground, create no entry, 
    // but indicate everything is OK
    if (!e || !u) {
        return std::make_tuple(nullptr, true);
    }

    // Create entry
    auto [ptr, ok] = sparsityMap_.insert(e, u);
    if (!ok) {
        s.set(Status::MatrixCreate, "Failed to create matrix entry.");
    }

    return std::make_tuple(ptr, ok);
}

StateIndex Circuit::allocateStates(LocalStateIndex n) {
    auto retval = statesCount_;
    statesCount_ += n;
    return retval;
}

bool Circuit::nodeOrdering(Status& s) {
    // Build node order vector
    nodeOrder.clear();
    for(auto& it : nodeMap) {
        nodeOrder.push_back(it.second);
    }
    
    // Sort nodes
    // Primary criterion: ground nodes before non-ground nodes
    // Secondary criterion: node name (lexicographically)
    struct {
        bool operator()(Node* a, Node* b) { 
            struct CstrLess nameComparison;
            // a ground, b not ground
            if (a->checkFlags(Node::Flags::Ground) && !b->checkFlags(Node::Flags::Ground)) {
                return true;
            }
            // a not ground, b ground
            if (!a->checkFlags(Node::Flags::Ground) && b->checkFlags(Node::Flags::Ground)) {
                return false;
            }
            // Both either ground or not ground, compare names
            return nameComparison(a->name().c_str(), b->name().c_str()); 
        };
    } nodeComparison;
    std::sort(nodeOrder.begin(), nodeOrder.end(), nodeComparison);

    // Do we have at least one node
    if (nodeOrder.size()<=0) {
        s.set(Status::BadNode, "No nodes defined.");
        return false;
    }
    
    // Do we have at least one ground node
    if (!nodeOrder[0]->checkFlags(Node::Flags::Ground)) {
        s.set(Status::NoGround, "No ground node defined.");
        return false;
    }

    // Enumerate unknowns, ground nodes are 0, the rest are enumerated starting with 1
    UnknownIndex nodeNumber = 1;
    for(auto it=nodeOrder.begin(); it!=nodeOrder.end(); ++it) {
        if ((*it)->checkFlags(Node::Flags::Ground)) {
            (*it)->setUnknownIndex(0);
        } else {
            (*it)->setUnknownIndex(nodeNumber);
            nodeNumber++;
        }
    }

    return true;
}

bool Circuit::mapUnknowns(Status& s) {
    // Reset internal structures
    auto nNodes = nodeOrder.size(); 
    unknownToNodes.clear();

    // Check number of nodes
    if (nodeMap.size()>std::numeric_limits<NodeIndex>::max()) {
        throw std::length_error("Too many nodes.");
    }

    // Prepare unknownToNodes unordered_multimap
    for(NodeIndex i=0; i<nNodes; i++) {
        unknownToNodes.insert({nodeOrder[i]->unknownIndex(), i});
    }
    
    // Do node collapsing (uses unknownToNodes unordered_multimap)
    for(auto& dev : devices) {
        if (!dev.get()->collapseNodes(*this, s)) {
            return false;
        }
        
    }

    // Create a vector of booleans indicating that an unknown number is present. 
    // There are at most as many unknowns as there are nodes. 
    // Assume by default a number is not present. 
    std::vector<bool> unknownToRenumber(nNodes, false);

    // Fill the vector
    for(auto it : nodeMap) {
        unknownToRenumber[it.second->unknownIndex()] = true;
    }

    // Prepare a vector holding new unknown numbers coresponding to old ones
    std::vector<UnknownIndex> newUnknownIndex(nNodes, 0);
    UnknownIndex atUnknown = 0;
    for(UnknownIndex i=0; i<nNodes; i++) {
        if (unknownToRenumber[i]) {
            newUnknownIndex[i] = atUnknown;
            atUnknown++;
        }
    }

    // Assign new unknown indices to nodes
    for(auto it : nodeMap) {
        auto oldUnknownIndex = it.second->unknownIndex();
        it.second->setUnknownIndex(newUnknownIndex[oldUnknownIndex]);
    }

    // We don't need the map anymore
    unknownToNodes.clear();

    // Set number of unknowns (do not cound ground)
    unknownCountExcludingGround = atUnknown-1;

    // Create mapping from unknown to representative node
    // Go through all nodes, prefer nodes that are not internal
    // Prefer nodes with lower node index (scan odering in reverse)
    unknownToReprNode.resize(atUnknown, nullptr);
    for(auto it=nodeOrder.rbegin(); it!=nodeOrder.rend(); ++it) {
        auto node = *it;
        auto unknownIndex = node->unknownIndex();
        auto existingNode = unknownToReprNode[unknownIndex];
        if (!existingNode) {
            // No node assigned yet
            unknownToReprNode[unknownIndex] = node;
        } else if (
            existingNode->checkFlags(Node::Flags::InternalDeviceNode)==node->checkFlags(Node::Flags::InternalDeviceNode)
        ) {
            // Both nodes are either internal or not internal
            // Replace the existing node because we prefer lower node index
            unknownToReprNode[unknownIndex] = node;
        } else if (
            existingNode->checkFlags(Node::Flags::InternalDeviceNode) &&
            !node->checkFlags(Node::Flags::InternalDeviceNode)
        )  {
            // Existing node is internal, new node is not internal
            // Replace existing node because we prefer non-internal nodes
            unknownToReprNode[unknownIndex] = node;
        }
    }

    return true;
}

bool Circuit::buildSparsityAndStates(Status& s) {
    // Clear sparsity map
    sparsityMap_.clear();

    // Create Jacobian entries, allocate state vector slots
    statesCount_ = 0;
    for(auto& dev : devices) {
        if (!dev.get()->populateStructures(*this, s)) {
            return false;
        }
    }

    // Add diagonal entries
    for(decltype(unknownCountExcludingGround) i=1; i<=unknownCountExcludingGround; i++) {
        Node* node = unknownToReprNode[i];
        auto [ptr, ok] = createJacobianEntry(node, node, s);
        if (!ok) {
            return false;
        }
    }

    return true;
}

bool Circuit::enumerateSystem(Status& s) {
    // Enumerate Jacobian entries
    sparsityMap_.enumerate();
    
    // Sparsity pattern creation and Jacobian/states binding is done in analysis

    return true;
}

bool Circuit::bind(
    KluRealMatrix* matResistReal, KluComplexMatrix* matResistCx, Component compResist, 
    KluRealMatrix* matReactReal, KluComplexMatrix* matReactCx, Component compReact, 
    Status& s
) {
    // Sanity checks
    if (matResistReal && matResistCx) {
        s.set(Status::Internal, "Cannot bind resistive part twice.");
        return false;
    }
    if (matResistReal && compResist==Component::ImagPart) {
        s.set(Status::Internal, "Cannot bind resistive part to imaginary part of real matrix.");
        return false;
    }
    if (matReactReal && matReactCx) {
        s.set(Status::Internal, "Cannot bind reactive part twice.");
        return false;
    }
    if (matReactReal && compReact==Component::ImagPart) {
        s.set(Status::Internal, "Cannot bind reactive part to imaginary part of real matrix.");
        return false;
    }
    // Call bind() for all devices
    for(auto& dev : devices) {
        if (!dev.get()->bind(
            *this, 
            matResistReal, matResistCx, compResist,  
            matReactReal, matReactCx, compReact, 
            s
        )) {
            return false;
        }
    }
    return true;
}

bool Circuit::evalAndLoad(EvalAndLoadSetup& els, bool (*deviceSelector)(Device*), Status& s) {
    if (!els.initialize(s)) {
        return false;
    }
    els.clearFlags();
    els.clearBounds();
    for(auto& dev : devices) {
        if (!deviceSelector || deviceSelector(dev.get())) {
            if (!dev.get()->evalAndLoad(*this, els, s)) {
                return false;
            }
        }
    }
    return true;
}

void Circuit::updateEvalFlags(EvalAndLoadSetup& els, Flags mask) {
    if (els.abortRequested) {
        setFlags(Flags::Abort & mask);
    }
    if (els.finishRequested) {
        setFlags(Flags::Finish & mask);
    }
    if (els.stopRequested) {
        setFlags(Flags::Stop & mask);
    }
}

double Circuit::solutionTolerance(Node* node, double oldSolution) {
    auto& options = simOptions.core();
    auto i = node->unknownIndex();
    // Solution tolerance
    double tol = std::fabs(oldSolution)*options.reltol;
    // Absolute solution tolerance differs for potential and flow nodes
    if (node->maskedFlags(Node::Flags::NodeTypeMask)==Node::Flags::PotentialNode) {
        // Potential
        if (tol<options.vntol) {
            tol = options.vntol;
        }
    } else {
        // Flow
        if (tol<options.abstol) {
            tol = options.abstol;
        }
    }
    return tol;
}

double Circuit::residualTolerance(Node* node, double maxRes) {
    auto& options = simOptions.core();
    auto i = node->unknownIndex();
    // Residual tolerance (Designer's Guide to Spice and Spectre, chapter 2.2.2)
    double tol = std::fabs(maxRes)*options.reltol;
    // Absolute residual tolerance differs for potential and flow nodes
    if (node->maskedFlags(Node::Flags::NodeTypeMask)==Node::Flags::PotentialNode) {
        // Potential node residual is current
        if (tol<options.restol) {
            tol = options.restol;
        }
    } else {
        // Flow node residual is voltage
        tol = options.vnrestol;
        if (tol<options.vnrestol) {
            tol = options.vnrestol;
        }
    }
    return tol;
}

bool Circuit::storeDcSolution(Id name, Vector<double>& solution, Status& s) {
    auto [it, inserted] = dcSolutionRepository.insert({name, AnnotatedSolution()});
    auto& repo = it->second;
    repo.set(*this, solution);

    return true;
}

const AnnotatedSolution* Circuit::retrieveDcSolution(Id name, Status& s) const {
    auto it = dcSolutionRepository.find(name);
    if (it==dcSolutionRepository.cend()) {
        s.set(Status::NotFound, "DC solution with label '"+std::string(name)+"' not found.");
        return nullptr;
    }
    return &(it->second);
}


Id UnknownNameResolver::operator()(MatrixEntryIndex u) {
    return circuit.reprNode(u+1)->name();
}

void Circuit::dumpDevices(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    for(auto& it : devices) {
        os << pfx << it->name() << "\n";
    }
}

void Circuit::dumpModels(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    for(auto& it : modelMap) {
        os << pfx << it.first << " (device=" << it.second.get()->device()->name() << ")" << "\n";
    }
}

void Circuit::dumpVariables(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');

    os << pfx << "Circuit variables\n";
    variables.dump(indent+2, os);
}

void Circuit::dumpOptions(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    auto n = simOptions.parameterCount();
    for(decltype(n) i=0; i<n; i++) {
        auto name = simOptions.parameterName(i);
        Value v;
        auto ok = simOptions.getParameter(name, v);
        if (ok) {
            os << pfx << name << " = " << v << "\n";
        }
    }
}

static void dumpInstance(std::ostream& os, Instance* inst, int indent) {
    os << std::string(indent, ' ') << inst->name() << " (model=" << inst->model()->name() 
       << ", device=" << inst->model()->device()->name()<< ")" << "\n"; 
    if (inst->childModels()) {
        for(auto it=inst->childModels()->cbegin(); it!=inst->childModels()->cend(); ++it) {
            os << std::string(indent+2, ' ') << "model " << (*it)->name() 
               << " (device=" << (*it)->device()->name() << ")" << "\n"; 
        }
    }
    if (inst->childInstances()) {
        for(auto it=inst->childInstances()->cbegin(); it!=inst->childInstances()->cend(); ++it) {
            dumpInstance(os, *it, indent+2);
        }
    }
}

void Circuit::dumpHierarchy(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    for(auto toplevelInstance : toplevelInstances_) {
        dumpInstance(os, toplevelInstance, indent);
    }
}

void Circuit::dumpNodes(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    if (nodeMap.size()>0) {
        auto nnames = nodeOrder.size();
        for(decltype(nnames) i=0; i<nnames; i++) {
            os << pfx << std::to_string(i) << " : " << std::string(nodeOrder[i]->name()) << " ";
            if (nodeOrder[i]->maskedFlags(Node::Flags::NodeTypeMask) == Node::Flags::FlowNode) {
                os << "(flow) ";
            }
            if (nodeOrder[i]->checkFlags(Node::Flags::InternalDeviceNode)) {
                os << "(internal) ";
            }
            os << "\n";
        }
    }
}

void Circuit::dumpUnknowns(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    auto n = unknownCountExcludingGround;
    for(decltype(n) i=0; i<=n; i++) {
        os << pfx << i << " : " << reprNode(i)->name() << "\n";
    }
}

void Circuit::dumpSparsity(int indent, std::ostream& os) const {
    sparsityMap_.dump(indent, os);
}

void Circuit::dumpSolution(std::ostream& os, const double* solution, const char* prefix) const {
    for(NodeIndex i=0; i<nodeOrder.size(); i++) {
        if (i!=0) {
            os << std::endl;
        }
        os << prefix << nodeOrder[i]->name() << " : " << nodeValue(nodeOrder[i], solution);
    }
}

void Circuit::dumpSolution(std::ostream& os, const Complex* solution, const char* prefix) const {
    for(NodeIndex i=0; i<nodeOrder.size(); i++) {
        if (i!=0) {
            os << std::endl;
        }
        auto nv = nodeValue(nodeOrder[i], solution);
        os << prefix << nodeOrder[i]->name() << " : " << nv.real();
        if (nv.imag()>=0) {
            os << "+";
        }
        os << nv.imag();
        os << "i";
    }
}

}
