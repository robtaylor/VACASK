#include "osdimodel.h"
#include "osdiinstance.h"
#include "hierdevice.h"
#include "circuit.h"
#include "common.h"


namespace NAMESPACE {

// We will have only one instance of this class so the identifier is fixed
HierarchicalDevice::HierarchicalDevice(Id name, Status& s) : Device(name) {
    setFlags(Flags::IsValid);
}

HierarchicalDevice::~HierarchicalDevice() {
}

bool HierarchicalDevice::operator==(const Device& other) const {
    const HierarchicalDevice* hdOther = dynamic_cast<const HierarchicalDevice*>(&other);
    if (hdOther && hdOther == this) {
        return true;
    }
    return false;
}

void HierarchicalDevice::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "Hierarchical device: " << std::string(name()) << "\n";
}


HierarchicalModel::HierarchicalModel(HierarchicalDevice* dev, Id name, Instance* parentInstance, const PTSubcircuitDefinition& parsedSubcircuit, Status& s)
    : Model(dev, name, parentInstance, parsedSubcircuit) {
    // Prepare terminal map, check uniqueness
    if (!buildTerminalMap(s)) {
        return;
    }
    
    // Prepare parameter map, check uniqueness
    if (!buildParameterMap(s)) {
        return;
    }

    // Clear Flags::NeedsSetup because hierarchical models do not need to be set up
    clearFlags(Flags::NeedsSetup);

    setFlags(Flags::IsValid);
}

HierarchicalModel::~HierarchicalModel() {
}


bool HierarchicalModel::buildTerminalMap(Status& s) {
    // Build map, check for duplicates
    TerminalIndex i=0;
    terminalMap.clear();
    auto& parsedSubcircuit = static_cast<const PTSubcircuitDefinition&>(parsedModel);
    auto& terminals = parsedSubcircuit.terminals();
    TerminalIndex n = terminals.size();
    for(decltype(n) i=0; i<n; i++) {
        auto id = terminals[i].name();
        auto [itPrev, inserted] = terminalMap.insert({id, i});
        if (!inserted) {
            s.set(Status::BadTerminal, "Terminal '"+std::string(id)+"' is not unique.");
            s.extend(terminals[i].location());
            if (terminals[itPrev->second].location()) {
                s.extend("Terminal was first defined here");
                s.extend(terminals[itPrev->second].location());
            }
            return false;
        }
    }
    return true;
}

bool HierarchicalModel::buildParameterMap(Status& s) {
    // Check uniqueness 
    auto& parsedSubcircuit = static_cast<const PTSubcircuitDefinition&>(parsedModel);
    if (!parsedSubcircuit.parameters().verify(s)) {
        return false;
    }

    // Build parameter map, check for uniqueness across expressions, too
    parameterMap.clear();
    // Scan parameters
    auto& parameters = parsedSubcircuit.parameters().values();
    ParameterIndex n = parameters.size();
    for(decltype(n) i=0; i<n; i++) {
        auto id = parameters[i].name();
        auto [itPrev, inserted] = parameterMap.insert({id, i});
    }
    return true;
}

std::tuple<TerminalIndex, bool> HierarchicalModel::terminalIndex(Id nodeName) const {
    auto it = terminalMap.find(nodeName);
    if (it==terminalMap.end()) {
        // Not a terminal
        return std::make_tuple(0, false);
    }
    return std::make_tuple(it->second, true);
}

Id HierarchicalModel::terminalName(TerminalIndex ndx) const {
    return static_cast<const PTSubcircuitDefinition&>(parsedModel).terminals()[ndx].name(); 
}

std::tuple<ParameterIndex, bool> HierarchicalModel::parameterIndex(Id name) const {
    auto it = parameterMap.find(name);
    if (it==parameterMap.end()) {
        return std::make_tuple(0, false);
    }
    return std::make_tuple(it->second, true);
}

Id HierarchicalModel::parameterName(ParameterIndex ndx) const {
    auto& parsedSubcircuit = static_cast<const PTSubcircuitDefinition&>(parsedModel);
    if (ndx<parameterCount())
        return parsedSubcircuit.parameters().values().at(ndx).name();
    else
        return Id::none;
}

std::tuple<Value::Type,bool> HierarchicalModel::parameterType(ParameterIndex ndx, Status& s) const {
    auto& parsedSubcircuit = static_cast<const PTSubcircuitDefinition&>(parsedModel);
    if (ndx>=parameterCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return std::make_tuple(Value::Type::Int, false);
    }
    return std::make_tuple(parsedSubcircuit.parameters().values().at(ndx).val().type(), true);
}

bool HierarchicalModel::getParameter(ParameterIndex ndx, Value& v, Status& s) const {
    auto& parsedSubcircuit = static_cast<const PTSubcircuitDefinition&>(parsedModel);
    if (ndx>=parameterCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return false;
    }
    v = parsedSubcircuit.parameters().values().at(ndx).val();
    return true;
}

std::tuple<bool,bool> HierarchicalModel::setParameter(ParameterIndex ndx, const Value& v, Status& s) {
    s.set(Status::Unsupported, std::string("Model parameters cannot be set on a subcircuit definition."));
    return std::make_tuple(false, false);
}

void HierarchicalModel::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    auto& parsedSubcircuit = static_cast<const PTSubcircuitDefinition&>(parsedModel);
    os << pfx << "Hierarchical model " << std::string(name()) << " of device " << device()->name() << "\n";
    if (terminalCount()>0) {
        os << pfx << "  Terminals: ";
        for(int i=0; i<terminalCount(); i++) {
            os << std::string(parsedSubcircuit.terminals().at(i).name()) << " ";
        }
        os << "\n";
    }
    if (parameterCount()>0) {
        os << pfx << "  Parameter defaults:\n";
        for(int i=0; i<parameterCount(); i++) {
            Value v;
            getParameter(i, v);
            os << pfx << "    " << std::string(parameterName(i)) << " = " << v << " (" << v.typeName() << ")\n";
        }
    }
    if (parsedSubcircuit.parameters().expressionCount()>0) {
        os << pfx << "  Expressions:\n";
        const auto& ex = parsedSubcircuit.parameters().expressions();
        for(int i=0; i<ex.size(); i++) {
            os << pfx << "    " << std::string(ex.at(i).name()) << " = " << ex.at(i).rpn() << "\n";
        }
    }
}

Instance* HierarchicalModel::createInstance(Circuit& circuit, Instance* parentInstance, RpnEvaluator& evaluator, Context* externalContext, const PTInstance& parsedInstance, InstantiationData& idata, Status& s) {
    // Instance with empty name is the toplevel instance, set its name here
    auto instanceName = parsedInstance.name();
    // Default subcircuit instance name is __toplevel__
    if (!instanceName) {
        instanceName = Id("__toplevel__");
    }
    // Is parent given, but not hierarchical
    if (parentInstance && !parentInstance->model()->device()->isHierarchical()) {
        return nullptr;
    }
    // Translate
    if (parentInstance) {
        instanceName = parentInstance->translate(instanceName);
    }
    // Create instance
    auto* instance = new HierarchicalInstance(this, instanceName, parentInstance, parsedInstance, s);
    if (!instance->checkFlags(Instance::Flags::IsValid)) {
        s.extend(parsedInstance.location());
        delete instance;
        return nullptr;
    }

    // Add to instanceMap of circuit
    if (!circuit.add(instance, s)) {
        delete instance;
        return nullptr;
    }

    // Sanity check (connections should not exceed terminals)
    if (parsedInstance.connections().size()>terminalCount()) {
        s.set(Status::Range, "Too many terminals specified at instantiation.");
        s.extend(parsedInstance.connections().at(terminalCount()).location());
        return nullptr;
    }

    // Bind connected terminals
    auto& terms = parsedInstance.connections();
    TerminalIndex i=0;
    for(auto it=terms.cbegin(); it!=terms.cend(); ++it, i++) {
        // Translate, if needed
        auto nodeName = it->name();
        nodeName = parentInstance->translateNode(circuit, nodeName);
        auto node = circuit.getNode(nodeName, Node::Flags::PotentialNode, s);
        if (node == nullptr) {
            s.extend(std::string("Failed to obtain node '"+std::string(nodeName)+"'. from simulator"));
            s.extend(it->location());
            return nullptr;
        }
        if (!instance->bindTerminal(i, node, s)) {
            s.extend(std::string("Failed to bind terminal ")+std::to_string(i+1)+".");
            s.extend(it->location());
            return nullptr;
        }
    }

    // Set instance's parameters, use the evaluator whose latest context is the parent instance's context
    // Note that toplevel instance has no context
    auto [ok1, changed] = instance->setParameters(parsedInstance.parameters(), evaluator, s);
    if (!ok1) {
        return nullptr;
    }

    // Build all children (models, nodes, instances)

    // Prepare instance build data
    // At this point it is used for detecting recursion
    idata.addAncestor(instance); 

    // Enter context, add to context search path if the instance being created is a toplevel instance, rebuild context
    auto [ok2, contextMarker] = instance->enterContext(circuit, externalContext, parentInstance==nullptr, true, s);
    if (!ok2) {
        return nullptr;
    }
    
    // Build subhierarchy
    auto hierarchyOk = instance->buildHierarchy(circuit, evaluator, idata, s);

    // Leave context
    Instance::revertContext(circuit, contextMarker);
    
    // Restore instance building data
    idata.removeAncestor(instance);

    if (!hierarchyOk) {
        return nullptr;
    }
    
    return instance;
}


HierarchicalInstance::HierarchicalInstance(HierarchicalModel* mod, Id name, Instance* parentInstance, const PTInstance& parsedInstance, Status& s) 
    : Instance(mod, name, parentInstance, parsedInstance), connectedTerminalCount(0), parameters(nullptr) {
    // Terminals array, resize it
    connections.resize(mod->terminalCount());

    // By default all terminals are unconnected
    for(TerminalIndex i=0; i<mod->terminalCount(); i++) {
        connections[i] = nullptr;
    }

    // Parameters array
    parameters = new Value[mod->parameterCount()];

    // Copy defaults from subcircuit definition
    for(ParameterIndex i=0; i<mod->parameterCount(); i++) {
        mod->getParameter(i, parameters[i]);
    }

    // Clear Flags::NeedsSetup because hierarchical instances do not need to be set up
    clearFlags(Flags::NeedsSetup);

    setFlags(Flags::IsValid);
}

HierarchicalInstance::~HierarchicalInstance() {
    delete [] parameters;
}

Id HierarchicalInstance::translateNode(Circuit& circuit, Id nodeName) {
    // First check if node is a terminal
    // Terminals shadow global nodes
    // If node is a terminal, return corresponding connection
    auto [ndx, isTerminal] = model()->terminalIndex(nodeName); 
    if (isTerminal) {
        // Is it connected? 
        if (connectedTerminalCount>ndx) {
            // Yes, get the node name of the node connected to this terminal
            return connections[ndx]->name();
        }
        // If terminal with node's name is not connected, treat it as an internal node
    }

    // The rest is same as default behavior defined for Instance
    return Instance::translateNode(circuit, nodeName);
}

bool HierarchicalInstance::bindTerminal(TerminalIndex n, Node* node, Status& s) {
    if (n>=terminalCount()) {
        s.set(Status::Range, "Too many connections specified.");
        return false;
    }
    connections[n] = node;
    // Update connected terminal count, assume that there are no unconnected 
    // terminals between two connected terminals
    if (n+1>connectedTerminalCount) {
        connectedTerminalCount = n+1;
    }
    return true;
}

Node* HierarchicalInstance::terminal(TerminalIndex n, Status& s) const {
    if (n>=terminalCount()) {
        s.set(Status::Range, "Terminal not found.");
        return nullptr;
    }
    return connections[n];
}

bool HierarchicalInstance::unbindTerminals(Circuit& circuit, Status& s) {
    for(decltype(connectedTerminalCount) i=0; i<connectedTerminalCount; i++) {
        if (!circuit.releaseNode(connections[i], s)) {
            return false;
        }
    }
    return true;
}

std::tuple<ParameterIndex, bool> HierarchicalInstance::parameterIndex(Id name) const {
    return model()->parameterIndex(name);
}

Id HierarchicalInstance::parameterName(ParameterIndex ndx) const {
    if (ndx<parameterCount())
        return model()->parameterName(ndx);
    else
        return Id::none;
}

std::tuple<Value::Type,bool> HierarchicalInstance::parameterType(ParameterIndex ndx, Status& s) const {
    if (ndx>=parameterCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return std::make_tuple(Value::Type::Int, false);
    }
    return std::make_tuple(parameters[ndx].type(), true);
}

bool HierarchicalInstance::getParameter(ParameterIndex ndx, Value& v, Status& s) const {
    if (ndx>=parameterCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return false;
    }
    v = parameters[ndx];
    return true;
}

std::tuple<bool,bool> HierarchicalInstance::setParameter(ParameterIndex ndx, const Value& v, Status& s) {
    bool changed = parameters[ndx] == v;
    // Mark instance for parameter propagation
    if (changed) {
        setFlags(Instance::Flags::ParamsChanged);
    }
    parameters[ndx] = v;

    return std::make_tuple(true, changed);
}

std::tuple<bool, size_t> HierarchicalInstance::enterContext(Circuit& circuit, Context* externalContext, bool addToPath, bool rebuild, Status& s) {
    // Get evaluator
    auto& evaluator = circuit.paramEvaluator();

    // Get context stack
    ContextStack& cs = evaluator.contextStack();

    // Get stack marker for reverting to initial stack state before abort
    auto stackMarker = evaluator.contextMarker();

    // Enter new context
    // If it is external, add it to path because external context is used 
    // only for toplevel instances whose context is global for all 
    // of the hierarchy. 
    cs.enter(externalContext, addToPath);
    
    if (rebuild) {
        // Clear context before rebuilding
        cs.at().clear();

        // Get parsed subcircuit definition
        auto& parsedSubcircuit = static_cast<const PTSubcircuitDefinition&>(model()->parsedModel);
        
        // Load value parameters
        for(ParameterIndex i=0; i<parameterCount(); i++) {
            Value v;
            getParameter(i, v);
            if (!cs.insert(parameterName(i), v, s)) {
                s.extend(parsedSubcircuit.parameters().values().at(i).location());
                revertContext(circuit, stackMarker);
                return std::make_tuple(false, stackMarker);
            }
        }
        // Load expression parameters, compute them from subcircuit definition
        auto& ep = parsedSubcircuit.parameters().expressions();
        for(auto it=ep.cbegin(); it!=ep.cend(); ++it) {
            Value res;
            if (!evaluator.evaluate(it->rpn(), res, s)) {
                revertContext(circuit, stackMarker);
                return std::make_tuple(false, stackMarker);
            }
            if (!cs.insert(it->name(), std::move(res), s)) {
                s.extend(it->location());
                revertContext(circuit, stackMarker);
                return std::make_tuple(false, stackMarker);
            }
        }
    }

    return std::make_tuple(true, stackMarker);
}

bool HierarchicalInstance::propagateParameters(Circuit& circuit, RpnEvaluator& evaluator, Status& s) {
    // We already have an established context
    
    // Propagate parameters to submodels
    auto& parsedSubcircuit = static_cast<const PTSubcircuitDefinition&>(model()->parsedModel);
    auto nSubModels = parsedSubcircuit.models().size();
    for(decltype(nSubModels) i=0; i<nSubModels; i++) {
        auto subModelPtr = childModels_[i];
        auto& parsedSubmodel = parsedSubcircuit.models()[i];
        // Propagate only expressions
        auto [ok, changed] = subModelPtr->setParameters(parsedSubmodel.parameters().expressions(), evaluator, s);
        if (!ok) {
            return false;
        }
        // Mark child model for setup
        if (changed) {
            subModelPtr->setFlags(Model::Flags::NeedsSetup);
        }
    }
    // Propagate parameters to subinstances
    auto nSubInstances = parsedSubcircuit.instances().size();
    for(decltype(nSubInstances) i=0; i<nSubInstances; i++) { 
        auto subInstancePtr = childInstances_[i];
        auto& parsedSubinstance = parsedSubcircuit.instances()[i]; 
        auto [ok, changed] = subInstancePtr->setParameters(parsedSubinstance.parameters().expressions(), evaluator, s);
        if (!ok) {
            return false;
        }
        // Mark child instance for parameter propagation
        if (changed) {
            subInstancePtr->setFlags(Instance::Flags::ParamsChanged);
        }
    }

    // Parameter change has been propagated
    clearFlags(Instance::Flags::ParamsChanged);

    // No need to mark instance for setup - hierarchical instances don't need setup
    
    return true;
}

bool HierarchicalInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s) {
    // We already have an established context
    
    // Get parsed subcircuit
    auto& parsedSubcircuit = static_cast<const PTSubcircuitDefinition&>(model()->parsedModel);

    // Bind unconnected terminals to internal nodes
    auto& defTerms = parsedSubcircuit.terminals();
    auto i = connectedTerminalCount;
    for(auto it=defTerms.begin()+connectedTerminalCount; it!=defTerms.end(); ++it, i++) {
        auto nodeName = it->name();
        nodeName = parent()->translateNode(circuit, nodeName);
        auto node = circuit.getNode(nodeName, Node::Flags::PotentialNode, s);
        if (node == nullptr) {
            s.extend(std::string("Failed to obtain internal node '"+std::string(nodeName)+"'. from simulator"));
            s.extend(it->location());
            return false;
        }
        if (!bindTerminal(i, node, s)) { 
            s.extend(std::string("Failed to bind terminal ")+std::to_string(i+1)+".");
            s.extend(it->location());
            return false;
        }
    }

    // Create models
    for(auto it=parsedSubcircuit.models().cbegin(); it!=parsedSubcircuit.models().cend(); ++it) {
        // Find device
        // This returns a Device* pointer
        auto* dev = circuit.findDevice(it->device());
        if (!dev || dev->isHierarchical()) {
            // Device not found or device is a hierarchical device
            s.set(Status::NotFound, std::string("Device '")+std::string(it->device())+"' not found.");
            s.extend(it->location());
            return false;
        }

        // Create model (generic method)
        auto* subMod = dev->createModel(circuit, this, evaluator, *it, s);
        if (!subMod) {
            return false;
        }
    }
    
    // Create instances
    for(auto it=parsedSubcircuit.instances().cbegin(); it!=parsedSubcircuit.instances().cend(); ++it) {
        // Find master
        // 1) try local models (prefix by instance path)
        auto localMasterId = translate(it->masterName());
        auto* mod = circuit.findModel(localMasterId);

        // 2) try local definitions of hierarchical devices
        //    (prefix by definition path)
        if (!mod) {
            auto localDefId = idata.translateDefinition(it->masterName());
            mod = circuit.findModel(localDefId);
        }

        // 3) try global models/definitions of hierarchical devices
        if (!mod) {
            mod = circuit.findModel(it->masterName());
        }

        // We do not try device names (i.e. default models)
        // If one wants to have default models, one should create them explicitly, e.g. 
        // add a line of the following form to the netlist
        //   model <device name> <device name>
        if (!mod) {
            s.set(Status::NotFound, std::string("Master '")+std::string(it->masterName())+"' not found.");
            s.extend(it->location());
            return false;
        }

        // Check for hierarchical recursion
        if (idata.isAncestor(mod)) {
            // Recursion detected
            s.set(Status::HierarchicalRecursion, "Hierarchical recursion detected.");
            s.extend(it->location());
            return false;
        }

        // Create child instance
        Instance *subInst = mod->createInstance(circuit, this, evaluator, nullptr, *it, idata, s);
        if (!subInst) {
            return false;
        }
    }

    // At this point all parameters of the hierarchical instance 
    // we created are propagated down the hierarchy. 
    // We can clear the FParamsChanged flag. 
    clearFlags(Instance::Flags::ParamsChanged);

    return true;
}

bool HierarchicalInstance::deleteHierarchy(Circuit& circuit, Status& s) {
    // Delete instances
    for(auto it=childInstances_.begin(); it!=childInstances_.end(); ++it) {
        // Delete subhierarchy
        if (!(*it)->deleteHierarchy(circuit, s)) {
            return false;
        }
        // Unbind nodes
        if (!(*it)->unbindTerminals(circuit, s)) {
            return false;
        }
        // Delete instance
        circuit.remove(*it);
    }

    // Delete models
    for(auto it=childModels_.begin(); it!=childModels_.end(); ++it) {
        // Delete model
        circuit.remove(*it);
    }

    return true;
}

void HierarchicalInstance::addChild(Instance* child) {
    childInstances_.push_back(child);
}

void HierarchicalInstance::addChild(Model* child) {
    childModels_.push_back(child);
}


void HierarchicalInstance::dump(int indent, const Circuit& circuit, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "Hierarchical instance " << std::string(name()) << " of definition " << std::string(model()->name()) << "\n";
    if (terminalCount()>0) {
        os << pfx << "  Terminals: ";
        for(int i=0; i<terminalCount(); i++) {
            os << connections[i]->name() << " ";
        }
        os << "\n";
    }
    if (parameterCount()>0) {
        os << pfx << "  Parameters:\n";
        for(int i=0; i<model()->parameterCount(); i++) {
            os << pfx << "    " << std::string(model()->parameterName(i)) << " = " << parameters[i] << "\n";
        }
    }
}

}
