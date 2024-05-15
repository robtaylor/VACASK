#ifndef __DEVBUILTIN_DEFINED
#define __DEVBUILTIN_DEFINED

#include "devbase.h"
#include "circuit.h"
#include "common.h"

namespace NAMESPACE {

template<typename ModelParams, typename InstanceParams, typename InstanceData> class BuiltinDevice : public Device {
public:
    BuiltinDevice(Id name, const Loc& location=Loc::bad); 
    virtual ~BuiltinDevice();

    BuiltinDevice           (const BuiltinDevice&)  = delete;
    BuiltinDevice           (      BuiltinDevice&&) = default;
    BuiltinDevice& operator=(const BuiltinDevice&)  = delete;
    BuiltinDevice& operator=(      BuiltinDevice&&) = default;

    virtual bool operator==(const Device& other) const { return dynamic_cast<const BuiltinDevice<ModelParams, InstanceParams, InstanceData>*>(&other)!=nullptr; };
    virtual bool isHierarchical() const { return false; };
    virtual std::tuple<bool, bool, bool> setup(Circuit& circuit, bool force, Status& s=Status::ignore);
    virtual bool collapseNodes(Circuit& circuit, Status& s=Status::ignore);
    virtual bool populateStructures(Circuit& circuit, Status& s=Status::ignore);
    virtual bool preAnalysis(Circuit& circuit, Status& s=Status::ignore);
    virtual bool bind(
        Circuit& circuit, 
        KluRealMatrix* matResistReal, KluComplexMatrix* matResistCx, Component compResist, 
        KluRealMatrix* matReactReal, KluComplexMatrix* matReactCx, Component compReact, 
        Status& s=Status::ignore
    );
    virtual bool evalAndLoad(Circuit& circuit, EvalAndLoadSetup& els);
    virtual Model* createModel(Circuit& circuit, Instance* parentInstance, RpnEvaluator& evaluator, const PTModel& ptModel, Status& s=Status::ignore);
    virtual bool isSource() const { return false; };
    virtual bool isVoltageSource() const { return false; };
    virtual void dump(int indent, std::ostream& os) const; 

    // Extra flags set at construction
    static const Device::Flags extraFlags;

    // Terminal count, node names
    TerminalIndex terminalCount;
    std::vector<Id> nodeIds;

    // Build internals
    void defineInternals() {};
    void buildInternals();

    std::unordered_map<Id, TerminalIndex> nodeNameMap;
};


template<typename ModelParams, typename InstanceParams, typename InstanceData> class BuiltinModel : public Model {
public:
    using DeviceType = BuiltinDevice<ModelParams, InstanceParams, InstanceData>;

    BuiltinModel(
        DeviceType* device, Id name, Instance* parent, 
        const PTModel& parsedModel, Status& s=Status::ignore
    );
    virtual ~BuiltinModel();

    BuiltinModel           (const BuiltinModel&)  = delete;
    BuiltinModel           (      BuiltinModel&&) = default;
    BuiltinModel& operator=(const BuiltinModel&)  = delete;
    BuiltinModel& operator=(      BuiltinModel&&) = default;

    virtual ParameterIndex parameterCount() const { return params.parameterCount(); };
    virtual std::tuple<ParameterIndex, bool> parameterIndex(Id name) const { return params.parameterIndex(name); };
    virtual Id parameterName(ParameterIndex ndx) const { return params.parameterName(ndx); };
    virtual std::tuple<Value::Type,bool> parameterType(ParameterIndex ndx, Status& s=Status::ignore) const { return params.parameterType(ndx, s); };
    virtual bool getParameter(ParameterIndex ndx, Value& v, Status& s=Status::ignore) const { return params.getParameter(ndx, v, s); };
    virtual std::tuple<bool,bool> setParameter(ParameterIndex ndx, const Value& v, Status& s=Status::ignore) { 
        auto [ok, changed] = params.setParameter(ndx, v, s); 
        if (changed) {
            setFlags(Flags::NeedsSetup);
        }
        return std::make_tuple(ok, changed);
    };
    virtual std::tuple<bool, bool, bool> setup(Circuit& circuit, bool force, Status& s=Status::ignore);
    virtual Instance* createInstance(Circuit& circuit, Instance* parentInstance, RpnEvaluator& evaluator, Context* externalContext, const PTInstance& parsedInstance, InstantiationData& idata, Status& s=Status::ignore);
    virtual void dump(int indent, std::ostream& os) const;

    // Device access (as DeviceType)
    DeviceType* device() { return reinterpret_cast<DeviceType*>(device_); };
    const DeviceType* device() const { return reinterpret_cast<const DeviceType*>(device_); };
    
    // Wrapper for inlining setup functions
    std::tuple<bool, bool, bool> setupCore(Circuit& circuit, bool force, Status& s=Status::ignore);

    // Sets up this particular model
    // Define specialization in cpp file for each builtin device type
    std::tuple<bool, bool, bool> setupWorker(Circuit& circuit, Status& s=Status::ignore) { return std::make_tuple(true, false, false); }

    // Sets up this particular model
    // Define specialization in cpp file for each builtin device type
    bool preAnalysisWorker(Circuit& circuit, Status& s=Status::ignore) { return true; };

private:
    IStruct<ModelParams> params;
};


template<typename ModelParams, typename InstanceParams, typename InstanceData> class BuiltinInstance : public Instance {
public:
    using ModelType = BuiltinModel<ModelParams, InstanceParams, InstanceData>;

    friend ModelType::DeviceType;

    BuiltinInstance(Model* model, Id name, Instance* parent, const PTInstance& parsedInstance, Status &s=Status::ignore);
    virtual ~BuiltinInstance();

    BuiltinInstance           (const BuiltinInstance&)  = delete;
    BuiltinInstance           (      BuiltinInstance&&) = default;
    BuiltinInstance& operator=(const BuiltinInstance&)  = delete;
    BuiltinInstance& operator=(      BuiltinInstance&&) = default;

    virtual ParameterIndex parameterCount() const { return params.parameterCount(); };
    virtual std::tuple<ParameterIndex, bool> parameterIndex(Id name) const { return params.parameterIndex(name); };
    virtual Id parameterName(ParameterIndex ndx) const { return params.parameterName(ndx); };
    virtual std::tuple<Value::Type,bool> parameterType(ParameterIndex ndx, Status& s=Status::ignore) const { return params.parameterType(ndx, s); };
    virtual bool getParameter(ParameterIndex ndx, Value& v, Status& s=Status::ignore) const { return params.getParameter(ndx, v, s); };
    virtual std::tuple<bool,bool> setParameter(ParameterIndex ndx, const Value& v, Status& s=Status::ignore) { 
        auto [ok, changed] = params.setParameter(ndx, v, s); 
        if (changed) {
            setFlags(Flags::ParamsChanged);
        }
        return std::make_tuple(ok, changed);
    };
    virtual std::tuple<ParameterIndex, bool> principalParameterIndex() const { return std::make_tuple(0, false); };
    virtual TerminalIndex staticNodeCount() const { return model()->device()->nodeIds.size(); };
    virtual TerminalIndex terminalCount() const { return model()->device()->terminalCount; };
    virtual std::tuple<TerminalIndex, bool> nodeIndex(Id name) const {
        auto it = model()->device()->nodeNameMap.find(name);
        if (it!=model()->device()->nodeNameMap.end()) {
            return std::make_tuple(it->second, true);
        } else {
            return std::make_tuple(0, false);
        }
    };
    virtual Id nodeName(TerminalIndex ndx) const { return model()->device()->nodeIds[ndx]; };
    virtual bool bindTerminal(TerminalIndex n, Node* node, Status& s=Status::ignore);
    virtual Node* terminal(TerminalIndex n, Status& s=Status::ignore) const;
    virtual bool unbindTerminals(Circuit& circuit, Status& s=Status::ignore);
    virtual std::tuple<bool, bool> subhierarchyChanged(Circuit& circuit, Status& s=Status::ignore) { return std::make_tuple(true, false); }; 
    virtual bool propagateParameters(Circuit& circuit, RpnEvaluator& evaluator, Status& s=Status::ignore) { 
        if (checkFlags(Flags::ParamsChanged)) {
            clearFlags(Flags::ParamsChanged); 
            setFlags(Flags::NeedsSetup);
        }
        return true; 
    };
    virtual bool deleteHierarchy(Circuit& circuit, Status& s=Status::ignore) { return true; }; 
    virtual bool buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s=Status::ignore) { return true; };  
    virtual std::tuple<EquationIndex,EquationIndex> sourceExcitation(Circuit& circuit) const { return std::make_tuple(0, 0); };
    virtual std::tuple<UnknownIndex,UnknownIndex> sourceResponse(Circuit& circuit) const { return std::make_tuple(0, 0); };
    virtual double scaledUnityExcitation() const { return 1.0; };
    virtual double responseScalingFactor() const { return 1.0; };
    virtual ParameterIndex opvarCount() const { return data.parameterCount(); };
    virtual std::tuple<ParameterIndex, bool> opvarIndex(Id name) const { return data.parameterIndex(name); };
    virtual Id opvarName(ParameterIndex ndx) const { return data.parameterName(ndx); };
    virtual std::tuple<Value::Type,bool> opvarType(ParameterIndex ndx, Status& s=Status::ignore) const { return data.parameterType(ndx, s); };
    virtual bool getOpvar(ParameterIndex ndx, Value& v, Status& s=Status::ignore) const { return data.getParameter(ndx, v, s); };
    virtual std::tuple<bool, OutputSource> opvarOutputSource(ParameterIndex ndx, Status& s=Status::ignore) const { return std::make_tuple(false, OutputSource()); };
    virtual std::tuple<bool, bool, bool> setup(Circuit& circuit, bool force, Status& s=Status::ignore);
    virtual void dump(int indent, const Circuit& circuit, std::ostream& os) const;

    // Model access (as ModelType)
    ModelType* model() { return reinterpret_cast<ModelType*>(model_); };
    const ModelType* model() const { return reinterpret_cast<const ModelType*>(model_); };

    // Wrapper for inlining setup functions
    std::tuple<bool, bool, bool> setupCore(Circuit& circuit, bool force, Status& s=Status::ignore);

    // Sets up this particular instance
    // Define specialization in cpp file for each builtin device type
    std::tuple<bool, bool, bool> setupWorker(Circuit& circuit, Status& s=Status::ignore) { return std::make_tuple(true, false, false); }; 

    // Sets up this particular model
    // Define specialization in cpp file for each builtin device type
    bool preAnalysisWorker(Circuit& circuit, Status& s=Status::ignore) { return true; };
    
    bool collapseNodesCore(Circuit& circuit, Status& s) { return true; };
    bool populateStructuresCore(Circuit& circuit, Status& s=Status::ignore);
    bool bindCore(
        Circuit& circuit, 
        KluRealMatrix* matResistReal, KluComplexMatrix* matResistCx, Component compResist, 
        KluRealMatrix* matReactReal, KluComplexMatrix* matReactCx, Component compReact, 
        Status& s=Status::ignore
    );
    bool evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els);
    
protected:
    // Get Jacobian entry pointer
    void jacEntryPtr(double*& destination, EquationIndex e, UnknownIndex u, KluRealMatrix* realMat, KluComplexMatrix* cxMat, Component cxComp);

    // Create internal nodes for unconnected terminals
    bool createNodesForUnconnectedTerminals(Circuit& circuit, Status& s=Status::ignore);

    // Check if all terminals are connected
    bool verifyTerminalsConnected(Status& s=Status::ignore);

    std::vector<Node*> nodes_;
    TerminalIndex connectedTerminalCount;
    StateIndex statesStartIndex;
    IStruct<InstanceParams> params;
    IStruct<InstanceData> data;
};

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
BuiltinDevice<ModelParams, InstanceParams, InstanceData>::BuiltinDevice(Id name, const Loc& location) 
    : Device(name, location) {
    setFlags(extraFlags);
    defineInternals();
    buildInternals();
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
void BuiltinDevice<ModelParams, InstanceParams, InstanceData>::buildInternals() {
    TerminalIndex i=0;
    for(auto id : nodeIds) {
        nodeNameMap.insert({id, i});
        i++;
    }
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
BuiltinDevice<ModelParams, InstanceParams, InstanceData>::~BuiltinDevice() {
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
std::tuple<bool, bool, bool> BuiltinDevice<ModelParams, InstanceParams, InstanceData>::setup(
    Circuit& circuit, bool force, Status& s
) {
    using ModelType = BuiltinModel<ModelParams, InstanceParams, InstanceData>;
    bool unknownsChanged = false;
    bool sparsityChanged = false;
    for(auto model : models()) {
        // This will set up model and all its instances
        auto [ok, tmpUnknowns, tmpSparsity] = static_cast<ModelType*>(model)->setupCore(circuit, force, s);
        unknownsChanged |= tmpUnknowns;
        sparsityChanged |= tmpSparsity;
        if (!ok) {
            return std::make_tuple(false, unknownsChanged, sparsityChanged);
        }
        
    }
    return std::tuple(true, unknownsChanged, sparsityChanged);
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
bool BuiltinDevice<ModelParams, InstanceParams, InstanceData>::preAnalysis(
    Circuit& circuit, Status& s
) {
    using ModelType = BuiltinModel<ModelParams, InstanceParams, InstanceData>;
    using InstanceType = BuiltinInstance<ModelParams, InstanceParams, InstanceData>;
    for(auto model : models()) {
        if (!static_cast<ModelType*>(model)->preAnalysisWorker(circuit, s)) {
            return false;
        }
        for(auto instance : model->instances()) {
            if (!static_cast<InstanceType*>(instance)->preAnalysisWorker(circuit, s)) {
                return false;
            }
        }
    
    }
    return true;
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
bool BuiltinDevice<ModelParams, InstanceParams, InstanceData>::collapseNodes(
    Circuit& circuit, Status& s
) {
    using InstanceType = BuiltinInstance<ModelParams, InstanceParams, InstanceData>;
    for(auto model : models()) {
        for(auto instance : model->instances()) {
            if (!static_cast<InstanceType*>(instance)->collapseNodesCore(circuit, s)) {
                return false;
            } 
        }
    }
    return true;
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
bool BuiltinDevice<ModelParams, InstanceParams, InstanceData>::populateStructures(
    Circuit& circuit, Status& s
) {
    using InstanceType = BuiltinInstance<ModelParams, InstanceParams, InstanceData>;
    for(auto model : models()) {
        for(auto instance : model->instances()) {
            if (!static_cast<InstanceType*>(instance)->populateStructuresCore(circuit, s)) {
                return false;
            }
        }
    }
    return true;
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
bool BuiltinDevice<ModelParams, InstanceParams, InstanceData>::bind(
    Circuit& circuit, 
    KluRealMatrix* matResistReal, KluComplexMatrix* matResistCx, Component compResist, 
    KluRealMatrix* matReactReal, KluComplexMatrix* matReactCx, Component compReact, 
    Status& s
) {
    using InstanceType = BuiltinInstance<ModelParams, InstanceParams, InstanceData>;
    // Call bind() for all instances
    for(auto model : models()) {
        for(auto instance : model->instances()) {
            if (!static_cast<InstanceType*>(instance)->bindCore(
                circuit, 
                matResistReal, matResistCx, compResist, 
                matReactReal, matReactCx, compReact, 
                s
            )) {
                return false;
            }
        }
    }
    return true;
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
bool BuiltinDevice<ModelParams, InstanceParams, InstanceData>::evalAndLoad(
    Circuit& circuit, EvalAndLoadSetup& els
) {
    using InstanceType = BuiltinInstance<ModelParams, InstanceParams, InstanceData>;
    for(auto model : models()) {
        if (model->instanceCount()==0) {
            continue;
        }
        for(auto instance : model->instances()) {
            if (!static_cast<InstanceType*>(instance)->evalAndLoadCore(circuit, els)) {
                return false;
            }
        }
    }
    
    return true;
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
Model* BuiltinDevice<ModelParams, InstanceParams, InstanceData>::createModel(
    Circuit& circuit, Instance* parentInstance, RpnEvaluator& evaluator, const PTModel& parsedModel, Status& s
) {
    using ModelType = BuiltinModel<ModelParams, InstanceParams, InstanceData>;

    auto name = parsedModel.name();

    // If we have a hierarchical parent translate name
    if (parentInstance) {
        name = parentInstance->translate(name);
    }

    // Create model
    auto* model = new ModelType(this, name, parentInstance, parsedModel, s);
    if (!model->checkFlags(Model::Flags::IsValid)) {
        s.extend(parsedModel.location());
        delete model;
        return nullptr;
    }

    // Add to modelMap of circuit
    if (!circuit.add(model, s)) {
        delete model;
        return nullptr;
    }

    // Set model's parameters, use the evaluator whose latest context is the parent instance's context
    auto [ok, changed] = model->setParameters(parsedModel.parameters(), evaluator, s);
    if (!ok) {
        return nullptr;
    }

    return model;
}


template<typename ModelParams, typename InstanceParams, typename InstanceData> 
BuiltinModel<ModelParams, InstanceParams, InstanceData>::BuiltinModel(
   DeviceType* device, Id name, Instance* parent, 
   const PTModel& parsedModel, Status& s
) : Model(device, name, parent, parsedModel) {
   Model::setFlags(Model::Flags::IsValid); 
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
BuiltinModel<ModelParams, InstanceParams, InstanceData>::~BuiltinModel() {
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
std::tuple<bool, bool, bool> BuiltinModel<ModelParams, InstanceParams, InstanceData>::setup(
    Circuit& circuit, bool force, Status& s
) {
    return setupCore(circuit, force, s);
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
std::tuple<bool, bool, bool> BuiltinModel<ModelParams, InstanceParams, InstanceData>::setupCore(
    Circuit& circuit, bool force, Status& s
) {
    using InstanceType = BuiltinInstance<ModelParams, InstanceParams, InstanceData>;

    // Do we need to set up all instances? 
    bool forceAllInstances = false;
    bool unknownsChanged = false;
    bool sparsityChanged = false;
    if (force || checkFlags(Flags::NeedsSetup)) {
        auto [ok, tmpUnknowns, tmpSparsity] = setupWorker(circuit, s);
        unknownsChanged |= tmpUnknowns;
        sparsityChanged |= tmpSparsity;
        if (!ok) {
            // The problem is big enough to abort simulation
            return std::make_tuple(false, unknownsChanged, sparsityChanged);
        }

        clearFlags(Flags::NeedsSetup);
        
        // After model setup all instances have to be set up
        forceAllInstances = true;
    }
    
    for(auto it : instances()) {
        auto [ok, tmpUnknowns, tmpSparsity] = static_cast<InstanceType*>(it)->setupCore(circuit, forceAllInstances, s);
        unknownsChanged |= tmpUnknowns;
        sparsityChanged |= tmpSparsity;
        if (!ok) {
            return std::make_tuple(false, unknownsChanged, sparsityChanged);
        }
    }

    return std::make_tuple(true, unknownsChanged, sparsityChanged);
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
Instance* BuiltinModel<ModelParams, InstanceParams, InstanceData>::createInstance(
    Circuit& circuit, Instance* parentInstance, RpnEvaluator& evaluator, Context* externalContext, const PTInstance& parsedInstance, InstantiationData& idata, Status& s
) {
    using InstanceType = BuiltinInstance<ModelParams, InstanceParams, InstanceData>;
    
    auto name = parsedInstance.name();
    
    // If we have a hierarchical parent translate name
    if (parentInstance) {
        name = parentInstance->translate(name);
    }

    // Create instance
    auto* instance = new InstanceType(this, name, parentInstance, parsedInstance, s);
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
    if (parsedInstance.connections().size()>device()->terminalCount) {
        s.extend(parsedInstance.connections().at(device()->terminalCount).location());
        return nullptr;
    }

    // Bind terminals
    auto& terms = parsedInstance.connections();
    TerminalIndex i=0;
    for(auto it=terms.cbegin(); it!=terms.cend(); ++it, i++) {
        // Translate, if needed
        auto nodeName = it->name();
        nodeName = parentInstance->translateNode(circuit, nodeName);
        auto node = circuit.getNode(nodeName, Node::Flags::PotentialNode, s);
        if (node == nullptr) {
            s.extend(std::string("Failed to obtain node '"+std::string(nodeName)+"' from simulator."));
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
    auto [ok, changed] = instance->setParameters(parsedInstance.parameters(), evaluator, s);
    if (!ok) {
        return nullptr;
    }

    // Due to setting parameters the ParamsChanged flag is probably set indicating 
    // that parameter propagation is needed. Because we just created the instance 
    // propagation is not needed so clear the flag. 
    instance->clearFlags(Instance::Flags::ParamsChanged);

    // Build all children (models, nodes, instances)
    // Use the evaluator whose latest context is the parent instance's context
    if (!instance->buildHierarchy(circuit, evaluator, idata, s)) {
        return nullptr;
    }

    return instance;
}



template<typename ModelParams, typename InstanceParams, typename InstanceData> 
BuiltinInstance<ModelParams, InstanceParams, InstanceData>::BuiltinInstance(Model* model, Id name, Instance* parent, const PTInstance& parsedInstance, Status &s) 
    : Instance(model, name, parent, parsedInstance), connectedTerminalCount(0) {
    // Create nodes (terminals+internal nodes) list
    // Resize nodes vector
    auto n = static_cast<ModelType*>(model)->device()->nodeIds.size();
    nodes_.resize(n);
    
    // By default all terminals/nodes are unconnected
    for(TerminalIndex i=0; i<n; i++) {
        nodes_[i] = nullptr;
    }

    setFlags(Flags::IsValid);
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
BuiltinInstance<ModelParams, InstanceParams, InstanceData>::~BuiltinInstance() {
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
bool BuiltinInstance<ModelParams, InstanceParams, InstanceData>::bindTerminal(TerminalIndex n, Node* node, Status& s) {
    if (n>=model()->device()->terminalCount) {
        s.set(Status::Range, "Too many connections specified.");
        return false;
    }
    nodes_[n] = node;
    // Update connected terminal count, assume that there are no unconnected 
    // terminals between two connected terminals
    if (n+1>connectedTerminalCount) {
        connectedTerminalCount = n+1;
    }
    return true;
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
Node* BuiltinInstance<ModelParams, InstanceParams, InstanceData>::terminal(TerminalIndex n, Status& s) const {
    if (n>=model()->device()->terminalCount) {
        s.set(Status::Range, "Terminal not found.");
        return nullptr;
    }
    return nodes_[n];
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
bool BuiltinInstance<ModelParams, InstanceParams, InstanceData>::unbindTerminals(Circuit& circuit, Status& s) {
    for(decltype(connectedTerminalCount) i=0; i<connectedTerminalCount; i++) {
        if (!circuit.releaseNode(nodes_[i], s)) {
            return false;
        }
    }
    return true;
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
std::tuple<bool, bool, bool> BuiltinInstance<ModelParams, InstanceParams, InstanceData>::setup(Circuit& circuit, bool force, Status& s) {
    return setupCore(circuit, force, s);
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
std::tuple<bool, bool, bool> BuiltinInstance<ModelParams, InstanceParams, InstanceData>::setupCore(Circuit& circuit, bool force, Status& s) {
    auto [ok, unknownsChanged, sparsityChanged] = setupWorker(circuit, s);
    if (ok) {
        clearFlags(Flags::NeedsSetup); 
    }
    return std::make_tuple(ok, unknownsChanged, sparsityChanged);
}
    
template<typename ModelParams, typename InstanceParams, typename InstanceData> 
void BuiltinInstance<ModelParams, InstanceParams, InstanceData>::jacEntryPtr(double*& destination, EquationIndex e, UnknownIndex u, KluRealMatrix* realMat, KluComplexMatrix* cxMat, Component cxComp) {
    if (realMat) {
        // Set destination to point to teal matrix entry
        destination = realMat->elementPtr(e, u, cxComp);
    } else if (cxMat) {
        // Set destination to point to complex matrix entry
        destination = cxMat->elementPtr(e, u, cxComp); 
    }
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
bool BuiltinInstance<ModelParams, InstanceParams, InstanceData>::createNodesForUnconnectedTerminals(Circuit& circuit, Status& s) {
    for(TerminalIndex i=connectedTerminalCount; i<model()->device()->terminalCount(); i++) {
        // Create node name
        Id nodeName = translate(model()->device()->nodeName(i));
        
        // Create/get node
        auto node = circuit.getNode(nodeName, Node::Flags::PotentialNode, s);
        if (node==nullptr) {
            s.extend(std::string("Failed to obtain internal node '"+std::string(nodeName)+"' from simulator."));
            s.extend(location());
            return false;
        }
        
        // Bind node, do not use bindTerminal() because it increases connectedTerminalCount
        nodes_[i] = node;
    }
    return true;
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
bool BuiltinInstance<ModelParams, InstanceParams, InstanceData>::verifyTerminalsConnected(Status& s) {
    // If we require all terminals to be connected, do this
    if (connectedTerminalCount<model()->device()->terminalCount) {
        s.set(Status::Conflicting, "Instance has "+std::to_string(model()->device()->terminalCount)+" terminal(s) but only "+std::to_string(connectedTerminalCount)+" connection(s).");
        return false;
    }
    return true;
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
void BuiltinDevice<ModelParams, InstanceParams, InstanceData>::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "Builtin device " << name() << "\n";
    auto nNodes = nodeIds.size();
    auto nTerminals = terminalCount;
    os << pfx << "  Static nodes (terminals + static internals = "+std::to_string(nNodes)+", terminals = "+std::to_string(nTerminals)+")\n";
    if (nNodes>0) {
        for(decltype(nNodes) i=0; i<nNodes; i++) {
            os << pfx << "    " << std::string(nodeIds[i]) << "\n";
        }
    }
    auto mpar = Introspection<ModelParams>::count();
    if (mpar>0) {
        os << pfx << "  Model parameters:\n";
        for(ParameterIndex i=0; i<mpar; i++) {
            os << pfx << "    " << std::string(Introspection<ModelParams>::name(i)) << "\n";
        }
    }
    auto ipar = Introspection<InstanceParams>::count();
    if (ipar>0) {
        os << pfx << "  Instance parameters:\n";
        for(ParameterIndex i=0; i<ipar; i++) {
            os << pfx << "    " << std::string(Introspection<InstanceParams>::name(i)) << "\n";
        }
    }
    auto opvars = Introspection<InstanceData>::count();
    if (ipar>0) {
        os << pfx << "  Opvars:\n";
        for(ParameterIndex i=0; i<opvars; i++) {
            os << pfx << "    " << Introspection<InstanceData>::name(i) << "\n";
        }
    }
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
void BuiltinModel<ModelParams, InstanceParams, InstanceData>::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "Builtin device model " << std::string(name()) << " of device " << device()->name() << "\n";
    if (parameterCount()>0) {
        os << pfx << "  Parameters:\n";
        auto np = parameterCount();
        for(decltype(np) i=0; i<np; i++) {
            Value v;
            getParameter(i, v);
            os << pfx << "    " << std::string(parameterName(i)) << " = " << v << " (" << v.typeName() << ")\n";
        }
    }
}

template<typename ModelParams, typename InstanceParams, typename InstanceData> 
void BuiltinInstance<ModelParams, InstanceParams, InstanceData>::dump(int indent, const Circuit& circuit, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "Builtin device instance " << std::string(name()) << " of model " << model()->name() << "\n";
    if (terminalCount()>0) {
        os << pfx << "  Terminals: ";
        auto termCount = terminalCount();
        for(decltype(termCount) i=0; i<termCount; i++) {
            os << terminal(i)->name() << " ";
        }
        os << "\n";
    }
    if (parameterCount()>0) {
        os << pfx << "  Parameters:\n";
        auto np = parameterCount();
        for(decltype(np) i=0; i<np; i++) {
            Value v;
            getParameter(i, v);
            os << pfx << "    " << std::string(parameterName(i)) << " = " << v << " (" << v.typeName() << ")\n";
        }
    }
}

}

#endif

