#ifndef __HIERDEVICE_DEFINED
#define __HIERDEVICE_DEFINED

#include <vector>
#include "identifier.h"
#include "devbase.h"
#include "parseroutput.h"
#include "common.h"


namespace NAMESPACE {

class Circuit;

class HierarchicalModel;

class HierarchicalDevice : public Device {
public:
    friend class HierarchicalModel;

    HierarchicalDevice(Id name, Status& s=Status::ignore);
    virtual ~HierarchicalDevice();

    HierarchicalDevice           (const HierarchicalDevice&)  = delete;
    HierarchicalDevice           (      HierarchicalDevice&&) = default;
    HierarchicalDevice& operator=(const HierarchicalDevice&)  = delete;
    HierarchicalDevice& operator=(      HierarchicalDevice&&) = default;

    virtual bool operator==(const Device& other) const;
    virtual bool isHierarchical() const { return true; }; 
    // Generic model creation is not allowed, constructor must be used
    virtual Model* createModel(Circuit& circuit, Instance* parentInstance, RpnEvaluator& evaluator, const PTModel& parsedSubcircuit, Status& s=Status::ignore) {
        DBGCHECK(true, "Internal error, requested generic model creation from hierarchical device.");
        return nullptr;
    }
    virtual void dump(int indent, std::ostream& os) const;
};


class HierarchicalInstance;

class HierarchicalModel : public Model {
public:
    friend class HierarchicalInstance;

    // PTSubcircuitDefinition is derived from PTModel. 
    // Therefore we store a reference to PTModel in Model and static_cast it
    // to PTSubcircuitDefinition, when we need it. 
    // The model is added to the list of models of the given device so that when
    // a device is destroyed all the models in the list are destroyed, too. 
    // HierarchicalModel(HierarchicalDevice* dev, Id name, Status& s=Status::ignore);
    HierarchicalModel(HierarchicalDevice* dev, Id name, Instance* parentInstance, const PTSubcircuitDefinition& parsedSubcircuit, Status& s=Status::ignore);
    virtual ~HierarchicalModel();

    HierarchicalModel           (const HierarchicalModel&)  = delete;
    HierarchicalModel           (      HierarchicalModel&&) = default;
    HierarchicalModel& operator=(const HierarchicalModel&)  = delete;
    HierarchicalModel& operator=(      HierarchicalModel&&) = default;

    virtual ParameterIndex parameterCount() const { 
        return static_cast<const PTSubcircuitDefinition&>(parsedModel_).parameters().valueCount(); 
    };
    virtual std::tuple<ParameterIndex, bool> parameterIndex(Id name) const;
    virtual Id parameterName(ParameterIndex ndx) const;
    virtual std::tuple<Value::Type,bool> parameterType(ParameterIndex ndx, Status& s=Status::ignore) const;
    virtual bool getParameter(ParameterIndex ndx, Value& v, Status& s=Status::ignore) const;
    virtual std::tuple<bool,bool> setParameter(ParameterIndex ndx, const Value& v, Status& s=Status::ignore);
    virtual Instance* createInstance(Circuit& circuit, Instance* parentInstance, RpnEvaluator& evaluator, Context* externalContext, const PTInstance& ptInstance, InstantiationData& idata, Status& s=Status::ignore);
    virtual void dump(int indent, std::ostream& os) const;

    // Device access (as HierarchicalDevice)
    HierarchicalDevice* device() { return static_cast<HierarchicalDevice*>(device_); };
    const HierarchicalDevice* device() const { return static_cast<HierarchicalDevice*>(device_); };
    
    // Get terminal count
    TerminalIndex terminalCount() const { 
        return static_cast<const PTSubcircuitDefinition&>(parsedModel_).terminals().size(); 
    };
    
    // Find terminal by name
    // Return value: terminal index, found
    std::tuple<TerminalIndex, bool> terminalIndex(Id nodeName) const;

    // Terminal name
    Id terminalName(TerminalIndex ndx) const;

    // Expose parser table
    const PTSubcircuitDefinition& definition() const { return static_cast<const PTSubcircuitDefinition&>(parsedModel_); };
    
protected:
    // Build a map from terminal name to terminal index
    bool buildTerminalMap(Status& s=Status::ignore);

    // Build a map from parameter name to parameter index
    bool buildParameterMap(Status& s=Status::ignore);

private:
    std::unordered_map<Id,TerminalIndex> terminalMap;
    std::unordered_map<Id,ParameterIndex> parameterMap;
};


class HierarchicalInstance : public Instance {
public:
    friend class HierarchicalModel;

    // The instance is added to the list of instances of the given model so that when
    // a model is destroyed all the instances in the list are destroyed, too. 
    HierarchicalInstance(HierarchicalModel* mod, Id name, Instance* parentInstance, const PTInstance& parsedInstance, Status& s=Status::ignore);
    virtual ~HierarchicalInstance();

    HierarchicalInstance           (const HierarchicalInstance&)  = delete;
    HierarchicalInstance           (      HierarchicalInstance&&) = default;
    HierarchicalInstance& operator=(const HierarchicalInstance&)  = delete;
    HierarchicalInstance& operator=(      HierarchicalInstance&&) = default;

    virtual ParameterIndex parameterCount() const { return model()->parameterCount(); };
    virtual std::tuple<ParameterIndex, bool> parameterIndex(Id name) const;
    virtual Id parameterName(ParameterIndex ndx) const;
    virtual std::tuple<Value::Type,bool> parameterType(ParameterIndex ndx, Status& s=Status::ignore) const;
    virtual bool getParameter(ParameterIndex ndx, Value& v, Status& s=Status::ignore) const;
    virtual std::tuple<bool,bool> setParameter(ParameterIndex ndx, const Value& v, Status& s=Status::ignore) ;
    // Principal parameter is the first parameter if there are any parameters
    virtual std::tuple<ParameterIndex, bool> principalParameterIndex() const { return std::make_tuple(0, parameterCount()>0);};
    virtual TerminalIndex staticNodeCount() const { return model()->terminalCount(); };
    virtual TerminalIndex terminalCount() const { return model()->terminalCount(); };
    virtual std::tuple<TerminalIndex, bool> nodeIndex(Id name) const { return model()->terminalIndex(name); };
    virtual Id nodeName(TerminalIndex ndx) const { return model()->terminalName(ndx); };
    virtual bool bindTerminal(TerminalIndex n, Node* node, Status& s=Status::ignore);
    virtual Node* terminal(TerminalIndex n, Status& s=Status::ignore) const;
    virtual bool unbindTerminals(Circuit& circuit, Status& s=Status::ignore);
    virtual std::vector<Instance*>* childInstances() { return &childInstances_; };
    virtual std::vector<Model*>* childModels() { return &childModels_; };
    virtual Id translateNode(Circuit& circuit, Id node);
    virtual void dump(int indent, const Circuit& circuit, std::ostream& os) const;

    // Model access (as HierarchicalModel)
    HierarchicalModel* model() { return static_cast<HierarchicalModel*>(model_); };
    const HierarchicalModel* model() const { return static_cast<HierarchicalModel*>(model_); };
    
    // Add child instances and models
    virtual void addChild(Instance* child);
    virtual void addChild(Model* child);
    
    // Context handling, parameter propagation, hierarchy building
    virtual std::tuple<bool, size_t> enterContext(Circuit& circuit, Context* externalContext, bool addToPath, bool rebuild, Status& s=Status::ignore);
    virtual bool propagateParameters(Circuit& circuit, RpnEvaluator& evaluator, Status& s=Status::ignore);
    virtual bool deleteHierarchy(Circuit& circuit, Status& s=Status::ignore);
    virtual bool buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s=Status::ignore);
    virtual std::tuple<bool, bool> subhierarchyChanged(Circuit& circuit, Status& s=Status::ignore); 

protected:
    bool buildBlock(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, const PTBlock& block, Status& s=Status::ignore);
    std::tuple<bool, bool> recomputeBlockConditions(Circuit& circuit, Status& s=Status::ignore); 
    
private:
    std::vector<Node*> connections;
    TerminalIndex connectedTerminalCount;
    Value* parameters;
    std::vector<Instance*> childInstances_;
    std::vector<Model*> childModels_;
};

}

#endif
