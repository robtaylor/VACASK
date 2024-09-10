#ifndef __OSDIINSTANCE_DEFINED
#define __OSDIINSTANCE_DEFINED

#include <string>
#include "status.h"
#include "osdimodel.h"
#include "devbase.h"
#include "common.h"


namespace NAMESPACE {

class Circuit;

class OsdiInstance : public Instance {
public:
    friend class OsdiDevice;

    // The instance is added to the list of instances of the given model so that when
    // a model is destroyed all the instances in the list are destroyed, too. 
    OsdiInstance(OsdiModel* mod, Id name, Instance* parentInstance, const PTInstance& parsedInstance, Status& s=Status::ignore);
    virtual ~OsdiInstance();

    OsdiInstance           (const OsdiInstance&)  = delete;
    OsdiInstance           (      OsdiInstance&&) = default;
    OsdiInstance& operator=(const OsdiInstance&)  = delete;
    OsdiInstance& operator=(      OsdiInstance&&) = default;

    virtual ParameterIndex parameterCount() const;
    virtual std::tuple<ParameterIndex, bool> parameterIndex(Id name) const;
    virtual Id parameterName(ParameterIndex ndx) const;
    virtual std::tuple<Value::Type,bool> parameterType(ParameterIndex ndx, Status& s=Status::ignore) const;
    virtual bool getParameter(ParameterIndex ndx, Value& v, Status& s=Status::ignore) const;
    virtual std::tuple<bool,bool> setParameter(ParameterIndex ndx, const Value& v, Status& s=Status::ignore);
    virtual std::tuple<bool,bool> parameterGiven(ParameterIndex ndx, Status& s=Status::ignore) const;
    
    // First parameter is $mfactor, second is principal parameter (if it exists)
    virtual std::tuple<ParameterIndex, bool> principalParameterIndex() const { return std::make_tuple(1, parameterCount()>1); };
    virtual TerminalIndex staticNodeCount() const;
    virtual TerminalIndex terminalCount() const;
    virtual std::tuple<TerminalIndex, bool> nodeIndex(Id name) const;
    virtual Id nodeName(TerminalIndex ndx) const;
    virtual bool bindTerminal(TerminalIndex n, Node* node, Status& s=Status::ignore);
    virtual Node* terminal(TerminalIndex n, Status& s=Status::ignore) const;
    virtual bool unbindTerminals(Circuit& circuit, Status& s=Status::ignore);
    virtual bool propagateParameters(Circuit& circuit, RpnEvaluator& evaluator, Status& s=Status::ignore);
    virtual bool deleteHierarchy(Circuit& circuit, Status& s=Status::ignore);
    virtual bool buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, Status& s=Status::ignore);
    virtual std::tuple<EquationIndex,EquationIndex> sourceExcitation(Circuit& circuit) const;
    virtual std::tuple<UnknownIndex,UnknownIndex> sourceResponse(Circuit& circuit) const;
    virtual ParameterIndex opvarCount() const;
    virtual std::tuple<ParameterIndex, bool> opvarIndex(Id name) const;
    virtual Id opvarName(ParameterIndex ndx) const;
    virtual std::tuple<Value::Type,bool> opvarType(ParameterIndex ndx, Status& s=Status::ignore) const;
    virtual bool getOpvar(ParameterIndex ndx, Value& v, Status& s=Status::ignore) const; 
    virtual std::tuple<bool, OutputSource> opvarOutputSource(ParameterIndex ndx) const;
    virtual std::tuple<bool, bool, bool> setup(Circuit& circuit, bool force, DeviceRequests* devReq, Status& s=Status::ignore);
    virtual void dump(int indent, const Circuit& circuit, std::ostream& os) const;

    // Model access (as OsdiModel)
    OsdiModel* model() { return static_cast<OsdiModel*>(model_); };
    const OsdiModel* model() const { return static_cast<const OsdiModel*>(model_); };

    // Model core access
    void* core() { return core_; };
    const void* core() const { return core_; };

    // Noise API
    virtual ParameterIndex noiseSourceCount() const { return model()->device()->noiseSourceCount(); };
    virtual Id noiseSourceName(ParameterIndex ndx) const { return model()->device()->noiseSourceName(ndx); };
    virtual std::tuple<ParameterIndex, bool> noiseSourceIndex(Id name) const { return model()->device()->noiseSourceIndex(name); }
    virtual std::tuple<EquationIndex, EquationIndex> noiseExcitation(Circuit& cir, ParameterIndex ndx) const;
    virtual bool loadNoise(Circuit& circuit, double freq, double* noiseDensity);

    // Helpers for inlining in device, model, and instance virtual functions
    std::tuple<bool, bool, bool> setupCore(Circuit& circuit, OsdiSimParas& sp, double temp, bool force, DeviceRequests* devReq, Status& s=Status::ignore);
    bool collapseNodesCore(Circuit& circuit, Status& s);
    bool populateStructuresCore(Circuit& circuit, Status& s=Status::ignore);
    bool bindCore(
        Circuit& circuit, 
        KluMatrixAccess* matResist, Component compResist, const std::optional<MatrixEntryPosition>& mepResist, 
        KluMatrixAccess* matReact, Component compReact, const std::optional<MatrixEntryPosition>& mepReact, 
        Status& s=Status::ignore
    );
    bool bypassCheckCore(Circuit& circuit, EvalSetup& evalSetup);
    bool evalCore(Circuit& circuit, OsdiSimInfo& simInfo, EvalSetup& evalSetup);
    bool loadCore(Circuit& circuit, LoadSetup& loadSetup);
    bool convergedCore(Circuit& circuit, ConvSetup& convSetup);
    
protected:
    OsdiFile::OsdiCollapsedNodesIndex collapsedNodesPatternSize() const { return model()->device()->collapsedNodesPatternSize(); };
    bool* collapsedNodesPattern() const { return getDataPtr<bool*>(core_, model()->device()->descriptor()->collapsed_offset); };
    OsdiFile::OsdiNodeIndex* nodeMappingArray() { 
        return getDataPtr<OsdiFile::OsdiNodeIndex*>(core_, model()->device()->descriptor()->node_mapping_offset);
    };
    double** resistiveJacobianPointers() { return getDataPtr<double**>(core_, model()->device()->descriptor()->jacobian_ptr_resist_offset); }; 
    double** reactiveJacobianPointer(OsdiFile::JacobianEntryIndex ndx) { 
        auto& entry = model()->device()->descriptor()->jacobian_entries[ndx];
        if (entry.react_ptr_off==OsdiFile::noJacobianEntry) {
            return nullptr;
        } else {
            return getDataPtr<double**>(core_, entry.react_ptr_off);
        }
    };
    // State index table of instance
    // Contains uint32_t offsets based on a pointer passed to eval()
    // in simInfo.prev_state and simInfo.next_state arrays of doubles
    OsdiFile::OsdiStateIndex* stateIndexTable() {
        return getDataPtr<OsdiFile::OsdiStateIndex*>(core_, model()->device()->descriptor()->state_idx_off);
    };
    
private:
    void* core_;
    std::vector<Node*> nodes_;
    TerminalIndex connectedTerminalCount;
    GlobalStorageIndex offsStates;
    GlobalStorageIndex offsDeviceStates;
};

}

#endif
