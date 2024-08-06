#ifndef __OSDIDEVICE_DEFINED
#define __OSDIDEVICE_DEFINED

#include <vector>
#include <unordered_map>
#include <type_traits>
#include "osdi.h"
#include "status.h"
#include "osdifile.h"
#include "value.h"
#include "devbase.h"
#include "common.h"


namespace NAMESPACE {

template<class T> T getDataPtr(void *objPtr, size_t offset) {
    return reinterpret_cast<T>( ((char*)objPtr)+offset);
}

class Circuit;

class OsdiModel;
class OsdiInstance;

class OsdiDevice : public Device {
public:
    friend class OsdiModel;

    enum class Error {
        OK,
        Range,
        NotFound, 
        NotModel, 
        NotInstance, 
        BadType
    };

    OsdiDevice(OsdiFile* of, int descriptorIndex, Id asName=Id::none, Loc location=Loc::bad, Status& s=Status::ignore);
    virtual ~OsdiDevice();

    OsdiDevice           (const OsdiDevice&)  = delete;
    OsdiDevice           (      OsdiDevice&&) = default;
    OsdiDevice& operator=(const OsdiDevice&)  = delete;
    OsdiDevice& operator=(      OsdiDevice&&) = default;

    virtual bool operator==(const Device& other) const; 
    Id name() const { return name_; };
    virtual std::tuple<bool, bool, bool> setup(Circuit& circuit, bool force, Status& s=Status::ignore);
    virtual bool collapseNodes(Circuit& circuit, Status& s=Status::ignore);
    virtual bool populateStructures(Circuit& circuit, Status& s=Status::ignore);
    virtual bool bind(
        Circuit& circuit, 
        KluMatrixAccess* matResist, Component compResist, 
        KluMatrixAccess* matReact, Component compReact, 
        Status& s=Status::ignore
    );
    virtual bool evalAndLoad(Circuit& circuit, EvalSetup* evalSetup, LoadSetup* loadSetup);
    virtual Model* createModel(Circuit& circuit, Instance* parentInstance, RpnEvaluator& evaluator, const PTModel& parsedModel, Status& s=Status::ignore);
    virtual void dump(int indent, std::ostream& os) const;

    // Osdi descriptor
    OsdiDescriptor* descriptor() { return descriptor_; };
    const OsdiDescriptor* descriptor() const { return descriptor_; };

    // Osdi file and device index
    const OsdiFile* file() const { return osdiFile; };
    OsdiFile::OsdiDeviceIndex deviceIndex() const { return index_; };
    
    // Parameter kind
    bool isInstanceParameter(OsdiFile::OsdiParameterId osdiId) const {
        return (descriptor_->param_opvar[osdiId].flags & PARA_KIND_MASK)==PARA_KIND_INST;
    };
    bool isModelParameter(OsdiFile::OsdiParameterId osdiId) const {
        // All parameters that are not opvars can be model parameters
        return (descriptor_->param_opvar[osdiId].flags & PARA_KIND_MASK)!=PARA_KIND_OPVAR;
    };
    bool isOpvar(OsdiFile::OsdiParameterId osdiId) const {
        return (descriptor_->param_opvar[osdiId].flags & PARA_KIND_MASK)==PARA_KIND_OPVAR;
    };

    // Parameter type (in form of Value::Type)
    Value::Type parameterType(OsdiFile::OsdiParameterId osdiId) const {
        return file()->parameterType(index_, osdiId);
    };

    // Translators
    // Id (name) -> simulator id of parameter/opvar
    std::tuple<ParameterIndex, bool> instanceParameterIndex(Id name) const {
        return osdiFile->instanceParameterIndex(index_, name);
    };
    std::tuple<ParameterIndex, bool> modelParameterIndex(Id name) const {
        return osdiFile->modelParameterIndex(index_, name);
    };
    std::tuple<ParameterIndex, bool> opvarIndex(Id name) const {
        return osdiFile->opvarIndex(index_, name);
    };
    // Number of parameters/opvars
    ParameterIndex instanceParameterCount() const { 
        return osdiFile->instanceParameterCount(index_); 
    };
    ParameterIndex modelParameterCount() const { 
        return osdiFile->modelParameterCount(index_); 
    };
    ParameterIndex opvarCount() const { 
        return osdiFile->opvarCount(index_); 
    };
    // Number of osdi parameters+opvars
    OsdiFile::OsdiParameterId osdiIdCount() const {
        return osdiFile->osdiIdCount(index_);
    };
    // Simulator id of parameter/opvar -> primary name
    Id instanceParameterName(ParameterIndex ndx) const {
        return osdiFile->instanceParameterName(index_, ndx);
    };
    Id modelParameterName(ParameterIndex ndx) const {
        return osdiFile->modelParameterName(index_, ndx);
    };
    Id opvarName(ParameterIndex ndx) const {
        return osdiFile->opvarName(index_, ndx);
    };
    // Osdi parameter id
    std::tuple<OsdiFile::OsdiParameterId, bool> osdiParameterId(Id name) const {
        return osdiFile->osdiParameterId(index_, name);
    }; 
    OsdiFile::OsdiParameterId instanceOsdiParameterId(ParameterIndex ndx) const {
        return osdiFile->instanceOsdiParameterId(index_, ndx);
    }; 
    OsdiFile::OsdiParameterId modelOsdiParameterId(ParameterIndex ndx) const {
        return osdiFile->modelOsdiParameterId(index_, ndx);
    }; 
    OsdiFile::OsdiParameterId opvarOsdiParameterId(ParameterIndex ndx) const {
        return osdiFile->opvarOsdiParameterId(index_, ndx);
    }; 
    // Osdi parameter id -> primary name
    Id parameterName(OsdiFile::OsdiParameterId osdiId) const {
        return osdiFile->opvarName(index_, osdiId);
    };

    // Terminal and node api
    TerminalIndex staticNodeCount() const { 
        return osdiFile->staticNodeCount(index_); 
    };
    TerminalIndex terminalCount() const { 
        return osdiFile->terminalCount(index_); 
    };
    std::tuple<TerminalIndex, bool> nodeIndex(Id name) const { 
        return osdiFile->nodeIndex(index_, name); 
    };
    inline Id nodeName(TerminalIndex ndx) const { 
        return osdiFile->nodeName(index_, ndx); 
    };

    // Number of noise sources
    inline ParameterIndex noiseSourceCount() const { return osdiFile->noiseSourceCount(index_); };

    // Noise source name
    inline Id noiseSourceName(ParameterIndex ndx) const { return osdiFile->noiseSourceName(index_, ndx); }; 

    // Noise source index
    inline std::tuple<ParameterIndex, bool> noiseSourceIndex(Id name) const { return osdiFile->noiseSourceIndex(index_, name); };

    // Access to nonzero entry indices
    auto& nonzeroResistiveResiduals() { return osdiFile->nonzeroResistiveResiduals(index_); };
    auto& nonzeroReactiveResiduals() { return osdiFile->nonzeroReactiveResiduals(index_); };
    auto& nonzeroResistiveJacobianEntries() { return osdiFile->nonzeroResistiveJacobianEntries(index_); };
    auto& nonzeroReactiveJacobianEntries() { return osdiFile->nonzeroReactiveJacobianEntries(index_); };

    // Noise source node indices
    inline std::tuple<OsdiFile::OsdiNodeIndex, OsdiFile::OsdiNodeIndex> noiseExcitation(ParameterIndex ndx) const {
        return osdiFile->noiseExcitation(index_, ndx);
    };

    // Parameter cleanup (called in instance and model destructor)
    bool freeValues(void* coreMod, void* coreInst);
    
    // Parameter access
    bool readParameter(OsdiFile::OsdiParameterId osdiId, void* coreMod, void* coreInst, Value& v, Status& s=Status::ignore) const;
    std::tuple<bool,bool> writeParameter(OsdiFile::OsdiParameterId osdiId, void* coreMod, void* coreInst, const Value& v, Status& s=Status::ignore);
    template<typename T> const T* parameterPtr(OsdiFile::OsdiParameterId osdiId, const void* coreMod, const void* coreInst) const;
    // Return value: ok, parameter given
    std::tuple<bool, bool> parameterGiven(OsdiFile::OsdiParameterId osdiId, void* coreMod, void* coreInst, Status& s=Status::ignore) const;

    static std::tuple<size_t, size_t> simParasSizes();
    static void populate(OsdiSimParas& sp, const SimulatorOptions& opt, const SimulatorInternals& internals, double* dblArray, char** chrPtrArray);
    
    bool processInitInfo(Circuit& circuit, OsdiInitInfo& initInfo, const char* typeString, Id name, Status& s=Status::ignore) const;

    OsdiFile::OsdiCollapsedNodesIndex collapsedNodesPatternSize() const { return descriptor_->num_collapsible; };
    OsdiFile::JacobianEntryIndex jacobianEntriesCount() const { return descriptor_->num_jacobian_entries; }; 
    const OsdiJacobianEntry& jacobianEntry(OsdiFile::JacobianEntryIndex ndx) const {
        return descriptor_->jacobian_entries[ndx];
    };
    OsdiFile::OsdiStateCount internalStateCount() const { return descriptor_->num_states; };
    
private:
    static const char* simParamNames[];
    static const char* simStrParamNames[];

    OsdiFile* osdiFile;
    OsdiFile::OsdiDeviceIndex index_;
    OsdiDescriptor* descriptor_;
};


// Template implementation

// Can only read parameter via pointer, writing is not allowed
template<typename T> const T* OsdiDevice::parameterPtr(OsdiFile::OsdiParameterId osdiId, const void* coreMod, const void* coreInst) const {
    // Check index
    if (osdiId>=osdiIdCount()) {
        return nullptr;
    }

    // Check kind
    if (!coreInst && !isModelParameter(osdiId)) {
        return nullptr;
    }

    if (coreInst && !(isInstanceParameter(osdiId) || isOpvar(osdiId))) {
        return nullptr;
    }

    auto t = parameterType(osdiId);
    // Check compatibility with requested pointer
    bool bad = false;
    if constexpr(std::is_same<T,int>::value) {
        bad = (t!=Value::Type::Int);
    } else if constexpr(std::is_same<T,double>::value) {
        bad = (t!=Value::Type::Real);
    } else {
        bad = true;
    }
    if (bad) {
        return nullptr;
    }

    OsdiFile::OsdiFlags flags = ACCESS_FLAG_READ;
    if (coreInst) {
        flags |= ACCESS_FLAG_INSTANCE;
    }
    return (T*)(descriptor_->access(const_cast<void*>(coreInst), const_cast<void*>(coreMod), osdiId, flags));
}

}

#endif
