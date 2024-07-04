#ifndef __DEVSOURCE_DEFINED
#define __DEVSOURCE_DEFINED

#include "exportdef.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <memory>
#include "osdifile.h"
#include "osdi.h"
#include "dynload.h"
#include "osdicallback.h"
#include "status.h"
#include "identifier.h"
#include "parseroutput.h"
#include "common.h"


namespace NAMESPACE {

class OsdiDevice;

class OsdiFile {
public:
    typedef uint32_t OsdiVersionType;
    typedef uint32_t OsdiDescriptorSize;
    typedef uint32_t OsdiDeviceIndex;
    typedef uint32_t OsdiLimitFunctionCount;
    typedef uint32_t OsdiParameterId;
    typedef uint32_t OsdiNoiseId;
    typedef uint32_t OsdiAliasIndex;
    typedef uint32_t OsdiNodeIndex;
    typedef uint32_t OsdiNodeCount;
    typedef uint32_t OsdiCollapsedNodesIndex;
    typedef uint32_t JacobianEntryIndex;
    typedef uint32_t OsdiFlags;
    typedef uint32_t OsdiErrorIndex;
    typedef uint32_t OsdiStateIndex;
    typedef uint32_t OsdiStateCount;

    static const JacobianEntryIndex noJacobianEntry = std::numeric_limits<JacobianEntryIndex>::max();

    OsdiFile(void*, std::string, Status& s=Status::ignore);

    OsdiFile           (const OsdiFile&)  = delete;
    OsdiFile           (      OsdiFile&&) = default;
    OsdiFile& operator=(const OsdiFile&)  = delete;
    OsdiFile& operator=(      OsdiFile&& from) = default;

    ~OsdiFile();

    bool isValid() const { return valid; }; 
    
    inline const std::string& fileName() const { return file; };
    
    inline OsdiDeviceIndex deviceCount() { return descriptorCount; };
    std::tuple<OsdiDeviceIndex,bool> deviceIndex(Id name, Status& s=Status::ignore);
    Id deviceIdentifier(OsdiDeviceIndex index=0, Status& s=Status::ignore);
    OsdiDescriptor *deviceDescriptor(OsdiDeviceIndex index, Status& s=Status::ignore);
    OsdiDescriptor *deviceDescriptor(Id name, Status& s=Status::ignore);

    // Factory function for creating a device source object from library file (API level 1)
    static OsdiFile* open(std::string file, const Loc& location, Status& s=Status::ignore);

    // Factory function for creating a device object from parsed netlist (API level 1)
    OsdiDevice *createDevice(OsdiDeviceIndex index, Id asName, Loc location=Loc::bad, Status& s=Status::ignore);
    OsdiDevice *createDevice(Id name, Id asName, Loc location=Loc::bad, Status& s=Status::ignore);

    // Translators
    // Id (name) -> simulator id of parameter/opvar
    inline std::tuple<ParameterIndex, bool> instanceParameterIndex(OsdiDeviceIndex deviceIndex, Id name) const {
        auto it = paramOsdiIdTranslators[deviceIndex].find(name);
        if (it==paramOsdiIdTranslators[deviceIndex].end()) {
            return std::make_tuple(0, false);
        }
        // Must be instance parameter
        if ((descriptors[deviceIndex]->param_opvar[it->second].flags & PARA_KIND_MASK) != PARA_KIND_INST) {
            return std::make_tuple(0, false);
        }
        return std::make_tuple(osdiIdSimInstIdLists[deviceIndex][it->second], true);
    };
    inline std::tuple<ParameterIndex, bool> modelParameterIndex(OsdiDeviceIndex deviceIndex, Id name) const {
        auto it = paramOsdiIdTranslators[deviceIndex].find(name);
        if (it==paramOsdiIdTranslators[deviceIndex].end()) {
            return std::make_tuple(0, false);
        }
        // Must be instance or model parameter (both are valid model parameters)
        auto tmp = descriptors[deviceIndex]->param_opvar[it->second].flags & PARA_KIND_MASK;
        if (tmp != PARA_KIND_INST && tmp != PARA_KIND_MODEL) {
            return std::make_tuple(0, false);
        }
        return std::make_tuple(osdiIdSimModIdLists[deviceIndex][it->second], true);
    };
    inline std::tuple<ParameterIndex, bool> opvarIndex(OsdiDeviceIndex deviceIndex, Id name) const {
        auto it = paramOsdiIdTranslators[deviceIndex].find(name);
        if (it==paramOsdiIdTranslators[deviceIndex].end()) {
            return std::make_tuple(0, false);
        }
        // Must be opvar parameter
        if ((descriptors[deviceIndex]->param_opvar[it->second].flags & PARA_KIND_MASK) != PARA_KIND_OPVAR) {
            return std::make_tuple(0, false);
        }
        return std::make_tuple(osdiIdSimInstIdLists[deviceIndex][it->second], true);
    };
    // Number of parameters/opvars
    inline ParameterIndex instanceParameterCount(OsdiDeviceIndex deviceIndex) const { 
        return instanceParamOsdiIdLists[deviceIndex].size(); 
    };
    inline ParameterIndex modelParameterCount(OsdiDeviceIndex deviceIndex) const { 
        return modelParamOsdiIdLists[deviceIndex].size(); 
    };
    inline ParameterIndex opvarCount(OsdiDeviceIndex deviceIndex) const { 
        return opvarOsdiIdLists[deviceIndex].size(); 
    };
    // Number of osdi parameters+opvars
    inline OsdiParameterId osdiIdCount(OsdiDeviceIndex deviceIndex) const {
        return descriptors[deviceIndex]->num_params+descriptors[deviceIndex]->num_opvars;
    };
    // Simulator id of parameter/opvar -> primary name
    inline Id instanceParameterName(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const {
        auto osdiId = instanceParamOsdiIdLists[deviceIndex][ndx];
        return osdiIdPrimaryParamName[deviceIndex][osdiId];
    };
    inline Id modelParameterName(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const {
        auto osdiId = modelParamOsdiIdLists[deviceIndex][ndx];
        return osdiIdPrimaryParamName[deviceIndex][osdiId];
    };
    inline Id opvarName(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const {
        auto osdiId = opvarOsdiIdLists[deviceIndex][ndx];
        return osdiIdPrimaryParamName[deviceIndex][osdiId];
    };
    // Osdi parameter id
    inline std::tuple<OsdiParameterId, bool> osdiParameterId(OsdiDeviceIndex deviceIndex, Id name) const {
        auto it = paramOsdiIdTranslators[deviceIndex].find(name);
        if (it==paramOsdiIdTranslators[deviceIndex].end()) {
            return std::make_tuple(0, false);
        }
        return std::make_tuple(it->second, true);
    }; 
    inline OsdiParameterId instanceOsdiParameterId(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const {
        return instanceParamOsdiIdLists[deviceIndex][ndx];
    }; 
    inline OsdiParameterId modelOsdiParameterId(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const {
        return modelParamOsdiIdLists[deviceIndex][ndx];
    }; 
    inline OsdiParameterId opvarOsdiParameterId(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const {
        return opvarOsdiIdLists[deviceIndex][ndx];
    }; 
    // Osdi parameter id -> primary name
    inline Id parameterName(OsdiDeviceIndex deviceIndex, OsdiParameterId osdiId) const {
        return osdiIdPrimaryParamName[deviceIndex][osdiId];
    };
    // Osdi parameter id -> Value type
    inline Value::Type parameterType(OsdiDeviceIndex deviceIndex, OsdiParameterId osdiId) const {
        auto& param = descriptors[deviceIndex]->param_opvar[osdiId];
        switch (param.flags & PARA_TY_MASK) {
            case PARA_TY_INT:
                if (param.len>0) {
                    return Value::Type::IntVec;
                } else {
                    return Value::Type::Int;
                }
                break;
            case PARA_TY_REAL:
                if (param.len>0) {
                    return Value::Type::RealVec;
                } else {
                    return Value::Type::Real;
                }
                break;
            default: // PARA_TY_STR
                if (param.len>0) {
                    return Value::Type::StringVec;
                } else {
                    return Value::Type::String;
                }
                break;
        }
    };

    // List of osdi ids of instance parameters that may need to be freed 
    const std::vector<OsdiParameterId>& allocatedInstanceParameterIds(OsdiDeviceIndex deviceIndex) const {
        return instParAllocatedOsdiId[deviceIndex];
    };
    // List of osdi ids of model parameters that may need to be freed 
    const std::vector<OsdiParameterId>& allocatedModelParamemeterIds(OsdiDeviceIndex deviceIndex) const {
        return modParAllocatedOsdiId[deviceIndex];
    };

    // Number of nodes
    TerminalIndex staticNodeCount(OsdiDeviceIndex deviceIndex) const { return descriptors[deviceIndex]->num_nodes; };

    // Number of terminals
    TerminalIndex terminalCount(OsdiDeviceIndex deviceIndex) const { return descriptors[deviceIndex]->num_terminals; };

    // Node index
    std::tuple<TerminalIndex, bool> nodeIndex(OsdiDeviceIndex deviceIndex, Id name) const { 
        auto it = nodeMaps[deviceIndex].find(name);
        if (it!=nodeMaps[deviceIndex].end()) {
            return std::make_tuple(it->second, true);
        } else {
            return std::make_tuple(0, false);
        }
    };

    // Node name
    inline Id nodeName(OsdiDeviceIndex deviceIndex, TerminalIndex ndx) const { return nodeNameLists[deviceIndex][ndx]; };

    // Number of noise sources
    inline ParameterIndex noiseSourceCount(OsdiDeviceIndex deviceIndex) const { return descriptors[deviceIndex]->num_noise_src; };

    // Noise source name
    inline Id noiseSourceName(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const { return noiseSourceNames[deviceIndex][ndx]; }; 

    // Noise source index
    inline std::tuple<ParameterIndex, bool> noiseSourceIndex(OsdiDeviceIndex deviceIndex, Id name) const {
        auto it = noiseSourceNameTranslators[deviceIndex].find(name);
        if (it==noiseSourceNameTranslators[deviceIndex].end()) {
            return std::make_tuple(0, false);
        }
        return std::make_tuple(it->second, true);
    };

    // Noise source nodes
    inline std::tuple<OsdiNodeIndex, OsdiNodeIndex> noiseExcitation(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const {
        return std::make_tuple(
            descriptors[deviceIndex]->noise_sources[ndx].nodes.node_1, 
            descriptors[deviceIndex]->noise_sources[ndx].nodes.node_2
        );
    };

    // Instance node state count
    LocalStorageIndex nodeStateCount(OsdiDeviceIndex deviceIndex) { return instanceNodeStateCounts[deviceIndex]; };

private:
    void* handle;
    std::string file;
    bool valid;
    void* descriptorArray;
    std::vector<OsdiDescriptor*> descriptors;
    OsdiDeviceIndex descriptorCount;
    size_t descriptorSize;
    
    // Map from device name to device index
    std::unordered_map<Id,OsdiDeviceIndex> deviceNameToIndex;

    // Map from dynamic library handle to OsdiFile structure
    static std::unordered_map<void*,std::unique_ptr<OsdiFile>> registry;

    // Vector of vectors of Osdi parameter names
    std::vector<std::vector<Id>> osdiIdPrimaryParamName;

    // Vector of maps from parameter name to osdi parameter index
    std::vector<std::unordered_map<Id,OsdiParameterId>>paramOsdiIdTranslators;

    // Vector of vectors of model parameter ids
    std::vector<std::vector<OsdiParameterId>> modelParamOsdiIdLists;

    // Vector of vectors of instance parameter ids
    std::vector<std::vector<OsdiParameterId>> instanceParamOsdiIdLists;
    
    // Vector of vectors of opvar ids
    std::vector<std::vector<OsdiParameterId>> opvarOsdiIdLists;

    // Vector of vectors of simulator instance parameter/opvar ids corresponding to osdi ids
    std::vector<std::vector<ParameterIndex>> osdiIdSimInstIdLists;

    // Vector of vectors of simulator model parameter ids corresponding to osdi ids
    std::vector<std::vector<ParameterIndex>> osdiIdSimModIdLists;

    // Vector of vectors of noise source names
    std::vector<std::vector<Id>> noiseSourceNames;

    // Vector of maps from noise source name to noise source index
    std::vector<std::unordered_map<Id, ParameterIndex>> noiseSourceNameTranslators;

    // Vector of instance node state counts
    std::vector<LocalStorageIndex> instanceNodeStateCounts;

    // Vector of vectors of node identifiers
    std::vector<std::vector<Id>> nodeNameLists;

    // Vector of maps from node name to node index
    std::vector<std::unordered_map<Id, OsdiNodeIndex>> nodeMaps;

    // We need these to quickly free all parameters that are allocated (strings, vectors)
    // These are lists of osdi ids of parameters that need to be freed if they were given
    std::vector<std::vector<OsdiParameterId>> instParAllocatedOsdiId;
    std::vector<std::vector<OsdiParameterId>> modParAllocatedOsdiId;

    // Limit functions
    static const OsdiLimitFunction limitFunctionTable[];
};

}

#endif
