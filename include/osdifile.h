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
#include "natures.h"
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
    typedef uint32_t OsdiNatureIndex;
    typedef uint32_t OsdiDisciplineIndex;

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

    // Get the vector of loaded files
    static const std::vector<OsdiFile*>& files() { return fileOrder; };

    // Translators
    // Id (name) -> simulator id of parameter/output variable
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
    inline std::tuple<ParameterIndex, bool> outvarIndex(OsdiDeviceIndex deviceIndex, Id name) const {
        auto it = paramOsdiIdTranslators[deviceIndex].find(name);
        if (it==paramOsdiIdTranslators[deviceIndex].end()) {
            return std::make_tuple(0, false);
        }
        // Must be an output variable
        if ((descriptors[deviceIndex]->param_opvar[it->second].flags & PARA_KIND_MASK) != PARA_KIND_OPVAR) {
            return std::make_tuple(0, false);
        }
        return std::make_tuple(osdiIdSimInstIdLists[deviceIndex][it->second], true);
    };
    // Number of parameters/output variables
    inline ParameterIndex instanceParameterCount(OsdiDeviceIndex deviceIndex) const { 
        return instanceParamOsdiIdLists[deviceIndex].size(); 
    };
    inline ParameterIndex modelParameterCount(OsdiDeviceIndex deviceIndex) const { 
        return modelParamOsdiIdLists[deviceIndex].size(); 
    };
    inline ParameterIndex outvarCount(OsdiDeviceIndex deviceIndex) const { 
        return outvarOsdiIdLists[deviceIndex].size(); 
    };
    // Number of osdi parameters+output variables
    inline OsdiParameterId osdiIdCount(OsdiDeviceIndex deviceIndex) const {
        return descriptors[deviceIndex]->num_params+descriptors[deviceIndex]->num_opvars;
    };
    // Simulator id of parameter/output variable -> primary name
    inline Id instanceParameterName(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const {
        auto osdiId = instanceParamOsdiIdLists[deviceIndex][ndx];
        return osdiIdPrimaryParamName[deviceIndex][osdiId];
    };
    inline Id modelParameterName(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const {
        auto osdiId = modelParamOsdiIdLists[deviceIndex][ndx];
        return osdiIdPrimaryParamName[deviceIndex][osdiId];
    };
    inline Id outvarName(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const {
        auto osdiId = outvarOsdiIdLists[deviceIndex][ndx];
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
    inline OsdiParameterId outvarOsdiParameterId(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const {
        return outvarOsdiIdLists[deviceIndex][ndx];
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

    // Number of unique noise sources
    inline ParameterIndex uniqueNoiseSourceCount(OsdiDeviceIndex deviceIndex) const { return noiseSourceNameTranslators[deviceIndex].size(); };

    // Noise source name
    inline Id noiseSourceName(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const { return noiseSourceNames[deviceIndex][ndx]; }; 

    // Unique noise source index from noise source index
    inline size_t uniqueNoiseSourceIndex(OsdiDeviceIndex deviceIndex, ParameterIndex ndx) const { return uniqueNoiseSourceIndices[deviceIndex][ndx]; }; 

    // Uniquie noise source index from noise source name
    inline std::tuple<ParameterIndex, bool> uniqueNoiseSourceIndex(OsdiDeviceIndex deviceIndex, Id name) const {
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

    // Does the device allow bypass
    bool allowsBypass(OsdiDeviceIndex deviceIndex) const { return allowsBypass_[deviceIndex]; };

    // Access to nonzero entry indices
    auto& nonzeroResistiveResiduals(OsdiDeviceIndex deviceIndex) { return nonzeroResistiveResNdx[deviceIndex]; };
    auto& nonzeroReactiveResiduals(OsdiDeviceIndex deviceIndex) { return nonzeroReactiveResNdx[deviceIndex]; };
    auto& nonzeroResistiveJacobianEntries(OsdiDeviceIndex deviceIndex) { return nonzeroResistiveJacNdx[deviceIndex]; };
    auto& nonzeroReactiveJacobianEntries(OsdiDeviceIndex deviceIndex) { return nonzeroReactiveJacNdx[deviceIndex]; };

    // Access to disciplines and natures
    const OsdiNature* nature(OsdiNatureIndex index) const { if (natures && index!=UINT32_MAX && index<naturesCount) return natures+index; else return nullptr; };
    const OsdiDiscipline* discipline(OsdiDisciplineIndex index) const { if (disciplines && index!=UINT32_MAX && index<disciplinesCount) return disciplines+index; else return nullptr; };

    // OsdiNatureRef to NatureTolerance and idt NatureTolerance
    std::tuple<NatureTolerance, NatureTolerance> natrefTolerances(OsdiNatureRef& natref) const {
        switch (natref.ref_type) {
            case NATREF_NATURE:
                return std::make_tuple(
                    NatureTolerance(
                        natureId[natref.index], 
                        natureAbstol[natref.index]
                    ), 
                    NatureTolerance(
                        natureIdtId[natref.index], 
                        natureIdtAbstol[natref.index]
                    )
                );
            case NATREF_DISCIPLINE_FLOW:
                return std::make_tuple(
                    NatureTolerance(
                        disciplineFlowId[natref.index], 
                        disciplineFlowAbstol[natref.index]
                    ), 
                    NatureTolerance(
                        disciplineFlowIdtId[natref.index], 
                        disciplineFlowIdtAbstol[natref.index]
                    )
                );
            case NATREF_DISCIPLINE_POTENTIAL:
                return std::make_tuple(
                    NatureTolerance(
                        disciplinePotentialId[natref.index], 
                        disciplinePotentialAbstol[natref.index]
                    ), 
                    NatureTolerance(
                        disciplinePotentialIdtId[natref.index], 
                        disciplinePotentialIdtAbstol[natref.index]
                    )
                );
        }
        return std::make_tuple(
            NatureTolerance(
                NatureRegistry::noNature, 
                std::numeric_limits<double>::infinity()
            ),
            NatureTolerance(
                NatureRegistry::noNature, 
                std::numeric_limits<double>::infinity()
            )
        );
    }

    // Absolute tolerances as NatureTolerance (unknown, unknown idt, residual, residual idt)
    std::tuple<NatureTolerance, NatureTolerance, NatureTolerance, NatureTolerance> tolerances(OsdiDeviceIndex deviceIndex, TerminalIndex node) const {
        auto desc = descriptors[deviceIndex];
        auto& un = desc->unknown_nature[node];
        auto& rn = desc->residual_nature[node];
        return std::tuple_cat(natrefTolerances(un), natrefTolerances(rn));
    }
    
    // Dump
    virtual void dump(int indent, std::ostream& os) const;

private:
    void* handle;
    std::string file;
    bool valid;
    void* descriptorArray;
    std::vector<std::string> namesArray; // (translated) names
    std::vector<OsdiDescriptor*> descriptors;
    OsdiDeviceIndex descriptorCount;
    size_t descriptorSize;
    OsdiNature* natures;
    OsdiDiscipline* disciplines;
    OsdiAttribute* attributes;
    std::vector<NatureId> natureId;
    std::vector<NatureId> disciplineFlowId;
    std::vector<NatureId> disciplinePotentialId;
    OsdiNatureIndex naturesCount;
    OsdiDisciplineIndex disciplinesCount;

    // Map from device name to device index
    std::unordered_map<Id,OsdiDeviceIndex> deviceNameToIndex;

    // Map from dynamic library handle to OsdiFile structure
    // This structure owns the OsdiFile object via its pointer
    static std::unordered_map<void*,std::unique_ptr<OsdiFile>> registry;

    // Map from path to OsdiFile structure
    static std::unordered_map<std::string,OsdiFile*> nameToFile;

    // Order in which files were opened
    static std::vector<OsdiFile*> fileOrder;

    // Vector of vectors of Osdi parameter names
    std::vector<std::vector<Id>> osdiIdPrimaryParamName;

    // Vector of maps from parameter name to osdi parameter index
    std::vector<std::unordered_map<Id,OsdiParameterId>>paramOsdiIdTranslators;

    // Vector of vectors of model parameter ids
    std::vector<std::vector<OsdiParameterId>> modelParamOsdiIdLists;

    // Vector of vectors of instance parameter ids
    std::vector<std::vector<OsdiParameterId>> instanceParamOsdiIdLists;
    
    // Vector of vectors of output variable ids
    std::vector<std::vector<OsdiParameterId>> outvarOsdiIdLists;

    // Vector of vectors of simulator instance parameter/output variable ids corresponding to osdi ids
    std::vector<std::vector<ParameterIndex>> osdiIdSimInstIdLists;

    // Vector of vectors of simulator model parameter ids corresponding to osdi ids
    std::vector<std::vector<ParameterIndex>> osdiIdSimModIdLists;

    // Vector of vectors of noise source names
    // Same name can appear multiple times in a single device
    std::vector<std::vector<Id>> noiseSourceNames;
    
    // Vector of vectors of unique noise source indices, one for each noise source
    std::vector<std::vector<size_t>> uniqueNoiseSourceIndices;

    // Vector of maps from noise source name to unique noise source index
    std::vector<std::unordered_map<Id, ParameterIndex>> noiseSourceNameTranslators;

    // Vector of vectors of node identifiers
    std::vector<std::vector<Id>> nodeNameLists;

    // Vector of maps from node name to node index
    std::vector<std::unordered_map<Id, OsdiNodeIndex>> nodeMaps;

    // We need these to quickly free all parameters that are allocated (strings, vectors)
    // These are lists of osdi ids of parameters that need to be freed if they were given
    std::vector<std::vector<OsdiParameterId>> instParAllocatedOsdiId;
    std::vector<std::vector<OsdiParameterId>> modParAllocatedOsdiId;

    // These lists are used for allocating
    // Vector of vectors of nonzero resistive Jacobian entry indices into jacobian_entries
    std::vector<std::vector<OsdiNodeIndex>> nonzeroResistiveJacNdx;

    // Vector of vectors of nonzero reactive Jacobian entry indices into jacobian_entries
    std::vector<std::vector<OsdiNodeIndex>> nonzeroReactiveJacNdx;

    // Vector of vectors of nonzero resistive residual entry indices into nodes
    std::vector<std::vector<OsdiNodeIndex>> nonzeroResistiveResNdx;

    // Vector of vectors of nonzero resistive residual entry indices into nodes
    std::vector<std::vector<OsdiNodeIndex>> nonzeroReactiveResNdx;

    // Vector of allows bypass flags
    std::vector<bool> allowsBypass_;

    // Limit functions
    static const OsdiLimitFunction limitFunctionTable[];

    // Process natures and disciplines 
    // - collect abstol values for each nature and its idt nature
    void processNaturesAndDisciplines();

    std::vector<double> natureAbstol;
    std::vector<NatureId> natureIdtId;
    std::vector<double> natureIdtAbstol;

    std::vector<double> disciplineFlowAbstol;
    std::vector<NatureId> disciplineFlowIdtId;
    std::vector<double> disciplineFlowIdtAbstol;
    std::vector<double> disciplinePotentialAbstol;
    std::vector<NatureId> disciplinePotentialIdtId;
    std::vector<double> disciplinePotentialIdtAbstol;
};

}

#endif
