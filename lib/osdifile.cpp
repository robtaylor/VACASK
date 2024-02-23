#include <string.h>
#include "dynload.h"
#include "osdicallback.h"
#include "limitfunctions.h"
#include "osdi.h"
#include "osdifile.h"
#include "osdidevice.h"
#include "identifier.h"
#include "devbase.h"
#include "simulator.h"
#include "common.h"


namespace NAMESPACE {

extern std::ostream *osdiOutStreamPtr;
extern std::ostream *osdiErrStreamPtr;

extern const OsdiLimitFunction* osdiLimitFunctionTable;

std::unordered_map<void*,std::unique_ptr<OsdiFile>> OsdiFile::registry;

static void toLowercase(std::string& s) { 
    for(auto it=s.begin(); it!=s.end(); ++it) {
        *it=std::tolower(*it); 
    }
}

OsdiFile::OsdiFile(void* handle_, std::string file_, Status& s) 
    : handle(handle_), file(file_), valid(false) {
    // Descriptors table
    descriptors = (OsdiDescriptor*)dynamicLibrarySymbol(handle, "OSDI_DESCRIPTORS");
    
    // Descriptor count
    auto ptr = (OsdiDeviceIndex*)dynamicLibrarySymbol(handle, "OSDI_NUM_DESCRIPTORS");
    if (ptr) {
        descriptorCount = *ptr;
    } else  {
        descriptorCount = 0;
        descriptors = nullptr;
    }
    if (!descriptors)  {
        s.set(Status::Empty, "OSDI file contains no models");
        return;
    }

    // Version
    auto major = ((OsdiVersionType*)dynamicLibrarySymbol(handle, "OSDI_VERSION_MAJOR"));
    auto minor = ((OsdiVersionType*)dynamicLibrarySymbol(handle, "OSDI_VERSION_MINOR"));
    if (!major || !minor) {
        s.set(Status::BadVersion, "Failed to retrieve OSDI interface version.");
        descriptors = nullptr;
        return;
    } else {
        if (*major!=0 || *minor!=3) {
            s.set(Status::BadVersion, "Wrong OSDI interface version. Only v0.3 is supported.");
            descriptors = nullptr;
            return;
        }
    }

    // Limit function table
    OsdiLimFunction* lft = (OsdiLimFunction*)dynamicLibrarySymbol(handle, "OSDI_LIM_TABLE");
    auto nlft = (OsdiLimitFunctionCount*)dynamicLibrarySymbol(handle, "OSDI_LIM_TABLE_LEN");
    if (lft && nlft && *nlft>0) {
        auto n = *nlft;
        decltype(n) nArgs, j;
        for(decltype(n) i=0; i<n; i++) {
            nArgs = -1;
            for(j=0; OsdiFile::limitFunctionTable[j].name; j++) {
                if (!strcmp(lft[i].name, OsdiFile::limitFunctionTable[j].name)) {
                    nArgs = OsdiFile::limitFunctionTable[j].nArgs;
                    if (nArgs == lft[i].num_args) {
                        lft[i].func_ptr = OsdiFile::limitFunctionTable[j].ptr;
                        break;
                    }
                }
            }
            if (nArgs==-1) {
                s.set(
                    Status::Unknown, 
                    std::string("Unknown limit function '")+lft[i].name+"'."
                );
                descriptors = nullptr;
            } else if (!OsdiFile::limitFunctionTable[j].ptr) {
                s.set(
                    Status::BadArguments,                     
                    std::string("Unexpected number of arguments (")+
                        std::to_string(nArgs)+", expected "+std::to_string(lft[i].num_args)+
                        ") for limit function '"+lft[i].name+"'."
                    
                );
                descriptors = nullptr;
            }
        }
    }

    // Stop on error
    if (!descriptors) {
        return;
    }
    
    // Logger
    void** funcPtr = ((void**)dynamicLibrarySymbol(handle, "osdi_log"));
    if (funcPtr) {
        *funcPtr = (void *)osdiLogMessage;
    }

    // Device index
    for(OsdiDeviceIndex i=0; i<descriptorCount; i++) {
        deviceNameToIndex[Id(descriptors[i].name)] = i;
    }

    // Build parameter name to id translators and lists of instance parameter ids
    paramOsdiIdTranslators.resize(descriptorCount);
    instanceParamOsdiIdLists.resize(descriptorCount);
    modelParamOsdiIdLists.resize(descriptorCount);
    opvarOsdiIdLists.resize(descriptorCount); 
    osdiIdSimInstIdLists.resize(descriptorCount);
    osdiIdSimModIdLists.resize(descriptorCount);
    instanceNodeStateCounts.resize(descriptorCount);
    instParGivenIndex.resize(descriptorCount);
    modParGivenIndex.resize(descriptorCount);
    osdiIdPrimaryParamName.resize(descriptorCount);
    noiseSourceNames.resize(descriptorCount);
    noiseSourceNameTranslators.resize(descriptorCount);
    nodeNameLists.resize(descriptorCount);
    nodeMaps.resize(descriptorCount);
    for(int i=0; i<descriptorCount; i++) {
        OsdiDescriptor* desc = descriptors+i;
        auto& translator = paramOsdiIdTranslators[i];
        osdiIdSimInstIdLists[i].resize(desc->num_params+desc->num_opvars);
        osdiIdSimModIdLists[i].resize(desc->num_params+desc->num_opvars);
        // osdiIdPrimaryParamName[i].resize(desc->num_params+desc->num_opvars);
        OsdiParameterId instParGivenCount = 0;
        OsdiParameterId modParGivenCount = 0;
        for(OsdiParameterId j=0; j<desc->num_params+desc->num_opvars; j++) {
            auto& paramInfo = desc->param_opvar[j];
            // Add lowercase name to list of primary parameter names
            std::string tmp = descriptors[i].param_opvar[j].name[0];
            toLowercase(tmp);
            osdiIdPrimaryParamName[i].push_back(tmp);
            // Add name and aliases to id translator
            for(int k=0; k<paramInfo.num_alias+1; k++) {
                std::string tmp = paramInfo.name[k];
                // Convert parameter name to lowercase
                toLowercase(tmp);
                translator[Id(tmp)] = j;
                // if (tmp=="type") {
                //     int a=1;
                // }
            }
            // Add id to instance/model/opvar osdi ID lists
            if ((paramInfo.flags & PARA_KIND_MASK) == PARA_KIND_INST) {
                // Instance and model parameter
                instanceParamOsdiIdLists[i].push_back(j);
                modelParamOsdiIdLists[i].push_back(j);
                // Add to osdi id to simulator id translator list
                osdiIdSimInstIdLists[i][j] = instanceParamOsdiIdLists[i].size()-1;
                osdiIdSimModIdLists[i][j] = modelParamOsdiIdLists[i].size()-1;
                // Param given flag index, only for strings and vectors
                if ((paramInfo.flags & PARA_TY_MASK)==PARA_TY_STR || paramInfo.len>0) {
                    instParGivenIndex[i].insert({j, instParGivenCount});
                    modParGivenIndex[i].insert({j, modParGivenCount});
                    instParGivenCount++;
                    modParGivenCount++;
                }
            } else if ((paramInfo.flags & PARA_KIND_MASK) == PARA_KIND_OPVAR) {
                // Opvar
                opvarOsdiIdLists[i].push_back(j);
                // Add to osdi id to simulator id translator list
                osdiIdSimInstIdLists[i][j] = opvarOsdiIdLists[i].size()-1;
            } else {
                // Model parameter
                modelParamOsdiIdLists[i].push_back(j);
                // Add to osdi id to simulator id translator list
                osdiIdSimModIdLists[i][j] = modelParamOsdiIdLists[i].size()-1;
                // Param given flag index, only for strings and vectors
                if ((paramInfo.flags & PARA_TY_MASK)==PARA_TY_STR || paramInfo.len>0) {
                    modParGivenIndex[i].insert({j, modParGivenCount});
                    modParGivenCount++;
                }
            }
            
        }
        for(OsdiNoiseId j=0; j<desc->num_noise_src; j++) {
            std::string tmp = desc->noise_sources[j].name;
            // Add lowercase noise source name to the list of noise source names and the translator
            toLowercase(tmp);
            noiseSourceNames[i].push_back(tmp);
            noiseSourceNameTranslators[i][Id(tmp)];
        }
        auto& nnList = nodeNameLists[i];
        auto& nodeMap = nodeMaps[i];
        for(decltype(desc->num_nodes) j=0; j<desc->num_nodes; j++) {
            std::string tmp = desc->nodes[j].name;
            toLowercase(tmp);
            Id id = tmp;
            nnList.push_back(id);
            nodeMap.insert({id, j});
        }
        // Count the number of states to allocate for instance
        instanceNodeStateCounts[i] = 0;
        for (decltype(desc->num_nodes) ni = 0; ni < desc->num_nodes; ni++) {
            // Two extra states (q=charge, c=current) per each node with reactive residual
            if (desc->nodes[ni].react_residual_off != UINT32_MAX) {
                instanceNodeStateCounts[i] += 2;
            }
        }
    }

    valid = true;
}

OsdiFile::~OsdiFile() {
    closeDynamicLibrary(handle);
}

Id OsdiFile::deviceIdentifier(OsdiDeviceIndex index, Status& s) {
    if (index>=descriptorCount || descriptors==nullptr) {
        s.set(Status::NotFound, 
            std::string("Device index ")+std::to_string(index)+" out of bounds [0.."+std::to_string(descriptorCount-1)+"]."
        );
        return Id::none;
    }
    return Id(descriptors[index].name);
}

std::tuple<OsdiFile::OsdiDeviceIndex,bool> OsdiFile::deviceIndex(Id name, Status& s) {
    auto it = deviceNameToIndex.find(name);
    if (it==deviceNameToIndex.end()) {
        s.set(Status::NotFound, 
            std::string("Device '")+std::string(name)+"' not found."
        );
        return std::make_tuple(0, false);
    }
    return std::make_tuple(it->second, true);
}

OsdiDescriptor* OsdiFile::deviceDescriptor(OsdiDeviceIndex index, Status& s) {
    if (index>=descriptorCount || descriptors==nullptr) {
        s.set(Status::NotFound, 
            std::string("Device index ")+std::to_string(index)+" out of bounds [0.."+std::to_string(descriptorCount-1)+"]."
        );
        return nullptr;
    }
    return descriptors+index;
}

OsdiDescriptor* OsdiFile::deviceDescriptor(Id name, Status& s) {
    auto [index, found] = deviceIndex(name, s);
    if (!found) {
        return nullptr;
    }
    return deviceDescriptor(index, s);
}

// Under Linux reference counting is used for opening/closing dll handles
OsdiFile* OsdiFile::open(std::string file, const Loc& location, Status& s) {
    // Open
    if (Simulator::fileDebug()) {
        Simulator::dbg() << "Loading file '" << file << "'.\n";
    }
    void* handle = openDynamicLibrary(file.c_str());
    if (!handle) {
        s.set(
            Status::NotFound, 
            std::string("File '")+file+"': "+dynamicLibraryError()
        );
        s.extend(location);
        return nullptr;
    }
    
    // Is it already loaded? 
    auto existing = registry.find(handle);
    if (existing!=registry.end()) {
        // Already loaded, close to decrease handle reference count
        closeDynamicLibrary(handle);
        return existing->second.get();
    }
    
    // Create object
    auto osdiFile = new OsdiFile(handle, file, s);
    // osdiFile now owns the handle
    if (!osdiFile->isValid()) {
        s.extend(location);
        delete osdiFile;
        // No need to close handle because the desctructor did that
        return nullptr;
    }

    // Use a unique_ptr so that the OsdiFile object will get deleted when the map entry is deleted
    auto [it, inserted] = registry.insert(std::make_pair(handle, std::unique_ptr<OsdiFile>(osdiFile)));
    if (inserted) {
        return it->second.get();
    } else {
        return nullptr;
    }
}

OsdiDevice *OsdiFile::createDevice(OsdiDeviceIndex index, Id asName, Loc location, Status& s) {
    if (index>=descriptorCount) {
        s.set(Status::Range, 
            std::string("Device decriptor index (")+std::to_string(index)+") out of range."
        );
        s.extend(location);
        return nullptr;
    }
    auto device = new OsdiDevice(this, index, asName, location, s);
    return device;
}


OsdiDevice *OsdiFile::createDevice(Id name, Id asName, const Loc location, Status& s) {
    auto [index, found] = deviceIndex(name, s);
    if (!found) {
        s.set(Status::NotFound, 
            std::string("Device decriptor '")+std::string(name)+"' not found."
        );
        s.extend(location);
        return nullptr;
    }
    auto device = new OsdiDevice(this, index, asName, location, s);
    return device;
}

}
