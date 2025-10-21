#include "natures.h"
#include "common.h"


namespace NAMESPACE {

std::unordered_map<NatureId,std::unordered_set<std::string>> NatureRegistry::natureFilesMap;

NatureId NatureRegistry::natureId(std::string file, std::string name) {
    return natureIdInternal(file, name);
}

NatureId NatureRegistry::natureIdInternal(std::string file, std::string name, bool staticId){
    Id id;
    if (staticId) {
        id = Id::createStatic(name.c_str());
    } else {
        id = Id(name);
    }
    auto [it, inserted] = NatureRegistry::natureFilesMap.try_emplace(id);
    it->second.insert(file);
    return id;
}

const NatureId NatureRegistry::noNature = NatureRegistry::natureIdInternal("", "", true);
const NatureId NatureRegistry::spiceVoltage = NatureRegistry::natureIdInternal("", ".voltage", true);
const NatureId NatureRegistry::spiceCurrent = NatureRegistry::natureIdInternal("", ".current", true);
const NatureId NatureRegistry::spiceFlux = NatureRegistry::natureIdInternal("", ".flux", true);
const NatureId NatureRegistry::spiceCharge = NatureRegistry::natureIdInternal("", ".charge", true);


const NaturesSubset::NatureIndex NaturesSubset::noNature     = 0; 
const NaturesSubset::NatureIndex NaturesSubset::spiceVoltage = 1;
const NaturesSubset::NatureIndex NaturesSubset::spiceCurrent = 2;
const NaturesSubset::NatureIndex NaturesSubset::spiceFlux    = 3;
const NaturesSubset::NatureIndex NaturesSubset::spiceCharge  = 4;

NaturesSubset::NaturesSubset() {
    reset();
}

void NaturesSubset::reset() {
    // Clear structures
    idToIndex.clear();
    indexToId.clear();

    // Register no nature and SPICE natures. 
    // Order is important because indices are defined statically. 
    registerNature(NatureRegistry::noNature);     // 0
    registerNature(NatureRegistry::spiceVoltage); // 1
    registerNature(NatureRegistry::spiceCurrent); // 2
    registerNature(NatureRegistry::spiceFlux);    // 3
    registerNature(NatureRegistry::spiceCharge);  // 4
}

NaturesSubset::NatureIndex NaturesSubset::registerNature(NatureId id) {
    if (idToIndex.contains(id)) {
        return idToIndex[id];
    } else {
        auto ndx = idToIndex.size();
        idToIndex[id] = ndx;
        indexToId.push_back(id);
        return ndx;
    }
}

}
