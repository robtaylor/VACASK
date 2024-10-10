#include "core.h"
#include "common.h"

namespace NAMESPACE {

AnalysisCore::AnalysisCore(OutputDescriptorResolver& parentResolver, Circuit& circuit) 
    : parentResolver(parentResolver), circuit(circuit), savesCount(0) {
}

size_t AnalysisCore::stateStorageSize() const { 
    return coreStates.size(); 
}

size_t AnalysisCore::allocateStateStorage(size_t n) { 
    auto nOld = coreStates.size();
    coreStates.resize(nOld+n); 
    
    // Initially no state is coherent nor valid
    for(decltype(nOld) i=nOld; i<nOld+n; i++) {
        coreStates[i].coherent = false;
        coreStates[i].valid = false;
    }
    
    return nOld;
};

void AnalysisCore::deallocateStateStorage(size_t n) { 
    auto nOld = coreStates.size();
    if (n==0 || n>nOld) {
        n = nOld;
    }
    coreStates.resize(nOld-n);
}

void AnalysisCore::makeStateIncoherent(size_t ndx) { 
    coreStates.at(ndx).coherent = false;
}

void AnalysisCore::clearOutputDescriptors() {
    outputDescriptors.clear();
    outputDescriptorIndices.clear();
    savesCount = 0;
}

bool AnalysisCore::addOutputDescriptor(const OutputDescriptor& descr) {
    // Avoids duplicates
    auto [it, inserted] = outputDescriptorIndices.insert({descr.name, outputDescriptors.size()});
    if (inserted) {
        outputDescriptors.push_back(descr);
    }
    return inserted;
}

bool AnalysisCore::addOutputDescriptor(OutputDescriptor&& descr) {
    // Avoids duplicates
    auto [it, inserted] = outputDescriptorIndices.insert({descr.name, outputDescriptors.size()});
    if (inserted) {
        outputDescriptors.push_back(std::move(descr));
    }
    return inserted;
}

bool AnalysisCore::addAllUnknowns(const PTSave& save) {
    clearError();
    if (save.objName() || save.subName()) {
        // No parameters should be passed for default
        lastError = Error::Arguments;
        errorExpectedArgCount = 0;
        return false;
    }
    // Go through all variables, skip index 0 (corresponds to ground node potential),
    // add descriptors for representative nodes
    auto n = circuit.unknownCount();
    for (decltype(n) i = 1; i <= n; i++) {
        // Representative node
        auto node = circuit.reprNode(i);
        addOutputDescriptor(OutputDescriptor(OutdSolComponent, node->name(), node->name()));
    }
    savesCount++;
    return true;
}
    
bool AnalysisCore::addAllNodes(const PTSave& save) {
    clearError();
    if (save.objName() || save.subName()) {
        // No parameters should be passed for default
        lastError = Error::Arguments;
        errorExpectedArgCount = 0;
        return false;
    }
    // Go through all nodes 
    // add saves of corresponding variables
    auto nn = circuit.nodeCount();
    for (decltype(nn) i = 0; i < nn; i++) {
        // Node
        auto node = circuit.node(i);
        // Skip ground nodes
        if (node->checkFlags(Node::Flags::Ground)) {
            continue;
        }
        addOutputDescriptor(OutputDescriptor(OutdSolComponent, node->name(), node->name()));
    }
    savesCount++;
    return true;
}

bool AnalysisCore::addNode(const PTSave& save) {
    clearError();
    if (!save.objName() || save.subName()) {
        // One parameter should be passed
        lastError = Error::Arguments;
        errorExpectedArgCount = 1;
        return false;
    }
    // Create descriptor
    addOutputDescriptor(OutputDescriptor(OutdSolComponent, save.objName(), save.objName()));
    savesCount++;
    return true;
}

bool AnalysisCore::addFlow(const PTSave& save) {
    clearError();
    if (!save.objName() || save.subName()) {
        // One parameter should be passed
        lastError = Error::Arguments;
        errorExpectedArgCount = 1;
        return false;
    }
    // Node name <objName>:flow(br) 
    // Use br instead of branch because branch is a reserved word in Verilog-A
    auto id = Id(std::string(save.objName()) + ":flow(br)");
    // Create descriptor
    addOutputDescriptor(OutputDescriptor(OutdSolComponent, id, id));
    savesCount++;
    return true;
}

bool AnalysisCore::addInstanceOpvar(const PTSave& save) {
    clearError();
    if (!save.objName() || !save.subName()) {
        // Both parameters should be passed
        lastError = Error::Arguments;
        errorExpectedArgCount = 2;
        return false;
    }
    // Create descriptor
    auto name = std::string(save.objName()) + "." + std::string(save.subName());
    addOutputDescriptor(OutputDescriptor(OutdOpvar, name, save.objName(), save.subName()));
    savesCount++;
    return true;
}

bool AnalysisCore::addAllTfZin(const PTSave& save, std::unordered_map<Id,size_t>& nameMap) {
    clearError();
    if (save.objName() || save.subName()) {
        // No parameters should be passed for default
        lastError = Error::Arguments;
        errorExpectedArgCount = 0;
        return false;
    }
    // Go through all independent sources
    // add descriptors for TF and Zin
    auto ndev = circuit.deviceCount();
    for(decltype(ndev) idev=0; idev<ndev; idev++) {
        auto dev = circuit.device(idev);
        if (!dev->isSource()) {
            continue;
        }
        auto nmod = dev->modelCount();
        for(decltype(nmod) imod=0; imod<nmod; imod++) {
            auto mod = dev->model(imod);
            auto ninst = mod->instanceCount();
            for(decltype(ninst) iinst=0; iinst<ninst; iinst++) {
                auto inst = mod->instance(iinst);
                auto objName = inst->name();
                // Will insert only if entry does not exist
                auto ndx = nameMap.size();
                auto [it, inserted] = nameMap.insert({objName, ndx});
                // If entry is already present, get the index
                if (!inserted) {
                    ndx = it->second;
                }
                auto name1 = std::string("tf(") + std::string(objName) + ")";
                addOutputDescriptor(OutputDescriptor(OutdTf, name1, objName, ndx));
                auto name2 = std::string("zin(") + std::string(objName) + ")";
                addOutputDescriptor(OutputDescriptor(OutdZin, name2, objName, ndx));
            }
        }
    }
    savesCount++;
    return true;
}

bool AnalysisCore::addTf(const PTSave& save, std::unordered_map<Id,size_t>& nameMap) {
    clearError();
    if (!save.objName() || save.subName()) {
        // One parameter should be passed
        lastError = Error::Arguments;
        errorExpectedArgCount = 1;
        return false;
    }
    // Will insert only if entry does not exist
    auto ndx = nameMap.size();
    auto [it, inserted] = nameMap.insert({save.objName(), ndx});
    // If entry is already present, get the index
    if (!inserted) {
        ndx = it->second;
    }
    // Create descriptor
    auto name = std::string("tf(") + std::string(save.objName()) + ")";
    // Get index
    addOutputDescriptor(OutputDescriptor(OutdTf, name, save.objName(), ndx));
    savesCount++;
    return true;
}

bool AnalysisCore::addZin(const PTSave& save, std::unordered_map<Id,size_t>& nameMap) {
    clearError();
    if (!save.objName() || save.subName()) {
        // One parameter should be passed
        lastError = Error::Arguments;
        errorExpectedArgCount = 1;
        return false;
    }
    // Will insert only if entry does not exist
    auto ndx = nameMap.size();
    auto [it, inserted] = nameMap.insert({save.objName(), ndx});
    // If entry is already present, get the index
    if (!inserted) {
        ndx = it->second;
    }
    // Create descriptor
    auto name = std::string("zin(") + std::string(save.objName()) + ")";
    addOutputDescriptor(OutputDescriptor(OutdZin, name, save.objName(), ndx));
    savesCount++;
    return true;
}

bool AnalysisCore::addYin(const PTSave& save, std::unordered_map<Id,size_t>& nameMap) {
    clearError();
    if (!save.objName() || save.subName()) {
        // One parameter should be passed
        lastError = Error::Arguments;
        errorExpectedArgCount = 1;
        return false;
    }
    // Will insert only if entry does not exist
    auto ndx = nameMap.size();
    auto [it, inserted] = nameMap.insert({save.objName(), ndx});
    // If entry is already present, get the index
    if (!inserted) {
        ndx = it->second;
    }
    // Create descriptor
    auto name = std::string("yin(") + std::string(save.objName()) + ")";
    addOutputDescriptor(OutputDescriptor(OutdYin, name, save.objName(), ndx));
    savesCount++;
    return true;
}

bool AnalysisCore::addAllNoiseContribInst(const PTSave& save, bool details) {
    clearError();
    if (save.objName() || save.subName()) {
        // No parameters should be passed for default
        lastError = Error::Arguments;
        errorExpectedArgCount = 0;
        return false;
    }
    // Go through all instances
    // add descriptors for total instance contribution if instance has at least one noise source
    auto ndev = circuit.deviceCount();
    for(decltype(ndev) idev=0; idev<ndev; idev++) {
        auto dev = circuit.device(idev);
        auto nmod = dev->modelCount();
        for(decltype(nmod) imod=0; imod<nmod; imod++) {
            auto mod = dev->model(imod);
            auto ninst = mod->instanceCount();
            for(decltype(ninst) iinst=0; iinst<ninst; iinst++) {
                auto inst = mod->instance(iinst);
                auto objName = inst->name();
                if (inst->noiseSourceCount()<=0) {
                    continue;
                }
                auto name = std::string("n(")+std::string(objName)+")";
                addOutputDescriptor(OutputDescriptor(OutdNoiseContribInst, name, objName, Id()));
                // Go through all entries
                if (details) {
                    for(ParameterIndex i=0; i<inst->noiseSourceCount(); i++) {
                        auto contrib = inst->noiseSourceName(i);
                        auto name = std::string("n(")+std::string(objName)+","+std::string(contrib)+")";
                        addOutputDescriptor(OutputDescriptor(OutdNoiseContribInstPartial, name, objName, contrib));
                    }
                }
            }
        }
    }
    savesCount++;
    return true;
}

bool AnalysisCore::addNoiseContribInst(const PTSave& save, bool details) {
    clearError();
    if (details) {
        // Expect only one argument
        if (!save.objName() || save.subName()) {
            lastError = Error::Arguments;
            errorExpectedArgCount = 1;
            return false;
        }
        // Get instance
        auto instName = save.objName();
        auto inst = circuit.findInstance(instName);
        if (inst) {
            // Instance contribution
            auto name = std::string("n(") + std::string(instName) + ")";
            addOutputDescriptor(OutputDescriptor(OutdNoiseContribInst, name, save.objName()));
            savesCount++;
            // Contributions of sources
            auto nNoise = inst->noiseSourceCount();
            for(decltype(nNoise) i=0; i<nNoise; i++) {
                auto name = std::string("n(") + std::string(instName) + "," + std::string(inst->noiseSourceName(i)) + ")";
                addOutputDescriptor(OutputDescriptor(OutdNoiseContribInstPartial, name, instName, inst->noiseSourceName(i)));
                savesCount++;
            }
        }
    } else {
        // One or two arguments
        if (!save.objName() && !save.subName()) {
            lastError = Error::Arguments;
            errorExpectedArgCount = 1;
            return false;
        }
        if (!save.subName()) {
            // One argument
            auto objName = save.objName();
            auto name = std::string("n(") + std::string(save.objName()) + ")";
            addOutputDescriptor(OutputDescriptor(OutdNoiseContribInst, name, save.objName()));
        } else {
            // Two arguments 
            auto name = std::string("n(") + std::string(save.objName()) + "," + std::string(save.subName())+ ")";
            addOutputDescriptor(OutputDescriptor(OutdNoiseContribInst, name, save.objName(), save.subName()));
        }
        savesCount++;
    }
    return true;
}
    
bool AnalysisCore::addRealVarOutputSource(bool strict, Id name, const Vector<double>& solution) {
    clearError();
    // Solution vector component
    auto node = circuit.findNode(name);
    // Get unknown
    if (node) {
        auto unknown = node->unknownIndex();
        outputSources.emplace_back(&solution, node->unknownIndex());
        return true;
    } else if (strict) {
        lastError = Error::NodeNotFound;
        errorId = name;
        return false;
    } else {
        outputSources.emplace_back();
    }
    return true;
}

bool AnalysisCore::addRealVarOutputSource(bool strict, Id name, const VectorRepository<double>& solution) {
    clearError();
    // Solution vector component
    auto node = circuit.findNode(name);
    // Get unknown
    if (node) {
        auto unknown = node->unknownIndex();
        outputSources.emplace_back(&solution, node->unknownIndex());
        return true;
    } else if (strict) {
        lastError = Error::NodeNotFound;
        errorId = name;
        return false;
    } else {
        outputSources.emplace_back();
    }
    return true;
}

bool AnalysisCore::addComplexVarOutputSource(bool strict, Id name, const Vector<Complex>& solution) {
    clearError();
    // Solution vector component
    auto node = circuit.findNode(name);
    // Get unknown
    if (node) {
        auto unknown = node->unknownIndex();
        outputSources.emplace_back(&solution, node->unknownIndex());
        return true;
    } else if (strict) {
        lastError = Error::NodeNotFound;
        errorId = name;
        return false;
    } else {
        outputSources.emplace_back();
    }
    return true;
}

bool AnalysisCore::addComplexVarOutputSource(bool strict, Id name, const VectorRepository<Complex>& solution) {
    clearError();
    // Solution vector component
    auto node = circuit.findNode(name);
    // Get unknown
    if (node) {
        auto unknown = node->unknownIndex();
        outputSources.emplace_back(&solution, node->unknownIndex());
        return true;
    } else if (strict) {
        lastError = Error::NodeNotFound;
        errorId = name;
        return false;
    } else {
        outputSources.emplace_back();
    }
    return true;
}

bool AnalysisCore::addOpvarOutputSource(bool strict, Id instance, Id opvar) {
    clearError();
    auto inst = circuit.findInstance(instance);
    if (inst) {
        // Find opvar
        auto [ndx, found] = inst->opvarIndex(opvar);
        if (found) {
            auto [ok, osrc] = inst->opvarOutputSource(ndx);
            if (!ok && strict) {
                return false;
            }
            outputSources.push_back(std::move(osrc));
        } else if (strict) {
            lastError = Error::OpvarNotFound;
            errorId = instance;
            errorId2 = opvar;
            return false;
        }
    } else if (strict) {
        lastError = Error::InstanceNotFound;
        errorId = instance;
        return false;
    } else {
        outputSources.emplace_back();
    }
    return true;
}

std::tuple<bool, UnknownIndex, UnknownIndex> AnalysisCore::getOutput(Value& v) {
    clearError();
    if (v.type()==Value::Type::String) {
        // Output is a string
        Id id = v.val<String>();
        auto node = circuit.findNode(id);
        if (!node) {
            lastError = Error::NodeNotFound;
            errorId = id;
            return std::make_tuple(false, 0, 0);
        }
        return std::make_tuple(true, node->unknownIndex(), 0);
    } else if (v.type()==Value::Type::StringVec) {
        // Output is a string vector
        Node* nodep;
        Node* noden;
        Id idp, idn;
        UnknownIndex up, un;
        StringVector& sv = v.val<StringVector>();
        switch (sv.size()) {
            case 1:
                idp = sv[0];
                nodep = circuit.findNode(idp);
                if (!nodep) {
                    lastError = Error::NodeNotFound;
                    errorId = idp;
                    return std::make_tuple(false, 0, 0);
                }
                return std::make_tuple(true, nodep->unknownIndex(), 0);
            case 2:
                idp = sv[0];
                nodep = circuit.findNode(idp);
                if (!nodep) {
                    lastError = Error::NodeNotFound;
                    errorId = idp;
                    return std::make_tuple(false, 0, 0);
                }
                up = nodep->unknownIndex();
                
                idn = sv[1];
                noden = circuit.findNode(idn);
                if (!noden) {
                    lastError = Error::NodeNotFound;
                    errorId = idn;
                    return std::make_tuple(false, 0, 0);
                }
                un = noden->unknownIndex();
                return std::make_tuple(true, up, un);
            default:
                lastError = Error::OutputSpec;
                return std::make_tuple(false, 0, 0);  
        }
    } else {
        lastError = Error::OutputType;
        return std::make_tuple(false, 0, 0);  
    }
}

std::tuple<bool, Instance*> AnalysisCore::getInput(Id name) {
    clearError();
    auto inst = circuit.findInstance(name);
    if (!inst) {
        lastError = Error::InstanceNotFound;
        errorId = name;
        return std::make_tuple(false, nullptr);
    }
    if (!inst->model()->device()->isSource()) {
        lastError = Error::InstanceNotSource;
        errorId = name;
        return std::make_tuple(false, nullptr);
    }
    return std::make_tuple(true, inst);
}

bool AnalysisCore::formatError(Status& s) const {
    switch (lastError) {
        case Error::Arguments:
            switch (errorExpectedArgCount) {
                case 0: 
                    s.set(Status::Save, "Save directive does not accept arguments.");
                    break;
                case 1:
                    s.set(Status::Save, "Save directive requires one argument.");
                    break;
                default:
                    s.set(Status::Save, "Save directive requires "+std::to_string(errorExpectedArgCount)+" arguments.");
                    break;
            }
            return false;
        case Error::NodeNotFound:
            s.set(Status::Analysis, std::string("Node '")+std::string(errorId)+"' not found.");
            return false;
        case Error::OpvarNotFound:
            s.set(Status::Analysis, std::string("Opvar '")+std::string(errorId2)+"' of instance '"+std::string(errorId)+"' not found.");
            return false;
        case Error::InstanceNotFound:
            s.set(Status::Analysis, std::string("Instance '")+std::string(errorId)+"' not found.");
            return false;
        case Error::OutputSpec:
            s.set(Status::Analysis, "Output must must be a single node or a node pair.");
            return false;
        case Error::OutputType:
            s.set(Status::Analysis, "Output specification must be a string or a string vector.");
            return false;
        case Error::InstanceNotSource:
            s.set(Status::Analysis, std::string("Instance '")+std::string(errorId)+"' is not an independent source.");
            return false;
        case Error::Descriptor:
            s.set(Status::Analysis, std::string("Failed to add output descriptor for '")+std::string(errorId)+"'.");
            return false;
    }
    return true;
}

void AnalysisCore::dump(std::ostream& os) const {
    os << "  Output variables: " << "\n";
    for(auto& desc : outputDescriptors) {
        os << "    " << desc.name << "\n";
    }
}

}
