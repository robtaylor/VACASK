#include "core.h"
#include "common.h"

namespace NAMESPACE {

AnalysisCore::AnalysisCore(Analysis& analysis, Circuit& circuit) 
    : analysis(analysis), circuit(circuit), savesCount(0) {
}

bool AnalysisCore::addCoreOutputDescriptors(Status& s) { 
    return true; 
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

bool AnalysisCore::addAllUnknowns(const PTSave& save, bool verify, Status& s) {
    if (save.objName() || save.subName()) {
        // No parameters should be passed for default
        s.set(Status::Save, "Default save must not have arguments.");
        s.extend(save.location());
        return false;
    }
    // Go through all variables, skip index 0 (corresponds to ground node potential),
    // add descriptors for representative nodes
    auto n = circuit.unknownCount();
    for (decltype(n) i = 1; i <= n; i++) {
        // Representative node
        auto node = circuit.reprNode(i);
        // No need to verify it, just add it
        addOutputDescriptor(OutputDescriptor(OutdSolComponent, node->name(), node->name()));
    }
    savesCount++;
    return true;
}
    
bool AnalysisCore::addAllNodes(const PTSave& save, bool verify, Status& s) {
    if (save.objName() || save.subName()) {
        // No parameters should be passed for default
        s.set(Status::Save, "Default full save must not have arguments.");
        s.extend(save.location());
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
        // No need to verify it, just add it
        addOutputDescriptor(OutputDescriptor(OutdSolComponent, node->name(), node->name()));
    }
    savesCount++;
    return true;
}

bool AnalysisCore::addNode(const PTSave& save, bool verify, Status& s) {
    if (!save.objName() || save.subName()) {
        // One parameter should be passed
        s.set(Status::Save, "Node value save directive requires one argument.");
        s.extend(save.location());
        return false;
    }
    // Create descriptor
    addOutputDescriptor(OutputDescriptor(OutdSolComponent, save.objName(), save.objName()));
    savesCount++;
    return true;
}

bool AnalysisCore::addFlow(const PTSave& save, bool verify, Status& s) {
    if (!save.objName() || save.subName()) {
        // One parameter should be passed
        s.set(Status::Save, "Flow value save directive requires one argument.");
        s.extend(save.location());
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

bool AnalysisCore::addInstanceOpvar(const PTSave& save, bool verify, Status& s) {
    if (!save.objName() || !save.subName()) {
        // Both parameters should be passed
        s.set(Status::Save, "Opvar save directive requires two arguments.");
        s.extend(save.location());
        return false;
    }
    // Create descriptor
    auto name = std::string(save.objName()) + "." + std::string(save.subName());
    addOutputDescriptor(OutputDescriptor(OutdOpvar, name, save.objName(), save.subName()));
    savesCount++;
    return true;
}

bool AnalysisCore::addAllTfZin(const PTSave& save, bool verify, std::unordered_map<Id,size_t>& nameMap, Status& s) {
    if (save.objName() || save.subName()) {
        // No parameters should be passed for default
        s.set(Status::Save, "Default save must not have arguments.");
        s.extend(save.location());
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

bool AnalysisCore::addTf(const PTSave& save, bool verify, std::unordered_map<Id,size_t>& nameMap, Status& s) {
    if (!save.objName() || save.subName()) {
        // One parameter should be passed
        s.set(Status::Save, "TF save directive requires one argument.");
        s.extend(save.location());
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

bool AnalysisCore::addZin(const PTSave& save, bool verify, std::unordered_map<Id,size_t>& nameMap, Status& s) {
    if (!save.objName() || save.subName()) {
        // One parameter should be passed
        s.set(Status::Save, "Zin save directive requires one argument.");
        s.extend(save.location());
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

bool AnalysisCore::addYin(const PTSave& save, bool verify, std::unordered_map<Id,size_t>& nameMap, Status& s) {
    if (!save.objName() || save.subName()) {
        // One parameter should be passed
        s.set(Status::Save, "Yin save directive requires one argument.");
        s.extend(save.location());
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

bool AnalysisCore::addAllNoiseContribInst(const PTSave& save, bool verify, bool details, Status& s) {
    if (save.objName() || save.subName()) {
        // No parameters should be passed for default
        s.set(Status::Save, "Default save must not have arguments.");
        s.extend(save.location());
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

bool AnalysisCore::addNoiseContribInst(const PTSave& save, bool verify, bool details, Status& s) {
    if (details) {
        // Expect only one argument
        if (!save.objName() || save.subName()) {
            s.set(Status::Save, "Full instance noise save must have exactly one argument.");
            s.extend(save.location());
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
            s.set(Status::Save, "Instance noise save must have at least one argument.");
            s.extend(save.location());
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
    
bool AnalysisCore::addRealVarOutputSource(bool strict, Id name, const Vector<double>& solution, Status& s) {
    // Solution vector component
    auto node = circuit.findNode(name);
    // Get unknown
    if (node) {
        auto unknown = node->unknownIndex();
        outputSources.emplace_back(&solution, node->unknownIndex());
        return true;
    } else if (strict) {
        s.set(Status::NotFound, std::string("Node '")+std::string(name)+"' not found.");
        return false;
    } else {
        outputSources.emplace_back();
    }
    return true;
}

bool AnalysisCore::addRealVarOutputSource(bool strict, Id name, const VectorRepository<double>& solution, Status& s) {
    // Solution vector component
    auto node = circuit.findNode(name);
    // Get unknown
    if (node) {
        auto unknown = node->unknownIndex();
        outputSources.emplace_back(&solution, node->unknownIndex());
        return true;
    } else if (strict) {
        s.set(Status::NotFound, std::string("Node '")+std::string(name)+"' not found.");
        return false;
    } else {
        outputSources.emplace_back();
    }
    return true;
}

bool AnalysisCore::addComplexVarOutputSource(bool strict, Id name, const Vector<Complex>& solution, Status& s) {
    // Solution vector component
    auto node = circuit.findNode(name);
    // Get unknown
    if (node) {
        auto unknown = node->unknownIndex();
        outputSources.emplace_back(&solution, node->unknownIndex());
        return true;
    } else if (strict) {
        s.set(Status::NotFound, std::string("Node '")+std::string(name)+"' not found.");
        return false;
    } else {
        outputSources.emplace_back();
    }
    return true;
}

bool AnalysisCore::addComplexVarOutputSource(bool strict, Id name, const VectorRepository<Complex>& solution, Status& s) {
    // Solution vector component
    auto node = circuit.findNode(name);
    // Get unknown
    if (node) {
        auto unknown = node->unknownIndex();
        outputSources.emplace_back(&solution, node->unknownIndex());
        return true;
    } else if (strict) {
        s.set(Status::NotFound, std::string("Node '")+std::string(name)+"' not found.");
        return false;
    } else {
        outputSources.emplace_back();
    }
    return true;
}

bool AnalysisCore::addOpvarOutputSource(bool strict, Id instance, Id opvar, Status& s) {
    auto inst = circuit.findInstance(instance);
    if (inst) {
        // Find opvar
        auto [ndx, found] = inst->opvarIndex(opvar);
        if (found) {
            Status tmps;
            auto [ok, osrc] = inst->opvarOutputSource(ndx, tmps);
            if (!ok && strict) {
                s = std::move(tmps);
                return false;
            }
            outputSources.push_back(std::move(osrc));
        } else if (strict) {
            s.set(Status::NotFound, std::string("Opvar '")+std::string(opvar)+"' of instance '"+std::string(instance)+"' not found.");
            return false;
        }
    } else if (strict) {
        s.set(Status::NotFound, std::string("Instance '")+std::string(instance)+"' not found.");
        return false;
    } else {
        outputSources.emplace_back();
    }
    return true;
}

std::tuple<bool, UnknownIndex, UnknownIndex> AnalysisCore::getOutput(Value& v, Status& s) {
    if (v.type()==Value::Type::String) {
        // Output is a string
        Id id = v.val<String>();
        auto node = circuit.findNode(id);
        if (!node) {
            s.set(Status::NotFound, std::string("Output node '")+std::string(id)+"' not found.");
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
                    s.set(Status::NotFound, std::string("Output node '")+std::string(idp)+"' not found.");
                    return std::make_tuple(false, 0, 0);
                }
                return std::make_tuple(true, nodep->unknownIndex(), 0);
            case 2:
                idp = sv[0];
                nodep = circuit.findNode(idp);
                if (!nodep) {
                    s.set(Status::NotFound, std::string("Output positive node '")+std::string(idp)+"' not found.");
                    return std::make_tuple(false, 0, 0);
                }
                up = nodep->unknownIndex();
                
                idn = sv[1];
                noden = circuit.findNode(idn);
                if (!noden) {
                    s.set(Status::NotFound, std::string("Output negative node '")+std::string(idn)+"' not found.");
                    return std::make_tuple(false, 0, 0);
                }
                un = noden->unknownIndex();
                return std::make_tuple(true, up, un);
            default:
                s.set(Status::BadArguments, "Output must must be a single node or a node pair.");
                return std::make_tuple(false, 0, 0);  
        }
    } else {
        s.set(Status::BadArguments, "Output must be a string or a string vector.");
        return std::make_tuple(false, 0, 0);  
    }
}

std::tuple<bool, Instance*> AnalysisCore::getInput(Id name, Status& s) {
    auto inst = circuit.findInstance(name);
    if (!inst) {
        s.set(Status::NotFound, std::string("Instance '")+std::string(name)+"' not found.");
        return std::make_tuple(false, nullptr);
    }
    if (!inst->model()->device()->isSource()) {
        s.set(Status::BadArguments, std::string("Instance '")+std::string(name)+"' is not an independent source.");
        return std::make_tuple(false, nullptr);
    }
    return std::make_tuple(true, inst);
}

void AnalysisCore::dump(std::ostream& os) const {
    os << "  Output variables: " << std::endl;
    for(auto& desc : outputDescriptors) {
        os << "    " << desc.name << std::endl;
    }
}

}
