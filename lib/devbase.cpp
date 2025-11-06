#include "devbase.h"
#include "circuit.h"
#include "common.h"

namespace NAMESPACE {

Device::Device(Id name, const Loc& location) 
    : name_(name), loc(location), tovh(0), novh(0) {
}

Device::~Device() {
}

bool Device::addModel(Model* model) { 
    if (model->device()!=this) {
        return false;
    }
    models_.push_back(model); 
    return true;
};
    

InstantiationData::InstantiationData() {
    // Start with no ancestor models
}

InstantiationData::InstantiationData(Instance* inst) {
    // For building ancestorPathStack we need the sequence of instances from toplevel to inst
    std::vector<Instance*> instPath;
    
    // Start at given instance, collect all ancestor models up to top level instance
    for(;inst; inst=inst->parent()) {
        ancestorModels_.insert(inst->model());
        instPath.push_back(inst);
    }

    // Build ancestor path stack (master names) from toplevel down to inst
    std::string ancestorPath;
    // Reverse iterator over instPath vector because toplevel instance is last
    for (auto it = instPath.rbegin(); it != instPath.rend(); ++it) {
        if (ancestorPath.size()>0) {
            ancestorPath += ":";
        }
        ancestorPath += std::string((*it)->model()->name());
        ancestorPathStack_.push_back(ancestorPath);
    }
}

bool InstantiationData::addAncestor(Instance* inst) {
    // Add ancestor model, return true if the model was inserted (was not in the set before)
    if (inst) {
        auto [it, inserted] = ancestorModels_.insert(inst->model());
        if (!inserted) {
            return false;
        }
        if (inst->parent()) {
            // Not toplevel instance
            auto& base = ancestorPathStack_.back();
            if (base.size()>0) {
                ancestorPathStack_.push_back(base+":"+std::string(inst->model()->name()));
            } else {
                ancestorPathStack_.push_back(inst->model()->name());
            }
        } else {
            // Toplevel instance, ancestor path is empty for all its children
            ancestorPathStack_.push_back("");
        }
    }
    return true;
}

bool InstantiationData::removeAncestor(Instance* inst) {
    // Remove ancestor from set, return true if the ancestor was found and removed
    if (inst) {
        auto numRemoved = ancestorModels_.erase(inst->model());
        if (numRemoved=0) {
            return false;
        }
        ancestorPathStack_.pop_back();
    }
    return true;
}

Id InstantiationData::translateDefinition(Id definitionName) { 
    if (ancestorPathStack_.size()>0) {
        return definitionName; 
    } else { 
        return ancestorPathStack_.back()+":"+std::string(definitionName); 
    }
}

bool InstantiationData::isAncestor(Model* model) {
    return ancestorModels_.contains(model);
}


Model::Model(Device* device, Id name, Instance* parent, const PTModel& parsedModel) 
    : name_(name), device_(device), parent_(parent), parsedModel_(parsedModel) {
    if (parent) {
        parent->addChild(this);
    }
    setFlags(Flags::NeedsSetup);
}

Model::~Model() {
}

bool Model::addInstance(Instance* instance) { 
    if (instance->model()!=this) {
        return false;
    }
    instances_.push_back(instance); 
    device_->instanceCount_++;

    return true;
};

// TODO: maybe find a faster way for checking this instead of linear search
bool Model::parameterIsFree(Id name) {
    auto& params = parsedModel_.parameters();
    for(auto& it : params.expressions()) {
        if (it.name()==name) {
            return false;
        }
    }
    return true;
}

Instance::Instance(Model* model, Id name, Instance* parent, const PTInstance& parsedInstance) 
    : name_(name), model_(model), parent_(parent), parsedInstance_(parsedInstance) {
    if (parent) {
        parent->addChild(this);
    }
    setFlags(Flags::NeedsSetup);
}

Instance::~Instance() {
}

// TODO: maybe find a faster way for checking this instead of linear search
bool Instance::parameterIsFree(Id name) {
    auto& params = parsedInstance_.parameters();
    for(auto& it : params.expressions()) {
        if (it.name()==name) {
            return false;
        }
    }
    return true;
}

Id Instance::translate(Id child) { 
    if (!parent_) {
        // Toplevel instance has an empty path, does not prepend it to its children
        return child;
    } else {
        // Other instances prepend their name to their children
        return Id(std::string(name_)+":"+std::string(child));
    }
}

Id Instance::translateNode(Circuit& cir, Id nodeName) {
    // Global nodes are not translated (this includes ground nodes, as those are global, too)
    if (cir.isGlobalNode(nodeName)) {
        return nodeName;
    }

    // Otherwise translate in the same way as subinstances/submodels are translated
    return translate(nodeName);
}

Id Instance::translatePeer(Id peer) {
    auto p = parent();
    if (p) {
        // We have a parent, use it to translate peer
        return p->translate(peer);
    } else {
        // No parent, somebody is trying to translate a peer of toplevel subcircuit
        // No translation
        return peer;
    }
}

std::tuple<Value::Type,bool> Instance::opvarType(Id name, Status& s) const {
    auto [ndx, found] = opvarIndex(name);
    if (!found) {
        s.set(Status::NotFound, std::string("Opvar '")+std::string(name)+"' not found.");
        return std::make_tuple(Value::Type::Int, false);
    }
    return opvarType(ndx, s);
}

bool Instance::getOpvar(Id name, Value& v, Status& s) const {
    auto [ndx, found] = opvarIndex(name);
    if (!found) {
        s.set(Status::NotFound, std::string("Opvar '")+std::string(name)+"' not found.");
        return false;
    }
    return getOpvar(ndx, v, s);
}

Node* Instance::getInternalNode(Circuit& circuit, const std::string& name, Node::Flags flags, Status& s) {
    Id nodeName = translate(name);
    auto node = circuit.getNode(nodeName, flags, s);
    if (node==nullptr) {
        s.extend(std::string("Failed to obtain internal node '"+std::string(nodeName)+"' from simulator."));
        s.extend(location());
        return nullptr;
    }
    node->setFlags(Node::Flags::InternalDeviceNode);
    return node;
}

std::tuple<bool, size_t> Instance::enterContext(Circuit& circuit, Context* externalContext, bool addToPath, bool rebuild, Status& s) { 
    return std::make_tuple(true, circuit.paramEvaluator().contextMarker()); 
}

bool Instance::revertContext(Circuit& circuit, size_t contextMarker) { 
    return circuit.paramEvaluator().revertContext(contextMarker); 
}

// Construct end iterator
Instance::HierarchicalIterator::HierarchicalIterator() 
    : forcePeer(false) {
}

// Construct begin iterator
Instance::HierarchicalIterator::HierarchicalIterator(pointer instance) 
    : forcePeer(false) {
    dummy.push_back(instance);
    stack.push_back(std::make_pair(&dummy, 0));
}

// *it
Instance::HierarchicalIterator::reference Instance::HierarchicalIterator::operator*() const { 
    auto [insts, pos] = stack.back();
    return *((*insts)[pos]);
};

// it->
Instance::HierarchicalIterator::pointer Instance::HierarchicalIterator::operator->() {
    auto [insts, pos] = stack.back();
    return (*insts)[pos];
};

// ++it
Instance::HierarchicalIterator& Instance::HierarchicalIterator::operator++() { 
    auto position = stack.back(); 
    auto instance = position.first->at(position.second);
    
    // Calling operator++() means we want to move to next instance. 
    // We have two cases depending on the instance pointed to by the iterator
    // 1) instance is hierarchical and has children
    //    go down (start exploring children) 
    //    i.e. push childInstances vector on stack and set counter to 0
    // 2) instance is not hierarchical or has no children
    //    move to a per instance (increase counter), 
    //    when peers are exhausted, move up (pop stack)
    
    // If forcePeer is true we must move to a peer even 
    // if there are further children to process. 
    if (
        instance->model()->device()->isHierarchical() && 
        instance->childInstances() && 
        instance->childInstances()->size()>0 &&
        !forcePeer 
    ) {
        // At hierarchical instance with children, go down
        stack.emplace_back(instance->childInstances(), 0);
    } else {
        // Not a hierarchical instance, or a hierarchical instance without children
        // Go to peer
        while (stack.size()>0) {
            // Repeat until a peer is found
            stack.back().second++;
            if (stack.back().second>=stack.back().first->size()) {
                // No more peers, go up
                stack.pop_back();
            } else {
                // Have peer
                break;
            }
            // No more instances
        }
    }
    // Clear forcePeer flag
    forcePeer = false;
    return *this;
}

Instance::HierarchicalIterator Instance::HierarchicalIterator::operator++(int) { 
    auto tmp = *this; 
    ++(*this); 
    return tmp; 
}

bool operator==(const Instance::HierarchicalIterator& a, const Instance::HierarchicalIterator& b) { 
    return a.stack == b.stack; 
}

bool operator!=(const Instance::HierarchicalIterator& a, const Instance::HierarchicalIterator& b) { 
    return a.stack != b.stack; 
};

void Instance::HierarchicalIterator::stopDescent() {
    // Next increment should not go down the current instance
    // Instead it should move to next peer instance
    forcePeer = true;
}

}
