#include <cstring>
#include "osdimodel.h"
#include "osdiinstance.h"
#include "circuit.h"
#include "libplatform.h"
#include "common.h"


namespace NAMESPACE {

OsdiModel::OsdiModel(OsdiDevice *device, Id name, Instance* parentInstance, const PTModel& parsedModel, Status& s) 
    : Model(device, name, parentInstance, parsedModel), core_(nullptr) {
    core_ = alignedAlloc(sizeof(max_align_t), device->descriptor()->model_size);
    memset(core_, 0, device->descriptor()->model_size);
    
    setFlags(Flags::IsValid);
}

OsdiModel::~OsdiModel() {
    // Free allocated values (strings and vectors)
    device()->freeValues(core_, nullptr);

    alignedFree(core_);
}

Instance* OsdiModel::createInstance(Circuit& circuit, Instance* parentInstance, RpnEvaluator& evaluator, Context* externalContext, const PTInstance& parsedInstance, InstantiationData& idata, Status& s) {
    auto name = parsedInstance.name();
    
    // If we have a hierarchical parent translate name
    if (parentInstance) {
        name = parentInstance->translate(name);
    }

    // Create instance
    auto* instance = new OsdiInstance(this, name, parentInstance, parsedInstance, s);
    if (!instance->checkFlags(Instance::Flags::IsValid)) {
        s.extend(parsedInstance.location());
        delete instance;
        return nullptr;
    }

    // Add to instanceMap of circuit
    if (!circuit.add(instance, s)) {
        delete instance;
        return nullptr;
    }

    // Sanity check (connections should not exceed terminals)
    auto nterm = device()->terminalCount();
    if (parsedInstance.connections().size()>nterm) {
        s.set(Status::Range, "Too many terminals specified at instantiation.");
        s.extend(parsedInstance.connections().at(nterm).location());
        return nullptr;
    }

    // Bind terminals
    auto& terms = parsedInstance.connections();
    TerminalIndex i=0;
    for(auto it=terms.cbegin(); it!=terms.cend(); ++it, i++) {
        // Translate, if needed
        auto nodeName = it->name();
        nodeName = parentInstance->translateNode(circuit, nodeName);
        auto node = circuit.getNode(nodeName, Node::Flags::PotentialNode, s);
        if (node == nullptr) {
            s.extend(std::string("Failed to obtain node '"+std::string(nodeName)+"' from simulator."));
            s.extend(it->location());
            return nullptr;
        }
        if (!instance->bindTerminal(i, node, s)) {
            s.extend(std::string("Failed to bind terminal ")+std::to_string(i+1)+".");
            s.extend(it->location());
            return nullptr;
        }
    }

    // Set instance's parameters, use the evaluator whose latest context is the parent instance's context
    auto [ok, changed] = instance->setParameters(parsedInstance.parameters(), evaluator, s);
    if (!ok) {
        return nullptr;
    }

    // Due to setting parameters the ParamsChanged flag is probably set indicating 
    // that parameter propagation is needed. Because we just created the instance 
    // propagation is not needed so clear the flag. 
    instance->clearFlags(Instance::Flags::ParamsChanged);

    // Build all children (models, nodes, instances)
    // Use the evaluator whose latest context is the parent instance's context
    // In our case we build only internal nodes
    if (!instance->buildHierarchy(circuit, evaluator, s)) {
        return nullptr;
    }
    
    return instance;
}

ParameterIndex OsdiModel::parameterCount() const {
    return device()->modelParameterCount();
}

std::tuple<ParameterIndex, bool> OsdiModel::parameterIndex(Id name) const {
    auto [index, found] = device()->modelParameterIndex(name);
    if (!found) {
        return std::make_tuple(0, false);
    }
    return std::make_tuple(index, true);
}

Id OsdiModel::parameterName(ParameterIndex ndx) const {
    if (ndx<parameterCount())
        return device()->modelParameterName(ndx);
    else
        return Id::none;
}

std::tuple<Value::Type,bool> OsdiModel::parameterType(ParameterIndex ndx, Status& s) const {
    if (ndx>=device()->modelParameterCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return std::make_tuple(Value::Type::Int, false);
    }
    auto osdiId = device()->modelOsdiParameterId(ndx);
    return std::make_tuple(device()->parameterType(osdiId), true);
}

bool OsdiModel::getParameter(ParameterIndex ndx, Value& v, Status& s) const {
    if (ndx>=device()->modelParameterCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return false;
    }
    auto osdiId = device()->modelOsdiParameterId(ndx);
    return device()->readParameter(osdiId, core_, nullptr, v);
}

std::tuple<bool,bool> OsdiModel::setParameter(ParameterIndex ndx, const Value& v, Status& s) {
    if (ndx>=device()->modelParameterCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return std::make_tuple(false, false);
    }
    auto osdiId = device()->modelOsdiParameterId(ndx);

    // Write
    auto [ok, changed] = device()->writeParameter(osdiId, core_, nullptr, v, s);

    // Mark model for setup
    if (changed) {
        setFlags(Model::Flags::NeedsSetup);
    }

    return std::make_tuple(ok, changed);
}

std::tuple<bool,bool> OsdiModel::parameterGiven(ParameterIndex ndx, Status& s) const {
    if (ndx>=device()->modelParameterCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return std::make_tuple(false, false);
    }
    auto osdiId = device()->modelOsdiParameterId(ndx);
    return device()->parameterGiven(osdiId, core_, nullptr, s);
}

// Set up this model (virtual method)
std::tuple<bool, bool, bool> OsdiModel::setup(Circuit& circuit, bool force, DeviceRequests* devReq, Status& s) {
    OsdiSimParas sp;
    auto& opt = circuit.simulatorOptions().core();
    auto& internals = circuit.simulatorInternals(); 

    // Allocate tables on stack
    auto [ndbl, nchrptr ] = device()->simParasSizes();
    double dblArray[ndbl];
    char* chrPtrArray[nchrptr];
    
    OsdiDevice::populateSimParas(sp, opt, internals, dblArray, chrPtrArray);
    // Verilog-A $temperature is in K, convert the value given by options (in C)
    auto retval = setupCore(circuit, sp, opt.temp+273.15, force, devReq, s);
    return retval;
}

// Set up this model (method used by friends for inlining)
std::tuple<bool, bool, bool> OsdiModel::setupCore(Circuit& circuit, OsdiSimParas& sp, double temp, bool force, DeviceRequests* devReq, Status& s) {
    auto handle = OsdiCallbackHandle {
        .kind = 1, 
        .name = const_cast<char*>(name().c_str())
    };
    OsdiInitInfo initInfo;
    
    // Do we need to set up all instances? 
    bool forceAllInstances = false;
    if (force || checkFlags(Flags::NeedsSetup)) {
        device()->descriptor()->setup_model((void*)&handle, core(), &sp, &initInfo);
        if (!device()->processInitInfo(circuit, initInfo, "Model", name(), devReq, s)) {
            // The problem is big enough to abort simulation
            return std::make_tuple(false, false, false);
        }

        clearFlags(Flags::NeedsSetup);
        
        // After model setup all instances have to be set up
        forceAllInstances = true;
    }
    
    bool unknownsChanged = false;
    bool sparsityChanged = false;
    for(auto it : instances()) {
        auto [ok, tmpUnknowns, tmpSparsity] = static_cast<OsdiInstance*>(it)->setupCore(circuit, sp, temp, forceAllInstances, devReq, s);
        unknownsChanged |= tmpUnknowns;
        sparsityChanged |= tmpSparsity;
        if (!ok) {
            return std::make_tuple(false, unknownsChanged, sparsityChanged);
        }
    }

    return std::make_tuple(true, unknownsChanged, sparsityChanged);
}

void OsdiModel::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "Osdi device model " << std::string(name()) << " of device " << device()->name() << "\n";
    if (parameterCount()>0) {
        os << pfx << "  Parameters:\n";
        auto np = parameterCount();
        for(decltype(np) i=0; i<np; i++) {
            Value v;
            getParameter(i, v);
            auto [ok, given] = parameterGiven(i);
            os << pfx << "    " << std::string(parameterName(i)) << " = " << v << " (" << v.typeName() << ")";
            if (!given) {
                os << " (not given)";
            }
            os << "\n";
        }
    }
}

}
