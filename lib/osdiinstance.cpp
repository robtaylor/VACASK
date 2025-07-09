#include <cstring>
#include "osdiinstance.h"
#include "circuit.h"
#include "simulator.h"
#include "libplatform.h"
#include "common.h"


namespace NAMESPACE {

// TODO: check unconnected terminals handling, i.e. BJT with unconnected bulk
//       how are unconnected terminals handled

OsdiInstance::OsdiInstance(OsdiModel* model, Id name, Instance* parentInstance, const PTInstance& parsedInstance, Status &s) 
    : Instance(model, name, parentInstance, parsedInstance), core_(nullptr), connectedTerminalCount(0) {
    core_ = alignedAlloc(sizeof(std::max_align_t), model->device()->descriptor()->instance_size);
    memset(core_, 0, model->device()->descriptor()->instance_size);

    // Create nodes (terminals+internal nodes) list
    // Resize nodes vector
    nodes_.resize(staticNodeCount());
    
    // By default all terminals/nodes are unconnected
    auto n = staticNodeCount();
    for(TerminalIndex i=0; i<n; i++) {
        nodes_[i] = nullptr;
    }

    // Internal nodes are created after setup

    setFlags(Flags::IsValid);
}

OsdiInstance::~OsdiInstance() {
    // Free allocated values (strings and vectors)
    model()->device()->freeValues(model()->core(), core_);
    
    alignedFree(core_);
}

TerminalIndex OsdiInstance::staticNodeCount() const {
    return model()->device()->staticNodeCount();
}

TerminalIndex OsdiInstance::terminalCount() const {
    return model()->device()->terminalCount();
}

std::tuple<TerminalIndex, bool> OsdiInstance::nodeIndex(Id name) const {
    return model()->device()->nodeIndex(name);
}

Id OsdiInstance::nodeName(TerminalIndex ndx) const {
    return model()->device()->nodeName(ndx);
}

bool OsdiInstance::bindTerminal(TerminalIndex n, Node* node, Status& s) {
    if (n>=terminalCount()) {
        s.set(Status::Range, "Too many connections specified.");
        return false;
    }
    nodes_[n] = node;
    // Update connected terminal count, assume that there are no unconnected 
    // terminals between two connected terminals
    if (n+1>connectedTerminalCount) {
        connectedTerminalCount = n+1;
    }
    return true;
}

Node* OsdiInstance::terminal(TerminalIndex n, Status& s) const {
    if (n>=terminalCount()) {
        s.set(Status::Range, "Terminal not found.");
        return nullptr;
    }
    return nodes_[n];
}

bool OsdiInstance::unbindTerminals(Circuit& circuit, Status& s) {
    for(decltype(connectedTerminalCount) i=0; i<connectedTerminalCount; i++) {
        if (!circuit.releaseNode(nodes_[i], s)) {
            return false;
        }
    }
    return true;
}

ParameterIndex OsdiInstance::parameterCount() const {
    return model()->device()->instanceParameterCount();
}

std::tuple<ParameterIndex, bool> OsdiInstance::parameterIndex(Id name) const {
    auto [index, found] = model()->device()->instanceParameterIndex(name);
    if (!found) {
        return std::make_tuple(0, false);
    }
    return std::make_tuple(index, true);
}

Id OsdiInstance::parameterName(ParameterIndex ndx) const {
    if (ndx<parameterCount())
        return model()->device()->instanceParameterName(ndx);
    else
        return Id::none;
}

std::tuple<Value::Type,bool> OsdiInstance::parameterType(ParameterIndex ndx, Status& s) const {
    if (ndx>=model()->device()->instanceParameterCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return std::make_tuple(Value::Type::Int, false);
    }
    auto osdiId = model()->device()->instanceOsdiParameterId(ndx);
    return std::make_tuple(model()->device()->parameterType(osdiId), true);
}

bool OsdiInstance::getParameter(ParameterIndex ndx, Value& v, Status& s) const {
    if (ndx>=model()->device()->instanceParameterCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return false;
    }
    auto osdiId = model()->device()->instanceOsdiParameterId(ndx);
    return model()->device()->readParameter(osdiId, nullptr, core_, v);
}

std::tuple<bool,bool> OsdiInstance::setParameter(ParameterIndex ndx, const Value& v, Status& s) {
    if (ndx>=model()->device()->instanceParameterCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return std::make_tuple(false, false);
    }
    auto osdiId = model()->device()->instanceOsdiParameterId(ndx);
    
    // Write
    auto [ok, changed] = model()->device()->writeParameter(osdiId, model()->core(), core_, v, s);
    
    // Mark instance for parameter propagation
    if (changed) {
        setFlags(Instance::Flags::ParamsChanged);
    }

    return std::make_tuple(ok, changed);
}

std::tuple<bool,bool> OsdiInstance::parameterGiven(ParameterIndex ndx, Status& s) const {
    if (ndx>=model()->device()->instanceParameterCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return std::make_tuple(false, false);
    }
    auto osdiId = model()->device()->instanceOsdiParameterId(ndx);
    return model()->device()->parameterGiven(osdiId, nullptr, core_, s);
}

ParameterIndex OsdiInstance::opvarCount() const {
    return model()->device()->opvarCount();
}

std::tuple<ParameterIndex, bool> OsdiInstance::opvarIndex(Id name) const {
    return model()->device()->opvarIndex(name);
}

Id OsdiInstance::opvarName(ParameterIndex ndx) const {
    return model()->device()->opvarName(ndx);
}

std::tuple<Value::Type,bool> OsdiInstance::opvarType(ParameterIndex ndx, Status& s) const {
    if (ndx>=model()->device()->opvarCount()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return std::make_tuple(Value::Type::Int, false);
    }
    auto osdiId = model()->device()->opvarOsdiParameterId(ndx);
    return std::make_tuple(model()->device()->parameterType(osdiId), true);
}

bool OsdiInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const {
    if (ndx>=model()->device()->opvarCount()) {
        s.set(Status::Range, std::string("Opvar index id=")+std::to_string(ndx)+" out of range.");
        return false;
    }
    auto osdiId = model()->device()->instanceOsdiParameterId(ndx);
    return model()->device()->readParameter(osdiId, const_cast<void*>(model()->core()), core_, v, s);
}

std::tuple<bool, OutputSource> OsdiInstance::opvarOutputSource(ParameterIndex ndx) const {
    if (ndx>=model()->device()->opvarCount()) {
        return std::make_tuple(false, OutputSource());
    }
    auto [t, ok] = opvarType(ndx);
    if (!ok) {
        return std::make_tuple(false, OutputSource());
    }

    auto osdiId = model()->device()->opvarOsdiParameterId(ndx);
    switch (t) {
        case Value::Type::Int: {
            const int* ptr = model()->device()->parameterPtr<int>(osdiId, model()->core(), core());
            if (ptr) {
                return std::make_tuple(true, OutputSource(ptr));
            }
            break;
        }
        case Value::Type::Real: {
            const double* ptr = model()->device()->parameterPtr<double>(osdiId, model()->core(), core());
            if (ptr) {
                return std::make_tuple(true, OutputSource(ptr));
            }
            break;
        }
    }
    return std::make_tuple(false, OutputSource());
}

std::tuple<EquationIndex,EquationIndex> OsdiInstance::sourceExcitation(Circuit& circuit) const {
    if (model()->device()->isSource()) {
        if (model()->device()->isVoltageSource()) {
            // Excitation is the voltage, specified in the flow node equation on RHS
            // +excitation is positive in the sense of residual
            // Branch flow node (+excitation), 0 (-excitation)
            return std::make_tuple(nodes_[2]->unknownIndex(), 0);
        } else {
            // Excitation is the current specified on KCL equations of pos node (with +) and neg node (with -)
            // +excitation is positive in the sense of residual (outflowing currents represent positive residual)
            // Pos node (+excitation), neg node (-excitation)
            return std::make_tuple(nodes_[0]->unknownIndex(), nodes_[1]->unknownIndex());
        }
    }
    return std::make_tuple(0, 0);
}

std::tuple<UnknownIndex,UnknownIndex> OsdiInstance::sourceResponse(Circuit& circuit) const {
    if (model()->device()->isSource()) {
        if (model()->device()->isVoltageSource()) {
            // Response is -flow, obtained as 0 - flow node value
            // 0 (+response), flow node (-response)
            return std::make_tuple(0, nodes_[2]->unknownIndex());
        } else {
            // Response is the voltage between pos node and neg node
            // Pos node (+response), neg node (-response)
            return std::make_tuple(nodes_[0]->unknownIndex(), nodes_[1]->unknownIndex());
        }
    }
    return std::make_tuple(0, 0);
}

std::tuple<EquationIndex, EquationIndex> OsdiInstance::noiseExcitation(Circuit& cir, ParameterIndex ndx) const {
    auto [n1, n2] = model()->device()->noiseExcitation(ndx);
    // If n1, n2 is UINT32_MAX it refers to global ground
    auto e1 = n1==UINT32_MAX ? 0 : nodes_[n1]->unknownIndex();
    auto e2 = n2==UINT32_MAX ? 0 : nodes_[n2]->unknownIndex();
    return std::make_tuple(e1, e2);
}

bool OsdiInstance::loadNoise(Circuit& circuit, double freq, double* noiseDensity) { 
    model()->device()->descriptor()->load_noise(core(), model()->core(), freq, noiseDensity);
    return true;
}

std::tuple<bool, bool, bool> OsdiInstance::setup(Circuit& circuit, bool force, DeviceRequests* devReq, Status& s) {
    OsdiSimParas sp;
    const auto& opt = circuit.simulatorOptions().core();
    auto& internals = circuit.simulatorInternals();

    // Allocate tables on stack
    auto [ndbl, nchrptr ] = model()->device()->simParasSizes();
    double dblArray[ndbl];
    char* chrPtrArray[nchrptr];
    
    OsdiDevice::populateSimParas(sp, opt, internals, dblArray, chrPtrArray);
    // Verilog-A $temperature is in K, convert the value given by options (in C)
    auto retval = setupCore(circuit, sp, opt.temp+273.15, force, devReq, s);
    return retval;
}

std::tuple<bool, bool, bool> OsdiInstance::setupCore(Circuit& circuit, OsdiSimParas& sp, double temp, bool force, DeviceRequests* devReq, Status& s) {
    auto handle = OsdiCallbackHandle {
        .kind = 1, 
        .name = const_cast<char*>(name().c_str())
    };
    OsdiInitInfo initInfo;

    // Assume unknwons did not change
    bool unknownsChanged = false;

    // If there is no setup history, assume unknowns changed
    if (!checkFlags(Instance::Flags::HasSetupHistory)) {
        unknownsChanged = true;
    }

    if (force || checkFlags(Flags::NeedsSetup)) {
        // setup_instance() API functionm writes only true values to collapsed nodes pattern
        // Therefore we must clear the pattern before calling setup_instance(). 
        // As we go, we store the old patern for comparison. 

        // Can we check for change in unknowns
        bool canCheckUnknownsChange = checkFlags(Instance::Flags::HasSetupHistory);
        
        // If we are not going to check, set cpSize to 1 so that we do not allocate unneccessary space
        OsdiFile::OsdiCollapsedNodesIndex cpSize;
        if (canCheckUnknownsChange) {
            cpSize = collapsedNodesPatternSize();
        } else {
            // We don't want to allocate stack memory for nothing
            cpSize=1;
        }
        
        // Get collapsed nodes pattern array
        auto cpat = collapsedNodesPattern();

        // Allocate storage for old pattern on stack
        bool cpOldArray[cpSize];

        // Copy pattern to storage so we can check for changes. 
        // Clear it as we go so that setup_instance() can rewrite it
        if (canCheckUnknownsChange) {
            for(decltype(cpSize) i=0; i<cpSize; i++) {
                cpOldArray[i] = cpat[i];
                cpat[i] = false;
            }
        }
        
        // Call setup_instance(), interpret status
        model()->device()->descriptor()->setup_instance((void*)&handle, core(), model()->core(), temp, connectedTerminalCount, &sp, &initInfo);
        auto setupOk = model()->device()->processInitInfo(circuit, initInfo, "Instance", name(), devReq, s);
        if (!setupOk) {
            // The problem is big enough to abort simulation
            // Restore collapsed nodes pattern first
            if (canCheckUnknownsChange) {
                for(decltype(cpSize) i=0; i<cpSize; i++) {
                    cpat[i] = cpOldArray[i];
                }
            }
            return std::make_tuple(false, false, false);
        }

        // Setup succeeded 
        clearFlags(Flags::NeedsSetup);
        
        // If instance has setup history, compare new node collapse pattern to old pattern. 
        if (canCheckUnknownsChange) {
            // Assume pattern did not change
            unknownsChanged = false;
            for(decltype(cpSize) i=0; i<cpSize; i++) {
                if (cpOldArray[i] != cpat[i]) {
                    // Pattern changed
                    unknownsChanged = true;
                    break;
                }
            }
        }

        // After successful setup instance has setup history
        setFlags(Instance::Flags::HasSetupHistory);
    }

    // OSDI instances cannot change sparsity pattern
    return std::make_tuple(true, unknownsChanged, false);
}

bool OsdiInstance::propagateParameters(Circuit& circuit, RpnEvaluator& evaluator, Status& s) { 
    if (checkFlags(Flags::ParamsChanged)) {
        clearFlags(Flags::ParamsChanged);
        setFlags(Flags::NeedsSetup);
    }
    return true; 
};
    

bool OsdiInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, Status& s) {
    // Create internal nodes
    // Start at the node following the last conected terminal
    auto* descr = model()->device()->descriptor();
    auto n = staticNodeCount();
    for(TerminalIndex i=connectedTerminalCount; i<n; i++) {
        // Create node name
        Id internalNodeName = translate(nodeName(i)); 
        
        // Create/get node
        auto node = circuit.getNode(internalNodeName, descr->nodes[i].is_flow ? Node::Flags::FlowNode : Node::Flags::PotentialNode, s);
        if (node==nullptr) {
            s.extend(std::string("Failed to obtain internal node '"+std::string(internalNodeName)+"' from simulator."));
            s.extend(location());
            return false;
        }
        node->setFlags(Node::Flags::InternalDeviceNode);
        
        // Bind node
        nodes_[i] = node;
    }
    return true;
}

bool OsdiInstance::deleteHierarchy(Circuit& circuit, Status& s) {
    // Delete internal nodes
    auto n = staticNodeCount();
    for(TerminalIndex i=connectedTerminalCount; i<n; i++) {
        if (!circuit.releaseNode(nodes_[i], s)) {
            return false;
        }
        nodes_[i] = nullptr;
    }
    return true;
}

bool OsdiInstance::collapseNodesCore(Circuit& circuit, Status& s) {
    auto descr = model()->device()->descriptor();
    auto cpSize = collapsedNodesPatternSize();
    for(decltype(cpSize) i=0; i<cpSize; i++) {
        auto& nodePair = descr->collapsible[i];
        if (collapsedNodesPattern()[i]) {
            if (nodePair.node_2==UINT32_MAX) {
                // Collapse to ground
                if (!circuit.collapseNodes(nodes_[nodePair.node_1], nullptr, s)) {
                    return false;
                }
            } else {
                // Collapse two nodes
                if (!circuit.collapseNodes(nodes_[nodePair.node_1], nodes_[nodePair.node_2], s)) {
                    return false;
                }
            }
        }
    }
    return true;
}

bool OsdiInstance::populateStructuresCore(Circuit& circuit, Status& s) {
    auto descr = model()->device()->descriptor();
    auto numEntries = model()->device()->jacobianEntriesCount();
    for(decltype(numEntries) i=0; i<numEntries; i++) {
        auto& entry = model()->device()->jacobianEntry(i);
        auto ne = nodes_[entry.nodes.node_1];
        auto nu = nodes_[entry.nodes.node_2];
        
        EntryFlags f = EntryFlags::NoFlags;
        if (entry.flags & JACOBIAN_ENTRY_RESIST) {
            f = f | EntryFlags::Resistive;
        }
        if (entry.flags & JACOBIAN_ENTRY_REACT) {
            f = f | EntryFlags::Reactive;
        }

        if (auto [_, ok] = circuit.createJacobianEntry(ne, nu, f, s); !ok) {
            return false;
        }
    }

    // Reserve internal states, store indices
    auto internalStateCount = model()->device()->internalStateCount();
    auto stateIndices = stateIndexTable();

    // Allocate entries for internal states
    offsStates = circuit.allocateStates(internalStateCount);
    
    // Set internal state indices (local indices, uint32_t)
    // Because these indices are relative to the chunk allocated for the instance
    // internal state indices start with 0
    for(decltype(internalStateCount) i=0; i<internalStateCount; i++) {
        stateIndices[i] = i;
    }

    // Reserve reactive residual node states
    auto nodeStateCount = model()->device()->nonzeroReactiveResiduals().size()*2;
    circuit.allocateStates(nodeStateCount);
    // Total state vector chunk allocated for the instance 
    // has internalStateCount+nodeStateCount entries

    // Reserve storage for device convergence test
    auto deviceStateCount = descr->num_inputs + 
        // model()->device()->nonzeroResistiveJacobianEntries().size() + 
        // model()->device()->nonzeroReactiveJacobianEntries().size() + 
        model()->device()->nonzeroResistiveResiduals().size() + 
        model()->device()->nonzeroReactiveResiduals().size();
    offsDeviceStates = circuit.allocateDeviceStates(deviceStateCount);
    
    return true;
}

bool OsdiInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, const std::optional<MatrixEntryPosition>& mepResist, 
    KluMatrixAccess* matReact, Component compReact, const std::optional<MatrixEntryPosition>& mepReact, 
    Status& s
) {
    auto descr = model()->device()->descriptor();
    // Bind nodes
    auto nodeMapping = nodeMappingArray();
    auto numNodes = model()->device()->staticNodeCount();
    // 1-based unknown index, store it in nodeMapping
    for(decltype(numNodes) i=0; i<numNodes; i++) {
        nodeMapping[i] = nodes_[i]->unknownIndex();
    }

    // Bind Jacobian entries
    auto numEntries = model()->device()->jacobianEntriesCount();
    auto jacResistArray = resistiveJacobianPointers();
    for(decltype(numEntries) i=0; i<numEntries; i++) {
        // Position contains local terminal/node indices (0-based) 
        auto& entry = model()->device()->jacobianEntry(i);
        
        // Translate them to circuit nodes
        auto ne = nodes_[entry.nodes.node_1];
        auto nu = nodes_[entry.nodes.node_2];
        
        // Translate to equations/unknowns
        auto e = ne->unknownIndex();
        auto u = nu->unknownIndex();
        
        // Set resistive Jacobian element pointer
        if (matResist && !(jacResistArray[i] = matResist->valuePtr(MatrixEntryPosition(e, u), compResist, mepResist))) {
            s.set(Status::BadConversion, "Matrix is of incorrect type.");
            return false;
        }

        // Set reactive Jacobian element pointer
        if (matReact) {
            auto reactivePointer = reactiveJacobianPointer(i);
            if (reactivePointer) {
                // Set reactive Jacobian pointer
                if (!(*reactivePointer = matReact->valuePtr(MatrixEntryPosition(e, u), compReact, mepReact))) {
                    s.set(Status::BadConversion, "Matrix is of incorrect type.");
                    return false;
                }
            }
        }
    }
    return true;
}

void jacobianWriteSanityCheck(OsdiDescriptor* descriptor, void* model, void* instance, bool resist, bool react) {
    auto nnz = descriptor->num_jacobian_entries;
    auto nnzres = descriptor->num_resistive_jacobian_entries;
    auto nnzreact = descriptor->num_reactive_jacobian_entries;
    
    // Temporary pointer storage
    double* resPtrs[nnz];
    double* reacPtrs[nnz];
    
    // Extracted Jacobian contributions via load() 
    double resVals[nnz] = {0};
    double reacVals[nnz] = {0};

    // Extracted Jacobian values via write(), no initialization
    double reswVals[nnzres];
    double reacwVals[nnzreact];

    // Store pointers, redirect to our structures
    auto instResPtrs = getDataPtr<double**>(instance, descriptor->jacobian_ptr_resist_offset);
    for(decltype(nnz) i=0; i<nnz; i++) {
        auto& jac = descriptor->jacobian_entries[i];
        resPtrs[i] = instResPtrs[i];
        instResPtrs[i] = resVals+i;
        if (jac.react_ptr_off != UINT32_MAX) {
            auto reactPtr =  getDataPtr<double**>(instance, jac.react_ptr_off);
            reacPtrs[i] = *reactPtr;
            *reactPtr = reacVals+i;
        } else {
            reacPtrs[i] = nullptr;
        }
    }

    // Get Jacobian contributions
    if (resist) {
        descriptor->load_jacobian_resist(instance, model);
        descriptor->write_jacobian_array_resist(instance, model, reswVals);
    }
    if (react) {
        descriptor->load_jacobian_react(instance, model, 1.0);
        descriptor->write_jacobian_array_react(instance, model, reacwVals);
    }

    // Check values
    auto resPtr = reswVals;
    auto reactPtr = reacwVals;
    for(decltype(nnz) i=0; i<nnz; i++) {
        auto& jac = descriptor->jacobian_entries[i];
        if (jac.flags & JACOBIAN_ENTRY_RESIST) {
            // Have resistive entry
            if (resist) {
                if (resVals[i] != *resPtr) {
                    throw std::logic_error("Jacobian mismatch, resistive i="+std::to_string(i));
                }
                resPtr++;
            }
        }
        if (jac.flags & JACOBIAN_ENTRY_REACT) {
            // Have reactive entry
            if (react) {
                if (reacVals[i] != *reactPtr) {
                    throw std::logic_error("Jacobian mismatch, reactive i="+std::to_string(i));
                }
                reactPtr++;
            }
        }
    }

    // Restore pointers
    for(decltype(nnz) i=0; i<nnz; i++) {
        auto& jac = descriptor->jacobian_entries[i];
        instResPtrs[i] = resPtrs[i];
        if (jac.react_ptr_off != UINT32_MAX) {
            auto reactPtr =  getDataPtr<double**>(instance, jac.react_ptr_off);
            *reactPtr = reacPtrs[i];
        }
    }
}

void jacobianLoadWithOffsetSanityCheck(OsdiDescriptor* descriptor, void* model, void* instance, bool resist, bool react) {
    auto nnz = descriptor->num_jacobian_entries;
    
    // Temporary pointer storage
    double* resPtrs[nnz];
    double* reacPtrs[nnz];
    
    // Extracted Jacobian contributions via load() 
    double resVals[nnz] = {0};
    double reacVals[nnz] = {0};

    // Extracted Jacobian values via load_with_offset()
    // 4 slots for each Jacobian entry, initialize to 0
    double reswVals[nnz*4] = {0};
    double reacwVals[nnz*4] = {0};

    // Store pointers, redirect to our resVals, reacVals
    auto instResPtrs = getDataPtr<double**>(instance, descriptor->jacobian_ptr_resist_offset);
    for(decltype(nnz) i=0; i<nnz; i++) {
        auto& jac = descriptor->jacobian_entries[i];
        resPtrs[i] = instResPtrs[i];
        instResPtrs[i] = resVals+i;
        if (jac.react_ptr_off != UINT32_MAX) {
            auto reactPtr =  getDataPtr<double**>(instance, jac.react_ptr_off);
            reacPtrs[i] = *reactPtr;
            *reactPtr = reacVals+i;
        } else {
            reacPtrs[i] = nullptr;
        }
    }

    // Get Jacobian contributions
    if (resist) {
        descriptor->load_jacobian_resist(instance, model);
    }
    if (react) {
        descriptor->load_jacobian_react(instance, model, 1.0);
    }

    // Redirect to reswVals, reacwVals
    for(decltype(nnz) i=0; i<nnz; i++) {
        auto& jac = descriptor->jacobian_entries[i];
        instResPtrs[i] = reswVals+4*i;
        if (jac.react_ptr_off != UINT32_MAX) {
            auto reactPtr =  getDataPtr<double**>(instance, jac.react_ptr_off);
            *reactPtr = reacwVals+4*i;
        }
    }

    // Get Jacobian contributions with offset
    if (resist) {
        // Fill slot 2
        descriptor->load_jacobian_with_offset_resist(instance, model, 2);
    }
    if (react) {
        // Fill slot 1
        descriptor->load_jacobian_with_offset_react(instance, model, 1);
    }

    // Check values
    for(decltype(nnz) i=0; i<nnz; i++) {
        auto& jac = descriptor->jacobian_entries[i];
        // Have resistive entry
        if (resist) {
            if (resVals[i] != reswVals[4*i+2]) {
                throw std::logic_error("Jacobian mismatch, resistive i="+std::to_string(i));
            }
            if (reswVals[4*i]!=0 || reswVals[4*i+1]!=0 || reswVals[4*i+3]!=0) {
                throw std::logic_error("Write outside slot, resistive i="+std::to_string(i));
            }
        }
        // Have reactive entry
        if (react) {
            if (reacVals[i] != reacwVals[4*i+1]) {
                throw std::logic_error("Jacobian mismatch, reactive i="+std::to_string(i));
            }
            if (reacwVals[4*i]!=0 || reacwVals[4*i+2]!=0 || reacwVals[4*i+3]!=0) {
                throw std::logic_error("Write outside slot, reactive i="+std::to_string(i));
            }
        }
    }

    // Restore pointers
    for(decltype(nnz) i=0; i<nnz; i++) {
        auto& jac = descriptor->jacobian_entries[i];
        instResPtrs[i] = resPtrs[i];
        if (jac.react_ptr_off != UINT32_MAX) {
            auto reactPtr =  getDataPtr<double**>(instance, jac.react_ptr_off);
            *reactPtr = reacPtrs[i];
        }
    }
}

bool OsdiInstance::inputBypassCheckCore(Circuit& circuit, EvalSetup& evalSetup) {
    // Get descriptor 
    auto model_ = model();
    auto device = model_->device();
    auto descr = device->descriptor();

    // Get options
    auto& options = circuit.simulatorOptions().core();
    auto bypasstol = options.nr_bypasstol;

    // Get stored inputs
    auto input1 = evalSetup.deviceStates+offsDeviceStates;

    // Get inputs table
    auto tab = descr->inputs;
    auto nInputs = descr->num_inputs;
    bool bypass = true;

    // We must compare input change between stored input and old_solution. 
    // This is the last delta computed by the NR algorithm. 
    // because there is no no delta available between two consecutive NR 
    // runs we better compute this delta each time anew. 
    for(decltype(nInputs) i=0; i<nInputs; i++) {
        // Get column indices within instance
        auto col1 = tab[i].node_1;
        auto col2 = tab[i].node_2;

        // Get nodes
        auto node1 = col1==UINT32_MAX ? nullptr : nodes_[col1];
        auto node2 = col2==UINT32_MAX ? nullptr : nodes_[col2];

        // Both nodes are ground, nothing to do
        if (!node1 && !node2) {
            continue;
        }
        
        // Compute input
        double v1 = 0;
        double abstol = std::numeric_limits<double>::infinity();
        if (node1) {
            auto u1 = node1->unknownIndex();
            v1 = evalSetup.oldSolution[u1];
            abstol = std::min(abstol, node1->checkFlags(Node::Flags::FlowNode) ? options.abstol : options.vntol);
        }
        if (node2) {
            auto u2 = node2->unknownIndex();
            v1 -= evalSetup.oldSolution[u2];
            abstol = std::min(abstol, node2->checkFlags(Node::Flags::FlowNode) ? options.abstol : options.vntol);
        }
        
        // Get stored input
        auto v0 = input1[i];

        // Get delta
        auto delta = v1-v0;

        // Compute reference value and tolerance
        auto ref = std::abs(v0);
        auto tol = std::max(ref*options.reltol, abstol);

        // Check
        if (std::abs(delta)>tol*bypasstol) {
            // Outside tolerances, no bypass
            bypass = false;
            // Can stop checking 
            break;
        }
    }

    return bypass;
}

bool OsdiInstance::outputBypassCheckCore(Circuit& circuit, EvalSetup& evalSetup) {
    // If deviceStates is nullptr, we are done
    if (!evalSetup.deviceStates) {
        return false;
    }

    // If limiting was applied, we are not allowed to bypass. 
    // Clear HasDeviceHistory flag because next time we get here 
    // there will be no history available. 
    if (checkFlags(Flags::LimitingApplied)) {
        // Clear HasDeviceHistory flag, not converged
        clearFlags(Flags::HasDeviceHistory);
        return false;
    }

    // Get descriptor 
    auto model_ = model();
    auto device = model_->device();
    auto descr = device->descriptor();

    // Assume converged if requestHighPrecision==false. 
    // If requestHighPrecision==true, assume not converged. 
    // In the latter case we will only store history 
    // without checking for convergence.
    bool converged = !evalSetup.requestHighPrecision;

    // If no device history is available assume not converged. 
    // Only update history. 
    if (!checkFlags(Flags::HasDeviceHistory)) {
        converged = false;
    }

    // As soon as converged becomes false, convergence checks are skipped, 
    // history is updated regardless of what the value of converged is. 
    
    // Get options
    auto& options = circuit.simulatorOptions().core();
    auto convtol = options.nr_convtol;

    // Nonzero Jacobian indices
    // auto& Jrnz = device->nonzeroResistiveJacobianEntries();
    // Nonzero residual indices
    auto& fnz = device->nonzeroResistiveResiduals();

    // Nonzero Jacobian indices
    // auto& Jcnz = device->nonzeroReactiveJacobianEntries();
    // Nonzero residual indices
    auto& qnz = device->nonzeroReactiveResiduals();

    // Previous input values
    auto input1 = evalSetup.deviceStates+offsDeviceStates;
    
    // Previous values of Jacobian contributions
    // auto Jr1 = input1 + descr->num_inputs;
    // auto Jc1 = Jr1 + Jrnz.size();
    // Previous values of residual contributions
    auto f1 = input1 + descr->num_inputs; // Jr1 + Jrnz.size();
    auto q1 = f1 + fnz.size(); // Jc1 + Jcnz.size();
    
    //
    // Check resistive residual and resistive Jacobian
    // 

    // Allocate arrays with new values (on stack)
    // (Jacobian only because we need to read those from core instance)
    // auto nn = nodes_.size();
    
    // Index of entry in deviceStates subvectors
    size_t i;

    // Loop through resistive residuals, check convergence, write to array of previous values
    // double rres[nn];
    i = 0;
    for(auto resNdx : fnz) {
        // Get node
        auto node = nodes_[resNdx];
        // Get unknown index
        auto u = node->unknownIndex();

        // Get new residual
        // If we arrive here limiting was not applied, no need to add linearized residual
        auto resOff = descr->nodes[resNdx].resist_residual_off;
        auto res = *getDataPtr<double*>(core(), resOff);
        // rres[resNdx] = res;
            
        // Skip this if we already know instance is not converged
        // Skip if this residual contributes to ground node
        if (converged && u!=0) {
            // Get absolute tolerance (residual tolerance)
            // flow nodes -> voltage
            // potential nodes -> current
            auto abstol = node->checkFlags(Node::Flags::FlowNode) ? options.vntol : options.abstol;
            // Compute tolerance
            auto tol = std::max(std::abs(f1[i])*options.reltol, abstol);
            // Compare
            if (std::abs(f1[i]-res) > tol*convtol) {
                // Simulator::dbg() << name() << " : f tol violated @ " << node->name() << 
                //    " " << res << " " << (f1[i]-res) << " " << tol << "\n";
                converged = false;
            }
        }

        // Store in previous values array regardless if tolernce was checked or not
        f1[i] = res;
        
        i++; 
    }

    /*
    // Old solution
    auto x1 = convSetup.oldSolution;

    // Read resistive Jacobian entries
    double Jr2[Jrnz.size()];
    descr->write_jacobian_array_resist(core(), model_->core(), Jr2);
    
    // Loop through resistive Jacobian entries, check convergence, write to array of previous values
    i = 0;
    for(auto jacNdx : Jrnz) {
        // Skip this if we already know instance is not converged
        if (converged) {
            // Get Jacobian entry
            auto& jacEnt = descr->jacobian_entries[jacNdx];

            // Get row and column index within instance
            auto row = jacEnt.nodes.node_1;
            auto col = jacEnt.nodes.node_2;

            auto rowNode = nodes_[row];
            auto colNode = nodes_[col];

            // Get equation number and unknown index
            auto eq = rowNode->unknownIndex();
            auto u = colNode->unknownIndex();

            // Skip check if row or column correspond to ground node
            if (eq!=0 && u!=0) {
                // Get absolute tolerance (resistive residual tolerance of equation)
                auto abstol = rowNode->checkFlags(Node::Flags::FlowNode) ? options.vntol : options.abstol;

                // Reference value is the residual, we already computed it
                auto res = rres[row];

                // Compute tolerance
                auto tol = std::max(std::abs(res)*options.reltol, abstol);

                // Compare
                if (std::abs((Jr2[i]-Jr1[i])*x1[u])>tol*convtol) {
                    // Simulator::dbg() << name() << " : Jr tol violated @ " << 
                    //    "(" << nodes_[row]->name() << ", " << nodes_[col]->name() << ")\n";
                    converged = false;
                }
            }
        }

        // Store in previous values array
        Jr1[i] = Jr2[i];
        i++; 
    }
    */

    // 
    // Check reactive residual and Jacobian
    //

    // Loop through reactive residuals, check convergence, write to array of previous values
    if (evalSetup.checkReactiveConvergece) {
        // double rreac[nn];

        i = 0;
        for(auto resNdx : qnz) {
            // Get node
            auto node = nodes_[resNdx];
            // Get unknown index
            auto u = node->unknownIndex();
            
            // Get new residual
            // If we arrive here limiting was not applied, no need to add linearized residual
            auto resOff = descr->nodes[resNdx].react_residual_off;
            auto res = *getDataPtr<double*>(core(), resOff);
            // rreac[resNdx] = res;

            // Skip this if we already know instance is not converged
            // Skip if this residual contributes to ground node
            if (converged && u!=0) {
                // Get absolute tolerance (idt residual tolerance)
                // potential nodes -> charge
                // flow nodes -> flux
                auto abstol = node->checkFlags(Node::Flags::FlowNode) ? options.fluxtol : options.chgtol;
                // Compute tolerance
                auto tol = std::max(std::abs(q1[i])*options.reltol, abstol);
                // Compare
                if (std::abs(q1[i]-res) > tol*convtol) {
                    converged = false;
                }
            }

            // Store in previous values array regardless if tolernce was checked or not
            q1[i] = res;
            
            i++;
        }
        
        /*
        // Read reactive Jacobian entries
        double Jc2[Jcnz.size()];
        descr->write_jacobian_array_react(core(), model_->core(), Jc2);
        
        // Loop through resistive Jacobian entries, check convergence, write to array of previous values
        i = 0;
        for(auto jacNdx : Jcnz) {
            // Skip this if we already know instance is not converged
            if (converged) {
                // Get Jacobian entry
                auto& jacEnt = descr->jacobian_entries[jacNdx];

                // Get row and column index within instance
                auto row = jacEnt.nodes.node_1;
                auto col = jacEnt.nodes.node_2;

                auto rowNode = nodes_[row];
                auto colNode = nodes_[col];

                // Get equation number and unknown index
                auto eq = rowNode->unknownIndex();
                auto u = colNode->unknownIndex();
                
                // Skip check if row or column correspond to ground node
                if (eq!=0 && u!=0) {
                    // Get absolute tolerance (resistive residual idt tolerance of equation)
                    auto abstol = rowNode->checkFlags(Node::Flags::FlowNode) ? options.fluxtol : options.chgtol;

                    // Reference value is the residual, we already computed it
                    auto res = rreac[row];

                    // Compute tolerance
                    auto tol = std::max(std::abs(res)*options.reltol, abstol);

                    // Compare
                    if (std::abs((Jc2[i]-Jc1[i])*x1[u])>tol*convtol) {
                        converged = false;
                    }
                }
            }

            // Store in previous values array
            Jc1[i] = Jc2[i];
            i++;
        }
        */
    }

    //
    // Compute input values from old_solution, but only if outputs are converged. 
    // We need input history only if we are going to check the input change for 
    // convergence. 
    if (converged) {
        auto tab = descr->inputs;
        auto nInputs = descr->num_inputs;
        for(decltype(nInputs) i=0; i<nInputs; i++) {
            // Get column indices within instance
            auto col1 = tab[i].node_1;
            auto col2 = tab[i].node_2;

            // Get nodes
            auto node1 = col1==UINT32_MAX ? nullptr : nodes_[col1];
            auto node2 = col2==UINT32_MAX ? nullptr : nodes_[col2];
            
            // Compute input value
            double vin = 0;
            if (node1) {
                auto u1 = node1->unknownIndex();
                vin = evalSetup.oldSolution[u1];
            }
            if (node2) {
                auto u2 = node2->unknownIndex();
                vin -= evalSetup.oldSolution[u2];
            }
            
            // Store input value
            input1[i] = vin;
        }
    }

    // Count converged instances
    evalSetup.instancesConvergenceChecks++;
    if (converged) {
        evalSetup.convergedInstances++; 
    }

    // Device has (output) history now
    // Input history is updated only if outputs are converged. 
    setFlags(Flags::HasDeviceHistory);

    return converged;
}

bool OsdiInstance::evalCore(Circuit& circuit, OsdiSimInfo& simInfo, EvalSetup& evalSetup) {
    // Get descriptor 
    auto model_ = model();
    auto device = model_->device();
    auto descr = device->descriptor();

    Accounting::Timepoint td0;
    
    // Prepare callback handle
    OsdiCallbackHandle handle = OsdiCallbackHandle {
        .kind = 3, 
        .name = const_cast<char*>(name().c_str())
    };

    // Update simInfo - for now no need to do this because 
    // evalSetup->initializeLimiting does not change 
    // during OsdiDevice::evalAndLoad(). 
    // updateSimInfo() only updates the OsdiSimParas entry
    // that is read when the model requests $simparam("iniLim"). 
    // and is based on evalSetup->initializeLimiting. 
    // OsdiSimParas entry was already set in populateSimParas() which 
    // is called at the beginning of OsdiDevice::evalAndLoad(). 
    // OsdiDevice::updateSimInfo(simInfo, evalSetup);

    // Set beginning of state vector chunk belonging to this instance
    simInfo.prev_state = evalSetup.oldStates + offsStates;
    simInfo.next_state = evalSetup.newStates + offsStates;

    // Get bound step offset
    auto bsOffs = descr->bound_step_offset;

    // TODO: when timestep is very small elements are bypassed when they should not be bypassed. 
    //       This is probably because changes between two consecutive timepoints are 
    //       small and are mistakenly identified as small enough for instance convergence. 
    //       A similar problem happens in inputBypassCheckCore(), but this time with 
    //       the components of the circuit's solution. 
    //       Must do something about this. 
    //       This bug does not affect nr_acctbypass and sweep_innerbypass where the 
    //       change between two consecutive circuit solutions is 0. 
    
    // Inactive element bypass
    bool nr_bypass = circuit.simulatorOptions().core().nr_bypass;

    // Check if device is bypassable and high precision is not requested
    bool bypass = false;
    evalSetup.bypassableInstances++;

    if (!device->checkFlags(Device::Flags::Bypassable)) {
        // Not bypassable, nothing to do
    } else if (evalSetup.requestHighPrecision) {
        // If high precision is required, bypass is out of question
        bypass = false;
        // If device output is converged, this is a bypass opportunity that was not taken
        if (checkFlags(Flags::OutputConverged)) {
            evalSetup.bypassOpportunuties++;
        }
        // Not bypassed, output no longer converged
        clearFlags(Flags::Bypassed|Flags::OutputConverged);
    } else if (evalSetup.forceBypass) {
        // Forcing a bypass (continuation bypass)
        bypass = true;
        // This is a bypass opportunity that was taken
        evalSetup.bypassOpportunuties++;
        // Bypass regardless of converged flag, do not change converged flag
        setFlags(Flags::Bypassed);
    } else if (evalSetup.allowBypass) {
        // Bypass is allowed (inactive element bypass)
        if (checkFlags(Flags::OutputConverged)) {
            // Device output is converged, this is a bypass opportunity
            evalSetup.bypassOpportunuties++;
            // Converged, check if we can bypass
            if (inputBypassCheckCore(circuit, evalSetup)) {
                // Bypassing
                bypass = true;
                setFlags(Flags::Bypassed);
            } else {
                // Not bypassing, leave converged mode
                bypass = false;
                clearFlags(Flags::Bypassed|Flags::OutputConverged);
            }
        } else {
            // Output not converged, not bypassing 
            bypass = false;
            // This is not a bypass opportunity
            // Clear Bypassed flag
            clearFlags(Flags::Bypassed);
        }
    } else {
        // Bypass not allowed, clear output converged flag
        bypass = false;
        clearFlags(Flags::Bypassed|Flags::OutputConverged);
    }

    // Evaluation, skip it if bypass is true
    if (bypass) {
        evalSetup.bypassedInstances++;
        // The reactive residuals will be stored in the states even if instance is bypassed. 
        // The stored value will be taken from the last instance evaluation. 
        // But we must make sure we copy the states of the limiting functions to new states.  
        // These values are normally written by the device evaluation function which we now bypass. 
        auto isCount = device->internalStateCount();
        if (isCount>0) {
            auto nodeStateIndex = offsStates;
            for(decltype(isCount) i=0; i<isCount; i++) {
                evalSetup.newStates[nodeStateIndex] = evalSetup.oldStates[nodeStateIndex];
                nodeStateIndex++;
            }
        }
    } else {
        if constexpr(devacct) {
            td0 = Accounting::wclk();
        }
        
        // Core evaluation, need to call it always to compute bound step
        auto evalFlags = descr->eval(&handle, core(), model_->core(), &simInfo);
        
        // if constexpr(devacct) {
        //     device->tovh += Accounting::wclkDelta(td0);
        // }

        // Handle evalFlags
        if (evalFlags & EVAL_RET_FLAG_LIM) {
            // If some variable x is linearized to xl the Jacobian is computed at xl instead of x
            // i.e. limiting takes place and EVAL_RET_FLAG_LIM is set
            // We are may not stop the NR loop until this flag is gone for all instances. 
            evalSetup.limitingApplied = true;
            setFlags(Flags::LimitingApplied);
        } else {
            evalSetup.limitingApplied = false;
            clearFlags(Flags::LimitingApplied);
        }

        if (evalFlags & EVAL_RET_FLAG_FATAL) {
            // Fatal error occurred, must abort simulation. 
            evalSetup.requests.abort = true;
        } 
        if (evalFlags & EVAL_RET_FLAG_FINISH) {
            // $finish was called asking the simulator to finish simulation 
            // (exit gracefully) if the current iteration converged. 
            evalSetup.requests.finish = true;
        } 
        if (evalFlags & EVAL_RET_FLAG_STOP) {
            // $stop was called asking the simulator to pause the simulation
            // if the current iteration converged. 
            evalSetup.requests.stop = true;
        }
    }

    // if constexpr(devacct) {
    //     td0 = Accounting::wclk();
    // }
    
    auto nodeStateIndex = offsStates + device->internalStateCount();
    if (evalSetup.integCoeffs || evalSetup.storeReactiveState) {
        for(auto i : device->nonzeroReactiveResiduals()) {
            // Get unknown index
            auto u = nodes_[i]->unknownIndex();
            // For ground nodes, skip the rest
            if (u!=0) {
                // Get residual offset and residual
                auto resOff = descr->nodes[i].react_residual_off;
                auto contrib = *getDataPtr<double*>(core(), resOff);
                // Add linearized residual 
                // (note that the negative of the linearized residual is stored)
                if (checkFlags(Flags::LimitingApplied)) {
                    auto offsLim = descr->nodes[i].react_limit_rhs_off;
                    if (offsLim!=UINT32_MAX) {
                        contrib -= *getDataPtr<double*>(core(), offsLim);
                    }
                }
                // Store reactive residual contribution in states vector
                evalSetup.newStates[nodeStateIndex] = contrib;
                // Differentiate if requested
                if (evalSetup.integCoeffs) {
                    double flow = evalSetup.integCoeffs->differentiate(contrib, nodeStateIndex);
                    // Store in states vector
                    evalSetup.newStates[nodeStateIndex+1] = flow;
                }
            }
            // Go to next node (skip two state vector entries - charge and flow)
            nodeStateIndex += 2;
        }
    }

    /*
    // Residual debugging
    Simulator::dbg() << "  Instance " << name() << "\n";
    for(uint32_t i=0; i<descr->num_nodes; i++) {
        auto resOff = descr->nodes[i].resist_residual_off;
        auto resLimOff = descr->nodes[i].resist_limit_rhs_off;
        auto reactOff = descr->nodes[i].react_residual_off;
        auto reactLimOff = descr->nodes[i].react_limit_rhs_off;
        
        Simulator::dbg() << "    " << i << " " << descr->nodes[i].name << ": offsets " << "res " << resOff << " reslim " << resLimOff;
        Simulator::dbg() << " react " << reactOff << " reactlim " << reactLimOff << " : ";
        Simulator::dbg()  << "resistive=";
        if (resOff!=UINT32_MAX) {
            Simulator::dbg()  << *getDataPtr<double*>(core(), resOff);
        }
        if (resLimOff!=UINT32_MAX) {
            Simulator::dbg()  << " - " << *getDataPtr<double*>(core(), resLimOff);
        }
        Simulator::dbg()  << ", reactive=";
        if (reactOff!=UINT32_MAX) {
            Simulator::dbg()  << *getDataPtr<double*>(core(), reactOff);
        }
        if (reactLimOff!=UINT32_MAX) {
            Simulator::dbg()  << " - " << *getDataPtr<double*>(core(), reactLimOff);
        }
        Simulator::dbg() << "\n"; 
    }
    */
    
    if (evalSetup.computeBoundStep) {
        if (bsOffs!=UINT32_MAX) {
            evalSetup.setBoundStep(*getDataPtr<double*>(core(), bsOffs));
        }
    }

    // Convergence check of residuals, update history. 
    // The following must hold to do this:
    // - nr_bypass enabled
    // - device is bypassable
    // - device is not bypassed
    //   - evalSetup.forceBypass is disabled and
    //   - device is not bypassed due to inactive element bypass (Flag::Bypassed)
    if (nr_bypass && device->checkFlags(Device::Flags::Bypassable) && !evalSetup.forceBypass && !checkFlags(Flags::Bypassed)) {
        auto converged = outputBypassCheckCore(circuit, evalSetup);
        if (converged) {
            // Mark outputs as converged
            setFlags(Flags::OutputConverged);
        } else {
            // Mark outputs as not converged
            clearFlags(Flags::OutputConverged);
        }
    }

    return true;
}

bool OsdiInstance::loadCore(Circuit& circuit, LoadSetup& loadSetup) {
    // Get descriptor
    auto model_ = model();
    auto device = model_->device();
    auto descr = device->descriptor();

    Accounting::Timepoint td0;
    if constexpr(devacct) {
        td0 = Accounting::wclk();
    }

    /*
    // For development
    jacobianWriteSanityCheck(
        model()->device()->descriptor(), model()->core(), core(), 
        loadSetup.loadResistiveJacobian || loadSetup.loadTransientJacobian, 
        loadSetup.loadReactiveJacobian || loadSetup.loadTransientJacobian
    );
    */

    /*
    // For development
    jacobianLoadWithOffsetSanityCheck(
        model()->device()->descriptor(), model()->core(), core(), 
        loadSetup.loadResistiveJacobian || loadSetup.loadTransientJacobian, 
        loadSetup.loadReactiveJacobian || loadSetup.loadTransientJacobian
    );
    */
    
    // Loading
    
    // Load Jacobian computed with limiting (if it was computed)
    if (loadSetup.loadResistiveJacobian) {
        if (loadSetup.jacobianLoadOffset==0) {
            descr->load_jacobian_resist(core(), model_->core());
        } else {
            descr->load_jacobian_with_offset_resist(core(), model_->core(), loadSetup.jacobianLoadOffset);
        }
    }
    if (loadSetup.loadReactiveJacobian) {
        if (loadSetup.jacobianLoadOffset==0) {
            descr->load_jacobian_react(core(), model_->core(), loadSetup.reactiveJacobianFactor);
        } else {
            descr->load_jacobian_with_offset_react(core(), model_->core(), loadSetup.jacobianLoadOffset); 
        }
    }
    if (loadSetup.loadTransientJacobian) {
        if (loadSetup.jacobianLoadOffset==0) {
            descr->load_jacobian_tran(core(), model_->core(), loadSetup.integCoeffs->leadingCoeff());
        } else {
            // Not supported, needs support in OpenVAF
            // At this point it is not needed. 
            return false;
        }
    }

    // Without limiting the residual is 
    //   g(x)
    // With limiting it is linearized above xl so it is actually
    //   g(xl) + Jg(xl) (x-xl)
    // g(xl) part is loaded by 
    //   load_residual_resist()
    // Jg(xl) (xl-x) = -Jg(xl) (x-xl) is the linearized rhs residual part and is loaded by 
    //   load_lim_rhs_resist()
    // The second one needs to be subtracted from the first one to get the 
    // actual residual, i.e.
    //   actual_residual_when_limiting = residual - linearized_rhs_residual
    // To maintain cache locality load_residual_*() and load_lim_rhs_*() must be 
    // called close together for the same instance. 
    
    // Load residual
    if (loadSetup.resistiveResidual) {
        descr->load_residual_resist(core(), model_->core(), loadSetup.resistiveResidual);
    }
    if (loadSetup.reactiveResidual) {
        descr->load_residual_react(core(), model_->core(), loadSetup.reactiveResidual);
    }

    // Load limited residual only if limiting was applied
    if (checkFlags(Flags::LimitingApplied)) {
        if (loadSetup.linearizedResistiveRhsResidual) {
            descr->load_limit_rhs_resist(core(), model_->core(), loadSetup.linearizedResistiveRhsResidual);
        }
        if (loadSetup.linearizedReactiveRhsResidual) {
            descr->load_limit_rhs_react(core(), model_->core(), loadSetup.linearizedReactiveRhsResidual); 
        }
    }

    // Update maximal resistive residual contribution
    if (loadSetup.maxResistiveResidualContribution) {
        for(auto i : model()->device()->nonzeroResistiveResiduals()) {
        // for(uint32_t i=0; i<descr->num_nodes; i++) {
            // Get unknown index
            auto u = nodes_[i]->unknownIndex();
            // For ground nodes, skip the rest
            if (u!=0) {
                // Get residual offset and residual
                auto offs = descr->nodes[i].resist_residual_off;
                double contrib = *getDataPtr<double*>(core(), offs);
                // Add linearized residual
                // (note that the negative of the linearized residual is stored)
                if (checkFlags(Flags::LimitingApplied)) {
                    auto offsLim = descr->nodes[i].resist_limit_rhs_off;
                    if (offsLim!=UINT32_MAX) {
                        contrib -= *getDataPtr<double*>(core(), offsLim);
                    }
                }
                // Update
                contrib = std::abs(contrib);
                if (contrib > loadSetup.maxResistiveResidualContribution[u]) {
                    loadSetup.maxResistiveResidualContribution[u] = contrib;
                }
            }
        }
    }
    
    // Update maximal reactive residual contribution and its derivative
    auto nodeStateIndex = offsStates + device->internalStateCount();
    if (
        loadSetup.reactiveResidualDerivative || 
        loadSetup.maxReactiveResidualContribution || 
        loadSetup.maxReactiveResidualDerivativeContribution
    ) {
        for(auto i : device->nonzeroReactiveResiduals()) {
            // Get unknown index
            auto u = nodes_[i]->unknownIndex();
            // For ground nodes, skip the rest
            if (u!=0) {
                // No need to compute it, retrieve it from states
                // Also retrieve derivative wrt time
                auto contrib = loadSetup.newStates[nodeStateIndex];
                auto flow = loadSetup.newStates[nodeStateIndex+1];
                // Store in rhs
                if (loadSetup.reactiveResidualDerivative) {
                    loadSetup.reactiveResidualDerivative[u] += flow;
                }
                // Update max
                if (loadSetup.maxReactiveResidualContribution) {
                    contrib = std::abs(contrib);
                    if (contrib > loadSetup.maxReactiveResidualContribution[u]) {
                        loadSetup.maxReactiveResidualContribution[u] = contrib;
                    }
                }
                if (loadSetup.maxReactiveResidualDerivativeContribution) {
                    flow = std::abs(flow);
                    if (flow > loadSetup.maxReactiveResidualDerivativeContribution[u]) {
                        loadSetup.maxReactiveResidualDerivativeContribution[u] = flow;
                    }
                }
            }
            // Go to next node
            nodeStateIndex += 2;
        }
    }

    return true;
}

void OsdiInstance::dump(int indent, const Circuit& circuit, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "Osdi device instance " << std::string(name()) << " of model " << model()->name() << "\n";
    if (terminalCount()>0) {
        os << pfx << "  Terminals: ";
        auto termCount = terminalCount();
        for(decltype(termCount) i=0; i<termCount; i++) {
            os << terminal(i)->name() << " ";
        }
        os << "\n";
    }
    auto desc = model()->device()->descriptor();
    if (desc->num_collapsible>0) {
        auto pattern = collapsedNodesPattern();
        auto n = collapsedNodesPatternSize();
        bool have = false;
        for(decltype(n) i=0; i<n; i++) {
            if (pattern[i]) {
                have = true;
                break;
            }
        }
        if (have) {
            os << pfx << "  Collapsed node pairs:\n";
            for(ParameterIndex i=0; i<desc->num_collapsible; i++) {
                if  (!pattern[i]) {
                    continue;
                }
                auto n1 = desc->collapsible[i].node_1;
                auto n2 = desc->collapsible[i].node_2;
                os << pfx << "    " << desc->nodes[n1].name;
                if (n2!=UINT32_MAX) {
                    os << ", " << desc->nodes[n2].name;
                } else {
                    os << ", (ground)";
                }
                os << "\n";
            }
        }
    }
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
