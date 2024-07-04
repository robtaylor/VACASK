#include <numbers>
#include "devctlsrc.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

template<> int Introspection<DevCtlSourceModelParams>::setup() {
    return 0;
}
instantiateIntrospection(DevCtlSourceModelParams);

DevCtlSourceModelParams::DevCtlSourceModelParams() {
}


template<> int Introspection<DevVctlSourceInstanceParams>::setup() {
    registerNamedMember(mfactor, "$mfactor");
    registerMember(gain);
    return 0;
}
instantiateIntrospection(DevVctlSourceInstanceParams);

DevVctlSourceInstanceParams::DevVctlSourceInstanceParams() {
}


template<> int Introspection<DevCctlSourceInstanceParams>::setup() {
    registerNamedMember(mfactor, "$mfactor");
    registerMember(gain);
    registerMember(ctlinst);
    registerMember(ctlnode);
    return 0;
}
instantiateIntrospection(DevCctlSourceInstanceParams);

DevCctlSourceInstanceParams::DevCctlSourceInstanceParams() {
}


template<> int Introspection<DevVccsInstanceData>::setup() {
    registerMember(ctl);
    registerMember(v);
    registerMember(i);
    return 0;
}
instantiateIntrospection(DevVccsInstanceData);

DevVccsInstanceData::DevVccsInstanceData() {
}


template<> int Introspection<DevVcvsInstanceData>::setup() {
    registerMember(ctl);
    registerMember(v);
    registerMember(i);
    return 0;
}
instantiateIntrospection(DevVcvsInstanceData);

DevVcvsInstanceData::DevVcvsInstanceData() {
}


template<> int Introspection<DevCccsInstanceData>::setup() {
    registerMember(ctl);
    registerMember(v);
    registerMember(i);
    return 0;
}
instantiateIntrospection(DevCccsInstanceData);

DevCccsInstanceData::DevCccsInstanceData() {
}


template<> int Introspection<DevCcvsInstanceData>::setup() {
    registerMember(ctl);
    registerMember(v);
    registerMember(i);
    return 0;
}
instantiateIntrospection(DevCcvsInstanceData);

DevCcvsInstanceData::DevCcvsInstanceData() {
}



template<typename InstanceType> 
Instance* findPeerInstance(Circuit& circuit, InstanceType& inst, Id name, Status& s) { 
    auto peerInstanceName = inst.translatePeer(name);
    auto* peerInstance = circuit.findInstance(peerInstanceName);
    if (!peerInstance) {
        s.set(Status::NotFound, "Peer instance '"+std::string(peerInstanceName)+"' not found.");
        s.extend(inst.location());
        return nullptr;
    }
    return peerInstance;
}

template<typename InstanceType> 
Node* findControl(Circuit& circuit, InstanceType& inst, Id instanceName, Id internalNodeName, Status& s) { 
    auto peerInstance = findPeerInstance(circuit, inst, instanceName, s);
    if (!peerInstance) {
        return nullptr;
    }
    Id nodeName = peerInstance->translate(internalNodeName);
    auto* node = circuit.findNode(nodeName);
    if (!node) {
        s.set(Status::NotFound, "Controlling unknown '"+std::string(nodeName)+"' not found.");
        s.extend(inst.location());
        return nullptr;
    }
    return node;
}

template<typename InstanceData> static bool getCtlsrcOpvar(InstanceData& data, ParameterIndex ndx, Value& v, Status& s) { 
    switch (ndx) {
    case 0:
        v = data.core().ctl;
        break;
    case 1:
        v = data.core().v;
        break;
    case 2:
        v = data.core().i;
        break;
    default:
        s.set(Status::Range, std::string("Opvar index id=")+std::to_string(ndx)+" out of range.");
        return false;
    }
    return true;
}

template<typename InstanceData> static std::tuple<bool, OutputSource> ctlsrcOpvarOutputSource(InstanceData& data, ParameterIndex ndx) { 
    switch (ndx) {
    case 0:
        return std::make_tuple(true, OutputSource(&data.core().ctl));
    case 1: 
        return std::make_tuple(true, OutputSource(&data.core().v));
    case 2: 
        return std::make_tuple(true, OutputSource(&data.core().i));
    default: 
        return std::make_tuple(false, OutputSource());
    }
}



template<> void BuiltinVccs::defineInternals() {
    nodeIds = { "p", "n", "cp", "cn" };
    terminalCount = 4;
}

template<> const Device::Flags BuiltinVccs::extraFlags = static_cast<Device::Flags>(0);

static ParameterIndex principalVccs = std::get<0>(Introspection<DevCctlSourceInstanceParams>::index("gain"));
template<> std::tuple<ParameterIndex, bool> BuiltinVccsInstance::principalParameterIndex() const {
    return std::make_tuple(principalVccs, true); // gain
}

template<> bool BuiltinVccsInstance::deleteHierarchy(Circuit& circuit, Status& s) { 
    return true; 
} 

template<> bool BuiltinVccsInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s) { 
    // If we require all terminals to be connected, do this
    if (!verifyTerminalsConnected(s)) { 
        return false;
    }
    return true; 
};  

template<> bool BuiltinVccsInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const { 
    return getCtlsrcOpvar(data, ndx, v, s);
}

template<> std::tuple<bool, OutputSource> BuiltinVccsInstance::opvarOutputSource(ParameterIndex ndx) const { 
    return ctlsrcOpvarOutputSource(data, ndx);
}

template<> bool BuiltinVccsInstance::populateStructuresCore(Circuit& circuit, Status& s) {
    // Create Jacobian entries
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[0], nodes_[2], s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[0], nodes_[3], s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[1], nodes_[2], s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[1], nodes_[3], s); !ok) {
        return false;
    }
    // No states to reserve
    return true;
}

template<> bool BuiltinVccsInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, 
    KluMatrixAccess* matReact, Component compReact, 
    Status& s
) {
    auto& d = data.core();

    // Unknown indices
    d.uP = nodes_[0]->unknownIndex();
    d.uN = nodes_[1]->unknownIndex();
    d.uCp = nodes_[2]->unknownIndex();
    d.uCn = nodes_[3]->unknownIndex();

    // Resistive Jacobian entry pointer
    if (matResist) {
        jacEntryPtr(d.jacPCp, d.uP, d.uCp, matResist, compResist);
        jacEntryPtr(d.jacPCn, d.uP, d.uCn, matResist, compResist);
        jacEntryPtr(d.jacNCp, d.uN, d.uCp, matResist, compResist);
        jacEntryPtr(d.jacNCn, d.uN, d.uCn, matResist, compResist);
    }
    
    // No reactive Jacobian entries
    
    return true;
}

template<> bool BuiltinVccsInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els) {
    auto& p = params.core();
    auto& d = data.core();
    auto& internals = circuit.simulatorInternals();
    
    // Evaluate
    if (!els.skipEvaluation && els.evaluateResistiveResidual) {
        d.ctl = els.oldSolution[d.uCp] - els.oldSolution[d.uCn]; // controlling voltage
        d.i = p.gain*d.ctl;  // current of one parallel instance
        d.v = els.oldSolution[d.uP] - els.oldSolution[d.uN]; // voltage across instance
    }

    // Load resistive Jacobian, transient load is identical because there is no reactive component
    if (els.loadResistiveJacobian || els.loadTransientJacobian) {
        // KCL
        *(d.jacPCp) += p.mfactor*p.gain;
        *(d.jacPCn) += -p.mfactor*p.gain;
        *(d.jacNCp) += -p.mfactor*p.gain;
        *(d.jacNCn) += p.mfactor*p.gain;
    }

    // Load resistive residual
    if (els.resistiveResidual) {
        els.resistiveResidual[d.uP] += p.mfactor*d.i;
        els.resistiveResidual[d.uN] += -p.mfactor*d.i;
    }

    // No limiting, so nothing to load for limited residual

    // Maximal residual contribution
    if (els.maxResistiveResidualContribution) {
        auto flowContrib = std::abs(p.mfactor*d.i);
        if (els.maxResistiveResidualContribution[d.uP]<flowContrib) {
            els.maxResistiveResidualContribution[d.uP] = flowContrib;
        }
        if (els.maxResistiveResidualContribution[d.uN]<flowContrib) {
            els.maxResistiveResidualContribution[d.uN] = flowContrib;
        }
    }

    // No reactive component, reactive residual derivative wrt. time is zero

    return true;
}


template<> void BuiltinVcvs::defineInternals() {
    nodeIds = { "p", "n", "cp", "cn", "flow(br)" };
    terminalCount = 4;
}

template<> const Device::Flags BuiltinVcvs::extraFlags = static_cast<Device::Flags>(0);

static ParameterIndex principalVcvs = std::get<0>(Introspection<DevCctlSourceInstanceParams>::index("gain"));
template<> std::tuple<ParameterIndex, bool> BuiltinVcvsInstance::principalParameterIndex() const {
    return std::make_tuple(principalVcvs, true); // gain
}

template<> bool BuiltinVcvsInstance::deleteHierarchy(Circuit& circuit, Status& s) { 
    if (!circuit.releaseNode(nodes_[4], s)) {
        return false;
    }
    nodes_[4] = nullptr;
    return true; 
} 

template<> bool BuiltinVcvsInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s) { 
    // If we require all terminals to be connected, do this
    if (!verifyTerminalsConnected(s)) { 
        return false;
    }
    
    // Create internal static flow node
    auto node = getInternalNode(circuit, nodeName(4), Node::Flags::FlowNode, s);
    if (!node) {
        return false;
    }

    // Bind static flow node
    nodes_[4] = node;

    return true; 
};  

template<> bool BuiltinVcvsInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const { 
    return getCtlsrcOpvar(data, ndx, v, s);
}

template<> std::tuple<bool, OutputSource> BuiltinVcvsInstance::opvarOutputSource(ParameterIndex ndx) const { 
    return ctlsrcOpvarOutputSource(data, ndx);
}

template<> bool BuiltinVcvsInstance::populateStructuresCore(Circuit& circuit, Status& s) {
    // Create Jacobian entries
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[0], nodes_[4], s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[1], nodes_[4], s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[4], nodes_[0], s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[4], nodes_[1], s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[4], nodes_[2], s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[4], nodes_[3], s); !ok) {
        return false;
    }
    // No states to reserve
    return true;
}

template<> bool BuiltinVcvsInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, 
    KluMatrixAccess* matReact, Component compReact, 
    Status& s
) {
    auto& d = data.core();

    // Unknown indices
    d.uP = nodes_[0]->unknownIndex();
    d.uN = nodes_[1]->unknownIndex();
    d.uCp = nodes_[2]->unknownIndex();
    d.uCn = nodes_[3]->unknownIndex();
    d.uFlow = nodes_[4]->unknownIndex();

    // Resistive Jacobian entry pointers
    if (matResist) {
        jacEntryPtr(d.jacPFlow,  d.uP,    d.uFlow, matResist, compResist);
        jacEntryPtr(d.jacNFlow,  d.uN,    d.uFlow, matResist, compResist);
        jacEntryPtr(d.jacFlowP,  d.uFlow, d.uP,    matResist, compResist);
        jacEntryPtr(d.jacFlowN,  d.uFlow, d.uN,    matResist, compResist);
        jacEntryPtr(d.jacFlowCp, d.uFlow, d.uCp,   matResist, compResist);
        jacEntryPtr(d.jacFlowCn, d.uFlow, d.uCn,   matResist, compResist);
    }

    // No reactive Jacobian entries
    
    return true;
}

template<> bool BuiltinVcvsInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els) {
    auto& p = params.core();
    auto& d = data.core();
    auto& internals = circuit.simulatorInternals();
    
    // Evaluate
    if (!els.skipEvaluation && els.evaluateResistiveResidual) {
        d.ctl = els.oldSolution[d.uCp] - els.oldSolution[d.uCn]; // controlling voltage
        d.v = els.oldSolution[d.uP] - els.oldSolution[d.uN]; // coltage across instance
        d.i = els.oldSolution[d.uFlow]; // current of one parallel instance
    }

    // Load resistive Jacobian, transient load is identical because there is no reactive component
    if (els.loadResistiveJacobian || els.loadTransientJacobian) {
        // KCL
        *(d.jacPFlow) += p.mfactor;
        *(d.jacNFlow) += -p.mfactor;
        // Control
        *(d.jacFlowP) += -1;
        *(d.jacFlowN) += 1;
        *(d.jacFlowCp) += p.gain;
        *(d.jacFlowCn) += -p.gain;
    }

    // Load resistive residual
    if (els.resistiveResidual) {
        els.resistiveResidual[d.uP] += p.mfactor*d.i;
        els.resistiveResidual[d.uN] += -p.mfactor*d.i;
        els.resistiveResidual[d.uFlow] += -d.v + d.ctl*p.gain;
    }

    // No limiting, so nothing to load for limited residual

    // Maximal residual contribution
    if (els.maxResistiveResidualContribution) {
        auto valueContrib = std::abs(-d.v + d.ctl*p.gain);
        auto flowContrib = std::abs(p.mfactor*d.i);
        if (els.maxResistiveResidualContribution[d.uP]<flowContrib) {
            els.maxResistiveResidualContribution[d.uP] = flowContrib;
        }
        if (els.maxResistiveResidualContribution[d.uN]<flowContrib) {
            els.maxResistiveResidualContribution[d.uN] = flowContrib;
        }
        if (els.maxResistiveResidualContribution[d.uFlow]<valueContrib) {
            els.maxResistiveResidualContribution[d.uFlow] = valueContrib;
        }
    }

    // No reactive component, reactive residual derivative wrt. time is zero

    return true;
}


template<> void BuiltinCccs::defineInternals() {
    nodeIds = { "p", "n" };
    terminalCount = 2;
}

template<> const Device::Flags BuiltinCccs::extraFlags = static_cast<Device::Flags>(0);

static ParameterIndex principalCccs = std::get<0>(Introspection<DevCctlSourceInstanceParams>::index("gain"));
template<> std::tuple<ParameterIndex, bool> BuiltinCccsInstance::principalParameterIndex() const {
    return std::make_tuple(principalCccs, true); // gain
}

template<> bool BuiltinCccsInstance::deleteHierarchy(Circuit& circuit, Status& s) { 
    return true; 
} 

template<> bool BuiltinCccsInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s) { 
    // If we require all terminals to be connected, do this
    if (!verifyTerminalsConnected(s)) { 
        return false;
    }
    return true; 
};  

template<> bool BuiltinCccsInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const { 
    return getCtlsrcOpvar(data, ndx, v, s);
}

template<> std::tuple<bool, OutputSource> BuiltinCccsInstance::opvarOutputSource(ParameterIndex ndx) const { 
    return ctlsrcOpvarOutputSource(data, ndx);
}

template<> bool BuiltinCccsInstance::populateStructuresCore(Circuit& circuit, Status& s) {
    // Find controlling node
    auto ctlNode = findControl(circuit, *this, params.core().ctlinst, params.core().ctlnode, s);
    if (!ctlNode) {
        return false;
    }
    data.core().uCtl = ctlNode->unknownIndex();

    // Create Jacobian entries
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[0], ctlNode, s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[1], ctlNode, s); !ok) {
        return false;
    }
    // No states to reserve
    return true;
}

template<> bool BuiltinCccsInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, 
    KluMatrixAccess* matReact, Component compReact, 
    Status& s
) {
    auto& d = data.core();

    // Unknown indices
    d.uP = nodes_[0]->unknownIndex();
    d.uN = nodes_[1]->unknownIndex();
    
    // Resistive Jacobian entry pointers
    if (matResist) {
        jacEntryPtr(d.jacPCtl,  d.uP,    d.uCtl, matResist, compResist);
        jacEntryPtr(d.jacNCtl,  d.uN,    d.uCtl, matResist, compResist);
    }
    
    // No reactive Jacobian entries
    
    return true;
}

template<> bool BuiltinCccsInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els) {
    auto& p = params.core();
    auto& d = data.core();
    auto& internals = circuit.simulatorInternals();
    
    // Evaluate
    if (!els.skipEvaluation && els.evaluateResistiveResidual) {
        d.ctl = els.oldSolution[d.uCtl]; // controlling current (unknown)
        d.i = p.gain*d.ctl;  // current of one parallel instance
        d.v = els.oldSolution[d.uP] - els.oldSolution[d.uN]; // voltage across instance
    }

    // Load resistive Jacobian, transient load is identical because there is no reactive component
    if (els.loadResistiveJacobian || els.loadTransientJacobian) {
        // KCL
        *(d.jacPCtl) += p.mfactor*p.gain;
        *(d.jacNCtl) += -p.mfactor*p.gain;
    }

    // Load resistive residual
    if (els.resistiveResidual) {
        els.resistiveResidual[d.uP] += p.mfactor*d.i;
        els.resistiveResidual[d.uN] += -p.mfactor*d.i;
    }

    // No limiting, so nothing to load for limited residual

    // Maximal residual contribution
    if (els.maxResistiveResidualContribution) {
        auto flowContrib = std::abs(p.mfactor*d.i);
        if (els.maxResistiveResidualContribution[d.uP]<flowContrib) {
            els.maxResistiveResidualContribution[d.uP] = flowContrib;
        }
        if (els.maxResistiveResidualContribution[d.uN]<flowContrib) {
            els.maxResistiveResidualContribution[d.uN] = flowContrib;
        }
    }

    // No reactive component, reactive residual derivative wrt. time is zero

    return true;
}


template<> void BuiltinCcvs::defineInternals() {
    nodeIds = { "p", "n", "flow(br)" };
    terminalCount = 2;
}

template<> const Device::Flags BuiltinCcvs::extraFlags = static_cast<Device::Flags>(0);

static ParameterIndex principalCcvs = std::get<0>(Introspection<DevCctlSourceInstanceParams>::index("gain"));
template<> std::tuple<ParameterIndex, bool> BuiltinCcvsInstance::principalParameterIndex() const {
    return std::make_tuple(principalCcvs, true); // gain
}

template<> bool BuiltinCcvsInstance::deleteHierarchy(Circuit& circuit, Status& s) { 
    if (!circuit.releaseNode(nodes_[2], s)) {
        return false;
    }
    nodes_[2] = nullptr;
    return true; 
} 

template<> bool BuiltinCcvsInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s) { 
    // If we require all terminals to be connected, do this
    if (!verifyTerminalsConnected(s)) { 
        return false;
    }
    
    // Create internal static flow node
    auto node = getInternalNode(circuit, nodeName(2), Node::Flags::FlowNode, s);
    if (!node) {
        return false;
    }

    // Bind static flow node
    nodes_[2] = node;

    return true; 
};  

template<> bool BuiltinCcvsInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const { 
    return getCtlsrcOpvar(data, ndx, v, s); 
}

template<> std::tuple<bool, OutputSource> BuiltinCcvsInstance::opvarOutputSource(ParameterIndex ndx) const { 
    return ctlsrcOpvarOutputSource(data, ndx);
}

template<> bool BuiltinCcvsInstance::populateStructuresCore(Circuit& circuit, Status& s) {
    // Find controlling node
    auto ctlNode = findControl(circuit, *this, params.core().ctlinst, params.core().ctlnode, s);
    if (!ctlNode) {
        return false;
    }
    data.core().uCtl = ctlNode->unknownIndex();

    // Create Jacobian entries
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[0], nodes_[2], s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[1], nodes_[2], s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[2], nodes_[0], s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[2], nodes_[1], s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(nodes_[2], ctlNode, s); !ok) {
        return false;
    }
    // No states to reserve
    return true;
}

template<> bool BuiltinCcvsInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, 
    KluMatrixAccess* matReact, Component compReact, 
    Status& s
) {
    auto& d = data.core();

    // Unknown indices
    d.uP = nodes_[0]->unknownIndex();
    d.uN = nodes_[1]->unknownIndex();
    d.uFlow = nodes_[2]->unknownIndex();

    // Resistive Jacobian entry pointers
    if (matResist) {
        jacEntryPtr(d.jacPFlow,   d.uP,    d.uFlow, matResist, compResist);
        jacEntryPtr(d.jacNFlow,   d.uN,    d.uFlow, matResist, compResist);
        jacEntryPtr(d.jacFlowP,   d.uFlow, d.uP,    matResist, compResist);
        jacEntryPtr(d.jacFlowN,   d.uFlow, d.uN,    matResist, compResist);
        jacEntryPtr(d.jacFlowCtl, d.uFlow, d.uCtl,  matResist, compResist);
    }
    
    // No reactive Jacobian entries
    
    return true;
}

template<> bool BuiltinCcvsInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els) {
    auto& p = params.core();
    auto& d = data.core();
    auto& internals = circuit.simulatorInternals();
    
    // Evaluate
    if (!els.skipEvaluation && els.evaluateResistiveResidual) {
        d.ctl = els.oldSolution[d.uCtl]; // controlling current
        d.v = els.oldSolution[d.uP] - els.oldSolution[d.uN]; // voltage across instance 
        d.i = els.oldSolution[d.uFlow]; // current of one parallel instance
    }

    // Load resistive Jacobian, transient load is identical because there is no reactive component
    if (els.loadResistiveJacobian || els.loadTransientJacobian) {
        // KCL
        *(d.jacPFlow) += p.mfactor;
        *(d.jacNFlow) += -p.mfactor;
        // Control
        *(d.jacFlowP) += -1;
        *(d.jacFlowN) += 1;
        *(d.jacFlowCtl) += p.gain;
    }

    // Load resistive residual
    if (els.resistiveResidual) {
        els.resistiveResidual[d.uP] += p.mfactor*d.i;
        els.resistiveResidual[d.uN] += -p.mfactor*d.i;
        els.resistiveResidual[d.uFlow] += -d.v + d.ctl*p.gain;
    }

    // No limiting, so nothing to load for limited residual

    // Maximal residual contribution
    if (els.maxResistiveResidualContribution) {
        auto valueContrib = std::abs(-d.v + d.ctl*p.gain);
        auto flowContrib = std::abs(p.mfactor*d.i);
        if (els.maxResistiveResidualContribution[d.uP]<flowContrib) {
            els.maxResistiveResidualContribution[d.uP] = flowContrib;
        }
        if (els.maxResistiveResidualContribution[d.uN]<flowContrib) {
            els.maxResistiveResidualContribution[d.uN] = flowContrib;
        }
        if (els.maxResistiveResidualContribution[d.uFlow]<valueContrib) {
            els.maxResistiveResidualContribution[d.uFlow] = valueContrib;
        }
    }

    // No reactive component, reactive residual derivative wrt. time is zero

    return true;
}


template<> int Introspection<DevMutualInstanceParams>::setup() {
    registerMember(k);
    registerMember(ind1);
    registerMember(ind2);
    registerMember(ctlnode1);
    registerMember(ctlnode2);
    return 0;
}
instantiateIntrospection(DevMutualInstanceParams);

DevMutualInstanceParams::DevMutualInstanceParams() {
}

template<> int Introspection<DevMutualInstanceData>::setup() {
    registerMember(mutual);
    return 0;
}
instantiateIntrospection(DevMutualInstanceData);

DevMutualInstanceData::DevMutualInstanceData() {
}


template<> void BuiltinMutual::defineInternals() {
    nodeIds = {};
    terminalCount = 0;
}

template<> const Device::Flags BuiltinMutual::extraFlags = static_cast<Device::Flags>(0);

static ParameterIndex principalMutual = std::get<0>(Introspection<DevMutualInstanceParams>::index("k"));
template<> std::tuple<ParameterIndex, bool> BuiltinMutualInstance::principalParameterIndex() const {
    return std::make_tuple(principalMutual, true); // gain
}

template<> std::tuple<bool, bool, bool> BuiltinMutualInstance::setupWorker(Circuit& circuit, Status& s) {
    auto& p = params.core();
    auto& d = data.core();

    // Check k
    if (p.k>1 || p.k<0) {
        s.set(Status::BadArguments, "Parameter 'k' must be between 0 and 1.");
        return std::make_tuple(false, false, false);
    }

    // Find both inductors
    auto ind1 = findPeerInstance(circuit, *this, p.ind1, s);
    if (!ind1) {
        return std::make_tuple(false, false, false); 
    }
    auto ind2 = findPeerInstance(circuit, *this, p.ind2, s);
    if (!ind2) {
        return std::make_tuple(false, false, false); 
    }
    d.ind1 = ind1;
    d.ind2 = ind2;

    // Get internal control nodes
    auto ctlNode1 = findControl(circuit, *this, p.ind1, p.ctlnode1, s);
    if (!ctlNode1) {
        return std::make_tuple(false, false, false); 
    }
    auto ctlNode2 = findControl(circuit, *this, p.ind2, p.ctlnode2, s);
    if (!ctlNode2) {
        return std::make_tuple(false, false, false); 
    }
    if (!ctlNode1->checkFlags(Node::Flags::FlowNode)) {
        s.set(Status::NotFound, "Instance '"+std::string(ind1->name())+"' has no flow node to bind to.");
        s.extend(location());
        return std::make_tuple(false, false, false); 
    }
    if (!ctlNode2->checkFlags(Node::Flags::FlowNode)) {
        s.set(Status::NotFound, "Instance '"+std::string(ind1->name())+"' has no flow node to bind to.");
        s.extend(location());
        return std::make_tuple(false, false, false); 
    }

    // Check if controlling unknowns changed, store indices
    bool sparsityChanged;
    if (!d.ctlNode1 || !d.ctlNode2) {
        // No control nodes yet
        sparsityChanged = true;
    } else {
        // Check if unknown indices changed
        sparsityChanged = (
            ctlNode1->unknownIndex() != d.ctlNode1->unknownIndex() || 
            ctlNode2->unknownIndex() != d.ctlNode2->unknownIndex()
        );
    }

    // Store control nodes
    d.ctlNode1 = ctlNode1;
    d.ctlNode2 = ctlNode2;

    // Store controlling unknowns
    d.uFlow1 = ctlNode1->unknownIndex();
    d.uFlow2 = ctlNode2->unknownIndex();
    
    // No change in unknowns, possible change in sparsity
    return std::make_tuple(true, false, sparsityChanged);
}

template<> bool BuiltinMutualInstance::preAnalysisWorker(Circuit& circuit, Status& s) {
    // Compute mutual inductance
    auto& p = params.core();
    auto& d = data.core();

    Value v;
    
    // Get L1
    if (!d.ind1->getParameter("l", v, s)) {
        s.extend(location());
        return false;
    }
    if (v.type()!=Value::Type::Real) {
        s.set(Status::BadConversion, "Parameter 'l' of '"+std::string(p.ind1)+"' is not a real number.");
        s.extend(location());
        return false;
    }
    auto l1 = v.val<Real>();
    // Get m1
    if (!d.ind1->getParameter("$mfactor", v, s)) {
        s.extend(location());
        return false;
    }
    if (v.type()!=Value::Type::Real) {
        s.set(Status::BadConversion, "Parameter '$mfactor' of '"+std::string(p.ind1)+"' is not a real number.");
        s.extend(location());
        return false;
    }
    auto m1 = v.val<Real>();
    // Get L2
    if (!d.ind2->getParameter("l", v, s)) {
        s.extend(location());
        return false;
    }
    if (v.type()!=Value::Type::Real) {
        s.set(Status::BadConversion, "Parameter 'l' of '"+std::string(p.ind2)+"' is not a real number.");
        s.extend(location());
        return false;
    }
    auto l2 = v.val<Real>();
    // Get L2
    if (!d.ind2->getParameter("$mfactor", v, s)) {
        s.extend(location());
        return false;
    }
    // $mfactor is a real parameter in OSDI, therefore we also make it real
    if (v.type()!=Value::Type::Real) {
        s.set(Status::BadConversion, "Parameter '$mfactor' of '"+std::string(p.ind2)+"' is not a real number.");
        s.extend(location());
        return false;
    }
    auto m2 = v.val<Real>();
    if (m1<=0 || m2<=0) {
        s.set(Status::BadArguments, "Multipliers ('$mfactor') of inductors must be greater than 0.");
        s.extend(location());
        return false;
    }
    // Compute mutual inductance
    d.mutual = p.k*std::sqrt(l1*l2/(m1*m2));

    // Note: coupling inductors with $mfactor!=1 does not make sense. 
    // $mfactor != 1 means that we have $mfactor mutually uncoupled inductors in parallel
    // which effectively reduces inductance by factor $mfactor. 
    // If, however, we couple two inductors (e.g. L1 and L2) via mutual, we couple each parallel 
    // instace of L1 with each parallel instance of L2 and vice versa. This implies that 
    // suddenly individual parallel instances of L1 are now also mutually coupled which 
    // in turn goes against the initial assumption that they are uncoupled. 
    // We try to make some sense by using L1/$mfactor1 and L2/$mfator2, but that is 
    // not completely sane. 

    return true;
}

template<> bool BuiltinMutualInstance::populateStructuresCore(Circuit& circuit, Status& s) {
    auto& d = data.core();

    // Find controlling node
    // Create Jacobian entries
    if (auto [ptr, ok] = circuit.createJacobianEntry(d.ctlNode1, d.ctlNode2, s); !ok) {
        return false;
    }
    if (auto [ptr, ok] = circuit.createJacobianEntry(d.ctlNode2, d.ctlNode1, s); !ok) {
        return false;
    }

    // States for M di1/dt and M di2/dt
    data.core().offsStates = circuit.allocateStates(2*2);

    return true;
}

template<> bool BuiltinMutualInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, 
    KluMatrixAccess* matReact, Component compReact, 
    Status& s
) {
    auto& d = data.core();
    
    // No resistive Jacobian entries
    
    // Reactive Jacobian entry pointers
    if (matReact) {
        jacEntryPtr(d.jacReact12, d.uFlow1, d.uFlow2, matReact, compReact);
        jacEntryPtr(d.jacReact21, d.uFlow2, d.uFlow1, matReact, compReact);
    }
        
    return true;
}

template<> bool BuiltinMutualInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els) {
    auto& p = params.core();
    auto& d = data.core();
    auto& internals = circuit.simulatorInternals();

    // Evaluate
    if (!els.skipEvaluation && els.evaluateReactiveResidual) {
        d.res1 = d.mutual * els.oldSolution[d.uFlow2];
        d.res2 = d.mutual * els.oldSolution[d.uFlow1];
    }

    // Load reactive Jacobian, transient load is identical because there is no resistive component
    if (els.loadReactiveJacobian) {
        auto factor = els.reactiveJacobianFactor;
        // Extra equations
        *(d.jacReact12) += d.mutual*factor;
        *(d.jacReact21) += d.mutual*factor;
    }

    if (els.loadTransientJacobian) {
        auto factor = els.integCoeffs->leadingCoeff();
        // Extra equations
        *(d.jacReact12) += d.mutual*factor;
        *(d.jacReact21) += d.mutual*factor;
    }

    // Load reactive residual
    if (els.reactiveResidual) {
        els.resistiveResidual[d.uFlow1] += d.res1;
        els.resistiveResidual[d.uFlow2] += d.res2;
    }

    // No limiting, so nothing to load for limited residual

    // Maximal residual contribution
    if (els.maxReactiveResidualContribution) {
        if (els.maxReactiveResidualContribution[d.uFlow1]<d.res1) {
            els.maxReactiveResidualContribution[d.uFlow1] = d.res1;
        }
        if (els.maxReactiveResidualContribution[d.uFlow2]<d.res2) {
            els.maxReactiveResidualContribution[d.uFlow2] = d.res2;
        }
    }

    // Reactive residual derivative wrt. time
    if (
        els.storeTerminalReactiveResidualState ||
        els.storeTerminalReactiveResidualDerivativeState ||
        els.reactiveResidualDerivative ||
        els.maxReactiveResidualDerivativeContribution
    ) { 
        // Store residual state
        if (els.storeTerminalReactiveResidualState) {
            els.futureStates[d.offsStates] = d.res1; 
            els.futureStates[d.offsStates+2] = d.res2; 
        }
        
        // Compute residual derivative
        if (
            els.storeTerminalReactiveResidualDerivativeState ||
            els.reactiveResidualDerivative || 
            els.maxReactiveResidualDerivativeContribution
        ) {
            // Differentiate (compute flow)
            double res1dot = els.integCoeffs->differentiate(d.res1, d.offsStates);
            double res2dot = els.integCoeffs->differentiate(d.res2, d.offsStates+2);

            // Store flow in states vector
            if (els.storeTerminalReactiveResidualDerivativeState) {
                els.futureStates[d.offsStates+1] = res1dot;
                els.futureStates[d.offsStates+3] = res2dot;
            }

            // Add flow to vector
            if (els.reactiveResidualDerivative) {
                els.reactiveResidualDerivative[d.uFlow1] += res1dot;
                els.reactiveResidualDerivative[d.uFlow2] += res2dot;
            }

            // Update max residual contribution
            if (els.maxReactiveResidualDerivativeContribution) {
                auto contrib1 = std::abs(res1dot);
                if (contrib1 > els.maxReactiveResidualDerivativeContribution[d.uFlow1]) {
                    els.maxReactiveResidualDerivativeContribution[d.uFlow1] = contrib1;
                }
                auto contrib2 = std::abs(res2dot);
                if (contrib2 > els.maxReactiveResidualDerivativeContribution[d.uFlow2]) {
                    els.maxReactiveResidualDerivativeContribution[d.uFlow2] = contrib2;
                }
            }
        }
    }

    return true;
}


}
