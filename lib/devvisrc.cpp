#include <numbers>
#include "devvisrc.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

static const double PI = std::numbers::pi;

template<> int Introspection<DevSourceModelParams>::setup() {
    return 0;
}
instantiateIntrospection(DevSourceModelParams);


DevSourceModelParams::DevSourceModelParams() {
}

template<> int Introspection<DevSourceInstanceParams>::setup() {
    registerNamedMember(mfactor, "$mfactor");
    registerMember(type);
    registerMember(delay);
    registerMember(dc);
    registerMember(val0);
    registerMember(val1);
    registerMember(period);
    registerMember(rise);
    registerMember(fall);
    registerMember(width);
    registerMember(sinedc);
    registerMember(ampl);
    registerMember(freq);
    registerMember(sinephase);
    registerMember(mag);
    registerMember(phase);
    return 0;
}
instantiateIntrospection(DevSourceInstanceParams);

DevSourceInstanceParams::DevSourceInstanceParams() {
    mfactor = 1;

    type = "dc";
    delay = 0.0;
    
    dc = 0.0;

    val0 = 0.0;
    val1 = 1.0;
    period = 0.0;
    rise = 1e-9;
    fall = 0.0;
    width = 0.0;

    sinedc = 0.0;
    ampl = 1.0;
    freq = 1e3;
    sinephase = 0.0;

    mag = 0.0;
    phase = 0.0;
}

static ParameterIndex principal = std::get<0>(Introspection<DevSourceInstanceParams>::index("dc"));

template<> int Introspection<DevVSourceInstanceData>::setup() {
    registerMember(v);
    registerMember(i);
    return 0;
}
instantiateIntrospection(DevVSourceInstanceData);

DevVSourceInstanceData::DevVSourceInstanceData() {
}

template<> int Introspection<DevISourceInstanceData>::setup() {
    registerMember(v);
    registerMember(i);
    return 0;
}
instantiateIntrospection(DevISourceInstanceData);

DevISourceInstanceData::DevISourceInstanceData() {
}

static Id typeSine = Id::createStatic("sine");
static Id typePulse = Id::createStatic("pulse");
static Id typeDc = Id::createStatic("dc");

template<typename InstanceParams, typename InstanceData> 
std::tuple<bool, bool, bool> sourceSetup(InstanceParams& params, InstanceData& data, Loc loc, Circuit& circuit, Status& s) { 
    auto& p = params.core();
    auto& d = data.core();

    // Store quick access code, check parameters
    if (p.type == typeSine) {
        d.typeCode = IndependentSourceType::Sine;
        if (p.freq<=0) {
            s.set(Status::BadArguments, "Frequency of sinusoidal transient must be greater than 0.");
            s.extend(loc);
            return std::make_tuple(false, false, false);
        }
    } else if (p.type == typePulse) {
        d.typeCode = IndependentSourceType::Pulse;
        if (p.rise<=0) {
            s.set(Status::BadArguments, "Rise time of pulse transient must be grater than 0.");
            s.extend(loc);
            return std::make_tuple(false, false, false);
        }
        
        // fall<=0 generates a step

        // width<=0 generates
        // - a step if fall<=0
        // - a triangle if fall>0

        // period<=0 generates a single pulse

        // If greater than 0, period must be greater than rise+fall+width
        if (p.period>0 && p.period<=p.rise+p.fall+p.width) {
            s.set(Status::BadArguments, "Period of pulse transient must be greater than rise+fall+width.");
            s.extend(loc);
            return std::make_tuple(false, false, false);
        }
    } else if (p.type == typeDc) {
        d.typeCode = IndependentSourceType::Dc;
    } else {
        s.set(Status::BadArguments, "Unknown transient waveform type.");
        s.extend(loc);
        return std::make_tuple(false, false, false);
    }

    return std::make_tuple(true, false, false); 
}

template<typename InstanceParams, typename InstanceData> 
std::tuple<double, double> sourceCompute(const InstanceParams& params, InstanceData& data, double time) {
    double val;
    double nextBreak = std::numeric_limits<double>::infinity();

    switch (data.typeCode) {
    case IndependentSourceType::Dc:
        val = params.dc;
        break;
    case IndependentSourceType::Sine:
        val = params.sinedc+params.ampl*std::sin(2*PI*params.freq*(time-params.delay)+params.phase*PI/180);
        break;
    case IndependentSourceType::Pulse:
        if (time<params.delay) {
            // Before waveform starts
            val = params.val0;
            nextBreak = params.delay;
        } else {
            // Waveform started, see where we are
            // Time since start of this repetition
            double basetime = params.delay;
            // Start of current period
            if (params.period<=0) {
                // Not periodic
                basetime = params.delay;
            } else {
                // Periodic
                basetime = params.delay + std::floor((time-params.delay)/params.period)*params.period;
            }
            // Significant time points
            auto t0 = basetime; // start of rise
            auto t1 = t0 + params.rise; // end of rise, start of top
            auto t2 = t1 + params.width; // start of fall
            auto t3 = t2 + params.fall; // end of fall
            auto t4 = t0 + params.period; // end of period
            // Relative time since start of period
            double reltime = time - basetime;
            if (time>=t0 && time<t1-timeRelativeTolerance*t1) {
                // Rising
                val = (time-t0)/params.rise*(params.val1-params.val0)+params.val0;
                nextBreak = t1;
            } else if (time>=t1 && time<t2-timeRelativeTolerance*t2) {
                // On top
                val = params.val1;
                // Set next break only if width>0
                if (params.width>0) {
                    nextBreak = t2;
                }
            } else if (params.fall<=0) {
                // Fall not set, stay on top, no breakpoint
                val = params.val1;
            } else if (time>=t2 && time<t3-timeRelativeTolerance*t3) {
                // Falling
                val = (time-t2)/params.fall*(params.val0-params.val1)+params.val1;
                nextBreak = t3;
            } else {
                // After fall, back on base level
                val = params.val0;
                // Beakpoint only if period is set
                if (params.period>0) {
                    nextBreak = t4;
                }
            } 
        }
        break;
    }
    
    // Simulator::dbg() << "val=" << val << " next break=" << nextBreak << "\n";
    // Simulator::dbg() << " next break=" << nextBreak << "\n";
    
    return std::make_tuple(val, nextBreak);
}

template<> void BuiltinVSource::defineInternals() {
    nodeIds = { "p", "n", "flow(br)" };
    terminalCount = 2;
}

template<> void BuiltinISource::defineInternals() {
    nodeIds = { "p", "n" };
    terminalCount = 2;
}

template<> bool BuiltinVSource::isSource() const {
    return true;
}

template<> bool BuiltinISource::isSource() const {
    return true;
}

template<> bool BuiltinVSource::isVoltageSource() const {
    return true;
}


template<> const Device::Flags BuiltinVSource::extraFlags = 
    Device::Flags::GeneratesAC | Device::Flags::GeneratesDCIncremental;

template<> const Device::Flags BuiltinISource::extraFlags = 
    Device::Flags::GeneratesAC | Device::Flags::GeneratesDCIncremental;

template<> std::tuple<ParameterIndex, bool> BuiltinVSourceInstance::principalParameterIndex() const {
    return std::make_tuple(principal, true); 
}

template<> std::tuple<ParameterIndex, bool> BuiltinISourceInstance::principalParameterIndex() const {
    return std::make_tuple(principal, true); 
}

template<> bool BuiltinVSourceInstance::deleteHierarchy(Circuit& circuit, Status& s) { 
    if (!circuit.releaseNode(nodes_[2], s)) {
        return false;
    }
    nodes_[2] = nullptr;
    return true; 
} 

template<> bool BuiltinVSourceInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s) { 
    // If we want to leave unconnected terminals hanging, do this 
    // // Create internal nodes for unconnected terminals
    // createNodesForUnconnectedTerminals();
    
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

template<> bool BuiltinISourceInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s) { 
    // If we want to leave unconnected terminals hanging, do this 
    // // Create internal nodes for unconnected terminals
    // createNodesForUnconnectedTerminals();

    // If we require all terminals to be connected, do this
    if (!verifyTerminalsConnected(s)) { 
        return false;
    }
    return true;
}

template<> std::tuple<EquationIndex,EquationIndex> BuiltinVSourceInstance::sourceExcitation(Circuit& circuit) const { 
    return std::make_tuple(nodes_[2]->unknownIndex(), 0); 
}

template<> std::tuple<UnknownIndex,UnknownIndex> BuiltinVSourceInstance::sourceResponse(Circuit& circuit) const { 
    return std::make_tuple(0, nodes_[2]->unknownIndex()); 
}

template<> double BuiltinVSourceInstance::scaledUnityExcitation() const { 
    // mag=1 produces a voltage of 1V, regardless of $mfactor. 
    return 1.0; 
}

template<> double BuiltinVSourceInstance::responseScalingFactor() const { 
    // Because the computed response (branch current) is 1/$mfactor times 
    // the total current of all parallel instances combined, the 
    // scaling factor must be $mfactor. 
    return params.core().mfactor; 
}

template<> std::tuple<EquationIndex,EquationIndex> BuiltinISourceInstance::sourceExcitation(Circuit& circuit) const { 
    return std::make_tuple(nodes_[0]->unknownIndex(), nodes_[1]->unknownIndex()); 
}

template<> std::tuple<UnknownIndex,UnknownIndex> BuiltinISourceInstance::sourceResponse(Circuit& circuit) const { 
    return std::make_tuple(nodes_[1]->unknownIndex(), nodes_[0]->unknownIndex()); 
}

template<> double BuiltinISourceInstance::scaledUnityExcitation() const { 
    // mag=1 produces a current of $mfactor A. 
    return params.core().mfactor; 
}

template<> double BuiltinISourceInstance::responseScalingFactor() const { 
    // Because the computed response (branch voltage) is the same for 
    // all parallel instances the scaling factor must be 1. 
    return 1.0; 
}

template<> bool BuiltinVSourceInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const { 
    switch (ndx) {
    case 0:
        v = data.core().v;
        break;
    case 1:
        v = data.core().i;
        break;
    default:
        s.set(Status::Range, std::string("Opvar index id=")+std::to_string(ndx)+" out of range.");
        return false;
    }
    return true;
}

template<> bool BuiltinISourceInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const { 
    switch (ndx) {
    case 0:
        v = data.core().v;
        break;
    case 1:
        v = data.core().i;
        break;
    default:
        s.set(Status::Range, std::string("Opvar index id=")+std::to_string(ndx)+" out of range.");
        return false;
    }
    return true;
}

template<> std::tuple<bool, OutputSource> BuiltinVSourceInstance::opvarOutputSource(ParameterIndex ndx, Status& s) const { 
    switch (ndx) {
    case 0:
        return std::make_tuple(true, OutputSource(&data.core().v));
    case 1: 
        return std::make_tuple(true, OutputSource(&data.core().i));
    default: 
        return std::make_tuple(false, OutputSource());
    }
}

template<> std::tuple<bool, OutputSource> BuiltinISourceInstance::opvarOutputSource(ParameterIndex ndx, Status& s) const { 
    switch (ndx) {
    case 0:
        return std::make_tuple(true, OutputSource(&data.core().v));
    case 1: 
        return std::make_tuple(true, OutputSource(&data.core().i));
    default: 
        return std::make_tuple(false, OutputSource());
    }
}

template<> std::tuple<bool, bool, bool> BuiltinVSourceInstance::setupWorker(Circuit& circuit, Status& s) {
    clearFlags(Flags::NeedsSetup); 
    return sourceSetup(params, data, location(), circuit, s);
}; 

template<> std::tuple<bool, bool, bool> BuiltinISourceInstance::setupWorker(Circuit& circuit, Status& s) { 
    clearFlags(Flags::NeedsSetup);
    return sourceSetup(params, data, location(), circuit, s);
}; 

template<> bool BuiltinVSourceInstance::populateStructuresCore(Circuit& circuit, Status& s) {
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
    // No states to reserve
    return true;
}

template<> bool BuiltinISourceInstance::populateStructuresCore(Circuit& circuit, Status& s) {
    // No Jacobian entries
    // No states to reserve
    return true;
}

template<> bool BuiltinVSourceInstance::bindCore(
    Circuit& circuit, 
    KluRealMatrix* matResistReal, KluComplexMatrix* matResistCx, Component compResist, 
    KluRealMatrix* matReactReal, KluComplexMatrix* matReactCx, Component compReact, 
    Status& s
) {
    auto& d = data.core();

    // Unknown indices
    d.uP = nodes_[0]->unknownIndex();
    d.uN = nodes_[1]->unknownIndex();
    d.uFlow = nodes_[2]->unknownIndex();

    // Resistive Jacobian entry pointers
    jacEntryPtr(d.jacPFlow, d.uP, d.uFlow, matResistReal, matResistCx, compResist);
    jacEntryPtr(d.jacNFlow, d.uN, d.uFlow, matResistReal, matResistCx, compResist);
    jacEntryPtr(d.jacFlowP, d.uFlow, d.uP, matResistReal, matResistCx, compResist);
    jacEntryPtr(d.jacFlowN, d.uFlow, d.uN, matResistReal, matResistCx, compResist);

    // No reactive Jacobian entries
    
    return true;
}

template<> bool BuiltinISourceInstance::bindCore(
    Circuit& circuit, 
    KluRealMatrix* matResistReal, KluComplexMatrix* matResistCx, Component compResist, 
    KluRealMatrix* matReactReal, KluComplexMatrix* matReactCx, Component compReact, 
    Status& s
) {
    auto& d = data.core();

    // Unknown indices
    d.uP = nodes_[0]->unknownIndex();
    d.uN = nodes_[1]->unknownIndex();
    
    // No Jacobian entries
    
    return true;
}

template<> bool BuiltinVSourceInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els, Status& s) {
    auto& p = params.core();
    auto& d = data.core();
    auto& internals = circuit.simulatorInternals();
    auto sourceFactor = internals.sourcescalefactor;
    
    // Evaluate
    if (!els.skipEvaluation) {
        if (els.evaluateResistiveResidual || els.computeNextBreakpoint) {
            auto [val, nextBreakpoint] = sourceCompute(p, d, internals.time);
            if (els.evaluateResistiveResidual) {
                d.v = val; // mfactor does not affect voltage source value
                d.i = els.oldSolution[d.uFlow]; // flow across one parallel instance
            }
            if (els.computeNextBreakpoint) {
                els.setBreakPoint(nextBreakpoint, internals); 
            }
        }
    }

    // Load resistive Jacobian, transient load is identical because there is no reactive component
    if (els.loadResistiveJacobian || els.loadTransientJacobian) {
        // KCL
        *(d.jacPFlow) += p.mfactor;
        *(d.jacNFlow) += -p.mfactor;
        // Extra equation
        *(d.jacFlowP) += -1.0;
        *(d.jacFlowN) += 1.0;
    }

    // Load resistive residual
    if (els.resistiveResidual) {
        auto valueRes = -els.oldSolution[d.uP] + els.oldSolution[d.uN] + sourceFactor*d.v;
        auto flowRes = p.mfactor*d.i;
        els.resistiveResidual[d.uP] += flowRes;
        els.resistiveResidual[d.uN] += -flowRes;
        els.resistiveResidual[d.uFlow] += valueRes;
    }

    // No limiting, so nothing to load for limited residual

    // Maximal residual contribution
    if (els.maxResistiveResidualContribution) {
        auto valueContrib = std::abs(-els.oldSolution[d.uP] + els.oldSolution[d.uN] + sourceFactor*d.v);
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

    // DC increment residual
    if (els.dcIncrementResidual) { 
        els.dcIncrementResidual[d.uFlow] += p.mag;
    }

    // AC residual
    if (els.acResidual) {
        double re = p.mag*std::cos(p.phase*PI/180);
        double im = p.mag*std::sin(p.phase*PI/180);
        els.acResidual[d.uFlow] += Complex(re, im);
    }

    // Set maximal frequency
    if (els.computeMaxFreq) {
        if (d.typeCode==IndependentSourceType::Sine) {
            els.setMaxFreq(p.freq);
        }
    }

    return true;
}

template<> bool BuiltinISourceInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els, Status& s) {
    auto& p = params.core();
    auto& d = data.core();
    auto& internals = circuit.simulatorInternals();
    auto sourceFactor = internals.sourcescalefactor;
    
    // Evaluate
    if (!els.skipEvaluation) {
        if (els.evaluateResistiveResidual || els.computeNextBreakpoint) {
            auto [val, nextBreakpoint] = sourceCompute(p, d, internals.time);
            if (els.evaluateResistiveResidual) {
                d.i = val; // current of one parallel instance
                d.v = els.oldSolution[d.uP] - els.oldSolution[d.uN]; // voltage across instance
            }
            if (els.computeNextBreakpoint) {
                els.setBreakPoint(nextBreakpoint, internals); 
            }
        } 
    }

    // No Jacobian contribution

    // Load resistive residual
    if (els.resistiveResidual) {
        auto valueRes = sourceFactor*p.mfactor*d.i;
        els.resistiveResidual[d.uP] += valueRes;
        els.resistiveResidual[d.uN] += -valueRes;
    }

    // No limiting, so nothing to load for limited residual

    // Maximal residual contribution
    if (els.maxResistiveResidualContribution) {
        auto valueContrib = std::abs(sourceFactor*p.mfactor*d.i); 
        if (els.maxResistiveResidualContribution[d.uP]<valueContrib) {
            els.maxResistiveResidualContribution[d.uP] = valueContrib;
        }
        if (els.maxResistiveResidualContribution[d.uN]<valueContrib) {
            els.maxResistiveResidualContribution[d.uN] = valueContrib;
        }
    }

    // No reactive component, reactive residual derivative wrt. time is zero

    // DC increment residual
    if (els.dcIncrementResidual) { 
        els.dcIncrementResidual[d.uP] += p.mfactor*p.mag;
        els.dcIncrementResidual[d.uN] += -p.mfactor*p.mag;
    }

    // AC residual
    if (els.acResidual) {
        double re = p.mfactor*p.mag*std::cos(p.phase*PI/180);
        double im = p.mfactor*p.mag*std::sin(p.phase*PI/180);
        els.acResidual[d.uP] += Complex(re, im);
        els.acResidual[d.uN] += -Complex(re, im);
    }

    // Set maximal frequency
    if (els.computeMaxFreq) {
        if (d.typeCode==IndependentSourceType::Sine) {
            els.setMaxFreq(p.freq);
        }
    }

    return true;
}


}
