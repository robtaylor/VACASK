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
    registerMember(theta);

    registerMember(td2);
    registerMember(tau1);
    registerMember(tau2);
    
    registerMember(wave);
    registerMember(offset);
    registerMember(scale);
    registerMember(stretch);
    registerMember(pwlperiod);
    registerMember(twidth);
    registerMember(allbrkpts);
    registerMember(slopetol);
    registerMember(reltol);
    
    registerMember(modfreq);
    registerMember(modphase);
    registerMember(modindex);
    
    registerMember(mag);
    registerMember(phase);
    return 0;
}
instantiateIntrospection(DevSourceInstanceParams);

DevSourceInstanceParams::DevSourceInstanceParams() {
    mfactor = 1;

    type = "dc";
    delay = 0.0;
    
    // type="dc"
    dc = 0.0;

    // type="pulse"
    val0 = 0.0;
    val1 = 1.0;
    period = 0.0;
    rise = 1e-9;
    fall = 0.0;
    width = 0.0;

    // type="sine"
    sinedc = 0.0;
    ampl = 1.0;
    freq = 1e3;
    sinephase = 0.0;
    theta = 0.0;

    // type="exp"
    td2 = 0.0;
    tau1 = 0.0;
    tau2 = 0.0;

    // type="pwl"
    wave = RealVector();
    offset = 0.0;
    scale = 1.0;
    stretch = 1.0;
    pwlperiod = 0.0;
    twidth = 0.0;
    allbrkpts = "auto";
    slopetol = 0.0;
    reltol = 0.1;

    // type="am" or "fm"
    // sinedc, ampl, freq, sinephase
    modfreq = 1e3;
    modphase = 0;
    modindex = 0.5; 
    
    // small signal parameters
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
static Id typeExp = Id::createStatic("exp");
static Id typeDc = Id::createStatic("dc");
static Id typePwl = Id::createStatic("pwl");
static Id typeAm = Id::createStatic("am");
static Id typeFm = Id::createStatic("fm");

static Id valYes = Id::createStatic("yes");
static Id valNo = Id::createStatic("no");
static Id valAuto = Id::createStatic("auto");

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
    } else if (p.type == typePwl) {
        d.typeCode = IndependentSourceType::Pwl;
        // No points .. error
        if (p.wave.size()==0) {
            s.set(Status::BadArguments, "Pwl waveform needs at least one point.");
            s.extend(loc);
            return std::make_tuple(false, false, false);
        }
        
        // Need pairs of values
        if (p.wave.size()%2 == 1) {
            s.set(Status::BadArguments, "Pwl waveform needs an even number of values.");
            s.extend(loc);
            return std::make_tuple(false, false, false);
        }
        
        // Check time scale monotonicity, extract endpoints and length

        // If there are at least 2 point, check scale, stretch, pwlperiod, twidth

        // Check slopetol and reltol

        // Breakpoint handling mechanism
    } else if (p.type == typeExp) {
        d.typeCode = IndependentSourceType::Exp;
        if (p.td2<=0) {
            s.set(Status::BadArguments, "Parameter td2 of exponential transient must be greater than 0.");
            s.extend(loc);
            return std::make_tuple(false, false, false);
        }
    } else if (p.type == typeAm) {
        d.typeCode = IndependentSourceType::Am;
    } else if (p.type == typeFm) {
        d.typeCode = IndependentSourceType::Fm;
    } else if (p.type == typeDc) {
        d.typeCode = IndependentSourceType::Dc;
    } else {
        s.set(Status::BadArguments, "Unknown transient waveform type.");
        s.extend(loc);
        return std::make_tuple(false, false, false);
    }

    return std::make_tuple(true, false, false); 
}

// A device model should not rely on tolerances. 
// Its only job is to produce consistent reponses, i.e. 
// in this case t5 should match t1 in the next period. 
template<typename InstanceParams, typename InstanceData> 
std::tuple<double, double> sourceCompute(const InstanceParams& params, InstanceData& data, double time) {
    double val;
    double nextBreak = std::numeric_limits<double>::infinity();

    switch (data.typeCode) {
    case IndependentSourceType::Dc:
        val = params.dc;
        break;
    case IndependentSourceType::Sine:
        if (time<params.delay) {
            // For t < delay the value is equal to value at t=delay
            val = params.sinedc+params.ampl*std::sin(params.phase*PI/180);
            nextBreak = params.delay;
        } else {
            // For t >= delay start sine at given phase
            val = params.sinedc+
                  params.ampl
                    *std::sin(2*PI*params.freq*(time-params.delay)+params.phase*PI/180)
                    *std::exp(-params.theta*(time-params.delay));
        }
        break;
    case IndependentSourceType::Exp:
        if (time<params.delay) {
            // For t < delay the value is equal to value at t=delay
            val = params.val0;
            nextBreak = params.delay;
        } else {
            auto t1 = params.delay + params.td2; // start of fall
            if (time<t1) {
                // Rising exponential
                val = params.val0 
                      + (params.val1-params.val0)*(1-std::exp(-(time-params.delay)/params.tau1));
                nextBreak = t1;
            } else {
                // Falling exponential
                val = params.val0 
                      + (params.val1-params.val0)*(1-std::exp(-(time-params.delay)/params.tau1))
                      + (params.val0-params.val1)*(1-std::exp(-(time-t1)/params.tau2));
            }
        }
        break;
    case IndependentSourceType::Am:
        if (time<params.delay) {
            // For t < delay the value is equal to value at t=delay
            val = params.sinedc+params.ampl*std::sin(params.phase*PI/180)*(
                1+params.modindex*std::sin(params.modphase*PI/180)
            );
            nextBreak = params.delay;
        } else {
            // For t >= delay start sine at given phase
            val = params.sinedc+params.ampl*std::sin(2*PI*params.freq*(time-params.delay)+params.phase*PI/180)*(
                1+params.modindex*std::sin(2*PI*params.modfreq*(time-params.delay)+params.modphase*PI/180)
            );
        }
        break;
    case IndependentSourceType::Fm:
        if (time<params.delay) {
            // For t < delay the value is equal to value at t=delay
            val = params.sinedc+params.ampl*std::sin(
                params.phase*PI/180+
                params.modindex*std::sin(params.modphase*PI/180)
            );
            nextBreak = params.delay;
        } else {
            // For t >= delay start sine at given phase
            val = params.sinedc+params.ampl*std::sin(
                2*PI*params.freq*(time-params.delay)+params.phase*PI/180+
                params.modindex*std::sin(2*PI*params.modfreq*(time-params.delay)+params.modphase*PI/180)
            );
        }
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
            int64_t cycle = 0; 
            // Start of current period
            if (params.period<=0) {
                // Not periodic
                basetime = params.delay;
            } else {
                // Periodic
                cycle = static_cast<int64_t>(std::floor((time-params.delay)/params.period));
                basetime = params.delay + cycle*params.period;
            }
            // Significant time points
            auto t0 = basetime; // start of rise
            auto t1 = t0 + params.rise; // end of rise, start of top
            auto t2 = t1 + params.width; // start of fall
            auto t3 = t2 + params.fall; // end of fall
            auto t4 = t0 + params.period; // end of period
            auto t5 = params.delay + (cycle+1)*params.period + params.rise; // next period, end of rise
            // Relative time since start of period
            double reltime = time - basetime;
            // Simulator::dbg().setf(std::ios::scientific, std::ios::floatfield);
            if (time<t1) {
                // Rising
                val = (time-t0)/params.rise*(params.val1-params.val0)+params.val0;
                nextBreak = t1;
            } else if (time<t2) {
                // On top
                val = params.val1;
                // Set next break only if width>0
                if (params.width>0) {
                    nextBreak = t2;
                }
            } else if (params.fall<=0) {
                // Fall not set, stay on top, no breakpoint
                val = params.val1;
            } else if (time<t3) {
                // Falling
                val = (time-t2)/params.fall*(params.val0-params.val1)+params.val1;
                nextBreak = t3;
            } else if (time<t4) {
                // After fall, back on base level
                val = params.val0;
                // Beakpoint only if period is set
                if (params.period>0) {
                    nextBreak = t4;
                }
            } else {
                if (params.period>0) {
                    // Periodic waveform, rising flank of next period
                    // We may end up here for the first point of next period due to tolerances, 
                    // Compute rising flank with origin at t4
                    val = (time-t4)/params.rise*(params.val1-params.val0)+params.val0;
                    nextBreak = t5;
                } else {
                    // Not periodic, back on base level after fall
                    val = params.val0;
                }
            }
        }
        break;
    }
    
    // Simulator::dbg() << "val=" << val << " next break=" << nextBreak << "\n";
    // Simulator::dbg() << "t=" << time << " next break=" << nextBreak << "\n";
    
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

template<> std::tuple<bool, OutputSource> BuiltinVSourceInstance::opvarOutputSource(ParameterIndex ndx) const { 
    switch (ndx) {
    case 0:
        return std::make_tuple(true, OutputSource(&data.core().v));
    case 1: 
        return std::make_tuple(true, OutputSource(&data.core().i));
    default: 
        return std::make_tuple(false, OutputSource());
    }
}

template<> std::tuple<bool, OutputSource> BuiltinISourceInstance::opvarOutputSource(ParameterIndex ndx) const { 
    switch (ndx) {
    case 0:
        return std::make_tuple(true, OutputSource(&data.core().v));
    case 1: 
        return std::make_tuple(true, OutputSource(&data.core().i));
    default: 
        return std::make_tuple(false, OutputSource());
    }
}

template<> std::tuple<bool, bool, bool> BuiltinVSourceInstance::setupWorker(Circuit& circuit, DeviceRequests* devReq, Status& s) {
    clearFlags(Flags::NeedsSetup); 
    return sourceSetup(params, data, location(), circuit, s);
}; 

template<> std::tuple<bool, bool, bool> BuiltinISourceInstance::setupWorker(Circuit& circuit, DeviceRequests* devReq, Status& s) { 
    clearFlags(Flags::NeedsSetup);
    return sourceSetup(params, data, location(), circuit, s);
}; 

template<> bool BuiltinVSourceInstance::populateStructuresCore(Circuit& circuit, Status& s) {
    // Create Jacobian entries
    if (auto [_, ok] = circuit.createJacobianEntry(nodes_[0], nodes_[2], EntryFlags::Resistive, s); !ok) {
        return false;
    }
    if (auto [_, ok] = circuit.createJacobianEntry(nodes_[1], nodes_[2], EntryFlags::Resistive, s); !ok) {
        return false;
    }
    if (auto [_, ok] = circuit.createJacobianEntry(nodes_[2], nodes_[0], EntryFlags::Resistive, s); !ok) {
        return false;
    }
    if (auto [_, ok] = circuit.createJacobianEntry(nodes_[2], nodes_[1], EntryFlags::Resistive, s); !ok) {
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
    KluMatrixAccess* matResist, Component compResist, const std::optional<MatrixEntryPosition>& mepResist, 
    KluMatrixAccess* matReact, Component compReact, const std::optional<MatrixEntryPosition>& mepReact, 
    Status& s
) {
    auto& d = data.core();

    // Unknown indices
    d.uP = nodes_[0]->unknownIndex();
    d.uN = nodes_[1]->unknownIndex();
    d.uFlow = nodes_[2]->unknownIndex();

    // Resistive Jacobian entry pointers
    if (matResist) {
        jacEntryPtr(d.jacPFlow, d.uP, d.uFlow, matResist, compResist, mepResist);
        jacEntryPtr(d.jacNFlow, d.uN, d.uFlow, matResist, compResist, mepResist);
        jacEntryPtr(d.jacFlowP, d.uFlow, d.uP, matResist, compResist, mepResist);
        jacEntryPtr(d.jacFlowN, d.uFlow, d.uN, matResist, compResist, mepResist);
    }

    // No reactive Jacobian entries
    
    return true;
}

template<> bool BuiltinISourceInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, const std::optional<MatrixEntryPosition>& mepResist, 
    KluMatrixAccess* matReact, Component compReact, const std::optional<MatrixEntryPosition>& mepReact, 
    Status& s
) {
    auto& d = data.core();

    // Unknown indices
    d.uP = nodes_[0]->unknownIndex();
    d.uN = nodes_[1]->unknownIndex();
    
    // No Jacobian entries
    
    return true;
}

template<> bool BuiltinVSourceInstance::evalCore(Circuit& circuit, EvalSetup& evalSetup) {
    auto& p = params.core();
    auto& d = data.core();
    auto& internals = circuit.simulatorInternals();
    auto sourceFactor = internals.sourcescalefactor;
    
    // Evaluate, placeholder for bypass implementation
    auto [val, nextBreakpoint] = sourceCompute(p, d, evalSetup.time);
    if (true) {
        if (evalSetup.evaluateResistiveResidual) {
            if (evalSetup.evaluateResistiveResidual) {
                d.flowResidual = p.mfactor*evalSetup.oldSolution[d.uFlow];
                d.eqResidual = -evalSetup.oldSolution[d.uP] + evalSetup.oldSolution[d.uN] + sourceFactor*val;
            }
            // Opvars
            d.v = sourceFactor*val; // mfactor does not affect voltage source value
            d.i = evalSetup.oldSolution[d.uFlow]; // flow across one parallel instance
        }
    }

    // Next breakpoint
    if (evalSetup.computeNextBreakpoint) {
        evalSetup.setBreakPoint(nextBreakpoint, internals); 
    }

    // Set maximal frequency
    if (evalSetup.computeMaxFreq) {
        if (d.typeCode==IndependentSourceType::Sine) {
            evalSetup.setMaxFreq(p.freq);
        }
    }

    return true;
}

template<> bool BuiltinVSourceInstance::loadCore(Circuit& circuit, LoadSetup& loadSetup) {
    auto& p = params.core();
    auto& d = data.core();
    auto& internals = circuit.simulatorInternals();
    auto sourceFactor = internals.sourcescalefactor;
    
    // Load resistive Jacobian, transient load is identical because there is no reactive component
    if (loadSetup.loadResistiveJacobian || loadSetup.loadTransientJacobian) {
        // KCL
        *(d.jacPFlow+loadSetup.jacobianLoadOffset) += p.mfactor;
        *(d.jacNFlow+loadSetup.jacobianLoadOffset) += -p.mfactor;
        // Extra equation
        *(d.jacFlowP+loadSetup.jacobianLoadOffset) += -1.0;
        *(d.jacFlowN+loadSetup.jacobianLoadOffset) += 1.0;
    }

    // Load resistive residual
    if (loadSetup.resistiveResidual) {
        loadSetup.resistiveResidual[d.uP] += d.flowResidual;
        loadSetup.resistiveResidual[d.uN] += -d.flowResidual;
        loadSetup.resistiveResidual[d.uFlow] += d.eqResidual;
    }

    // No limiting, so nothing to load for limited residual

    // Maximal residual contribution
    if (loadSetup.maxResistiveResidualContribution) {
        auto flowContrib = std::abs(d.flowResidual);
        auto eqContrib = std::abs(d.eqResidual);
        if (loadSetup.maxResistiveResidualContribution[d.uP]<flowContrib) {
            loadSetup.maxResistiveResidualContribution[d.uP] = flowContrib;
        }
        if (loadSetup.maxResistiveResidualContribution[d.uN]<flowContrib) {
            loadSetup.maxResistiveResidualContribution[d.uN] = flowContrib;
        }
        if (loadSetup.maxResistiveResidualContribution[d.uFlow]<eqContrib) {
            loadSetup.maxResistiveResidualContribution[d.uFlow] = eqContrib;
        }
    }

    // No reactive component, reactive residual derivative wrt. time is zero

    // DC increment residual
    if (loadSetup.dcIncrementResidual) { 
        loadSetup.dcIncrementResidual[d.uFlow] += p.mag;
    }

    // AC residual
    if (loadSetup.acResidual) {
        double re = p.mag*std::cos(p.phase*PI/180);
        double im = p.mag*std::sin(p.phase*PI/180);
        loadSetup.acResidual[d.uFlow] += Complex(re, im);
    }

    return true;
}

template<> bool BuiltinISourceInstance::evalCore(Circuit& circuit, EvalSetup& evalSetup) {
    auto& p = params.core();
    auto& d = data.core();
    auto& internals = circuit.simulatorInternals();
    auto sourceFactor = internals.sourcescalefactor;
    
    // Evaluate, placeholder for bypass implementation
    auto [val, nextBreakpoint] = sourceCompute(p, d, evalSetup.time);
    if (true) {
        if (evalSetup.evaluateResistiveResidual) {
            if (evalSetup.evaluateResistiveResidual) {
                d.flowResidual = sourceFactor*p.mfactor*val;
                // Opvars
                d.i = sourceFactor*val; // current of one parallel instance
                d.v = evalSetup.oldSolution[d.uP] - evalSetup.oldSolution[d.uN]; // voltage across instance
            }  
        } 
    }

    // Next breakjpoint
    if (evalSetup.computeNextBreakpoint) {
        evalSetup.setBreakPoint(nextBreakpoint, internals); 
    }

    // Set maximal frequency
    if (evalSetup.computeMaxFreq) {
        if (d.typeCode==IndependentSourceType::Sine) {
            evalSetup.setMaxFreq(p.freq);
        }
    }

    return true;
}

template<> bool BuiltinISourceInstance::loadCore(Circuit& circuit, LoadSetup& loadSetup) {
    auto& p = params.core();
    auto& d = data.core();
    auto& internals = circuit.simulatorInternals();
    auto sourceFactor = internals.sourcescalefactor;
    
    // Load resistive residual
    if (loadSetup.resistiveResidual) {
        loadSetup.resistiveResidual[d.uP] += d.flowResidual;
        loadSetup.resistiveResidual[d.uN] += -d.flowResidual;
    }

    // No limiting, so nothing to load for limited residual

    // Maximal residual contribution
    if (loadSetup.maxResistiveResidualContribution) {
        auto flowContrib = std::abs(d.flowResidual); 
        if (loadSetup.maxResistiveResidualContribution[d.uP]<flowContrib) {
            loadSetup.maxResistiveResidualContribution[d.uP] = flowContrib;
        }
        if (loadSetup.maxResistiveResidualContribution[d.uN]<flowContrib) {
            loadSetup.maxResistiveResidualContribution[d.uN] = flowContrib;
        }
    }

    // No reactive component, reactive residual derivative wrt. time is zero

    // DC increment residual
    if (loadSetup.dcIncrementResidual) { 
        loadSetup.dcIncrementResidual[d.uP] += p.mfactor*p.mag;
        loadSetup.dcIncrementResidual[d.uN] += -p.mfactor*p.mag;
    }

    // AC residual
    if (loadSetup.acResidual) {
        double re = p.mfactor*p.mag*std::cos(p.phase*PI/180);
        double im = p.mfactor*p.mag*std::sin(p.phase*PI/180);
        loadSetup.acResidual[d.uP] += Complex(re, im);
        loadSetup.acResidual[d.uN] += -Complex(re, im);
    }

    return true;
}

}
