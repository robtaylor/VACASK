#ifndef __DEVVISRC_DEFINED
#define __DEVVISRC_DEFINED

#include "value.h"
#include "devbuiltin.h"
#include "common.h"


namespace NAMESPACE {

enum class IndependentSourceType : char { Dc, Sine, Pulse, Sffm, Exp, Pwl };

struct DevSourceModelParams {
    DevSourceModelParams();
};

struct DevSourceInstanceParams {
    // $mfactor is a real parameter in OSDI, therefore we also make it real
    Real mfactor; // Number of parallel instances 
    Id type;
    Real delay;
    // DC
    Real dc;
    // Pulse
    Real val0;
    Real val1;
    Real period;
    Real rise;
    Real fall;
    Real width;
    // Sine
    Real sinedc;
    Real ampl;
    Real freq;
    Real sinephase; // degrees
    // Exp
    Real td1;
    Real tau1;
    Real td2;
    Real tau2;
    Real expperiod;
    // Pwl
    RealVector wave;
    Real offset;
    Real scale;
    Real stretch;
    // AC, DC incremental
    Real mag;
    Real phase; // degrees (only for AC)

    DevSourceInstanceParams();
};


struct DevVSourceInstanceData {
    IndependentSourceType typeCode;

    UnknownIndex uP;
    UnknownIndex uN;
    UnknownIndex uFlow;
    double* jacFlowP;
    double* jacFlowN;
    double* jacPFlow;
    double* jacNFlow;

    double v;  // Voltage across instance
    double i;  // Current of one parallel instances

    DevVSourceInstanceData();
};

struct DevISourceInstanceData {
    IndependentSourceType typeCode;

    UnknownIndex uP;
    UnknownIndex uN;
    
    double i;  // Current of one parallel instance
    double v;  // Voltage across instance

    DevISourceInstanceData();
};

template<typename InstanceParams, typename InstanceData> 
std::tuple<double, double> sourceCompute(const InstanceParams& params, InstanceData& data, double time);


using BuiltinVSource = BuiltinDevice<DevSourceModelParams, DevSourceInstanceParams, DevVSourceInstanceData>;
using BuiltinVSourceModel = BuiltinModel<DevSourceModelParams, DevSourceInstanceParams, DevVSourceInstanceData>;
using BuiltinVSourceInstance = BuiltinInstance<DevSourceModelParams, DevSourceInstanceParams, DevVSourceInstanceData>;

using BuiltinISource = BuiltinDevice<DevSourceModelParams, DevSourceInstanceParams, DevISourceInstanceData>;
using BuiltinISourceModel = BuiltinModel<DevSourceModelParams, DevSourceInstanceParams, DevISourceInstanceData>;
using BuiltinISourceInstance = BuiltinInstance<DevSourceModelParams, DevSourceInstanceParams, DevISourceInstanceData>;


// Specializations of methods - need to write these prototypes otherwise 
// gcc will optimize them out in Release build

template<> void BuiltinVSource::defineInternals();
template<> void BuiltinISource::defineInternals();
template<> bool BuiltinVSource::isSource() const;
template<> bool BuiltinISource::isSource() const;
template<> bool BuiltinVSource::isVoltageSource() const;
template<> std::tuple<ParameterIndex, bool> BuiltinVSourceInstance::principalParameterIndex() const;
template<> std::tuple<ParameterIndex, bool> BuiltinISourceInstance::principalParameterIndex() const;
template<> bool BuiltinVSourceInstance::deleteHierarchy(Circuit& circuit, Status& s);
template<> bool BuiltinVSourceInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s);
template<> bool BuiltinISourceInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s);
template<> std::tuple<EquationIndex,EquationIndex> BuiltinVSourceInstance::sourceExcitation(Circuit& circuit) const;
template<> std::tuple<UnknownIndex,UnknownIndex> BuiltinVSourceInstance::sourceResponse(Circuit& circuit) const;
template<> double BuiltinVSourceInstance::scaledUnityExcitation() const;
template<> double BuiltinVSourceInstance::responseScalingFactor() const;
template<> std::tuple<EquationIndex,EquationIndex> BuiltinISourceInstance::sourceExcitation(Circuit& circuit) const;
template<> std::tuple<UnknownIndex,UnknownIndex> BuiltinISourceInstance::sourceResponse(Circuit& circuit) const;
template<> double BuiltinISourceInstance::scaledUnityExcitation() const;
template<> double BuiltinISourceInstance::responseScalingFactor() const;
template<> bool BuiltinVSourceInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const;
template<> bool BuiltinISourceInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const;
template<> std::tuple<bool, OutputSource> BuiltinVSourceInstance::opvarOutputSource(ParameterIndex ndx, Status& s) const;
template<> std::tuple<bool, OutputSource> BuiltinISourceInstance::opvarOutputSource(ParameterIndex ndx, Status& s) const;
template<> std::tuple<bool, bool, bool> BuiltinVSourceInstance::setupWorker(Circuit& circuit, Status& s);
template<> std::tuple<bool, bool, bool> BuiltinISourceInstance::setupWorker(Circuit& circuit, Status& s);
template<> bool BuiltinVSourceInstance::populateStructuresCore(Circuit& circuit, Status& s);
template<> bool BuiltinISourceInstance::populateStructuresCore(Circuit& circuit, Status& s);
template<> bool BuiltinVSourceInstance::bindCore(
    Circuit& circuit, 
    KluRealMatrix* matResistReal, KluComplexMatrix* matResistCx, Component compResist, 
    KluRealMatrix* matReactReal, KluComplexMatrix* matReactCx, Component compReact, 
    Status& s
);
template<> bool BuiltinISourceInstance::bindCore(
    Circuit& circuit, 
    KluRealMatrix* matResistReal, KluComplexMatrix* matResistCx, Component compResist, 
    KluRealMatrix* matReactReal, KluComplexMatrix* matReactCx, Component compReact, 
    Status& s
);
template<> bool BuiltinVSourceInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els);
template<> bool BuiltinISourceInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els);

}

#endif

