#ifndef __DEVVISRC_DEFINED
#define __DEVVISRC_DEFINED

#include "value.h"
#include "devbuiltin.h"
#include "common.h"


namespace NAMESPACE {

enum class IndependentSourceType : char { Dc, Sine, Pulse, Sffm, Exp, Pwl, Am, Fm };

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
    Real theta;
    // Exp
    // val0, val1, delay
    Real td2;
    Real tau1;
    Real tau2;
    // Pwl
    RealVector wave;
    Real offset;
    Real scale;
    Real stretch;
    Real pwlperiod;
    Real twidth;
    Id allbrkpts; // yes, no, auto
    Real slopetol;
    Real reltol;
    // AM, FM
    Real modfreq;
    Real modphase;
    Real modindex; 
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

    double flowResidual;
    double eqResidual;

    double v;  // Voltage across instance
    double i;  // Current of one parallel instances

    DevVSourceInstanceData();
};

struct DevISourceInstanceData {
    IndependentSourceType typeCode;

    UnknownIndex uP;
    UnknownIndex uN;
    
    double flowResidual;
    
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
template<> std::tuple<bool, OutputSource> BuiltinVSourceInstance::opvarOutputSource(ParameterIndex ndx) const;
template<> std::tuple<bool, OutputSource> BuiltinISourceInstance::opvarOutputSource(ParameterIndex ndx) const;
template<> std::tuple<bool, bool, bool> BuiltinVSourceInstance::setupCore(Circuit& circuit, CommonData& commons, DeviceRequests* devReq, Status& s);
template<> std::tuple<bool, bool, bool> BuiltinISourceInstance::setupCore(Circuit& circuit, CommonData& commons, DeviceRequests* devReq, Status& s);
template<> bool BuiltinVSourceInstance::setStaticTolerancesCore(Circuit& circuit, CommonData& commons, Status& s);
template<> bool BuiltinISourceInstance::setStaticTolerancesCore(Circuit& circuit, CommonData& commons, Status& s);
template<> bool BuiltinVSourceInstance::populateStructuresCore(Circuit& circuit, Status& s);
template<> bool BuiltinISourceInstance::populateStructuresCore(Circuit& circuit, Status& s);
template<> bool BuiltinVSourceInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, const std::optional<MatrixEntryPosition>& mepResist, 
    KluMatrixAccess* matReact, Component compReact, const std::optional<MatrixEntryPosition>& mepReact, 
    Status& s
);
template<> bool BuiltinISourceInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, const std::optional<MatrixEntryPosition>& mepResist, 
    KluMatrixAccess* matReact, Component compReact, const std::optional<MatrixEntryPosition>& mepReact, 
    Status& s
);

template<> bool BuiltinVSourceInstance::evalCore(Circuit& circuit, CommonData& commons, EvalSetup& evalSetup);    
template<> bool BuiltinISourceInstance::evalCore(Circuit& circuit, CommonData& commons, EvalSetup& evalSetup);    
template<> bool BuiltinVSourceInstance::loadCore(Circuit& circuit, CommonData& commons, LoadSetup& loadSetup);    
template<> bool BuiltinISourceInstance::loadCore(Circuit& circuit, CommonData& commons, LoadSetup& loadSetup);    

}

#endif

