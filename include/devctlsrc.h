#ifndef __DEVCTLSRC_DEFINED
#define __DEVCTLSRC_DEFINED

#include "value.h"
#include "devbuiltin.h"
#include "common.h"


namespace NAMESPACE {

struct DevCtlSourceModelParams {
    DevCtlSourceModelParams();
};

struct DevVctlSourceInstanceParams {
    // $mfactor is a real parameter in OSDI, therefore we also make it real
    Real mfactor {1}; // Number of parallel instances 
    Real gain {1.0};
    
    DevVctlSourceInstanceParams();
};

struct DevCctlSourceInstanceParams {
    // $mfactor is a real parameter in OSDI, therefore we also make it real
    Real mfactor {1}; // Number of parallel instances 
    Real gain {1.0};
    Id ctlinst {""}; // Controlling instance
    Id ctlnode {"flow(br)"}; // Internal node name within controlling instance
    
    DevCctlSourceInstanceParams();
};


struct DevVccsInstanceData {
    UnknownIndex uP;
    UnknownIndex uN;
    UnknownIndex uCp;
    UnknownIndex uCn;
    double* jacPCp;
    double* jacPCn;
    double* jacNCp;
    double* jacNCn;

    double flowResidual;
    
    double ctl; // Controlling value
    double v;   // Voltage across instance
    double i;   // Current of one parallel instances
    
    DevVccsInstanceData();
};


struct DevVcvsInstanceData {
    UnknownIndex uP;
    UnknownIndex uN;
    UnknownIndex uCp;
    UnknownIndex uCn;
    UnknownIndex uFlow;
    double* jacPFlow;
    double* jacNFlow;
    double* jacFlowP;
    double* jacFlowN;
    double* jacFlowCp;
    double* jacFlowCn;

    double flowResidual;
    double eqResidual;

    double ctl; // Controlling value
    double v;   // Voltage across instance
    double i;   // Current of one parallel instances
    
    DevVcvsInstanceData();
};


struct DevCccsInstanceData {
    UnknownIndex uP;
    UnknownIndex uN;
    UnknownIndex uCtl;
    double* jacPCtl;
    double* jacNCtl;

    double flowResidual;

    double ctl; // Controlling value
    double v;   // Voltage across instance
    double i;   // Current of one parallel instances
    
    DevCccsInstanceData();
};


struct DevCcvsInstanceData {
    UnknownIndex uP;
    UnknownIndex uN;
    UnknownIndex uFlow;
    UnknownIndex uCtl;
    double* jacPFlow;
    double* jacNFlow;
    double* jacFlowP;
    double* jacFlowN;
    double* jacFlowCtl;

    double flowResidual;
    double eqResidual;
    
    double ctl; // Controlling value
    double v;   // Voltage across instance
    double i;   // Current of one parallel instances
    
    DevCcvsInstanceData();
};

using BuiltinVccs = BuiltinDevice<DevCtlSourceModelParams, DevVctlSourceInstanceParams, DevVccsInstanceData>;
using BuiltinVccsModel = BuiltinModel<DevCtlSourceModelParams, DevVctlSourceInstanceParams, DevVccsInstanceData>;
using BuiltinVccsInstance = BuiltinInstance<DevCtlSourceModelParams, DevVctlSourceInstanceParams, DevVccsInstanceData>;

using BuiltinVcvs = BuiltinDevice<DevCtlSourceModelParams, DevVctlSourceInstanceParams, DevVcvsInstanceData>;
using BuiltinVcvsModel = BuiltinModel<DevCtlSourceModelParams, DevVctlSourceInstanceParams, DevVcvsInstanceData>;
using BuiltinVcvsInstance = BuiltinInstance<DevCtlSourceModelParams, DevVctlSourceInstanceParams, DevVcvsInstanceData>;

using BuiltinCccs = BuiltinDevice<DevCtlSourceModelParams, DevCctlSourceInstanceParams, DevCccsInstanceData>;
using BuiltinCccsModel = BuiltinModel<DevCtlSourceModelParams, DevCctlSourceInstanceParams, DevCccsInstanceData>;
using BuiltinCccsInstance = BuiltinInstance<DevCtlSourceModelParams, DevCctlSourceInstanceParams, DevCccsInstanceData>;

using BuiltinCcvs = BuiltinDevice<DevCtlSourceModelParams, DevCctlSourceInstanceParams, DevCcvsInstanceData>;
using BuiltinCcvsModel = BuiltinModel<DevCtlSourceModelParams, DevCctlSourceInstanceParams, DevCcvsInstanceData>;
using BuiltinCcvsInstance = BuiltinInstance<DevCtlSourceModelParams, DevCctlSourceInstanceParams, DevCcvsInstanceData>;

// Specializations of methods - need to write these prototypes otherwise 
// gcc will optimize them out in Release build

// VCCS
template<> void BuiltinVccs::defineInternals();
template<> std::tuple<ParameterIndex, bool> BuiltinVccsInstance::principalParameterIndex() const;
template<> bool BuiltinVccsInstance::deleteHierarchy(Circuit& circuit, Status& s);
template<> bool BuiltinVccsInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s);
template<> bool BuiltinVccsInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const;
template<> std::tuple<bool, OutputSource> BuiltinVccsInstance::opvarOutputSource(ParameterIndex ndx) const;
template<> bool BuiltinVccsInstance::populateStructuresCore(Circuit& circuit, Status& s);
template<> bool BuiltinVccsInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, 
    KluMatrixAccess* matReact, Component compReact, 
    Status& s
);
template<> bool BuiltinVccsInstance::evalCore(Circuit& circuit, EvalSetup& els);
template<> bool BuiltinVccsInstance::loadCore(Circuit& circuit, LoadSetup& els);
template<> bool BuiltinVccsInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els);

// VCVS
template<> void BuiltinVcvs::defineInternals();
template<> std::tuple<ParameterIndex, bool> BuiltinVcvsInstance::principalParameterIndex() const;
template<> bool BuiltinVcvsInstance::deleteHierarchy(Circuit& circuit, Status& s);
template<> bool BuiltinVcvsInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s);
template<> bool BuiltinVcvsInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const;
template<> std::tuple<bool, OutputSource> BuiltinVcvsInstance::opvarOutputSource(ParameterIndex ndx) const;
template<> bool BuiltinVcvsInstance::populateStructuresCore(Circuit& circuit, Status& s);
template<> bool BuiltinVcvsInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, 
    KluMatrixAccess* matReact, Component compReact, 
    Status& s
);
template<> bool BuiltinVcvsInstance::evalCore(Circuit& circuit, EvalSetup& els);
template<> bool BuiltinVcvsInstance::loadCore(Circuit& circuit, LoadSetup& els);
template<> bool BuiltinVcvsInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els);

// CCCS
template<> void BuiltinCccs::defineInternals();
template<> std::tuple<ParameterIndex, bool> BuiltinCccsInstance::principalParameterIndex() const;
template<> bool BuiltinCccsInstance::deleteHierarchy(Circuit& circuit, Status& s);
template<> bool BuiltinCccsInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s);
template<> bool BuiltinCccsInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const;
template<> std::tuple<bool, OutputSource> BuiltinCccsInstance::opvarOutputSource(ParameterIndex ndx) const;
template<> bool BuiltinCccsInstance::populateStructuresCore(Circuit& circuit, Status& s);
template<> bool BuiltinCccsInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, 
    KluMatrixAccess* matReact, Component compReact, 
    Status& s
);
template<> bool BuiltinCccsInstance::evalCore(Circuit& circuit, EvalSetup& els);
template<> bool BuiltinCccsInstance::loadCore(Circuit& circuit, LoadSetup& els);
template<> bool BuiltinCccsInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els);

// CCVS
template<> void BuiltinCcvs::defineInternals();
template<> std::tuple<ParameterIndex, bool> BuiltinCcvsInstance::principalParameterIndex() const; 
template<> bool BuiltinCcvsInstance::deleteHierarchy(Circuit& circuit, Status& s);
template<> bool BuiltinCcvsInstance::buildHierarchy(Circuit& circuit, RpnEvaluator& evaluator, InstantiationData& idata, Status& s);
template<> bool BuiltinCcvsInstance::getOpvar(ParameterIndex ndx, Value& v, Status& s) const;
template<> std::tuple<bool, OutputSource> BuiltinCcvsInstance::opvarOutputSource(ParameterIndex ndx) const;
template<> bool BuiltinCcvsInstance::populateStructuresCore(Circuit& circuit, Status& s);
template<> bool BuiltinCcvsInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, 
    KluMatrixAccess* matReact, Component compReact, 
    Status& s
);
template<> bool BuiltinCcvsInstance::evalCore(Circuit& circuit, EvalSetup& els);
template<> bool BuiltinCcvsInstance::loadCore(Circuit& circuit, LoadSetup& els);
template<> bool BuiltinCcvsInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els);


struct DevMutualInstanceParams {
    Real k {0.0};
    Id ind1 {""}; // inductor1
    Id ind2 {""}; // inductor2
    Id ctlnode1 {"flow(br)"}; 
    Id ctlnode2 {"flow(br)"}; 
    
    DevMutualInstanceParams();
};

struct DevMutualInstanceData {
    Instance* ind1;
    Instance* ind2;
    Node* ctlNode1 {nullptr};
    Node* ctlNode2 {nullptr};
    UnknownIndex uFlow1 {0};
    UnknownIndex uFlow2 {0};
    double mutual;
    double* jacReact12;
    double* jacReact21;
    double reacRes1;
    double reacRes2;

    GlobalStorageIndex offsStates;
    
    DevMutualInstanceData();
};

using BuiltinMutual = BuiltinDevice<DevCtlSourceModelParams, DevMutualInstanceParams, DevMutualInstanceData>;
using BuiltinMutualModel = BuiltinModel<DevCtlSourceModelParams, DevMutualInstanceParams, DevMutualInstanceData>;
using BuiltinMutualInstance = BuiltinInstance<DevCtlSourceModelParams, DevMutualInstanceParams, DevMutualInstanceData>;

// Specializations of methods - need to write these prototypes otherwise 
// gcc will optimize them out in Release build

template<> void BuiltinMutual::defineInternals();
template<> std::tuple<ParameterIndex, bool> BuiltinMutualInstance::principalParameterIndex() const; 
template<> std::tuple<bool, bool, bool> BuiltinMutualInstance::setupWorker(Circuit& circuit, Status& s);
template<> bool BuiltinMutualInstance::preAnalysisWorker(Circuit& circuit, Status& s);
template<> bool BuiltinMutualInstance::populateStructuresCore(Circuit& circuit, Status& s);
template<> bool BuiltinMutualInstance::bindCore(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, 
    KluMatrixAccess* matReact, Component compReact, 
    Status& s
);
template<> bool BuiltinMutualInstance::evalCore(Circuit& circuit, EvalSetup& evalSetup);
template<> bool BuiltinMutualInstance::loadCore(Circuit& circuit, LoadSetup& loadSetup);
template<> bool BuiltinMutualInstance::evalAndLoadCore(Circuit& circuit, EvalAndLoadSetup& els);

}

#endif

