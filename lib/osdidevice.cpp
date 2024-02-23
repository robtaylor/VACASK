#include <cstddef>
#include <cstdlib>
#include "osdidevice.h"
#include "osdimodel.h"
#include "osdiinstance.h"
#include "parseroutput.h"
#include "circuit.h"
#include "simulator.h"
#include "common.h"


namespace NAMESPACE {

OsdiDevice::OsdiDevice(OsdiFile* of, int descriptorIndex, Id asName, Loc location, Status &s) 
    : Device(asName ? asName : Id(of->deviceDescriptor(descriptorIndex)->name), location), osdiFile(of), index_(descriptorIndex) {
    descriptor_ = of->deviceDescriptor(descriptorIndex);
    setFlags(Flags::IsValid);
}

OsdiDevice::~OsdiDevice() {
}

bool OsdiDevice::operator==(const Device& other) const {
    const OsdiDevice* devOther = dynamic_cast<const OsdiDevice*>(&other);
    if (devOther && devOther->osdiFile==osdiFile && devOther->index_==index_) {
        return true;
    }
    return false;
}

Model* OsdiDevice::createModel(Circuit& circuit, Instance* parentInstance, RpnEvaluator& evaluator, const PTModel& parsedModel, Status& s) {
    auto name = parsedModel.name();

    // If we have a hierarchical parent translate name
    if (parentInstance) {
        name = parentInstance->translate(name);
    }

    // Create model
    auto* model = new OsdiModel(this, name, parentInstance, parsedModel, s);
    if (!model->checkFlags(Model::Flags::IsValid)) {
        s.extend(parsedModel.location());
        delete model;
        return nullptr;
    }

    // Add to modelMap of circuit
    if (!circuit.add(model, s)) {
        delete model;
        return nullptr;
    }

    // Set model's parameters, use the evaluator whose latest context is the parent instance's context
    auto [ok, changed] = model->setParameters(parsedModel.parameters(), evaluator, s);
    if (!ok) {
        return nullptr;
    }

    return model;
}

bool OsdiDevice::freeValues(std::vector<bool>& paramGiven, void* coreMod, void* coreInst) {
    OsdiFile::OsdiFlags flags = ACCESS_FLAG_SET;
    if (coreInst) {
        flags |= ACCESS_FLAG_INSTANCE;
    }
    
    if (coreInst) {
        for(auto [osdiId, givenIndex] : osdiFile->instanceParamemeterGivenMap(index_)) {
            if (!paramGiven[givenIndex]) {
                continue;
            }
            auto t = parameterType(osdiId);
            switch (t) {
                case Value::Type::String: {
                    auto ptr = (char**)(descriptor_->access(coreInst, coreMod, osdiId, flags));
                    free(*ptr);
                    break;
                }
            }
        }
    } else {
        for(auto [osdiId, givenIndex] : osdiFile->modelParamemeterGivenMap(index_)) {
            if (!paramGiven[givenIndex]) {
                continue;
            }
            auto t = parameterType(osdiId);
            switch (t) {
                case Value::Type::String: {
                    auto ptr = (char**)(descriptor_->access(nullptr, coreMod, osdiId, flags));
                    free(*ptr);
                    break;
                }
            }
        }
    }
    return true;
}

std::tuple<bool,bool> OsdiDevice::writeParameter(OsdiFile::OsdiParameterId osdiId, void* coreMod, void* coreInst, const Value& v, bool wasGiven, Status& s) {
    // Check index
    if (osdiId>=osdiIdCount()) {
        s.set(Status::Range, std::string("OSDI parameter id=")+std::to_string(osdiId)+" out of range.");
        return std::make_tuple(false, false);
    }

    // Check kind
    if (!coreInst && !isModelParameter(osdiId)) {
        s.set(Status::NotFound, std::string("OSDI parameter id=")+std::to_string(osdiId)+" is not a model parameter.");
        return std::make_tuple(false, false);
    }

    if (coreInst && !isInstanceParameter(osdiId)) {
        s.set(Status::NotFound, std::string("OSDI parameter id=")+std::to_string(osdiId)+" is not an instance parameter.");
        return std::make_tuple(false, false);
    }

    // Get type
    auto t = parameterType(osdiId);
    if (t!=Value::Type::Int && t!=Value::Type::Real && t!=Value::Type::String) {
        s.set(Status::Unsupported, std::string("OSDI parameter id=")+std::to_string(osdiId)+" has unsupported type ("+Value::typeCodeToName(t)+").");
        return std::make_tuple(false, false);
    }

    // Convert value
    const Value* vwrite = &v;
    Value vconv;
    if (v.type()!=t) {
        vconv = v;
        if (!vconv.convertInPlace(t, s)) {
            s.extend(std::string("Value conversion failed for OSDI parameter id=")+std::to_string(osdiId)+".");
            return std::make_tuple(false, false);
        }
        vwrite = &vconv;
    }

    // Write
    OsdiFile::OsdiFlags flags = ACCESS_FLAG_SET;
    if (coreInst) {
        flags |= ACCESS_FLAG_INSTANCE;
    }
    bool changed = false;
    switch (t) {
        case Value::Type::Int: {
            auto ptr = (int*)(descriptor_->access(coreInst, coreMod, osdiId, flags));
            changed = *ptr != vwrite->val<const Int>();
            *ptr = vwrite->val<const Int>();
            break;
        }
        case Value::Type::Real: {
            auto ptr = (double*)(descriptor_->access(coreInst, coreMod, osdiId, flags));
            changed = *ptr != vwrite->val<const Real>();
            *ptr = vwrite->val<const Real>();
            break;
        }
        case Value::Type::String: {
            auto ptr = (char**)(descriptor_->access(coreInst, coreMod, osdiId, flags));
            auto& src = vwrite->val<const String>();
            // Check for change
            changed = *ptr && (src != *ptr);
            // Free old value and allocate new
            if (wasGiven) {
                free(*ptr);
            }
            *ptr = (char*)malloc(src.size()+1);
            // Copy data
            size_t i;
            for(i=0; i<src.size(); i++) {
                (*ptr)[i] = src.c_str()[i];
            }
            (*ptr)[i] = 0;
            break;
        }
    }
    
    return std::make_tuple(true, changed);
}

bool OsdiDevice::readParameter(OsdiFile::OsdiParameterId osdiId, void* coreMod, void* coreInst, Value& v, Status& s) const {
    // Check index
    if (osdiId>=osdiIdCount()) {
        s.set(Status::Range, std::string("OSDI parameter id=")+std::to_string(osdiId)+" out of range.");
        return false;
    }

    // Check kind
    if (!coreInst && !isModelParameter(osdiId)) {
        s.set(Status::NotFound, std::string("OSDI parameter id=")+std::to_string(osdiId)+" is not a model parameter.");
        return false;
    }

    if (coreInst && !isInstanceParameter(osdiId)) {
        s.set(Status::NotFound, std::string("OSDI parameter id=")+std::to_string(osdiId)+" is not an instance parameter.");
        return false;
    }

    // Get type
    auto t = parameterType(osdiId);
    if (t!=Value::Type::Int && t!=Value::Type::Real && t!=Value::Type::String) {
        s.set(Status::Unsupported, std::string("OSDI parameter id=")+std::to_string(osdiId)+" has unsupported type ("+Value::typeCodeToName(t)+").");
        return false;
    }

    // Read
    OsdiFile::OsdiFlags flags = ACCESS_FLAG_READ;
    if (coreInst) {
        flags |= ACCESS_FLAG_INSTANCE;
    }
    switch (t) {
        case Value::Type::Int: {
            auto ptr = (int*)(descriptor_->access(coreInst, coreMod, osdiId, flags));
            v = Int(*ptr);
            break;
        }
        case Value::Type::Real: {
            auto ptr = (double*)(descriptor_->access(coreInst, coreMod, osdiId, flags));
            v = Real(*ptr);
            break;
        }
        case Value::Type::String: {
            auto ptr = (char**)(descriptor_->access(coreInst, coreMod, osdiId, flags));
            // Should not be null, but anyway
            if (*ptr) {
                v = std::move(std::string(*ptr));
            } else {
                v = "";
            }
            break;
        }
    }
    
    return true;
}

std::tuple<bool, bool, bool> OsdiDevice::setup(Circuit& circuit, bool force, Status& s) {
    bool unknownsChanged = false;
    bool sparsityChanged = false;
    const auto& opt = circuit.simulatorOptions().core();
    const auto& internals = circuit.simulatorInternals();
    OsdiSimParas sp;
    populate(sp, opt, internals);
    for(auto model : models()) {
        // Verilog-A $temperature is in K, convert the value given by options (in C)
        auto [ok, tmpUnknowns, tmpSparsity] = static_cast<OsdiModel*>(model)->setupCore(circuit, sp, opt.temp+273.15, force, s);
        unknownsChanged |= tmpUnknowns;
        sparsityChanged |= tmpSparsity;
        if (!ok) {
            depopulate(sp);
            return std::make_tuple(false, unknownsChanged, sparsityChanged);
        }
        
    }
    depopulate(sp);
    return std::tuple(true, unknownsChanged, sparsityChanged);
}

bool OsdiDevice::collapseNodes(Circuit& circuit, Status& s) {
    for(auto model : models()) {
        for(auto instance : model->instances()) {
            if (!static_cast<OsdiInstance*>(instance)->collapseNodesCore(circuit, s)) {
                return false;
            }
            
        }
    }
    return true;
}

bool OsdiDevice::populateStructures(Circuit& circuit, Status& s) {
    for(auto model : models()) {
        for(auto instance : model->instances()) {
            // std::cout << std::string(name()) << " " << std::string(model->name()) << " " << std::string(instance->name()) << std::endl;
            if (!static_cast<OsdiInstance*>(instance)->populateStructuresCore(circuit, s)) {
                return false;
            }
        }
    }
    return true;
}

bool OsdiDevice::bind(
    Circuit& circuit, 
    KluRealMatrix* matResistReal, KluComplexMatrix* matResistCx, Component compResist, 
    KluRealMatrix* matReactReal, KluComplexMatrix* matReactCx, Component compReact, 
    Status& s
) {
    // Call bind() for all instances
    for(auto model : models()) {
        for(auto instance : model->instances()) {
            if (!static_cast<OsdiInstance*>(instance)->bindCore(
                circuit, 
                matResistReal, matResistCx, compResist, 
                matReactReal, matReactCx, compReact, 
                s
            )) {
                return false;
            }
        }
    }
    return true;
}

bool OsdiDevice::evalAndLoad(Circuit& circuit, EvalAndLoadSetup& els, Status& s) {
    auto& opt = circuit.simulatorOptions().core();
    auto& internals = circuit.simulatorInternals();
    OsdiSimInfo simInfo;

    populate(simInfo.paras, opt, internals);
    simInfo.abstime = internals.time;
    simInfo.prev_solve = els.oldSolution;

    simInfo.flags = 0;
    if (els.evaluateResistiveJacobian) {
        simInfo.flags |= CALC_RESIST_JACOBIAN;
    }
    if (els.evaluateReactiveJacobian) {
        simInfo.flags |= CALC_REACT_JACOBIAN;
    }
    if (els.evaluateResistiveResidual) {
        simInfo.flags |= CALC_RESIST_RESIDUAL; 
    }
    if (els.evaluateReactiveResidual) {
        simInfo.flags |= CALC_REACT_RESIDUAL; 
    }
    if (els.evaluateLinearizedResistiveRhsResidual) {
        simInfo.flags |= CALC_RESIST_LIM_RHS; 
    }
    if (els.evaluateLinearizedReactiveRhsResidual) {
        simInfo.flags |= CALC_REACT_LIM_RHS; 
    }
    if (els.evaluateNoise) {
        simInfo.flags |= CALC_NOISE; 
    }
    if (els.evaluateOpvars) {
        simInfo.flags |= CALC_OP; 
    }
    
    if (els.enableLimiting) {
        simInfo.flags |= ENABLE_LIM; 
    }
    if (els.initializeLimiting) {
        simInfo.flags |= INIT_LIM;
    }

    if (els.staticAnalysis) {
        simInfo.flags |= ANALYSIS_STATIC;
    }
    if (els.dcAnalysis) {
        simInfo.flags |= ANALYSIS_DC;
    }
    if (els.acAnalysis) {
        simInfo.flags |= ANALYSIS_AC;
    }
    if (els.tranAnalysis) {
        simInfo.flags |= ANALYSIS_TRAN;
    }
    if (els.noiseAnalysis) {
        simInfo.flags |= ANALYSIS_NOISE;
    }
    if (els.nodesetEnabled) {
        simInfo.flags |= ANALYSIS_NODESET;
    }
    if (els.icEnabled) {
        simInfo.flags |= ANALYSIS_IC;
    }

    for(auto model : models()) {
        for(auto instance : model->instances()) {
            if (!static_cast<OsdiInstance*>(instance)->evalAndLoadCore(circuit, simInfo, els, s)) {
                depopulate(simInfo.paras);
                return false;
            }
        }
    }
    depopulate(simInfo.paras);

    return true;
}

const char* OsdiDevice::simParamNames[] = {
    "gmin", // minimum conductance to place in parallel with nonlinear branches (simulator gmin)
    "gdev", // extra conductance in parallel with nonlinear branches during homotopy
    "tnom", 
    "minr", 
    "scale", 
    "iteration", 
    "simulatorVersion", 
    "simulatorSubversion", 
    "sourceScaleFactor",
    "initializeLimiting",  
    nullptr
    // TODO: According to VAMS LRM table 9-27 the following are missing
    //   imax
    //   imelt
    //   shrink
    //   timeUnit
    //   timePrecision
};

const char *OsdiDevice::simStrParamNames[] = { 
    "analysis_name",
    "analysis_type", 
    "cwd", 
    nullptr
    // TODO: According to VAMS LRM table 9-28 the following are missing
    //   module
    //   instance
    //   path
};

 void OsdiDevice::populate(OsdiSimParas& sp, const SimulatorOptions& opt, const SimulatorInternals& internals) {
    double* simParamValues = new double[sizeof(OsdiDevice::simParamNames)/sizeof(char*)];
    
    // Because most Verilog-A devices use only gmin, we set it to internals gmin+gdev
    // Those that implement this properly will actually use 2x gdev when gdev!=0
    simParamValues[0] = internals.gmin + internals.gdev; 
    simParamValues[1] = internals.gdev; 
    // $simparam(tnom) should return the tnom value given by options (in C)
    // No conversion needed. 
    simParamValues[2] = opt.tnom;
    simParamValues[3] = opt.minr;
    simParamValues[4] = opt.scale;
    simParamValues[5] = internals.iteration;
    simParamValues[6] = Simulator::majorVersion;
    simParamValues[7] = Simulator::minorVersion;
    simParamValues[8] = internals.sourcescalefactor;
    simParamValues[9] = internals.initalizeLimiting; 

    sp.names = const_cast<char**>(simParamNames); 
    sp.vals = simParamValues; 
    
    char** simStrParamValues = new char*[sizeof(OsdiDevice::simStrParamNames)/sizeof(char*)];
    
    simStrParamValues[0] = const_cast<char*>(internals.analysis_name.c_str());
    simStrParamValues[1] = const_cast<char*>(internals.analysis_type.c_str());
    simStrParamValues[2] = const_cast<char*>(Simulator::startupPath().c_str());

    sp.names_str = const_cast<char**>(simStrParamNames); 
    sp.vals_str = simStrParamValues;
}

void OsdiDevice::depopulate(OsdiSimParas& sp) {
    delete [] sp.vals;
    sp.vals = nullptr;
    
    delete [] sp.vals_str;
    sp.vals_str = nullptr;
}

bool OsdiDevice::processInitInfo(Circuit& circuit, OsdiInitInfo& initInfo, const char* typeString, Id name, Status& s) const {
    // If no flags that cause the simulation to abort are set and no error is recorded, we are done
    if (!(initInfo.flags & EVAL_RET_FLAG_FATAL) && initInfo.num_errors==0) {
        return true; 
    }
    
    // Prepare message prefix
    std::string pfx = std::string(typeString)+" '"+std::string(name)+"'";

    // Add error messages to status
    for (OsdiFile::OsdiErrorIndex i=0; i<initInfo.num_errors; i++) {
        OsdiInitError *err = &(initInfo.errors[i]);
        switch (err->code) {
            case INIT_ERR_OUT_OF_BOUNDS: {
                char *param = descriptor_->param_opvar[err->payload.parameter_id].name[0];
                circuit.setFlags(Circuit::Flags::Abort);
                s.extend(pfx+": parameter '"+std::string(param)+"' is out of bounds.");
                break;
        }
        default:
            circuit.setFlags(Circuit::Flags::Abort);
            s.extend(pfx+": unknown OSDI error code "+std::to_string(err->code)+".");
            break;
        }
    }

    // Must free list of errors, because we own it now
    // Use libc free() because it was allocated by libc malloc()
    if (initInfo.num_errors>0)
        free(initInfo.errors);

    // Handle flags
    if (initInfo.flags & EVAL_RET_FLAG_FINISH) {
        circuit.setFlags(Circuit::Flags::Finish);
    }
    if (initInfo.flags & EVAL_RET_FLAG_STOP) {
        circuit.setFlags(Circuit::Flags::Stop);
    }
    if (initInfo.flags & EVAL_RET_FLAG_FATAL) {
        circuit.setFlags(Circuit::Flags::Abort);
        s.extend(pfx+": Fatal error during setup. Aborting simulation.");
    }

    if (initInfo.num_errors>0 || circuit.checkFlags(Circuit::Flags::Abort)) {
        return false;
    } else {
        return true;
    }
}

void OsdiDevice::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "OSDI device " << std::string(name()) << " : " << file()->fileName() << " : " << index_ << "\n";
    if (descriptor_->num_nodes>0) {
        os << "  Nodes (terminals+internals=" << descriptor_->num_nodes << ", terminals=" << descriptor_->num_terminals << "):\n";
        for(OsdiFile::OsdiNodeIndex i=0; i<descriptor_->num_nodes; i++) {
            os << pfx << "    ";
            os << descriptor_->nodes[i].name;
            if (descriptor_->nodes[i].is_flow) {
                os << " (flow)";
            }
            os << ", units \"" << descriptor_->nodes[i].units << "\"";
            os << ", residual units \"" << descriptor_->nodes[i].residual_units << "\"\n";
        }
    }
    os << pfx << "  Number of internal states: " << descriptor_->num_states << "\n";
    if (modelParameterCount()>0) {
        os << pfx << "  Model parameters:\n";
        for(ParameterIndex i=0; i<modelParameterCount(); i++) {
            os << pfx << "    " << std::string(modelParameterName(i)) << " id=" << modelOsdiParameterId(i);
            auto& p = descriptor_->param_opvar[modelOsdiParameterId(i)];
            if (p.num_alias>0) {
                os << "\n" << pfx << "      Aliases: ";
                for(OsdiFile::OsdiAliasIndex j=0; j<p.num_alias; j++) {
                    os << p.name[1+j] << " ";
                }
            }
            os << "\n";
        }
    }
    if (instanceParameterCount()>0) {
        os << pfx << "  Instance parameters:\n";
        for(ParameterIndex i=0; i<instanceParameterCount(); i++) {
            os << pfx << "    " << std::string(instanceParameterName(i)) << " id=" << instanceOsdiParameterId(i);
            auto& p = descriptor_->param_opvar[instanceOsdiParameterId(i)];
            if (p.num_alias>0) {
                os << "\n" << pfx << "      Aliases: ";
                for(OsdiFile::OsdiAliasIndex j=0; j<p.num_alias; j++) {
                    os << p.name[1+j] << " ";
                }
            }
            os << "\n";
        }
    }
    if (opvarCount()>0) {
        os << pfx << "  Opvars:\n";
        for(ParameterIndex i=0; i<opvarCount(); i++) {
            os << pfx << "    " << std::string(opvarName(i)) << " id=" << opvarOsdiParameterId(i);
            auto& p = descriptor_->param_opvar[opvarOsdiParameterId(i)];
            if (p.num_alias>0) {
                os << "\n" << pfx << "      Aliases: ";
                for(OsdiFile::OsdiAliasIndex j=0; j<p.num_alias; j++) {
                    os << p.name[1+j] << " ";
                }
            }
            os << "\n";
        }
    }
    if (noiseSourceCount()>0) {
        os << pfx << "  Noise contributions:\n";
        for(ParameterIndex i=0; i<noiseSourceCount(); i++) {
            os << pfx << "    " << std::string(noiseSourceName(i)) << "\n";
        }
    }
}

}