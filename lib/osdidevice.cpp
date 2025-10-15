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
    : Device(asName ? asName : Id(of->deviceIdentifier(descriptorIndex)), location), osdiFile(of), index_(descriptorIndex) {
    descriptor_ = of->deviceDescriptor(descriptorIndex);
    setFlags(Flags::IsValid);
    if (osdiFile->allowsBypass(index_)) {
        setFlags(Flags::Bypassable);
    }}

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

bool OsdiDevice::freeValues(void* coreMod, void* coreInst) {
    OsdiFile::OsdiFlags flags = ACCESS_FLAG_SET;
    if (coreInst) {
        flags |= ACCESS_FLAG_INSTANCE;
    }
    
    if (coreInst) {
        // Free allocated parameters in instance structure
        for(auto osdiId : osdiFile->allocatedInstanceParameterIds(index_)) {
            if (auto [ok, given] = parameterGiven(osdiId, coreMod, coreInst); !ok || !given) {
                continue;
            }
            auto t = parameterType(osdiId);
            // Only need to free strings, vectors are preallocated in instance core structure
            switch (t) {
                case Value::Type::String: {
                    // TODO: handle string vectors
                    auto ptr = (char**)(descriptor_->access(coreInst, coreMod, osdiId, flags));
                    free(*ptr);
                    break;
                }
            }
        }
    } else {
        // Free allocated parameters in model structure
        for(auto osdiId : osdiFile->allocatedModelParamemeterIds(index_)) {
            if (auto [ok, given] = parameterGiven(osdiId, coreMod, coreInst); !ok || !given) {
                continue;
            }
            auto t = parameterType(osdiId);
            // Only need to free strings, vectors are preallocated in model core structure
            switch (t) {
                case Value::Type::String: {
                    // TODO: handle string vectors
                    auto ptr = (char**)(descriptor_->access(nullptr, coreMod, osdiId, flags));
                    free(*ptr);
                    break;
                }
            }
        }
    }
    return true;
}

std::tuple<bool,bool> OsdiDevice::writeParameter(OsdiFile::OsdiParameterId osdiId, void* coreMod, void* coreInst, const Value& v, Status& s) {
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

    if (coreInst && isOpvar(osdiId)) {
        s.set(Status::NotFound, std::string("OSDI parameter id=")+std::to_string(osdiId)+" is an opvar and cannot be written.");
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
            auto [ok, wasGiven] = parameterGiven(osdiId, coreMod, coreInst, s);
            if (!ok) {
                return std::make_tuple(false, false);
            }
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

std::tuple<bool, bool> OsdiDevice::parameterGiven(OsdiFile::OsdiParameterId osdiId, void* coreMod, void* coreInst, Status& s) const {
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

    if (coreInst && isOpvar(osdiId)) {
        s.set(Status::NotFound, std::string("OSDI parameter id=")+std::to_string(osdiId)+" is an opvar and cannot be given.");
        return std::make_tuple(false, false);
    }

    // Get given flag
    bool flag;
    if (coreInst) {
        // Instance
        flag = (bool)(descriptor_->given_flag_instance(coreInst, osdiId));
    } else {
        // Model
        flag = (bool)(descriptor_->given_flag_model(coreMod, osdiId));
    }

    return std::make_tuple(true, flag);
}

std::tuple<bool, bool, bool> OsdiDevice::setup(Circuit& circuit, CommonData& commons, bool force, DeviceRequests* devReq, Status& s) {
    bool unknownsChanged = false;
    bool sparsityChanged = false;
    const auto& opt = circuit.simulatorOptions().core();
    OsdiSimParas sp;
    
    // Allocate tables on stack
    auto [ndbl, nchrptr ] = simParasSizes();
    double dblArray[ndbl];
    char* chrPtrArray[nchrptr];
    
    populateSimParas(sp, opt, commons, dblArray, chrPtrArray);
    for(auto model : models()) {
        // Verilog-A $temperature is in K, convert the value given by options (in C)
        auto [ok, tmpUnknowns, tmpSparsity] = static_cast<OsdiModel*>(model)->setupCore(circuit, sp, opt.temp+273.15, force, devReq, s);
        unknownsChanged |= tmpUnknowns;
        sparsityChanged |= tmpSparsity;
        if (!ok) {
            return std::make_tuple(false, unknownsChanged, sparsityChanged);
        }
        
    }
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
            // std::cout << std::string(name()) << " " << std::string(model->name()) << " " << std::string(instance->name()) << "\n";
            if (!static_cast<OsdiInstance*>(instance)->populateStructuresCore(circuit, s)) {
                return false;
            }
        }
    }
    return true;
}

bool OsdiDevice::bind(
    Circuit& circuit, 
    KluMatrixAccess* matResist, Component compResist, const std::optional<MatrixEntryPosition>& mepResist, 
    KluMatrixAccess* matReact, Component compReact, const std::optional<MatrixEntryPosition>& mepReact, 
    Status& s
) {
    // Call bind() for all instances
    for(auto model : models()) {
        for(auto instance : model->instances()) {
            if (!static_cast<OsdiInstance*>(instance)->bindCore(
                circuit, 
                matResist, compResist, mepResist, 
                matReact, compReact, mepReact, 
                s
            )) {
                return false;
            }
        }
    }
    return true;
}

bool OsdiDevice::evalAndLoad(Circuit& circuit, CommonData& commons, EvalSetup* evalSetup, LoadSetup* loadSetup) {
    auto& opt = circuit.simulatorOptions().core();
    OsdiSimInfo simInfo;

    // Allocate tables on stack
    // Need to do this in the top level block because otherwise 
    // it gets freed as soon as we leave the block. 
    auto [ndbl, nchrptr ] = simParasSizes();
    double dblArray[ndbl];
    char* chrPtrArray[nchrptr];
    
    if (evalSetup) {
        populateSimParas(simInfo.paras, opt, commons, dblArray, chrPtrArray);
        simInfo.abstime = evalSetup->time;
        simInfo.prev_solve = evalSetup->oldSolution;

        simInfo.flags = 0;
        if (evalSetup->evaluateResistiveJacobian) {
            simInfo.flags |= CALC_RESIST_JACOBIAN;
        }
        if (evalSetup->evaluateReactiveJacobian) {
            simInfo.flags |= CALC_REACT_JACOBIAN;
        }
        if (evalSetup->evaluateResistiveResidual) {
            simInfo.flags |= CALC_RESIST_RESIDUAL; 
        }
        if (evalSetup->evaluateReactiveResidual) {
            simInfo.flags |= CALC_REACT_RESIDUAL; 
        }
        if (evalSetup->evaluateLinearizedResistiveRhsResidual) {
            simInfo.flags |= CALC_RESIST_LIM_RHS; 
        }
        if (evalSetup->evaluateLinearizedReactiveRhsResidual) {
            simInfo.flags |= CALC_REACT_LIM_RHS; 
        }
        if (evalSetup->evaluateNoise) {
            simInfo.flags |= CALC_NOISE; 
        }
        if (evalSetup->evaluateOpvars) {
            simInfo.flags |= CALC_OP; 
        }
        
        if (evalSetup->enableLimiting) {
            simInfo.flags |= ENABLE_LIM; 
        }
        if (evalSetup->initializeLimiting) {
            simInfo.flags |= INIT_LIM;
        }

        if (evalSetup->staticAnalysis) {
            simInfo.flags |= ANALYSIS_STATIC;
        }
        if (evalSetup->dcAnalysis) {
            simInfo.flags |= ANALYSIS_DC;
        }
        if (evalSetup->acAnalysis) {
            simInfo.flags |= ANALYSIS_AC;
        }
        if (evalSetup->tranAnalysis) {
            simInfo.flags |= ANALYSIS_TRAN;
        }
        if (evalSetup->noiseAnalysis) {
            simInfo.flags |= ANALYSIS_NOISE;
        }
        if (evalSetup->nodesetEnabled) {
            simInfo.flags |= ANALYSIS_NODESET;
        }
        if (evalSetup->icEnabled) {
            simInfo.flags |= ANALYSIS_IC;
        }
    }

    if constexpr(devacct) {
        novh++;
    }
    
    for(auto model : models()) {
        if (model->instanceCount()==0) {
            continue;
        }
        for(auto instance : model->instances()) {
            if (evalSetup && !static_cast<OsdiInstance*>(instance)->evalCore(circuit, commons, simInfo, *evalSetup)) {
                return false;
            }
            if (loadSetup) {
                auto lst = static_cast<OsdiInstance*>(instance)->loadCore(circuit, commons, *loadSetup);
                if (!lst) {
                    return false;
                }
            }
        }
    }
    
    return true;
}

const char* OsdiDevice::simParamNames[] = {
    "iniLim",  
    "gmin", // minimum conductance to place in parallel with nonlinear branches (simulator gmin)
    "gdev", // extra conductance in parallel with nonlinear branches during homotopy
    "tnom", 
    "minr", 
    "scale", 
    "iteration", 
    "simulatorVersion", 
    "simulatorSubversion", 
    "sourceScaleFactor",
    "reltol", 
    "vntol", 
    "abstol", 
    "chgtol", 
    "fluxtol", 
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

// Will allocate double and char* arrays on stack to make it faster
// We need sizes for that
std::tuple<size_t, size_t> OsdiDevice::simParasSizes() {
    return std::make_tuple(
        sizeof(OsdiDevice::simParamNames)/sizeof(char*), 
        sizeof(OsdiDevice::simStrParamNames)/sizeof(char*)
    );
}

void OsdiDevice::populateSimParas(OsdiSimParas& sp, const SimulatorOptions& opt, const CommonData& internals, double* dblArray, char** chrPtrArray) {
    // dblArray and chrPtrArray should be allocated on stack to save time
    // simParasSizes() reports the reuired size of these two arrays
    double* simParamValues = dblArray; 
    
    simParamValues[0] = 0; // iniLim
    // Because most Verilog-A devices use only gmin, we set it to internals gmin+gdev
    // Those that implement this properly will actually use 2x gdev when gdev!=0
    simParamValues[1] = internals.gmin + internals.gdev; 
    simParamValues[2] = internals.gdev; 
    // $simparam(tnom) should return the tnom value given by options (in C)
    // No conversion needed. 
    simParamValues[3] = opt.tnom;
    simParamValues[4] = opt.minr;
    simParamValues[5] = opt.scale;
    simParamValues[6] = internals.iteration;
    simParamValues[7] = Simulator::majorVersion;
    simParamValues[8] = Simulator::minorVersion;
    simParamValues[9] = internals.sourcescalefactor;
    
    simParamValues[10] = opt.reltol;
    simParamValues[11] = opt.vntol;
    simParamValues[12] = opt.abstol;
    simParamValues[13] = opt.chgtol;
    simParamValues[14] = opt.fluxtol;

    sp.names = const_cast<char**>(simParamNames); 
    sp.vals = simParamValues; 

    char** simStrParamValues = chrPtrArray; 
    
    simStrParamValues[0] = const_cast<char*>(internals.analysis_name.c_str());
    simStrParamValues[1] = const_cast<char*>(internals.analysis_type.c_str());
    simStrParamValues[2] = const_cast<char*>(Simulator::startupPath().c_str());

    sp.names_str = const_cast<char**>(simStrParamNames); 
    sp.vals_str = simStrParamValues;
}

void OsdiDevice::updateSimInfo(OsdiSimInfo& simInfo, EvalSetup& evalSetup) {
    simInfo.paras.vals[0] = evalSetup.initializeLimiting ? 1 : 0;
}

bool OsdiDevice::processInitInfo(Circuit& circuit, OsdiInitInfo& initInfo, const char* typeString, Id name, DeviceRequests* devReq, Status& s) const {
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
                s.extend(pfx+": parameter '"+std::string(param)+"' is out of bounds.");
                break;
            }
            default:
                s.extend(pfx+": unknown OSDI error code "+std::to_string(err->code)+".");
                break;
        }
    }
    if (devReq && initInfo.num_errors>0) {
        if (devReq) {
            devReq->abort = true;
        }
    }

    // Must free list of errors, because we own it now
    // Use libc free() because it was allocated by libc malloc()
    if (initInfo.num_errors>0)
        free(initInfo.errors);

    // Handle flags
    if (devReq) {
        if (initInfo.flags & EVAL_RET_FLAG_FINISH) {
            devReq->finish = true;
        }
        if (initInfo.flags & EVAL_RET_FLAG_STOP) {
            devReq->stop = true;
        }
        if (initInfo.flags & EVAL_RET_FLAG_FATAL) {
            devReq->abort = true;
        }
    }
    bool error = initInfo.num_errors>0 || (initInfo.flags & EVAL_RET_FLAG_FATAL);

    if (initInfo.flags & EVAL_RET_FLAG_FATAL) {
        s.extend(pfx+": Fatal error during setup. Aborting simulation.");
    }

    if (error) {
        return false;
    } else {
        return true;
    }
}

void OsdiDevice::dump(int indent, std::ostream& os) const {
    // Nature formatter
    auto natref = [this](OsdiNatureRef ref) -> const std::string {
        switch (ref.ref_type) {
            case NATREF_NATURE: 
                if (auto nat = file()->nature(ref.index)) {
                    return nat->name;
                }
                break;
            case NATREF_DISCIPLINE_FLOW: 
            case NATREF_DISCIPLINE_POTENTIAL: 
                if (auto disc = file()->discipline(ref.index)) {
                    auto dname = disc->name;
                    if (ref.ref_type==NATREF_DISCIPLINE_FLOW) {
                        return std::string(dname)+".flow"; 
                    } else {
                        return std::string(dname)+".potential";
                    }
                }
                break;
        }
        return "<none>";
    };
    
    auto tolToStr = [](double tol, int width) -> std::string {
        std::stringstream ss;
        ss << std::left << std::setw(width);
        ss.str(""); 
        if (std::isinf(tol)) {
            ss << "<none>";
        } else {
            ss << std::defaultfloat << tol;
        }
        return ss.str(); 
    };

    std::string pfx = std::string(indent, ' ');
    os << pfx << "OSDI device " << std::string(name()) << " : " << file()->fileName() << " : " << index_ << "\n";
    if (descriptor_->num_nodes>0) {
        os << "  Nodes (terminals+internals=" << descriptor_->num_nodes << ", terminals=" << descriptor_->num_terminals << "):\n";
        for(OsdiFile::OsdiNodeIndex i=0; i<descriptor_->num_nodes; i++) {
            os << pfx << "    ";
            os << nodeName(i);
            if (descriptor_->nodes[i].is_flow) {
                os << " (flow)";
            }
            auto un = descriptor_->unknown_nature[i];
            auto rn = descriptor_->residual_nature[i];
            os << ", " << natref(un) << " [" << descriptor_->nodes[i].units << "]";
            os << ", residual " << natref(rn) << " [" << descriptor_->nodes[i].residual_units << "]";
            if (descriptor_->nodes[i].resist_residual_off != UINT32_MAX) {
                os << ", resistive";
            }
            if (descriptor_->nodes[i].resist_limit_rhs_off != UINT32_MAX) {
                os << "+limiting";
            }
            if (descriptor_->nodes[i].react_residual_off != UINT32_MAX) {
                os << ", reactive";
            }
            if (descriptor_->nodes[i].react_limit_rhs_off != UINT32_MAX) {
                os << "+limiting";
            }
            
            os << "\n";
        }
        os << "  Absolute tolerances in Verilog-A mode:\n";
        for(OsdiFile::OsdiNodeIndex i=0; i<descriptor_->num_nodes; i++) {
            os << pfx << "    ";
            os << nodeName(i) << ":";
            auto [u, idtu, r, idtr] = file()->tolerances(index_, i);
            os << " abstol=" << tolToStr(u, 8) << " idt_abstol=" << tolToStr(idtu, 8);
            os << " res_abstol=" << tolToStr(r, 8) << " res_idt_abstol=" << tolToStr(idtr, 8);
            os << "\n";
        }
    }
    if (descriptor_->num_collapsible>0) {
        os << pfx << "  Collapsible node pairs:\n";
        for(ParameterIndex i=0; i<descriptor_->num_collapsible; i++) {
            auto n1 = descriptor_->collapsible[i].node_1;
            auto n2 = descriptor_->collapsible[i].node_2;
            os << pfx << "    " << descriptor_->nodes[n1].name;
            if (n2!=UINT32_MAX) {
                os << ", " << descriptor_->nodes[n2].name;
            } else {
                os << ", (ground)";
            }
            os << "\n";
        }
    }
    if (descriptor_->num_inputs>0) {
        os << pfx << "  Model inputs:\n";
        for(ParameterIndex i=0; i<descriptor_->num_inputs; i++) {
            auto n1 = descriptor_->inputs[i].node_1;
            auto n2 = descriptor_->inputs[i].node_2;
            os << pfx << "    ";
            if (n1!=UINT32_MAX) {
                os << nodeName(n1);
            } else {
                os << "(ground)";
            }
            if (n2!=UINT32_MAX) {
                os << " - " << nodeName(n2);
            } else {
                os << " - (ground)";
            }
            os << "\n";
        }
    }
    if (descriptor_->num_jacobian_entries>0) {
        os << pfx << "  Jacobian entries (";
        os << descriptor_->num_jacobian_entries << " nonzeros, ";
        os << descriptor_->num_resistive_jacobian_entries << " resistive, ";
        os << descriptor_->num_reactive_jacobian_entries << " reactive):\n";
        for(ParameterIndex i=0; i<descriptor_->num_jacobian_entries; i++) {
            auto& jac = descriptor_->jacobian_entries[i];
            auto n1 = jac.nodes.node_1;
            auto n2 = jac.nodes.node_2;
            os << pfx << "    (" << nodeName(n1) << ", " << nodeName(n2) << "): ";
            bool comma = false;
            if (jac.flags & JACOBIAN_ENTRY_RESIST) {
                comma = true;
                os << "resistive";
                if (jac.flags & JACOBIAN_ENTRY_RESIST_CONST) {
                    os << " constant";
                }
            } 
            if (jac.flags & JACOBIAN_ENTRY_REACT) {
                if (comma) {
                    os << ", ";
                }
                os << "reactive";
                if (jac.flags & JACOBIAN_ENTRY_REACT_CONST) {
                    os << " constant";
                }
            } 
            os << "\n";
        }
    }
    os << pfx << "  Number of internal states: " << descriptor_->num_states << "\n";
    if (modelParameterCount()>0) {
        os << pfx << "  Model parameters:\n";
        for(ParameterIndex i=0; i<modelParameterCount(); i++) {
            os << pfx << "    " << "id=" << modelOsdiParameterId(i) << ": " << std::string(modelParameterName(i));
            auto& p = descriptor_->param_opvar[modelOsdiParameterId(i)];
            os << ": \"" << p.description << "\"";
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
            os << pfx << "   " << " id=" << instanceOsdiParameterId(i) << ": " << std::string(instanceParameterName(i));
            auto& p = descriptor_->param_opvar[instanceOsdiParameterId(i)];
            os << ": \"" << p.description << "\"";
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
            os << pfx << "   " << " id=" << opvarOsdiParameterId(i) << ": " << std::string(opvarName(i));
            auto& p = descriptor_->param_opvar[opvarOsdiParameterId(i)];
            os << ": \"" << p.description << "\"";
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
            auto [n1, n2] = noiseExcitation(i);
            os << pfx << "    " << std::string(noiseSourceName(i)) << " : "
               << ((n1!=UINT32_MAX) ? nodeName(n1) : "(ground)") << " - "
               << ((n2!=UINT32_MAX) ? nodeName(n2) : "(ground)") << "\n";
        }
    }
    os << pfx << "  Model structure size   : " << descriptor_->model_size << " bytes\n";
    os << pfx << "  Instance structure size: " << descriptor_->instance_size << " bytes\n";
}

}