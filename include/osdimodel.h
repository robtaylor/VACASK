#ifndef __OSDIMODEL_DEFINED
#define __OSDIMODEL_DEFINED

#include <string>
#include <stdexcept>
#include "status.h"
#include "osdidevice.h"
#include "value.h"
#include "parseroutput.h"
#include "devbase.h"
#include "common.h"


namespace NAMESPACE {

class Circuit;

class OsdiInstance;

class OsdiModel : public Model {
public:
    friend class OsdiDevice;
    friend class OsdiInstance;

    // The model is added to the list of models of the given device so that when
    // a device is destroyed all the models in the list are destroyed, too. 
    OsdiModel(OsdiDevice* device, Id name, Instance* parentInstance, const PTModel& parsedModel, Status& s=Status::ignore);
    virtual ~OsdiModel();

    OsdiModel           (const OsdiModel&)  = delete;
    OsdiModel           (      OsdiModel&&) = default;
    OsdiModel& operator=(const OsdiModel&)  = delete;
    OsdiModel& operator=(      OsdiModel&&) = default;

    virtual ParameterIndex parameterCount() const;
    virtual std::tuple<ParameterIndex, bool> parameterIndex(Id name) const;
    virtual Id parameterName(ParameterIndex ndx) const;
    virtual std::tuple<Value::Type,bool> parameterType(ParameterIndex ndx, Status& s=Status::ignore) const;
    virtual bool getParameter(ParameterIndex ndx, Value& v, Status& s=Status::ignore) const;
    virtual std::tuple<bool, bool> setParameter(ParameterIndex ndx, const Value& v, Status& s=Status::ignore);
    virtual std::tuple<bool,bool> parameterGiven(ParameterIndex ndx, Status& s=Status::ignore) const;

    virtual std::tuple<bool, bool, bool> setup(Circuit& circuit, CommonData& commons, bool force, DeviceRequests* devReq, Status& s=Status::ignore);
    virtual Instance* createInstance(Circuit& circuit, Instance* parentInstance, RpnEvaluator& evaluator, Context* externalContext, const PTInstance& parsedInstance, InstantiationData& idata, Status& s=Status::ignore);
    virtual void dump(int indent, std::ostream& os) const;

    // Device access (as OsdiDevice)
    OsdiDevice* device() { return static_cast<OsdiDevice*>(device_); };
    const OsdiDevice* device() const { return static_cast<const OsdiDevice*>(device_); };

    // Model core access
    void* core() { return core_; };
    const void* core() const { return core_; };
    
    // Helpers for inlining in model and instance virtual functions
    std::tuple<bool, bool, bool> setupWrapper(Circuit& circuit, OsdiSimParas& sp, double temp,  bool force, DeviceRequests* devReq, Status& s=Status::ignore);
    
private:
    void* core_;
};

}

#endif
