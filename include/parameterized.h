#ifndef __PARAMETERIZED_DEFINED
#define __PARAMETERIZED_DEFINED

#include <stdint.h>
#include <tuple>
#include "rpneval.h"
#include "parseroutput.h"
#include "identifier.h"
#include "value.h"
#include "status.h"
#include "introspection.h"
#include "common.h"


namespace NAMESPACE {

// Parameterized object with parameter interface
class Parameterized {
public:
    // typedef ParameterIndex Index;

    Parameterized();

    // No data members, no need to restrict copy/move

    virtual ParameterIndex parameterCount() const = 0;
    virtual std::tuple<ParameterIndex, bool> parameterIndex(Id name) const = 0;
    virtual Id parameterName(ParameterIndex ndx) const = 0;

    virtual std::tuple<Value::Type,bool> parameterType(ParameterIndex ndx, Status& s=Status::ignore) const = 0;
    virtual std::tuple<Value::Type,bool> parameterType(Id name, Status& s=Status::ignore) const;

    virtual bool getParameter(ParameterIndex ndx, Value& v, Status& s=Status::ignore) const = 0;
    virtual std::tuple<bool,bool> setParameter(ParameterIndex ndx, const Value& v, Status& s=Status::ignore) = 0;
    
    virtual bool getParameter(Id name, Value& v, Status& s=Status::ignore) const;
    virtual std::tuple<bool,bool> setParameter(Id name, const Value& v, Status& s=Status::ignore); 

    virtual std::tuple<bool,bool> parameterGiven(ParameterIndex ndx, Status& s=Status::ignore);
    virtual std::tuple<bool,bool> parameterGiven(Id name, Status& s=Status::ignore);
    
    // Set parameters from parsed netlist
    std::tuple<bool,bool> setParameters(const std::vector<PTParameterValue>& params, Status& s=Status::ignore);
    std::tuple<bool,bool> setParameters(const std::vector<PTParameterExpression>& params, RpnEvaluator& eval, Status& s=Status::ignore);
    std::tuple<bool,bool> setParameters(const PTParameters& params, RpnEvaluator& eval, Status& s=Status::ignore);
    std::tuple<bool,bool> setParameters(const PTParameterMap& params, RpnEvaluator& eval, Status& s=Status::ignore);

    void dump(std::ostream& os, const char* prefix="") const;
};


// Parameterized introspectible struct
// Introspection must be set up when core struct (T) is defined
template<typename CoreT> class IStruct : public Parameterized {
public:
    typedef CoreT Core;

    IStruct();
    IStruct(CoreT& other);

    IStruct           (const IStruct&)  = default;
    IStruct           (      IStruct&&) = default;
    IStruct& operator=(const IStruct&)  = default;
    IStruct& operator=(      IStruct&&) = default;


    virtual ParameterIndex parameterCount() const;
    virtual std::tuple<ParameterIndex, bool> parameterIndex(Id name) const;
    virtual Id parameterName(ParameterIndex ndx) const;

    virtual std::tuple<Value::Type,bool> parameterType(ParameterIndex ndx, Status& s=Status::ignore) const;

    virtual bool getParameter(ParameterIndex ndx, Value& v, Status& s=Status::ignore) const;
    virtual std::tuple<bool,bool> setParameter(ParameterIndex ndx, const Value& v, Status& s=Status::ignore);
    virtual std::tuple<bool,bool> parameterGiven(ParameterIndex ndx, Status& s=Status::ignore);

    virtual bool getParameter(Id name, Value& v, Status& s=Status::ignore) const;
    virtual std::tuple<bool,bool> setParameter(Id name, const Value& v, Status& s=Status::ignore); 
    virtual std::tuple<bool,bool> parameterGiven(Id name, Status& s=Status::ignore);
        
    const CoreT& core() const { return core_; };
    CoreT& core() { return core_; };

private:
    CoreT core_;
};

// Tempate implementation
template<typename CoreT> IStruct<CoreT>::IStruct() {
}

template<typename CoreT> IStruct<CoreT>::IStruct(CoreT& other) 
    : core_(other) {
}

template<typename CoreT> ParameterIndex IStruct<CoreT>::parameterCount() const {
    ParameterIndex count = Introspection<CoreT>::count(); 
    return count; 
}

template<typename CoreT> std::tuple<ParameterIndex, bool> IStruct<CoreT>::parameterIndex(Id name) const {
    auto [ndx, ok] = Introspection<CoreT>::index(name); 
    return std::make_tuple(static_cast<ParameterIndex>(ndx), ok);
}

template<typename CoreT> Id IStruct<CoreT>::parameterName(ParameterIndex ndx) const {
    return Introspection<CoreT>::name(ndx);
}

template<typename CoreT> std::tuple<Value::Type,bool> IStruct<CoreT>::parameterType(ParameterIndex ndx, Status& s) const {
    return Introspection<CoreT>::type(ndx);
}

template<typename CoreT> bool IStruct<CoreT>::getParameter(ParameterIndex ndx, Value& v, Status& s) const {
    return Introspection<CoreT>::get(core_, ndx, v, s);
}

template<typename CoreT> std::tuple<bool,bool> IStruct<CoreT>::setParameter(ParameterIndex ndx, const Value& v, Status& s) {
    return Introspection<CoreT>::set(core_, ndx, v, s);
}

template<typename CoreT> std::tuple<bool,bool> IStruct<CoreT>::parameterGiven(ParameterIndex ndx, Status& s) {
    return Introspection<CoreT>::given(core_, ndx, s);
}

template<typename CoreT> bool IStruct<CoreT>::getParameter(Id name, Value& v, Status& s) const {
    return Introspection<CoreT>::get(core_, name, v, s);
}

template<typename CoreT> std::tuple<bool,bool> IStruct<CoreT>::setParameter(Id name, const Value& v, Status& s) {
    return Introspection<CoreT>::set(core_, name, v, s);
}

template<typename CoreT> std::tuple<bool,bool> IStruct<CoreT>::parameterGiven(Id name, Status& s) {
    return Introspection<CoreT>::given(core_, name, s);
}


}

#endif
