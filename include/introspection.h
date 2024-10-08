#ifndef __INTROSPECTION_DEFINED
#define __INTROSPECTION_DEFINED

#include "value.h"
#include "identifier.h"
#include "status.h"
#include <string>
#include <unordered_map>
#include "common.h"


namespace NAMESPACE {

// Ids are exposed as string values
typedef struct StructMember {
    enum class Type {
        Int, Real, String, Id, IntVec, RealVec, StringVec, IdVec, Value, ValueVec
    };
    Id name;
    Type type;
    size_t offset;
    
    StructMember(Id n, Type t, size_t o) : name(n), type(t), offset(o) {};

    template<typename T> static Type typeCode() { return Type::Int; };
    static Value::Type valueType(Type t);
} StructMember;

template<> inline StructMember::Type StructMember::typeCode<Int>() { return StructMember::Type::Int; }; 
template<> inline StructMember::Type StructMember::typeCode<Real>() { return StructMember::Type::Real; }; 
template<> inline StructMember::Type StructMember::typeCode<String>() { return StructMember::Type::String; }; 
template<> inline StructMember::Type StructMember::typeCode<Id>() { return StructMember::Type::Id; }; 
template<> inline StructMember::Type StructMember::typeCode<IntVector>() { return StructMember::Type::IntVec; }; 
template<> inline StructMember::Type StructMember::typeCode<RealVector>() { return StructMember::Type::RealVec; }; 
template<> inline StructMember::Type StructMember::typeCode<StringVector>() { return StructMember::Type::StringVec; }; 
template<> inline StructMember::Type StructMember::typeCode<std::vector<Id>>() { return StructMember::Type::IdVec; }; 
template<> inline StructMember::Type StructMember::typeCode<Value>() { return StructMember::Type::Value; }; 
template<> inline StructMember::Type StructMember::typeCode<ValueVector>() { return StructMember::Type::ValueVec; }; 


template<typename T> class Introspection {
public:
    using IntrocpectedStruct = T;

    static size_t count();
    static std::tuple<size_t, bool> index(Id name);
    static Id name(size_t ndx);
    static std::tuple<Value::Type,bool> type(size_t ndx);

    static std::tuple<bool,bool> set(T& st, size_t ndx, const Value& v, Status& s=Status::ignore);
    static std::tuple<bool,bool> set(T& st, Id name, const Value& v, Status& s=Status::ignore);
    static bool get(const T& st, size_t ndx, Value& v, Status& s=Status::ignore);
    static bool get(const T& st, Id name, Value& v, Status& s=Status::ignore);
    static std::tuple<bool,bool> given(const T& st, size_t ndx, Status& s=Status::ignore);
    static std::tuple<bool,bool> given(const T& st, Id name, Status& s=Status::ignore);
    
    static void registerMember_(Id name, StructMember::Type t, size_t offset);
    static int setup();

    template<typename Tsub> static Tsub& memberRef(T& st, size_t ndx);
    template<typename Tsub> static const Tsub& memberRef(const T& st, size_t ndx);

    static bool memberEqual(const T& st1, const T& st2, size_t ndx);

private:
    static std::vector<StructMember> members;
    static std::unordered_map<Id,size_t> memberIndex;
    static int isRegistered;
};

template<typename T> std::vector<StructMember> Introspection<T>::members;
template<typename T> std::unordered_map<Id,size_t> Introspection<T>::memberIndex;

// Macro for registering a struct member
#define registerMember(m) \
    registerMember_(Id::createStatic(#m), StructMember::typeCode<decltype(IntrocpectedStruct::m)>(), offsetof(IntrocpectedStruct, m))
#define registerNamedMember(m, name) \
    registerMember_(Id::createStatic(name), StructMember::typeCode<decltype(IntrocpectedStruct::m)>(), offsetof(IntrocpectedStruct, m))

// Macro for instantiating an introspection (call after setup() specialization)
#define instantiateIntrospection(S) \
    template class Introspection<S>; \
    template<> int Introspection<S>::isRegistered = Introspection<S>::setup(); 


// Template implementation

template<typename T> size_t Introspection<T>::count() { 
    return members.size(); 
}

template<typename T> std::tuple<size_t, bool> Introspection<T>::index(Id name) {
    auto it = memberIndex.find(name);
    if (it==memberIndex.end()) {
        return std::make_tuple(0, false);
    }
    return std::make_tuple(it->second, true);
}

template<typename T> Id Introspection<T>::name(size_t ndx) { 
    if (ndx<count()) 
        return members[ndx].name; 
    else 
        return Id::none; 
}

template<typename T> std::tuple<Value::Type,bool> Introspection<T>::type(size_t ndx) {
    if (ndx<count()) 
        return std::make_tuple(StructMember::valueType(members[ndx].type), true); 
    else 
        return std::make_tuple(Value::Type::Int, false); 
}


template<typename T> std::tuple<bool,bool> Introspection<T>::set(T& st, size_t ndx, const Value& v, Status& s) {
    bool changed = false;

    // Check range
    if (ndx>=count()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return std::make_tuple(false, changed);
    }

    // Get type and offset
    auto t = members[ndx].type;
    auto offs = members[ndx].offset;

    // Compute target Value type
    auto vt = StructMember::valueType(t);
    
    // Value type does not need conversion, just assign it
    if (t==StructMember::Type::Value) { 
        auto& ref = Introspection<T>::memberRef<Value>(st, ndx);
        auto changed = ref != v;
        ref = v;
        return std::make_tuple(true, changed);
    } else if (t==StructMember::Type::ValueVec && v.type()==Value::Type::ValueVec) { 
        auto& ref = Introspection<T>::memberRef<ValueVector>(st, ndx);
        auto changed = ref != v;
        ref = v.val<const ValueVector>();
        return std::make_tuple(true, changed);
    }

    // Convert
    const Value* vwrite = &v;
    Value vconv;
    if (v.type()!=vt) {
        vconv = v;
        if (!vconv.convertInPlace(vt, s)) {
            s.extend("Value conversion failed for parameter id='"+std::to_string(ndx)+"'.");
            return std::make_tuple(false, changed);
        }
        vwrite = &vconv;
    }

    // Set
    switch (t) {
        case StructMember::Type::Int: {
            auto& ref = Introspection<T>::memberRef<Int>(st, ndx);
            changed = ref != vwrite->val<const Int>();
            ref = vwrite->val<const Int>();
            break;
        }
        case StructMember::Type::Real: {
            auto& ref = Introspection<T>::memberRef<Real>(st, ndx);
            changed = ref != vwrite->val<const Real>();
            ref = vwrite->val<const Real>();
            break;
        }
        case StructMember::Type::String: {
            auto& ref = Introspection<T>::memberRef<String>(st, ndx);
            changed = ref != vwrite->val<const String>();
            ref = vwrite->val<const String>();
            break;
        }
        case StructMember::Type::Id: {
            auto& ref = Introspection<T>::memberRef<Id>(st, ndx);
            changed = ref != vwrite->val<const String>();
            ref = vwrite->val<const String>();
            break;
        }
        case StructMember::Type::IntVec: {
            auto& ref = Introspection<T>::memberRef<IntVector>(st, ndx);
            changed = ref != vwrite->val<const IntVector>();
            ref = vwrite->val<const IntVector>();
            break;
        }
        case StructMember::Type::RealVec: {
            auto& ref = Introspection<T>::memberRef<RealVector>(st, ndx);
            changed = ref != vwrite->val<const RealVector>();
            ref = vwrite->val<const RealVector>();
            break;
        }
        case StructMember::Type::StringVec: {
            auto& ref = Introspection<T>::memberRef<StringVector>(st, ndx);
            changed = ref != vwrite->val<const StringVector>();
            ref = vwrite->val<const StringVector>();
            break;
        }
        case StructMember::Type::IdVec: {
            auto& ref = Introspection<T>::memberRef<std::vector<Id>>(st, ndx);
            auto& srcRef = vwrite->val<const StringVector>();
            auto refLen = ref.size();
            auto srcLen = srcRef.size();
            if (refLen!=srcLen) {
                ref.resize(srcLen);
                changed = true;
            }
            for(decltype(srcLen) i=0; i<srcLen; i++) {
                if (!changed) {
                    changed |= ref[i] != srcRef[i];
                }
                ref[i] = srcRef[i];
            }
            break;
        }
    }
    
    return std::make_tuple(true, changed);
}

template<typename T> std::tuple<bool,bool> Introspection<T>::set(T& st, Id name, const Value& v, Status& s) {
    auto [ndx, found] = Introspection<T>::index(name);
    // Find
    if (!found) {
        s.set(Status::NotFound, std::string("Parameter '")+std::string(name)+"' not found.");
        return std::make_tuple(false, false);
    }

    return Introspection::set(st, ndx, v, s);
}

template<typename T> bool Introspection<T>::get(const T& st, size_t ndx, Value& v, Status& s) {
    // Check range
    if (ndx>=count()) {
        s.set(Status::Range, std::string("Parameter index id=")+std::to_string(ndx)+" out of range.");
        return false;
    }

    // Get type and offset
    auto t = members[ndx].type;
    auto offs = members[ndx].offset;

    // Compute source Value type
    auto vt = StructMember::valueType(t);
    
    // Get
    switch (t) {
        case StructMember::Type::Int:
            v = Introspection<T>::memberRef<Int>(st, ndx);
            break;
        case StructMember::Type::Real:
            v = Introspection<T>::memberRef<Real>(st, ndx);
            break;
        case StructMember::Type::String:
            v = Introspection<T>::memberRef<String>(st, ndx);
            break;
        case StructMember::Type::Id:
            v = std::string(Introspection<T>::memberRef<Id>(st, ndx));
            break;
        case StructMember::Type::IntVec:
            v = Introspection<T>::memberRef<IntVector>(st, ndx);
            break;
        case StructMember::Type::RealVec:
            v = Introspection<T>::memberRef<RealVector>(st, ndx);
            break;
        case StructMember::Type::StringVec:
            v = Introspection<T>::memberRef<StringVector>(st, ndx);
            break;
        case StructMember::Type::IdVec:
            { 
                auto& srcRef = Introspection<T>::memberRef<std::vector<Id>>(st, ndx);;
                std::vector<std::string> sv;
                for(auto id : srcRef) {
                    sv.push_back(id);
                }
                v = std::move(sv);
            }
            break;
        case StructMember::Type::Value:
            v = Introspection<T>::memberRef<Value>(st, ndx);
            break;
        case StructMember::Type::ValueVec:
            v = Introspection<T>::memberRef<ValueVector>(st, ndx);
            break;
    }

    return true;
}

template<typename T> bool Introspection<T>::get(const T& st, Id name, Value& v, Status& s) {
    auto [ndx, found] = Introspection<T>::index(name);
    // Find
    if (!found) {
        s.set(Status::NotFound, std::string("Parameter '")+std::string(name)+"' not found.");
        return false;
    }

    return Introspection::get(st, ndx, v, s);
}

template<typename T> std::tuple<bool,bool> Introspection<T>::given(const T& st, size_t ndx, Status& s) {
    // Assume all struct values are always given
    return std::make_tuple(true, true);
}

template<typename T> std::tuple<bool,bool> Introspection<T>::given(const T& st, Id name, Status& s) {
    auto [ndx, found] = Introspection<T>::index(name);
    // Find
    if (!found) {
        s.set(Status::NotFound, std::string("Parameter '")+std::string(name)+"' not found.");
        return std::make_tuple(false, false);;
    }

    return Introspection::given(st, ndx, s);
}

template<typename T> void Introspection<T>::registerMember_(Id name, StructMember::Type t, size_t offset) {
    members.push_back(std::move(StructMember(name, t, offset)));
    memberIndex.insert({name, members.size()-1});
}

template<typename T> template<typename Tsub> Tsub& Introspection<T>::memberRef(T& st, size_t ndx) {
    auto offs = members[ndx].offset;
    return *reinterpret_cast<Tsub*>((char*)(&st)+offs);
}
    
template<typename T> template<typename Tsub> const Tsub& Introspection<T>::memberRef(const T& st, size_t ndx) {
    auto offs = members[ndx].offset;
    return *reinterpret_cast<const Tsub*>((char*)(&st)+offs);
}

template<typename T> bool Introspection<T>::memberEqual(const T& st1, const T& st2, size_t ndx) {
    auto t = members[ndx].type;
    switch (t) {
        case StructMember::Type::Int:
            return Introspection<T>::memberRef<Int>(st1, ndx) == Introspection<T>::memberRef<Int>(st2, ndx);
        case StructMember::Type::Real:
            return Introspection<T>::memberRef<Real>(st1, ndx) == Introspection<T>::memberRef<Real>(st2, ndx);
        case StructMember::Type::String:
            return Introspection<T>::memberRef<String>(st1, ndx) == Introspection<T>::memberRef<String>(st2, ndx);
        case StructMember::Type::Id:
            return Introspection<T>::memberRef<Id>(st1, ndx) == Introspection<T>::memberRef<Id>(st2, ndx);
        case StructMember::Type::IntVec:
            return Introspection<T>::memberRef<IntVector>(st1, ndx) == Introspection<T>::memberRef<IntVector>(st2, ndx);
        case StructMember::Type::RealVec:
            return Introspection<T>::memberRef<RealVector>(st1, ndx) == Introspection<T>::memberRef<RealVector>(st2, ndx);
        case StructMember::Type::StringVec:
            return Introspection<T>::memberRef<StringVector>(st1, ndx) == Introspection<T>::memberRef<StringVector>(st2, ndx);
        case StructMember::Type::IdVec:
            return Introspection<T>::memberRef<std::vector<Id>>(st1, ndx) == Introspection<T>::memberRef<std::vector<Id>>(st2, ndx);
        case StructMember::Type::Value:
            return Introspection<T>::memberRef<Value>(st1, ndx) == Introspection<T>::memberRef<Value>(st2, ndx);
        case StructMember::Type::ValueVec:
            return Introspection<T>::memberRef<ValueVector>(st1, ndx) == Introspection<T>::memberRef<ValueVector>(st2, ndx);
    }

    return false;
}

}

#endif
