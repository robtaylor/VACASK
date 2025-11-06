#include <sstream>
#include "value.h"
#include "common.h"


namespace NAMESPACE {

std::string Value::typeCodeToName(Value::Type t) {
    switch (t) {
        case Type::Int: return "integer";
        case Type::Real: return "real";
        case Type::String: return "string";
        case Type::IntVec: return "integer vector";
        case Type::RealVec: return "real vector";
        case Type::StringVec: return "string vector";
        case Type::ValueVec: return "value vector";
    }
    return "";
}

template <typename Tl, typename Tr> bool valueEqual(const Value& l, const Value& r) {
    if constexpr(Value::IsVectorType<Tl>::value) {
        // Vector - vector
        auto& vl = l.val<const Tl>();
        auto& vr = r.val<const Tr>();
        auto nl = vl.size();
        auto nr = vr.size();
        if (nl!=nr) {
            return false;
        }
        for(decltype(nl) i=0; i<nl; i++) {
            if (vl[i]!=vr[i]) {
                return false;
            }
        }
    } else {
        // Scalar - scalar
        auto& vl = l.val<const Tl>();
        auto& vr = r.val<const Tr>();
        if (vl!=vr) {
            return false;
        }
    }
    return true;
}

bool operator==(const Value& l, const Value& r) {
    switch (type_pair(l.type_, r.type_)) {
        case type_pair(Value::Type::Int,    Value::Type::Int):    return valueEqual<Int, Int>(l, r);
        case type_pair(Value::Type::Int,    Value::Type::Real):   return valueEqual<Int, Real>(l, r);
        case type_pair(Value::Type::Real,   Value::Type::Int):    return valueEqual<Real, Int>(l, r);
        case type_pair(Value::Type::Real,   Value::Type::Real):   return valueEqual<Real, Real>(l, r);
        case type_pair(Value::Type::String, Value::Type::String): return valueEqual<String, String>(l, r);

        case type_pair(Value::Type::IntVec,    Value::Type::IntVec):    return valueEqual<IntVector, IntVector>(l, r);
        case type_pair(Value::Type::IntVec,    Value::Type::RealVec):   return valueEqual<IntVector, RealVector>(l, r);
        case type_pair(Value::Type::RealVec,   Value::Type::IntVec):    return valueEqual<RealVector, IntVector>(l, r);
        case type_pair(Value::Type::RealVec,   Value::Type::RealVec):   return valueEqual<RealVector, RealVector>(l, r);
        case type_pair(Value::Type::StringVec, Value::Type::StringVec): return valueEqual<StringVector, StringVector>(l, r);
        case type_pair(Value::Type::ValueVec,  Value::Type::ValueVec):  return valueEqual<ValueVector, ValueVector>(l, r);
    }
    return false;
}

bool operator!=(const Value& l, const Value& r) {
    return !(l==r);
}

template<typename Tfrom, typename Tto> bool convertValue(Value& src, Value& dest) {
    if constexpr(Value::IsVectorType<Tfrom>::value) {
        // vector -> vector
        auto& vfrom = src.val<Tfrom>();
        auto& vto = dest.val<Tto>();
        auto n = vfrom.size();
        vto.resize(n);
        for(decltype(n) i=0; i<n; i++) {
            vto[i] = Value::ScalarType<Tto>(vfrom[i]);
        }
    } else {
        // scalar -> scalar
        dest.val<Tto>() = Tto(src.val<Tfrom>());
    }
    return true;
}

bool Value::convert(Type to, Value& dest, Status &s) {
    if (type_==to) {
        dest = *this;
        return true;
    }

    // Prepare destination
    if (dest.type_ != to) {
        dest.~Value();
        switch (to) {
            case Type::String:
                dest.sVal = new String();
                break;
            case Type::IntVec:
                dest.iVec = new IntVector();
                break;
            case Type::RealVec:
                dest.rVec = new RealVector();
                break;
            case Type::StringVec:
                dest.sVec = new StringVector();
                break;
            case Type::ValueVec:
                dest.vVec = new ValueVector();
                break;
        }
        dest.type_ = to;
    }

    // Conversion of empty vector to any other vector is possible
    if (isVector() && size()==0 && dest.isVector()) {
        // Nothing to do, destination is already of correct type and length 0
        switch (to) {
            case Type::IntVec: 
            case Type::RealVec: 
            case Type::StringVec: 
            case Type::ValueVec: 
                return true;
            default: 
                return false;
        }
    }

    switch (type_pair(type_, to)) {
        case type_pair(Type::Int,  Type::Real): return convertValue<Int, Real>(*this, dest);
        case type_pair(Type::Real, Type::Int): return convertValue<Real, Int>(*this, dest);
        case type_pair(Type::IntVec,  Type::RealVec): return convertValue<IntVector, RealVector>(*this, dest);
        case type_pair(Type::RealVec, Type::IntVec): return convertValue<RealVector, IntVector>(*this, dest);
    }

    s.set(
        Status::BadConversion, 
        std::string("Cannot convert ")+Value::typeCodeToName(type_)+" into "+Value::typeCodeToName(to)+"."
    );
    return false;
}

template<typename Tfrom, typename Tto> bool convertValueInPlace(Value& v) {
    if constexpr(Value::IsVectorType<Tfrom>::value) {
        // vector -> vector
        Tfrom& old = v.val<Tfrom>();
        auto n = old.size();
        Tto newVec(n);
        for(decltype(n) i=0; i<n; i++) {
            newVec[i] = Value::ScalarType<Tto>(old[i]);
        }
        v = std::move(newVec);
    } else {
        // scalar -> scalar
        v = Tto(v.val<Tfrom>());
    }
    return true;
}

bool Value::convertInPlace(Type to, Status &s) {
    if (type_==to)
        return true;

    // Conversion of empty vector to any other vector is possible
    if (isVector() && size()==0 && (to & ValueType::VectorBit)!=0) {
        switch (to) {
            case Type::IntVec: *this = std::move(IntVector(0)); return true;
            case Type::RealVec: *this = std::move(RealVector(0)); return true;
            case Type::StringVec: *this = std::move(StringVector(0)); return true;
            case Type::ValueVec: *this = std::move(ValueVector(0)); return true;
            default: return false;
        }
    }

    switch (type_pair(type_, to)) {
        case type_pair(Type::Int,     Type::Real): return convertValueInPlace<Int, Real>(*this);
        case type_pair(Type::Real,    Type::Int): return convertValueInPlace<Real, Int>(*this);
        case type_pair(Type::IntVec,  Type::RealVec): return convertValueInPlace<IntVector, RealVector>(*this);
        case type_pair(Type::RealVec, Type::IntVec): return convertValueInPlace<RealVector, IntVector>(*this);
    }

    s.set(
        Status::BadConversion, 
        std::string("Cannot convert ")+Value::typeCodeToName(type_)+" into "+Value::typeCodeToName(to)+"."
    );
    return false;
}

bool Value::getScalar(Value& v, size_t ndx, Status& s) const {
    if (!isVector()) {
        s.set(Status::BadConversion, "Cannot extract components from a scalar.");
        return false;
    }
    if (ndx<0 || ndx>=size()) {
        s.set(Status::Range, "Vector index out of range.");
        return false;
    }
    switch (type_) {
        case Type::IntVec:
            v = (*iVec)[ndx];
            break;
        case Type::RealVec:
            v = (*rVec)[ndx];
            break;
        case Type::StringVec:
            v = (*sVec)[ndx];
            break;
        case Type::ValueVec:
            v = (*vVec)[ndx];
            break;
    }
    return true;
}

std::string Value::str() const {
    std::stringstream s;
    s << *this;
    return s.str();
}

std::ostream& operator<<(std::ostream& os, const Value& obj) {
    switch(obj.type_) {
        case Value::Type::Int: os << obj.iVal; break;
        case Value::Type::Real: os << obj.rVal; break;
        case Value::Type::String: os << "\"" << *(obj.sVal) << "\""; break;
        case Value::Type::IntVec: 
            os << "[";
            for(auto it=obj.iVec->cbegin(); it!=obj.iVec->cend(); ++it) { 
                if (it!=obj.iVec->cbegin()) {
                    os << ", ";
                }
                os << *it;
            }
            os << "]";
            break;
        case Value::Type::RealVec: 
            os << "[";
            for(auto it=obj.rVec->cbegin(); it!=obj.rVec->cend(); ++it) { 
                if (it!=obj.rVec->cbegin()) {
                    os << ", ";
                }
                os << *it;
            }
            os << "]";
            break;
        case Value::Type::StringVec: 
            os << "[";
            for(auto it=obj.sVec->cbegin(); it!=obj.sVec->cend(); ++it) { 
                if (it!=obj.sVec->cbegin()) {
                    os << ", ";
                }
                os << "\"" << *it << "\""; 
            }
            os << "]";
            break;
        case Value::Type::ValueVec: 
            os << "[";
            for(auto it=obj.vVec->cbegin(); it!=obj.vVec->cend(); ++it) { 
                if (it!=obj.vVec->cbegin()) {
                    os << ", ";
                }
                os << "\"" << *it << "\""; 
            }
            os << "]";
            break;
    }
    return os;
}

}
