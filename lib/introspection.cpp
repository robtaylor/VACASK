#include "introspection.h"
#include "common.h"

namespace NAMESPACE {

Value::Type StructMember::valueType(Type t) {
    switch (t) {
        case Type::Int: return Value::Type::Int;
        case Type::Real: return Value::Type::Real;
        case Type::String: return Value::Type::String;
        case Type::Id: return Value::Type::String;
        case Type::IntVec: return Value::Type::IntVec;
        case Type::RealVec: return Value::Type::RealVec;
        case Type::StringVec: return Value::Type::StringVec;
        case Type::IdVec: return Value::Type::StringVec;
        default: return Value::Type::Value;
    }
}

}
