#ifndef __RPNBUILTIN_DEFINED
#define __RPNBUILTIN_DEFINED

#include "value.h"
#include "status.h"
#include "rpnexpr.h"
#include "rpnstack.h"
#include "rpnfunctor.h"
#include <cmath>
#include <type_traits>
#include "common.h"


namespace NAMESPACE {

// General RPN builtin function/operator
typedef bool (*RpnBuiltinFunc)(RpnStack& stack, Rpn::Arity argc, Status& s);

// Unary operator / single argument function
// Applied in the same way to all vector components
template<typename Tin, typename F> bool rpnMathFunc1(Value* arg, Status& s=Status::ignore) {
    // Functor
    F functor;
    
    // Scalar agument type
    using Tsin = Value::ScalarType<Tin>;
    
    // Scalar result type
    using Tsout = decltype(functor(Tsin()));

    // Actual result type (vectorize if neccessary)
    using Tout = Value::VectorType<Tsout,Value::IsVectorType<Tin>::value>;
    
    // Check if we need temporary storage, prepare it
    Value *dest;
    Value tmpVal;
    if (Value::typeCode<Tout>()==arg->type()) {
        // Reuse argv as destination
        dest = arg;
    } else {
        // Temporary destination
        dest = &tmpVal;
        if (arg->isVector()) {
            // Initialize vector, set value type
            tmpVal = std::move(Tout(arg->size()));
        }
        // For scalar result the type is set at result assignment
    }

    // Compute    
    if constexpr(Value::IsVectorType<Tin>::value) {
        // Vector -> vector
        for(size_t i=0; i<arg->size(); i++) {
            auto& v = arg->val<Tin>()[i];
            auto& dv = dest->val<Tin>()[i];
            if (!functor.ok(v, s)) {
                return false;
            }
            dv = functor(v);
        }
    } else {
        // Scalar -> scalar
        auto& v = arg->val<Tin>();
        if (!functor.ok(v, s)) {
            return false;
        }
        *dest = functor(v);
    }

    // Swap with arg
    swap(*arg, *dest);

    return true;
}

// Binary operator / two argument function
// Applied in the same way to all vector components
template<typename Tin1, typename Tin2, typename F> bool rpnMathFunc2(Value* arg1, Value* arg2, Status& s=Status::ignore) {
    // Functor
    F functor;

    // Scalar agument type
    using Tsin1 = Value::ScalarType<Tin1>;
    using Tsin2 = Value::ScalarType<Tin2>;
    
    // Scalar result type
    using Tsout = decltype(functor(Tsin1(), Tsin2()));

    // Actual result type (vectorize if neccessary)
    using Tout = Value::VectorType<Tsout,Value::IsVectorType<Tin1>::value||Value::IsVectorType<Tin2>::value>;

    // Get vector size, check vector compatibility
    size_t n = 0;
    if constexpr(Value::IsVectorType<Tin1>::value && Value::IsVectorType<Tin2>::value) {
        if (arg1->size()!=arg2->size()) {
            s.set(Status::BadArguments, "Vector size mismatch.");
            return false;
        }
        n = arg1->size();
    } else if constexpr(Value::IsVectorType<Tin1>::value) {
        n = arg1->size();
    } else if constexpr(Value::IsVectorType<Tin2>::value) {
        n = arg2->size();
    }

    // Check if we need temporary storage, prepare it
    Value *dest;
    Value tmpVal;
    if (Value::typeCode<Tout>()==arg1->type()) {
        // Reuse argv[0] as destination
        dest = arg1;
    } else if (Value::typeCode<Tout>()==arg2->type()) {
        // Reuse argv[0] as destination
        dest = arg2;
    } else {
        // Temporary destination
        dest = &tmpVal;
        if constexpr(Value::IsVectorType<Tout>::value) {
            // Initialize vector, set value type
            tmpVal = std::move(Tout(n));
        }
        // For scalar result the type is set at result assignment
    }

    // Compute
    if constexpr(!Value::IsVectorType<Tin1>::value && !Value::IsVectorType<Tin2>::value) {
        // Scalar - scalar -> scalar
        auto& v1 = arg1->val<Tin1>();
        auto& v2 = arg2->val<Tin2>();
        if (!functor.ok(v1, v2, s)) {
            return false;
        }
        *dest = functor(v1, v2);
    } else if constexpr(Value::IsVectorType<Tin1>::value && !Value::IsVectorType<Tin2>::value) {
        // Vector - scalar -> vector
        for(size_t i=0; i<n; i++) {
            auto& v1 = arg1->val<Tin1>()[i];
            auto& v2 = arg2->val<Tin2>();
            if (!functor.ok(v1, v2, s)) {
                return false;
            }
            dest->val<Tout>()[i] = functor(v1, v2);
        }
    } else if constexpr(!Value::IsVectorType<Tin1>::value && Value::IsVectorType<Tin2>::value) {
        // Scalar - vector -> vector
         for(size_t i=0; i<n; i++) {
            auto& v1 = arg1->val<Tin1>();
            auto& v2 = arg2->val<Tin2>()[i];
            if (!functor.ok(v1, v2, s)) {
                return false;
            }
            dest->val<Tout>()[i] = functor(v1, v2);
        }
    } else {
        // Vector - vector -> vector
        for(size_t i=0; i<n; i++) {
            auto& v1 = arg1->val<Tin1>()[i];
            auto& v2 = arg2->val<Tin2>()[i];
            if (!functor.ok(v1, v2, s)) {
                return false;
            }
            dest->val<Tout>()[i] = functor(v1, v2);
        }
    }

    // Swap with arg1
    swap(*arg1, *dest);

    return true;
}

// Function with one argument that results in an aggregate scalar value
template<typename Tin, typename F> bool rpnAggregateFunc1(Value* arg, Status& s=Status::ignore) {
    // Functor
    F functor;
    
    // Check if we need temporary storage, prepare it
    Value *dest;
    Value tmpVal;

    // Result type
    using Tout = decltype(functor(Tin()));

    if (Value::typeCode<Tout>()==arg->type()) {
        // Reuse argv as destination
        dest = arg;
    } else {
        // Temporary destination
        dest = &tmpVal;
    }

    // Compute
    auto& v = arg->val<Tin>();
    if (!functor.ok(v, s)) {
        return false;
    }
    *dest = std::move(functor(v));

    // Swap with arg
    swap(*arg, *dest);

    return true;
}

// Math function with 1 argument, applied component-wise
template<typename F> inline bool mathFuncComp1(RpnStack& stack, Rpn::Arity argc, Status& s) { 
    Value* vp = stack.get(); 
    DBGCHECK(!vp, "Internal error. Attempt to get value from empty stack."); 
    bool coreSt = false;
    switch (vp->type()) { 
        case Value::Type::Int:     coreSt = rpnMathFunc1<Int, F>(vp, s); break; 
        case Value::Type::Real:    coreSt = rpnMathFunc1<Real, F>(vp, s); break; 
        case Value::Type::IntVec:  coreSt = rpnMathFunc1<IntVector, F>(vp, s); break; 
        case Value::Type::RealVec: coreSt = rpnMathFunc1<RealVector, F>(vp, s); break; 
        default: 
            s.set( 
                Status::Unsupported, 
                std::string("Operation not supported for ")+vp->typeName()+"." 
            ); 
    } 
    // No need to pop anything, the result replaced the argument
    return coreSt;
}

// Math function with 2 arguments, applied component-wise
template<typename F> inline bool mathFuncComp2(RpnStack& stack, Rpn::Arity argc, Status& s) { 
    Value* v1p = stack.get(1); 
    DBGCHECK(!v1p, "Internal error. Attempt to get value from empty stack."); 
    Value* v2p = stack.get(); 
    DBGCHECK(!v2p, "Internal error. Attempt to get value from empty stack."); 
    bool coreSt = false;
    switch (type_pair(v1p->type(), v2p->type())) { 
        case type_pair(Value::Type::Int, Value::Type::Int):     coreSt = rpnMathFunc2<Int, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Int, Value::Type::Real):    coreSt = rpnMathFunc2<Int, Real, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Int, Value::Type::IntVec):  coreSt = rpnMathFunc2<Int, IntVector, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Int, Value::Type::RealVec): coreSt = rpnMathFunc2<Int, RealVector, F>(v1p, v2p, s); break; 
        
        case type_pair(Value::Type::Real, Value::Type::Int):     coreSt = rpnMathFunc2<Real, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Real, Value::Type::Real):    coreSt = rpnMathFunc2<Real, Real, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Real, Value::Type::IntVec):  coreSt = rpnMathFunc2<Real, IntVector, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Real, Value::Type::RealVec): coreSt = rpnMathFunc2<Real, RealVector, F>(v1p, v2p, s); break; 
        
        case type_pair(Value::Type::IntVec, Value::Type::Int):     coreSt = rpnMathFunc2<IntVector, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::IntVec, Value::Type::Real):    coreSt = rpnMathFunc2<IntVector, Real, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::IntVec, Value::Type::IntVec):  coreSt = rpnMathFunc2<IntVector, IntVector, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::IntVec, Value::Type::RealVec): coreSt = rpnMathFunc2<IntVector, RealVector, F>(v1p, v2p, s); break; 
        
        case type_pair(Value::Type::RealVec, Value::Type::Int):     coreSt = rpnMathFunc2<RealVector, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::RealVec, Value::Type::Real):    coreSt = rpnMathFunc2<RealVector, Real, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::RealVec, Value::Type::IntVec):  coreSt = rpnMathFunc2<RealVector, IntVector, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::RealVec, Value::Type::RealVec): coreSt = rpnMathFunc2<RealVector, RealVector, F>(v1p, v2p, s); break; 
        
        default: 
            s.set( 
                Status::Unsupported, 
                std::string("Operation not supported for ")+v1p->typeName()+" and "+v2p->typeName()+"." 
            ); 
    } 
    // Result replaces first argument, pop second argument
    stack.pop();
    return coreSt;
}

// Selector, 2 arguments, second argument must be an Int or an Int vector
bool mathFuncSelector2(RpnStack& stack, Rpn::Arity argc, Status& s);

// Aggregates scalars/vectors into a numeric value
template<typename F> inline bool mathAggregateFunc1(RpnStack& stack, Rpn::Arity argc, Status& s) { 
    Value* vp = stack.get(); 
    DBGCHECK(!vp, "Internal error. Attempt to get value from empty stack."); 
    bool coreSt = false;
    switch (vp->type()) { 
        case Value::Type::Int:     coreSt = rpnAggregateFunc1<Int, F>(vp, s); break; 
        case Value::Type::Real:    coreSt = rpnAggregateFunc1<Real, F>(vp, s); break; 
        case Value::Type::String:    coreSt = rpnAggregateFunc1<String, F>(vp, s); break; 
        
        case Value::Type::IntVec:     coreSt = rpnAggregateFunc1<IntVector, F>(vp, s); break; 
        case Value::Type::RealVec:    coreSt = rpnAggregateFunc1<RealVector, F>(vp, s); break; 
        case Value::Type::StringVec:    coreSt = rpnAggregateFunc1<StringVector, F>(vp, s); break; 
        
        default: 
            s.set( 
                Status::Unsupported, 
                std::string("Function does not support ")+vp->typeName()+"." 
            ); 
    } 
    // No need to pop anything, the result replaced the argument
    return coreSt;
}

// Aggregates numerical scalars/vectors into a numeric value
template<typename F> inline bool mathAggregateNumFunc1(RpnStack& stack, Rpn::Arity argc, Status& s) { 
    Value* vp = stack.get(); 
    DBGCHECK(!vp, "Internal error. Attempt to get value from empty stack."); 
    bool coreSt = false;
    switch (vp->type()) { 
        case Value::Type::Int:     coreSt = rpnAggregateFunc1<Int, F>(vp, s); break; 
        case Value::Type::Real:    coreSt = rpnAggregateFunc1<Real, F>(vp, s); break; 
        
        case Value::Type::IntVec:     coreSt = rpnAggregateFunc1<IntVector, F>(vp, s); break; 
        case Value::Type::RealVec:    coreSt = rpnAggregateFunc1<RealVector, F>(vp, s); break; 
        
        default: 
            s.set( 
                Status::Unsupported, 
                std::string("Function does not support ")+vp->typeName()+"." 
            ); 
    } 
    // No need to pop anything, the result replaced the argument
    return coreSt;
}

// Relational operator with 2 arguments, applied component-wise
template<typename F> inline bool mathRelOpComp2(RpnStack& stack, Rpn::Arity argc, Status& s) { 
    Value* v1p = stack.get(1); 
    DBGCHECK(!v1p, "Internal error. Attempt to get value from empty stack."); 
    Value* v2p = stack.get(); 
    DBGCHECK(!v2p, "Internal error. Attempt to get value from empty stack."); 
    bool coreSt = false;
    switch (type_pair(v1p->type(), v2p->type())) { 
        case type_pair(Value::Type::Int, Value::Type::Int):     coreSt = rpnMathFunc2<Int, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Int, Value::Type::Real):    coreSt = rpnMathFunc2<Int, Real, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Int, Value::Type::IntVec):  coreSt = rpnMathFunc2<Int, IntVector, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Int, Value::Type::RealVec): coreSt = rpnMathFunc2<Int, RealVector, F>(v1p, v2p, s); break; 
        
        case type_pair(Value::Type::Real, Value::Type::Int):     coreSt = rpnMathFunc2<Real, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Real, Value::Type::Real):    coreSt = rpnMathFunc2<Real, Real, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Real, Value::Type::IntVec):  coreSt = rpnMathFunc2<Real, IntVector, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Real, Value::Type::RealVec): coreSt = rpnMathFunc2<Real, RealVector, F>(v1p, v2p, s); break; 
        
        case type_pair(Value::Type::IntVec, Value::Type::Int):     coreSt = rpnMathFunc2<IntVector, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::IntVec, Value::Type::Real):    coreSt = rpnMathFunc2<IntVector, Real, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::IntVec, Value::Type::IntVec):  coreSt = rpnMathFunc2<IntVector, IntVector, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::IntVec, Value::Type::RealVec): coreSt = rpnMathFunc2<IntVector, RealVector, F>(v1p, v2p, s); break; 
        
        case type_pair(Value::Type::RealVec, Value::Type::Int):     coreSt = rpnMathFunc2<RealVector, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::RealVec, Value::Type::Real):    coreSt = rpnMathFunc2<RealVector, Real, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::RealVec, Value::Type::IntVec):  coreSt = rpnMathFunc2<RealVector, IntVector, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::RealVec, Value::Type::RealVec): coreSt = rpnMathFunc2<RealVector, RealVector, F>(v1p, v2p, s); break; 
        
        case type_pair(Value::Type::String, Value::Type::String):     coreSt = rpnMathFunc2<String, String, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::String, Value::Type::StringVec):     coreSt = rpnMathFunc2<String, StringVector, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::StringVec, Value::Type::String):     coreSt = rpnMathFunc2<StringVector, String, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::StringVec, Value::Type::StringVec):     coreSt = rpnMathFunc2<StringVector, StringVector, F>(v1p, v2p, s); break; 

        default: 
            s.set( 
                Status::Unsupported, 
                std::string("Operation not supported for ")+v1p->typeName()+" and "+v2p->typeName()+"." 
            ); 
    } 
    // Result replaces first argument, pop second argument
    stack.pop();
    return coreSt;
}

// Bitwise operator with 1 argument, applied component-wise
template<typename F> inline bool mathBitOpComp1(RpnStack& stack, Rpn::Arity argc, Status& s) { 
    Value* vp = stack.get(); 
    DBGCHECK(!vp, "Internal error. Attempt to get value from empty stack."); 
    bool coreSt = false;
    switch (vp->type()) { 
        case Value::Type::Int:     coreSt = rpnMathFunc1<Int, F>(vp, s); break; 
        case Value::Type::IntVec:  coreSt = rpnMathFunc1<IntVector, F>(vp, s); break; 
        default: 
            s.set( 
                Status::Unsupported, 
                std::string("Operation not supported for ")+vp->typeName()+"." 
            ); 
    } 
    // No need to pop anything, the result replaced the argument
    return coreSt;
}

// Bitwise operator with 2 arguments, applied component-wise
template<typename F> inline bool mathBitOpComp2(RpnStack& stack, Rpn::Arity argc, Status& s) { 
    Value* v1p = stack.get(1); 
    DBGCHECK(!v1p, "Internal error. Attempt to get value from empty stack."); 
    Value* v2p = stack.get(); 
    DBGCHECK(!v2p, "Internal error. Attempt to get value from empty stack."); 
    bool coreSt = false;
    switch (type_pair(v1p->type(), v2p->type())) { 
        case type_pair(Value::Type::Int, Value::Type::Int):     coreSt = rpnMathFunc2<Int, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::IntVec, Value::Type::Int):     coreSt = rpnMathFunc2<IntVector, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Int, Value::Type::IntVec):     coreSt = rpnMathFunc2<Int, IntVector, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::IntVec, Value::Type::IntVec):     coreSt = rpnMathFunc2<IntVector, IntVector, F>(v1p, v2p, s); break; 
        
        default: 
            s.set( 
                Status::Unsupported, 
                std::string("Operation not supported for ")+v1p->typeName()+" and "+v2p->typeName()+"." 
            ); 
    } 
    // Result replaces first argument, pop second argument
    stack.pop();
    return coreSt;
}

// Logical operator with 1 argument, applied only to scalars
template<typename F> inline bool mathLogicOp1(RpnStack& stack, Rpn::Arity argc, Status& s) { 
    Value* vp = stack.get(); 
    DBGCHECK(!vp, "Internal error. Attempt to get value from empty stack."); 
    bool coreSt = false;
    switch (vp->type()) { 
        case Value::Type::Int:     coreSt = rpnMathFunc1<Int, F>(vp, s); break; 
        case Value::Type::Real:     coreSt = rpnMathFunc1<Real, F>(vp, s); break; 
        case Value::Type::String:     coreSt = rpnMathFunc1<String, F>(vp, s); break; 
        
        default: 
            s.set( 
                Status::Unsupported, 
                std::string("Operation not supported for ")+vp->typeName()+". Use all() or any() for vectors." 
            ); 
    } 
    // No need to pop anything, the result replaced the argument
    return coreSt;
}

// Relational operator with 2 arguments, applied only to scalars
template<typename F> inline bool mathLogicOp2(RpnStack& stack, Rpn::Arity argc, Status& s) { 
    Value* v1p = stack.get(1); 
    DBGCHECK(!v1p, "Internal error. Attempt to get value from empty stack."); 
    Value* v2p = stack.get(); 
    DBGCHECK(!v2p, "Internal error. Attempt to get value from empty stack."); 
    bool coreSt = false;
    switch (type_pair(v1p->type(), v2p->type())) { 
        case type_pair(Value::Type::Int, Value::Type::Int):       coreSt = rpnMathFunc2<Int, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Int, Value::Type::Real):      coreSt = rpnMathFunc2<Int, Real, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Int, Value::Type::String):    coreSt = rpnMathFunc2<Int, String, F>(v1p, v2p, s); break; 
        
        case type_pair(Value::Type::Real, Value::Type::Int):       coreSt = rpnMathFunc2<Real, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Real, Value::Type::Real):      coreSt = rpnMathFunc2<Real, Real, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::Real, Value::Type::String):    coreSt = rpnMathFunc2<Real, String, F>(v1p, v2p, s); break; 
        
        case type_pair(Value::Type::String, Value::Type::Int):       coreSt = rpnMathFunc2<String, Int, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::String, Value::Type::Real):      coreSt = rpnMathFunc2<String, Real, F>(v1p, v2p, s); break; 
        case type_pair(Value::Type::String, Value::Type::String):    coreSt = rpnMathFunc2<String, String, F>(v1p, v2p, s); break; 
        
        default: 
            s.set( 
                Status::Unsupported, 
                std::string("Operation not supported for ")+v1p->typeName()+" and "+v2p->typeName()+". Use all() or any() for vectors." 
            ); 
    } 
    // Result replaces first argument, pop second argument
    stack.pop();
    return coreSt;
}

// Type checking function
template<Value::Type typeCode> inline bool scalarTypeCheck(RpnStack& stack, Rpn::Arity argc, Status& s) {
    Value* vp = stack.get(); 
    DBGCHECK(!vp, "Internal error. Attempt to get value from empty stack."); 
    // Get type code, remove vector bit
    *vp = Int(vp->scalarType()==typeCode);
    return true;
}

// Vector checking function
bool vectorCheck(RpnStack& stack, Rpn::Arity argc, Status& s);

// Vector checking function
bool listCheck(RpnStack& stack, Rpn::Arity argc, Status& s);

// Conversion function
template<Value::Type typeCode> inline bool typeConversion(RpnStack& stack, Rpn::Arity argc, Status& s) {
    Value* vp = stack.get(); 
    DBGCHECK(!vp, "Internal error. Attempt to get value from empty stack."); 
    // Need conversion?
    if (vp->scalarType()==typeCode) {
        return true;
    }
    // Destination type code
    Value::Type destType;
    if (vp->isVector()) {
        destType = Value::Type(typeCode | Value::Type::VectorBit);
    }
    // Convert
    return vp->convertInPlace(destType, s);
}

// Vector/list length function
bool len(RpnStack& stack, Rpn::Arity argc, Status& s);

// Construct a vector with n components that have the same value
bool vectorBuild(RpnStack& stack, Rpn::Arity argc, Status& s);

// Construct a vector holding a range of values
bool vectorRange(RpnStack& stack, Rpn::Arity argc, Status& s);

// Pack scalars and vectors in a vector [ 1, 2, 3 ]
bool vectorPack(RpnStack& stack, Rpn::Arity argc, Status& s);

// Pack scalars, vectors, and lists in a list [ 1; 2; 3 ]
bool listPack(RpnStack& stack, Rpn::Arity argc, Status& s);

// Concatenate lists [ a : b : c ]
bool listMerge(RpnStack& stack, Rpn::Arity argc, Status& s);

// Min and max wrapper (1 or 2 arguments)
bool minWrapper(RpnStack& stack, Rpn::Arity argc, Status& s);
bool maxWrapper(RpnStack& stack, Rpn::Arity argc, Status& s);

}

#endif
