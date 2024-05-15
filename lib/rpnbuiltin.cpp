#include "rpnbuiltin.h"
#include <cmath>
#include "common.h"


namespace NAMESPACE {

bool vectorCheck(RpnStack& stack, Rpn::Arity argc, Status& s) {
    Value* vp = stack.get(); 
    DBGCHECK(!vp, "Internal error. Attempt to get value from empty stack."); 
    // Is it a vector type, but not a list
    *vp = Int(vp->isVector() && vp->type()!=Value::Type::ValueVec);
    return true;
}

bool listCheck(RpnStack& stack, Rpn::Arity argc, Status& s) {
    Value* vp = stack.get(); 
    DBGCHECK(!vp, "Internal error. Attempt to get value from empty stack."); 
    // Is it a list
    *vp = Int(vp->type()==Value::Type::ValueVec);
    return true;
}

bool len(RpnStack& stack, Rpn::Arity argc, Status& s) {
    Value* vp = stack.get(); 
    DBGCHECK(!vp, "Internal error. Attempt to get value from empty stack."); 
    if (!vp->isVector()) {
        s.set(Status::BadArguments, "Argument is not a vector/list."); 
        return false;
    } else {
        *vp = Int(vp->size());
    }
    return true;
}

// Variable argument count, possible counts: 1, 2
bool vectorBuild(RpnStack& stack, Rpn::Arity argc, Status& s) {
    // Number of values
    Value* np = stack.get(argc-1); 
    DBGCHECK(!np, "Internal error. Attempt to get value from empty stack."); 
    if (np->type()!=Value::Type::Int) {
        s.set(Status::BadArguments, "First argument must be an integer."); 
        return false; 
    }
    auto n = np->val<Int>();
    if (n<0) {
        s.set(Status::BadArguments, "First argument must be nonnegative."); 
        return false; 
    }
    
    // Value to put in vector
    Value dflVal = Int(0);
    Value *vp = &dflVal;
    if (argc>1) {
        vp = stack.get(argc-2); 
        DBGCHECK(!vp, "Internal error. Attempt to get value from empty stack."); 
        DBGCHECK(vp->isVector(), "Second argument must be a scalar."); 
    }

    // Construct vector
    Value res;
    switch (vp->type()) {
        case Value::Type::Int:
            res = std::move(IntVector(n, vp->val<Int>()));
            break;
        case Value::Type::Real:
            res = std::move(RealVector(n, vp->val<Real>()));
            break;
        case Value::Type::String:
            res = std::move(StringVector(n, vp->val<String>()));
            break;
    }

    swap(*np, res);

    if (argc>1) {
        stack.pop();
    }

    return true;
}

// Variable argument count, possible counts: 1, 2, 3
bool vectorRange(RpnStack& stack, Rpn::Arity argc, Status& s) {
    DBGCHECK(stack.size()<argc, "Internal error. Attempt to get value from empty stack."); 
    // Collect arguments
    auto t = Value::Type::Int;
    Value* args[3];
    args[0] = stack.get(argc-1);
    if (argc>1) {
        args[1] = stack.get(argc-2);
    }
    if (argc>2) {
        args[2] = stack.get(argc-3);
    }
    
    // Check arguments, get maximal type
    for(Rpn::Arity i=0; i<argc; i++) {
        if (args[i]->isVector() || args[i]->scalarType()==Value::Type::String) {
            s.set(Status::BadArguments, std::string("Argument ")+std::to_string(i+1)+" must be a scalar numeric value."); 
            return false; 
        }
        if (args[i]->type()==Value::Type::Real) {
            t = Value::Type::Real;
        }
    }

    // Convert types
    for(Rpn::Arity i=0; i<argc; i++) {
        args[i]->convertInPlace(t);
    }

    // Check range, build vector
    size_t n;
    Value res;
    if (t==Value::Type::Int) {
        Int v0, v1, dv;
        v0 = 0;
        v1 = args[0]->val<Int>();
        dv = 1;
        if (argc>1) {
            v0 = v1;
            v1 = args[1]->val<Int>();
        } 
        if (argc>2) {
            dv = args[2]->val<Int>();
        }
        if (dv==0) {
            s.set(Status::BadArguments, std::string("Zero step is not allowed.")); 
            return false; 
        } 
        auto nd = std::round((v1-v0)/dv);
        if (nd<0) {
            nd=0;
        }
        n = size_t(nd);
        auto vec = IntVector(n);
        for(size_t i=0; i<n; i++) {
            vec[i] = v0 + i*dv;
        }
        res = std::move(vec);
    } else {
        Real v0, v1, dv;
        v0 = 0;
        v1 = args[0]->val<Real>();
        dv = 1;
        if (argc>1) {
            v0 = v1;
            v1 = args[1]->val<Real>();
        } 
        if (argc>2) {
            dv = args[2]->val<Real>();
        }
        if (dv==0) {
            s.set(Status::BadArguments, std::string("Zero step is not allowed.")); 
            return false; 
        } 
        auto nd = std::round((v1-v0)/dv);
        if (nd<0) {
            nd=0;
        }
        n = size_t(nd);
        auto vec = RealVector(n);
        for(size_t i=0; i<n; i++) {
            vec[i] = v0 + i*dv;
        }
        res = std::move(vec);
    }

    swap(*(args[0]), res);

    // Pop all but the first argument
    stack.pop(argc-1);
    return true;
}

// Core of vectorPack preprocessing
bool vectorPackPreprocess(Value& comp, Value::Type& t, size_t& j, bool& first, Status& s) {
    if (comp.type()==Value::Type::Value) {
        // Scalar Value type (list entry type) is prohibited
        s.set(Status::Unsupported, "Scalar Value type is not supported by vector pack."); 
        return false; 
    } else if (comp.type()==Value::Type::ValueVec) {
        // Recursively check all components of list
        auto& valVec = comp.val<ValueVector>();
        auto n = valVec.size();
        for(decltype(n) i=0; i<n; i++) {
            if (!vectorPackPreprocess(valVec[i], t, j, first, s)) {
                s.extend("  component "+std::to_string(i)+" of ");
                return false;
            }
        }
        return true;
    } else {
        // Ordinary types
        if (first) {
            // Initialize result type
            t = comp.scalarType();
            first = false;
        } else {
            // Check compatibility, expand type
            if (comp.isVector() && comp.size()==0) {
                // Empty vectors are compatible with everything
            } else if (t==comp.scalarType()) {
                // Same type is compatible
            } else if (t==Value::Type::Int && comp.scalarType()==Value::Type::Real) {
                // Compatibility with Int, expand Int
                t = Value::Type::Real;
            } else if (t==Value::Type::Real && comp.scalarType()==Value::Type::Int) {
                // Compatibility with Real, nothing to do as Real is the largest type
            } else {
                s.set(Status::BadArguments, "Vector component type incompatibility at "); 
                return false; 
            }
        }

        // Vector or scalar
        if (comp.isVector()) {
            j += comp.size();
        } else {
            j++;
        }
    }

    return true;
}

// Ordinary type vector packing core
template<typename Tout, typename T> void ordinaryVectorPackCore(Value& res, Value& comp, size_t& j) {
    if constexpr(Value::IsVectorType<T>::value) {
        // Add a vector
        for(size_t i=0; i<comp.size(); i++) {
            res.val<Tout>()[j] = comp.val<T>()[i];
            j++;
        }
    } else {
        // Add a scalar
        res.val<Tout>()[j] = comp.val<T>();
        j++;
    }
}

// Core of vectorPack
void vectorPackCore(Value& res, Value& comp, size_t& j) {
    if (comp.type()==Value::Type::ValueVec) {
        // Loop through values in list, recursively call vectorPackCore()
        auto& valVec = comp.val<ValueVector>();
        for(auto& it : valVec) {
            vectorPackCore(res, it, j);
        }
    } else {
        // Packing ordinary scalar or vector
        switch (type_pair(res.scalarType(), comp.type())) {
            case type_pair(Value::Type::Int, Value::Type::Int):     ordinaryVectorPackCore<IntVector,Int>(res, comp, j); break;
            case type_pair(Value::Type::Int, Value::Type::IntVec):  ordinaryVectorPackCore<IntVector,IntVector>(res, comp, j); break;
            case type_pair(Value::Type::Int, Value::Type::RealVec): ordinaryVectorPackCore<IntVector,RealVector>(res, comp, j); break;
            
            case type_pair(Value::Type::Real, Value::Type::Int):     ordinaryVectorPackCore<RealVector,Int>(res, comp, j); break;
            case type_pair(Value::Type::Real, Value::Type::Real):    ordinaryVectorPackCore<RealVector,Real>(res, comp, j); break;
            case type_pair(Value::Type::Real, Value::Type::IntVec):  ordinaryVectorPackCore<RealVector,IntVector>(res, comp, j); break;
            case type_pair(Value::Type::Real, Value::Type::RealVec): ordinaryVectorPackCore<RealVector,RealVector>(res, comp, j); break;
            
            case type_pair(Value::Type::String, Value::Type::String):    ordinaryVectorPackCore<StringVector,String>(res, comp, j); break;
            case type_pair(Value::Type::String, Value::Type::StringVec): ordinaryVectorPackCore<StringVector,StringVector>(res, comp, j); break;
        }
    }
}

// Variable argument count, possible counts >=0
bool vectorPack(RpnStack& stack, Rpn::Arity argc, Status& s) {
    // Check argument count
    DBGCHECK(stack.size()<argc, "Internal error. Attempt to get value from empty stack."); 
    
    // Default scalar type is Int (for zero-length vectors), initial length is 0
    auto st = Value::Type::Int;
    size_t n = 0; 

    // Determine type, count entries
    // Empty vectors are Int vectors
    bool first = true;
    for(Rpn::Arity i=1; i<=argc; i++) {
        Value* vp = stack.get(argc-i);
        if (!vectorPackPreprocess(*vp, st, n, first, s)) {
            s.extend("  component "+std::to_string(i)+".");
        }
    }

    // Create vector
    Value res;
    switch (st) {
        case Value::Type::Int:
            res = std::move(IntVector(n));
            break;
        case Value::Type::Real:
            res = std::move(RealVector(n));
            break;
        case Value::Type::String:
            res = std::move(StringVector(n));
            break;
    }

    size_t j = 0;
    for(Rpn::Arity i=1; i<=argc; i++) {
        Value* vp = stack.get(argc-i);
        vectorPackCore(res, *vp, j);
    }

    if (argc>0) {
        // Swap with first argument
        swap(*stack.get(argc-1), res);
    } else {
        // Push if no arguments given
        stack.push(std::move(res));
    }

    // Pop all but the first argument, negative values of argc-1 pop nothing
    stack.pop(argc-1);
    return true;
}

bool listPack(RpnStack& stack, Rpn::Arity argc, Status& s) {
    // Create empty list
    Value res(std::move(ValueVector()));

    // Get Value vector
    auto& vVec = res.val<ValueVector>();

    // Go through arguments, move them to list entries
    for(Rpn::Arity i=1; i<=argc; i++) {
        Value* vp = stack.get(argc-i);
        vVec.push_back(std::move(*vp));
    }
    
    if (argc>0) {
        // Swap with first argument
        swap(*stack.get(argc-1), res);
    } else {
        // Push if no arguments given
        stack.push(std::move(res));
    }

    // Pop all but the first argument
    stack.pop(argc-1);
    return true;
}

bool listMerge(RpnStack& stack, Rpn::Arity argc, Status& s) {
    // Create empty list
    Value res(std::move(ValueVector()));

    // Get Value vector
    auto& vVec = res.val<ValueVector>();

    // Go through arguments
    for(Rpn::Arity i=1; i<=argc; i++) {
        Value* vp = stack.get(argc-i);
        
        if (vp->type()==Value::Type::ValueVec) {
            // List arguments are unpacked and moved to list entries
            for(auto &it : vp->val<ValueVector>()) {
                vVec.push_back(std::move(it));
            }
        } else {
            // Other arguments are moved to list entries
            vVec.push_back(std::move(*vp));
        }
        
        vVec.push_back(std::move(*vp));
    }
    
    if (argc>0) {
        // Swap with first argument
        swap(*stack.get(argc-1), res);
    } else {
        // Push if no arguments given
        stack.push(std::move(res));
    }

    // Pop all but the first argument
    stack.pop(argc-1);
    return true;
}

bool minWrapper(RpnStack& stack, Rpn::Arity argc, Status& s) {
    if (argc==2) {
        // Two arguments, component-wise
        return mathFuncComp2<FwMin>(stack, argc, s);
    } else {
        // Single argument, aggregation across vector/scalar
        return mathAggregateNumFunc1<FwMinAggregate>(stack, argc, s);
    }
}

bool maxWrapper(RpnStack& stack, Rpn::Arity argc, Status& s) {
    if (argc==2) {
        // Two arguments, component-wise
        return mathFuncComp2<FwMax>(stack, argc, s);
    } else {
        // Single argument, aggregation across vector/scalar
        return mathAggregateNumFunc1<FwMaxAggregate>(stack, argc, s);
    }
}

template<typename Tin, typename Tsel> bool rpnSelectorFunc2(Value* arg1, Value* arg2, Status& s=Status::ignore) {
    // Functor
    FwSelector functor;

    // Data vector
    auto& dat = arg1->val<Tin>();

    Value res;
    if constexpr(Value::IsVectorType<Tsel>::value) {
        // Result is a vector
        auto& sel = arg2->val<Tsel>();
        size_t n = sel.size();
        Tin resVec(n);
        for(size_t i=0; i<n; i++) {
            if (!functor.ok(dat, sel[i], s)) {
                return false;
            }
            resVec[i] = functor(dat, sel[i]);
        }
        res = std::move(resVec);
    } else {
        // Result is a scalar
        Tsel ndx = arg2->val<Tsel>();
        if (!functor.ok(dat, ndx, s)) {
            return false;
        }
        res = functor(dat, ndx);
    }

    // Swap with arg1
    swap(*arg1, res);

    return true;
}

bool mathFuncSelector2(RpnStack& stack, Rpn::Arity argc, Status& s) { 
    Value* v1p = stack.get(1); 
    DBGCHECK(!v1p, "Internal error. Attempt to get value from empty stack."); 
    Value* v2p = stack.get(); 
    DBGCHECK(!v2p, "Internal error. Attempt to get value from empty stack."); 
    // Check second value
    if (v2p->scalarType()!=Value::Type::Int) {
        s.set(Status::BadArguments, "Selector must be an integer or an integer vector."); 
        return false; 
    }
    bool coreSt = false;
    switch (type_pair(v1p->type(), v2p->type())) { 
        case type_pair(Value::Type::IntVec, Value::Type::Int):     coreSt = rpnSelectorFunc2<IntVector, Int>(v1p, v2p, s); break; 
        case type_pair(Value::Type::IntVec, Value::Type::IntVec):  coreSt = rpnSelectorFunc2<IntVector, IntVector>(v1p, v2p, s); break; 
        
        case type_pair(Value::Type::RealVec, Value::Type::Int):     coreSt = rpnSelectorFunc2<RealVector, Int>(v1p, v2p, s); break; 
        case type_pair(Value::Type::RealVec, Value::Type::IntVec):  coreSt = rpnSelectorFunc2<RealVector, IntVector>(v1p, v2p, s); break; 
        
        case type_pair(Value::Type::StringVec, Value::Type::Int):     coreSt = rpnSelectorFunc2<StringVector, Int>(v1p, v2p, s); break; 
        case type_pair(Value::Type::StringVec, Value::Type::IntVec):  coreSt = rpnSelectorFunc2<StringVector, IntVector>(v1p, v2p, s); break; 
        
        case type_pair(Value::Type::ValueVec, Value::Type::Int):     coreSt = rpnSelectorFunc2<ValueVector, Int>(v1p, v2p, s); break; 
        case type_pair(Value::Type::ValueVec, Value::Type::IntVec):  coreSt = rpnSelectorFunc2<ValueVector, IntVector>(v1p, v2p, s); break; 
        
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

}
