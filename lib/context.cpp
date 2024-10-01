#include "context.h"
#include "rpneval.h"
#include "common.h"
#include <cmath>
#include <numbers>


namespace NAMESPACE {

static const double PI = std::numbers::pi;

Context ContextStack::consts{ 
    { Id::createStatic("M_E"),         exp(1)                     },  // e or exp(1)
    { Id::createStatic("M_LOG2E"),     1/log(2)                   },  // log2(e)
    { Id::createStatic("M_LOG10E"),    1/log(10)                  },  // log10(e)
    { Id::createStatic("M_LN2"),       log(2)                     },  // ln(2)
    { Id::createStatic("M_LN10"),      log(10)                    },  // ln(10)
    { Id::createStatic("M_PI"),        PI                         },  // PI
    { Id::createStatic("M_TWO_PI"),    2*PI                       },  // 2 PI
    { Id::createStatic("M_PI_2"),      PI/2                       },  // PI/2
    { Id::createStatic("M_PI_4"),      PI/4                       },  // PI/4
    { Id::createStatic("M_1_PI"),      1/PI                       },  // 1/PI
    { Id::createStatic("M_2_PI"),      2/PI                       },  // 2/PI
    { Id::createStatic("M_2_SQRTPI"),  2/sqrt(PI)                 },  // 2/PI
    { Id::createStatic("M_SQRT2"),     sqrt(2)                    },  // 2
    { Id::createStatic("M_SQRT1_2"),   1/sqrt(2)                  },  // 1⁄2
    { Id::createStatic("M_DEGPERRAD"), 180/PI                     },  // Number of degrees per radian
    { Id::createStatic("P_Q"),         1.6021918e-19              },  // Charge of electron in coulombs
    { Id::createStatic("P_C"),         2.997924562e8              },  // Speed of light in vacuum in meters/sec
    { Id::createStatic("P_K"),         1.3806226e-23              },  // Boltzman’s constant in joules/Kelvin
    { Id::createStatic("P_H"),         6.6260755e-34              },  // Planck’s constant in joules times seco
    { Id::createStatic("P_EPS0"),      8.85418792394420013968e-12 },  // Permittivity of vacuum in farads/meter
    { Id::createStatic("P_U0"),        PI*4.0e-7                  },  // Permeability of vacuum in henrys/meter
    { Id::createStatic("P_CELSIUS0"),  273.15                     },  // Zero Celsius in Kelvin
};

ContextStack::Builtins ContextStack::builtins = {
    // One argument, math
    { Id::createStatic("sin"),      { 1, 1, true, mathFuncComp1<FwSin> } }, 
    { Id::createStatic("cos"),      { 1, 1, true, mathFuncComp1<FwCos> } }, 
    { Id::createStatic("tan"),      { 1, 1, true, mathFuncComp1<FwTan> } }, 
    { Id::createStatic("asin"),     { 1, 1, true, mathFuncComp1<FwAsin> } }, 
    { Id::createStatic("acos"),     { 1, 1, true, mathFuncComp1<FwAcos> } }, 
    { Id::createStatic("atan"),     { 1, 1, true, mathFuncComp1<FwAtan> } }, 
    { Id::createStatic("sinh"),     { 1, 1, true, mathFuncComp1<FwSinh> } }, 
    { Id::createStatic("cosh"),     { 1, 1, true, mathFuncComp1<FwCosh> } }, 
    { Id::createStatic("tanh"),     { 1, 1, true, mathFuncComp1<FwTanh> } }, 
    { Id::createStatic("asinh"),    { 1, 1, true, mathFuncComp1<FwAsinh> } },  
    { Id::createStatic("acosh"),    { 1, 1, true, mathFuncComp1<FwAcosh> } },  
    { Id::createStatic("atanh"),    { 1, 1, true, mathFuncComp1<FwAtanh> } }, 
    
    { Id::createStatic("log"),      { 1, 1, true, mathFuncComp1<FwLn> } }, 
    { Id::createStatic("ln"),       { 1, 1, true, mathFuncComp1<FwLn> } }, 
    { Id::createStatic("log10"),    { 1, 1, true, mathFuncComp1<FwLog10> } },  
    { Id::createStatic("exp"),      { 1, 1, true, mathFuncComp1<FwExp> } }, 
    { Id::createStatic("sqrt"),     { 1, 1, true, mathFuncComp1<FwSqrt> } }, 
    
    { Id::createStatic("abs"),      { 1, 1, true, mathFuncComp1<FwAbs> } }, 
    { Id::createStatic("integer"),  { 1, 1, true, mathFuncComp1<FwTrunc> } }, 
    { Id::createStatic("floor"),    { 1, 1, true, mathFuncComp1<FwFloor> } },  
    { Id::createStatic("ceil"),     { 1, 1, true, mathFuncComp1<FwCeil> } }, 
    { Id::createStatic("round"),    { 1, 1, true, mathFuncComp1<FwRound> } }, 
    { Id::createStatic("sgn"),      { 1, 1, true, mathFuncComp1<FwSgn> } }, 

    // Component-wise check for inf/nan
    { Id::createStatic("isinf"),    { 1, 1, true, mathFuncComp1<FwIsInf> } }, 
    { Id::createStatic("isnan"),    { 1, 1, true, mathFuncComp1<FwIsNan> } }, 
    { Id::createStatic("isfinite"), { 1, 1, true, mathFuncComp1<FwIsFinite> } }, 
    
    // Two arguments, math
    { Id::createStatic("pow"),      { 2, 2, true, mathFuncComp2<FwPower> } },
    { Id::createStatic("hypot"),    { 2, 2, true, mathFuncComp2<FwHypot> } }, 
    { Id::createStatic("atan2"),    { 2, 2, true, mathFuncComp2<FwAtan2> } },
    { Id::createStatic("sign"),     { 2, 2, true, mathFuncComp2<FwSign> } },
    { Id::createStatic("fmod"),     { 2, 2, true, mathFuncComp2<FwFmod> } },

    // One argument(aggregation) or two arguments (component-wise), math
    { Id::createStatic("min"),      { 1, 2, true, minWrapper } },
    { Id::createStatic("max"),      { 1, 2, true, maxWrapper } },
    
    // Boolean aggregation
    { Id::createStatic("all"),      { 1, 1, true, mathAggregateFunc1<FwAll> } },
    { Id::createStatic("any"),      { 1, 1, true, mathAggregateFunc1<FwAny> } },

    // Numerical aggregation
    { Id::createStatic("sum"),      { 1, 1, true, mathAggregateNumFunc1<FwSum> } },
    { Id::createStatic("prod"),     { 1, 1, true, mathAggregateNumFunc1<FwProd> } },
    
    // Type discovery
    // For vectors
    { Id::createStatic("isint"),    { 1, 1, true, scalarTypeCheck<Value::Type::Int> } }, 
    { Id::createStatic("isreal"),   { 1, 1, true, scalarTypeCheck<Value::Type::Real> } }, 
    { Id::createStatic("isstring"), { 1, 1, true, scalarTypeCheck<Value::Type::String> } }, 
    // Is vector/list
    { Id::createStatic("isvector"), { 1, 1, true, vectorCheck } }, 
    { Id::createStatic("islist"),   { 1, 1, true, listCheck } }, 
    
    // Type conversion
    { Id::createStatic("int"),      { 1, 1, true, typeConversion<Value::Type::Int> } }, 
    { Id::createStatic("real"),     { 1, 1, true, typeConversion<Value::Type::Real> } }, 
    // TODO: presently supports only type conversions built into Value class 
    //       in future add conversion of number to string and vice versa
    { Id::createStatic("string"),   { 1, 1, true, typeConversion<Value::Type::String> } }, 

    // Vector length
    { Id::createStatic("len"),      { 1, 1, true, len} }, 

    // Vector construction
    { Id::createStatic("vector"),   { 1, 2, true, vectorBuild } }, // length[, value]    type = type(value)
    { Id::createStatic("range"),    { 1, 3, true, vectorRange } }, // to | from, to | from, to, step,   type=maxtype(from, to, step)
    
    // String
    /*
    { Id::createStatic("join"),     { 1, 2, true, nullptr} }, // string vector[, separator]
    { Id::createStatic("split"),    { 1, 2, true, nullptr} }, // string, separator

    // Selection 
    { Id::createStatic("where"),    { 1, 3, true, nullptr} }, // where(s) -> indices of nonzeros, where(s, x, y) -> x if s nonzero, else y
    */
};


Context::Context()  {
}

Context::Context(std::initializer_list<std::pair<const Id, Value>> inl) 
    : data(inl) {
}

const Value* Context::get(Id name) const {
    auto it = data.find(name);
    if (it!=data.end()) {
        return &(it->second);
    } else {
        return nullptr;
    }
}

Value* Context::get(Id name) {
    auto it = data.find(name);
    if (it!=data.end()) {
        return &(it->second);
    } else {
        return nullptr;
    }
}

bool Context::insert(Id name, const Value& v) {
    auto it = data.find(name);
    if (it!=data.end()) {
        it->second = v;
        return false;
    } else {
        data.insert({name, v});
        return true;
    }
}

bool Context::insert(Id name, Value&& v) {
    auto it = data.find(name);
    if (it!=data.end()) {
        it->second = std::move(v);
        return false;
    } else {
        data.insert({name, std::move(v)});
        return true;
    }
}

std::tuple<bool, bool> Context::insertAndCheck(Id name, const Value& v) {
    bool changed = false;
    auto it = data.find(name);
    if (it!=data.end()) {
        changed = v != it->second;
        it->second = v;
        return std::make_tuple(false, changed);
    } else {
        data.insert({name, v});
        return std::make_tuple(true, changed);
    }
}

std::tuple<bool, bool> Context::insertAndCheck(Id name, Value&& v) {
    bool changed = false;
    auto it = data.find(name);
    if (it!=data.end()) {
        changed = v != it->second;
        it->second = std::move(v);
        return std::make_tuple(false, changed);
    } else {
        data.insert({name, std::move(v)});
        return std::make_tuple(true, changed);
    }
}

void Context::clear() {
    data.clear();
}


void Context::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    for(auto& it : data) {
        os << pfx << it.first << " = " << it.second << "\n";
    }
}

ContextStack::ContextStack() {
}

void ContextStack::clear() {
    stack.clear();
    searchPath.clear();
}

Context& ContextStack::at(size_t i) {
    auto& entry = stack.at(i);
    auto* cx = std::get_if<Context>(&entry);
    if (!cx) {
        return *std::get<Context*>(entry);
    }
    return *cx;
}

const Context& ContextStack::at(size_t i) const {
    auto& entry = stack.at(i);
    auto* cx = std::get_if<Context>(&entry);
    if (!cx) {
        return *std::get<Context*>(entry);
    }
    return *cx;
}

Context& ContextStack::at() {
    if (stack.size()==0) {
        throw std::out_of_range("No context available.");
    }
    return at(stack.size()-1);
}

const Context& ContextStack::at() const {
    if (stack.size()==0) {
        throw std::out_of_range("No context available.");
    }
    return at(stack.size()-1);
}

const Value* ContextStack::get(Id name) const {
    // Do we have any context
    if (stack.size()>0) {
        // Look in top context if it exists
        auto& cx = at();
        auto ptr = cx.get(name);
        if (ptr) {
            return ptr;
        }
        
        // Look in search path
        for(auto ci=searchPath.rbegin(); ci!=searchPath.rend(); ++ci) {
            // If context from search path is the top context, skip it
            if (*ci==stack.size()-1) {
                continue;
            }
            // Get context and search in it
            auto& cx = at(*ci);
            auto ptr = cx.get(name);
            if (ptr) {
                return ptr;
            }
        }
    }

    // Look in constants
    auto ptr = consts.get(name);
    if (ptr) {
        return ptr;
    }

    return nullptr;
}

const Builtin* ContextStack::getBuiltin(Id name) {
    auto it = builtins.find(name);
    if (it==builtins.end()) {
        return nullptr;
    }
    return &(it->second);
}

bool ContextStack::insert(Id name, const Value& v, Status& s) {
    if (isConstant(name)) {
        s.set(Status::Redefinition, "Redefinition of constants is not allowed.");
        return false;
    }
    auto& cx = at();
    cx.insert(name, v);
    return true;
}

bool ContextStack::insert(Id name, Value&& v, Status& s) {
    if (isConstant(name)) {
        s.set(Status::Redefinition, "Redefinition of constants is not allowed.");
        return false;
    }
    auto& cx = at();
    cx.insert(name, std::move(v));
    return true;
}

std::tuple<bool, bool> ContextStack::insertAndCheck(Id name, const Value& v, Status& s) {
    bool changed = false;
    if (isConstant(name)) {
        s.set(Status::Redefinition, "Redefinition of constants is not allowed.");
        return std::make_tuple(false, changed);
    }
    auto& cx = at();
    auto ptr = cx.get(name);
    if (ptr) {
        changed = v != *ptr;
    }
    cx.insert(name, v);
    return std::make_tuple(true, changed);
}

std::tuple<bool, bool> ContextStack::insertAndCheck(Id name, Value&& v, Status& s) {
    bool changed = false;
    if (isConstant(name)) {
        s.set(Status::Redefinition, "Redefinition of constants is not allowed.");
        return std::make_tuple(false, changed);
    }
    auto& cx = at();
    auto ptr = cx.get(name);
    if (ptr) {
        changed = v != *ptr;
    }
    cx.insert(name, std::move(v));
    return std::make_tuple(true, changed);
}

void ContextStack::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    for(int i=0; i<stack.size(); i++) {
        os << pfx << "Context level " << i << "\n";
        for(auto& it : stack) {
            if (std::holds_alternative<Context>(it)) {
                std::get<Context>(it).dump(indent+2, os);
            } else {
                std::get<Context*>(it)->dump(indent+2, os);
            }
        }
    }
}

}
