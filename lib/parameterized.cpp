#include "parameterized.h"
#include "common.h"

namespace NAMESPACE {

Parameterized::Parameterized() {
}

std::tuple<Value::Type,bool> Parameterized::parameterType(Id name, Status& s) const {
    auto [ndx, found] = parameterIndex(name);
    if (!found) {
        s.set(Status::NotFound, std::string("Parameter '")+std::string(name)+"' not found.");
        return std::make_tuple(Value::Type::Int, false);
    }
    return parameterType(ndx, s);
}

bool Parameterized::getParameter(Id name, Value& v, Status& s) const {
    auto [ndx, found] = parameterIndex(name);
    if (!found) {
        s.set(Status::NotFound, std::string("Parameter '")+std::string(name)+"' not found.");
        return false;
    }
    return getParameter(ndx, v, s);
}

std::tuple<bool,bool> Parameterized::setParameter(Id name, const Value& v, Status& s) {
    auto [ndx, found] = parameterIndex(name);
    if (!found) {
        s.set(Status::NotFound, std::string("Parameter '")+std::string(name)+"' not found.");
        return std::make_tuple(false, false);
    }
    return setParameter(ndx, v, s);
}

std::tuple<bool,bool> Parameterized::parameterGiven(ParameterIndex ndx, Status& s) {
    // By default parameter is always given (structure default value is used)
    return std::make_tuple(true, true);
}

std::tuple<bool,bool> Parameterized::parameterGiven(Id name, Status& s) {
    auto [ndx, found] = parameterIndex(name);
    if (!found) {
        s.set(Status::NotFound, std::string("Parameter '")+std::string(name)+"' not found.");
        return std::make_tuple(false, false);
    }
    return parameterGiven(ndx, s);
}

std::tuple<bool,bool> Parameterized::setParameters(const std::vector<PTParameterValue>& params, Status& s) {
    bool changed = false;
    for(auto it=params.cbegin(); it!=params.cend(); ++it) {
        auto [ok, ch] = setParameter(it->name(), it->val(), s);
        changed = changed | ch;
        if (!ok) {
            s.extend(it->location());
            return std::make_tuple(false, changed);
        }
    }
    return std::make_tuple(true, changed);
}

std::tuple<bool,bool> Parameterized::setParameters(const std::vector<PTParameterExpression>& params, RpnEvaluator& eval, Status& s) {
    // Assume the context is already set up in evaluator
    bool changed = false;
    for(auto it=params.cbegin(); it!=params.cend(); ++it) {
        Value res;
        if (!eval.evaluate(it->rpn(), res, s)) {
            return std::make_tuple(false, changed);
        }
        auto [ok, ch] = setParameter(it->name(), res, s);
        changed = changed | ch;
        if (!ok) {
            s.extend(it->location());
            return std::make_tuple(false, changed);
        }
    }
    return std::make_tuple(true, changed);
}

std::tuple<bool,bool> Parameterized::setParameters(const PTParameters& params, RpnEvaluator& eval, Status& s) {
    auto [ok1, changed] = setParameters(params.values(), s);
    if (!ok1) {
        return std::make_tuple(false, changed);
    }
    auto [ok2, ch] = setParameters(params.expressions(), eval, s);
    changed = changed | ch;
    if (!ok2) {
        return std::make_tuple(false, changed);
    }
    return std::make_tuple(true, changed);
}

std::tuple<bool,bool> Parameterized::setParameters(const PTParameterMap& params, RpnEvaluator& eval, Status& s) {
    // Go through parameter map, set values
    bool changed = false;
    for(auto& it : params) {
        if (std::holds_alternative<const PTParameterValue*>(it.second)) {
            auto* pv = std::get<const PTParameterValue*>(it.second); 
            auto [ok, ch] = setParameter(it.first, pv->val(), s);
            changed |= ch;
            if (!ok) {
                s.extend(pv->location());
                return std::make_tuple(false, changed);
            }
        } else {
            auto* pe = std::get<const PTParameterExpression*>(it.second); 
            Value res;
            if (!eval.evaluate(pe->rpn(), res, s)) {
                return std::make_tuple(false, changed);
            }
            auto [ok, ch] = setParameter(it.first, res, s);
            changed = changed | ch;
            if (!ok) {
                s.extend(pe->location());
                return std::make_tuple(false, changed);
            }
        }
    }
    return std::make_tuple(true, changed);
}

void Parameterized::dump(std::ostream& os, const char* prefix) const {
    for(ParameterIndex i=0; i<parameterCount(); i++) {
        if (i>0) {
            os << "\n";
        }
        Value v;
        getParameter(i, v);
        os << prefix << std::string(parameterName(i)) << " = " << v;
    }
}

}
