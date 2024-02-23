#include "circuit.h"
#include "common.h"


namespace NAMESPACE {

// Setting simulator options
// Checks for changes in options that affect hierarchy or mapping

// Set single option
std::tuple<bool, bool> Circuit::setOption(Id name, const Value& v, Status& s) {
    auto [ndx, found] = simOptions.parameterIndex(name);
    if (!found) {
        s.set(Status::NotFound, "Simulator option '"+std::string(name)+"' not found.");
        return std::make_tuple(false, false);
    }
    auto [ok, changed] = simOptions.setParameter(ndx, v, s);
    if (!ok) {
        return std::make_tuple(false, false);
    }
    if (changed) {
        ;
        if (SimulatorOptions::mappingAffectingOptions.contains(name)) {
            setFlags(Flags::MappingAffectingOptionsChanged);
        }
        if (SimulatorOptions::hierarchyAffectingOptions.contains(name)) {
            setFlags(Flags::HierarchyAffectingOptionsChanged);
        }
    }
    return std::make_tuple(ok, changed);
}

// Set all options
bool Circuit::setOptions(IStruct<SimulatorOptions>& opt) {
    bool changed = simOptions.core()!=opt.core();
    bool hierarchyOptChanged = hierarchyAffectingOptionsChanged(opt.core());
    bool mappingOptChanged = mappingAffectingOptionsChanged(opt.core());
    if (hierarchyOptChanged) {
        setFlags(Flags::HierarchyAffectingOptionsChanged);
    }
    if (mappingOptChanged) {
        setFlags(Flags::MappingAffectingOptionsChanged);
    }
    
    simOptions.core() = opt.core();
    
    return changed;
}

// Set options based on PTParameters
std::tuple<bool, bool> Circuit::setOptions(const PTParameters& params, Status& s) {
    bool changed = false;
    for(auto& it : params.values()) {
        auto [ok, ch] = setOption(it.name(), it.val(), s); 
        if (!ok) {
            s.extend(it.location());
            return std::make_tuple(false, changed);
        }
        changed = changed || ch;
    }
    auto& e = paramEvaluator();
    for(auto& it : params.expressions()) {
        Value v;
        if (!e.evaluate(it.rpn(), v, s)) {
            return std::make_tuple(false, changed);
        }
        auto [ok, ch] = setOption(it.name(), v, s); 
        if (!ok) {
            s.extend(it.location());
            return std::make_tuple(false, changed);
        }
        changed = changed || ch;
    }
    return std::make_tuple(true, changed);
}

// Check if hierarchy affecting options changed when setting all options
bool Circuit::hierarchyAffectingOptionsChanged(SimulatorOptions& opt) {
    return simOptions.core().optionsDiffer(SimulatorOptions::hierarchyAffectingOptions, opt);
}

// Check if mapping affecting options changed when setting all options
bool Circuit::mappingAffectingOptionsChanged(SimulatorOptions& opt) {
    return simOptions.core().optionsDiffer(SimulatorOptions::mappingAffectingOptions, opt);
}

}
