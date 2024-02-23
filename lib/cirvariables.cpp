#include "circuit.h"
#include "common.h"


namespace NAMESPACE {

// Update global context with global parameters
bool Circuit::updateGlobalContext(Status& s) {
    auto& cs = paramEvaluator_.contextStack();

    // Clear context stack of parameter evaluator
    cs.clear();

    // Insert variables context and add it to context search path
    cs.enter(&variables, true);

    return true;
}

const Value* Circuit::getVariable(Id name, Status& s) const {
    auto ptr = variables.get(name);
    if (!ptr) {
        s.set(Status::NotFound, "Variable '"+std::string(name)+"' not found.");
        return nullptr;
    }
    return ptr;
}

bool Circuit::setVariable(Id name, const Value& v, Status& s) {
    auto [dummy, changed] = variables.insertAndCheck(name, v);
    if (changed) {
        setFlags(Flags::VariablesChanged);
    }
    return true;
}

bool Circuit::clearVariables(Status& s) {
    variables.clear();
    setFlags(Flags::VariablesChanged);
    return true;
}

}
