#include <unordered_set>
#include "circuit.h"
#include "common.h"


namespace NAMESPACE {

// Update global context with variables
// Called in 
//   Circuit constructor 
//   Circuit::elaborateChanges() when the VariablesChanged flag is enabled
bool Circuit::updateGlobalContext(Status& s) {
    auto& cs = paramEvaluator_.contextStack();

    // Clear context stack of parameter evaluator
    cs.clear();

    // Insert variables context and add it to context search path
    cs.enter(&variables, true);

    return true;
}

// Used by ParameterSweeper when binding to a variable and storing the variable state
const Value* Circuit::getVariable(Id name, Status& s) const {
    auto ptr = variables.get(name);
    if (!ptr) {
        s.set(Status::NotFound, "Variable '"+std::string(name)+"' not found.");
        return nullptr;
    }
    return ptr;
}

// Used by 
//   the ParameterSweeper for setting a variable
//   the command interpreter when runing the var command and setting the PYTHON variable
std::tuple<bool, bool> Circuit::setVariable(Id name, const Value& v, Status& s) {
    auto [inserted, changed] = variables.insertAndCheck(name, v);
    // Newly inserted variable or a change to an existing variable sets the VariablesChnaged flag
    changed |= inserted;
    if (changed) {
        setFlags(Flags::VariablesChanged);
    }
    return std::make_tuple(true, changed);
}

// Used by command interpreter for clearing the variables in the clear command 
bool Circuit::clearVariables(Status& s) {
    variables.clear();
    setFlags(Flags::VariablesChanged);
    return true;
}

}
