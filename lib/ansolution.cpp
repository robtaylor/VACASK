#include <variant>
#include "ansolution.h"
#include "circuit.h"
#include "common.h"


namespace NAMESPACE {

AnnotatedSolution::AnnotatedSolution() {
}

void AnnotatedSolution::setNames(Circuit& circuit) {
    names_.clear();
    auto n = circuit.unknownCount();
    // Store all nodes (also ground). Later in Forces::set() we ignore it. 
    // In OperatingPointCore::runSolver() we check if the stored solution is 
    // consistent with the current circuit by comparing solution length 
    // (it should be equal to number_of_unknowns+1). 
    // If we stop storing ground, we also have to chage runSolver() to detect
    // consistency correctly (now it should be equal to number_of_unknowns). 
    // Also, we should change the code in OperatingPointCore::runSolver() 
    // that copies the annotated solution to solution vector in 
    // ordinary continue mode . 
    for(decltype(n) i=0; i<=n; i++) {
        auto node = circuit.reprNode(i);
        names_.push_back(node->name());
    }
}

}
