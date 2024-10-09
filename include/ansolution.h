#ifndef __ANSOLUTION_DEFINED
#define __ANSOLUTION_DEFINED

#include "ansupport.h"
#include "node.h"
#include "value.h"
#include "common.h"


namespace NAMESPACE {

class Circuit;

class AnnotatedSolution {
public:
    AnnotatedSolution();

    AnnotatedSolution           (const AnnotatedSolution&)  = delete;
    AnnotatedSolution           (      AnnotatedSolution&&) = default;
    AnnotatedSolution& operator=(const AnnotatedSolution&)  = delete;
    AnnotatedSolution& operator=(      AnnotatedSolution&&) = default;

    void setNames(Circuit& circuit);
    
    const Vector<double>& values() const { return values_; };
    Vector<double>& values() { return values_; };
    
    const std::vector<Id>& names() const { return names_; };
    std::vector<Id>& names() { return names_; };

    const Vector<double>& auxData() const { return auxData_; };
    Vector<double>& auxData() { return auxData_; };
    
private:
    // Solution vector
    // - dc: one component per unknown, index 0 is ground (bucket)
    // - hb: nt components per unknown (1 for DC, 2 for each frequency), 
    //       no bucket - index 0 is first unknown
    Vector<double> values_;

    // Names of unknowns for cross matching across slightly different circuits
    std::vector<Id> names_;

    // Vector of auxiliary data
    // - hb: list of frequencies including DC (first component)
    Vector<double> auxData_;
};

struct PreprocessedUserForces {
    std::vector<double> nodeValues;
    std::vector<Node*> nodes;
    std::vector<Id> nodeIds;
    std::vector<double> nodePairValues;
    std::vector<std::tuple<Node*, Node*>> nodePairs;
    std::vector<std::tuple<Id, Id>> nodeIdPairs;

    PreprocessedUserForces() {};

    PreprocessedUserForces           (const PreprocessedUserForces&)  = delete;
    PreprocessedUserForces           (      PreprocessedUserForces&&) = default;
    PreprocessedUserForces& operator=(const PreprocessedUserForces&)  = delete;
    PreprocessedUserForces& operator=(      PreprocessedUserForces&&) = default;

    // Preprocess user specified nodeset/ic parameter values, store node ptrs instead of string names. 
    // If syntax is bad, return error. Otherwise simply store the forced value, 
    // node and id (or node pair and ids). 
    // If node is not found, the corresponding force is ignored. 
    // Checks if all extradiagonal sparsity map entries are present. 
    // This is phase 1 of user forces processing. 
    // Return value: ok, needs to add entries to sparsity map
    std::tuple<bool, bool> set(Circuit& circuit, ValueVector& userForces, Status& s=Status::ignore);

    void clear() {
        nodeValues.clear();
        nodes.clear();
        nodeIds.clear();
        nodePairValues.clear();
        nodePairs.clear();
        nodeIdPairs.clear();
    };
};

class Forces {
public:
    enum class Error {
        OK, 
        ConflictNode, 
        ConflictDelta, 
    };

    Forces();

    Forces           (const Forces&)  = delete;
    Forces           (      Forces&&) = default;
    Forces& operator=(const Forces&)  = delete;
    Forces& operator=(      Forces&&) = default;

    // Clear error
    void clearError() { lastError = Error::OK; }; 

    // Format error, return false on error - this function is not cheap (works with strings)
    bool formatError(Status& s=Status::ignore) const; 

    // Resolve to actual forces, check for conflicts
    // This is phase 2 of user forces processing. 
    // Applies only to OP analysis (nodeset vector) and transient analysis (ic vector). 
    bool set(Circuit& circuit, const PreprocessedUserForces& userForces, bool uicMode, bool abortOnError);

    // Clear forces
    void clear();

    // Scale forces
    void resizeUnknownForces(size_t n);
    
    // Set a force on an unknown
    bool setForceOnUnknown(Node* node, double value);
    
    // Get forces for unknowns
    const Vector<double>& unknownValue() const { return unknownValue_; }; 
    const Vector<bool>& unknownForced() const { return unknownForced_; }; 

    // Get forces for unknown differences
    const Vector<double>& deltaValue() const { return deltaValue_; }; 
    const Vector<std::tuple<UnknownIndex, UnknownIndex>>& deltaIndices() const { return deltaIndices_; }; 

    void dump(Circuit& circuit, std::ostream& os) const;
    
private:
    Vector<double> unknownValue_;
    Vector<bool> unknownForced_;
    Vector<double> deltaValue_;
    Vector<std::tuple<UnknownIndex, UnknownIndex>> deltaIndices_;

    Error lastError;
    Node* errorNode1;
    Node* errorNode2;
};

}

#endif
