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

    void set(Circuit& circuit, Vector<double>& vec);
    
    const Vector<double>& values() const { return values_; };
    const std::vector<Id>& names() const { return names_; };

private:
    Vector<double> values_;
    std::vector<Id> names_;
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

    // Create from stored solution, check for conflicts
    bool set(Circuit& circuit, const AnnotatedSolution& solution, bool abortOnError);
    
    // Resolve to actual forces, check for conflicts
    // This is phase 2 of user forces processing. 
    bool set(Circuit& circuit, const PreprocessedUserForces& userForces, bool uicMode, bool abortOnError);

    // Clear forces
    void clear();
    
    // Get forces for unknowns
    const Vector<double>& unknownValue() const { return unknownValue_; }; 
    const Vector<bool>& unknownForced() const { return unknownForced_; }; 

    // Get forces for unknown differences
    const Vector<double>& deltaValue() const { return deltaValue_; }; 
    const Vector<std::tuple<UnknownIndex, UnknownIndex>>& deltaIndices() const { return deltaIndices_; }; 

    void dump(Circuit& circuit, std::ostream& os) const;
    
private:
    bool setForceOnUnknown(Circuit& circuit, Node* node, double value);
    
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
