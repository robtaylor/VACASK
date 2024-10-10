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


std::tuple<bool, bool> PreprocessedUserForces::set(Circuit& circuit, ValueVector& userForces, Status& s) {
    clear();

    // 0 -> 1 -> 2
    // 
    // 0 = have nothing, 
    // 1 = have node1, 
    // 2 = have node2, 
    // 3 = have number for single node (end)
    // 4 = have number for node pair (end)
    int state = 0; 
    Node* node1;
    Node* node2;
    Id id1, id2;
    double value; 
    size_t nsNdx = 0;
    bool haveAllEntries = true;
    for(auto& it : userForces) {
        switch (state) {
            case 0:
                node1 = node2 = nullptr;
                if (it.type()!=Value::Type::String) {
                    s.set(Status::BadArguments, "Expecting a string at position "+std::to_string(nsNdx)+".");
                    return std::make_tuple(false, false);
                }
                id1 = it.val<String>();
                node1 = circuit.findNode(id1);
                // Node not found is an error
                // if (!node1) {
                //     s.set(Status::BadArguments, "Cannot find node '"+std::string(id1)+"' during force preprocessing.");
                //     return std::make_tuple(false, false);
                // }
                state = 1;
                break;
            case 1:
                switch (it.type()) {
                    case Value::Type::Int:
                        value = it.val<Int>();
                        state = 3;
                        break;
                    case Value::Type::Real:
                        value = it.val<Real>();
                        state = 3;
                        break;
                    case Value::Type::String:
                        id2 = it.val<String>();
                        node2 = circuit.findNode(id2);
                        // Node not found is an error
                        // if (!node2) {
                        //     s.set(Status::BadArguments, "Cannot find node '"+std::string(id2)+"' during force preprocessing.");
                        //     return std::make_tuple(false, false);
                        // }
                        state = 2;
                        break;
                    default:
                        s.set(Status::BadArguments, "Expecting a string, an integer, or a real at position "+std::to_string(nsNdx)+".");
                        return std::make_tuple(false, false);
                }
                break;
            case 2:
                switch (it.type()) {
                    case Value::Type::Int:
                        value = it.val<Int>();
                        state = 4;
                        break;
                    case Value::Type::Real:
                        value = it.val<Real>();
                        state = 4;
                        break;
                    default:
                        s.set(Status::BadArguments, "Expecting an integer or a real at position "+std::to_string(nsNdx)+".");
                        return std::make_tuple(false, false);
                }
                break;
        }
        if (state==3) {
            // Have single node, ignore force if node is not found
            if (node1) {
                nodes.push_back(node1);
                nodeIds.push_back(id1);
                nodeValues.push_back(value);
            }
            state = 0;
        } else if (state==4) {
            // Have node pair, ignore force if node is not found
            if (node1 && node2) {
                nodePairs.push_back(std::make_tuple(node1, node2)); 
                nodeIdPairs.push_back(std::make_tuple(id1, id2));
                nodePairValues.push_back(value); 
            }
            
            // Check existence of extradiagonal entries, but only if both nodes were found. 
            // No need to check if haveAllEntries is already false. 
            if (node1 && node2 && haveAllEntries) {
                auto u1 = node1->unknownIndex();
                auto u2 = node2->unknownIndex();
                auto entry12 = circuit.sparsityMap().find(MatrixEntryPosition(u1, u2));
                auto entry21 = circuit.sparsityMap().find(MatrixEntryPosition(u2, u1));
                haveAllEntries = haveAllEntries && entry12 && entry21;
            }
            state = 0;
        }
        nsNdx++;
    }
    
    return std::make_tuple(true, !haveAllEntries);
}


Forces::Forces() {
}

bool Forces::setForceOnUnknown(Node* node, double value) {
    // Unknown
    auto u = node->unknownIndex();
    // Is it a ground node? 
    if (u==0) {
        // If yes, ignore the force. 
        return true;
    }
    // Is it conflicting with a previous nodeset
    if (unknownForced_[u] && unknownValue_[u]!=value) {
        lastError = Error::ConflictNode;
        errorNode1 = node;
        return false;
    }
    unknownValue_[u] = value;
    unknownForced_[u] = true;

    return true;
}

bool Forces::set(Circuit& circuit, const PreprocessedUserForces& preprocessed, bool uicMode, bool abortOnError) {
    // Clear forced values
    unknownValue_.clear();
    unknownForced_.clear();
    deltaValue_.clear();
    deltaIndices_.clear();

    // Number of unknowns
    auto n = circuit.unknownCount();

    // Make space for forces on unknowns, set them by default to 0
    unknownValue_.resize(n+1, 0.0);
    unknownForced_.resize(n+1, false);

    bool error = false;
    
    // Set forces on unknowns
    auto nNodeForces = preprocessed.nodes.size();
    for(decltype(nNodeForces) i=0; i<nNodeForces; i++) {
        // Check if node was found
        auto node = preprocessed.nodes[i];
        auto value = preprocessed.nodeValues[i];
        if (!setForceOnUnknown(node, value)) {
            error = true;
            if (abortOnError) {
                return false;
            }
        }
    }
    
    // Set delta forces of the form v(x,0) or v(0,x), check node pairs
    auto nDeltaForces = preprocessed.nodePairs.size();  
    for(decltype(nNodeForces) i=0; i<nDeltaForces; i++) {
        // Check if both nodes were found
        auto [node1, node2] = preprocessed.nodePairs[i];
        auto [id1, id2] = preprocessed.nodeIdPairs[i]; 
        
        // Get unknowns and value
        auto u1 = node1->unknownIndex();
        auto u2 = node2->unknownIndex();
        auto value = preprocessed.nodePairValues[i];

        // Check if both nodes are ground? 
        if (u1==0 && u2==0) {
            // If yes, ignore force
            continue;
        } else if (u1==0) {
            // Check if first node is ground, convert it to a force on an unknown
            // v(0,x) = value -> v(x)=-value
            if (!setForceOnUnknown(node2, -value)) {
                error = true;
                if (abortOnError) {
                    return false;
                }
            }
        } else if (u2==0) {
            // v(x,0) = value -> v(x)=value
            if (!setForceOnUnknown(node1, value)) {
                error = true;
                if (abortOnError) {
                    return false;
                }
            }
        }
    }

    // Set real delta forces
    for(decltype(nNodeForces) i=0; i<nDeltaForces; i++) {
        // Check if both nodes were found
        auto [node1, node2] = preprocessed.nodePairs[i];
        auto [id1, id2] = preprocessed.nodeIdPairs[i]; 
        
        // Get unknowns and value
        auto u1 = node1->unknownIndex();
        auto u2 = node2->unknownIndex();
        auto value = preprocessed.nodePairValues[i];

        // Check if both nodes are ground
        if (u1==0 && u2==0) {
            // Skip this node pair
            continue;
        } else if (u1==0) {
            // v(0,x) = value -> v(x)=-value, already handled
            continue;
        } else if (u2==0) {
            // v(x,0) = value -> v(x)=value, already handled
            continue;
        } else {
            // Actual delta force 
            // Is force already set on both nodes
            if (unknownForced_[u1] && unknownForced_[u2]) {
                // Both nodes are forced
                // Does delta force conflict with node forces
                if (unknownValue_[u1]-unknownValue_[u2]!=value) {
                    lastError = Error::ConflictDelta;
                    errorNode1 = node1;
                    errorNode2 = node2;
                    error = true;
                    if (abortOnError) {
                        return false;
                    }
                } else {
                    // Matches node forces, no need to add it, skip
                    continue;
                }
            } else {
                // At least one node is not forced yet
                if (uicMode) {
                    // As UIC forces, apply to nodes
                    // One of the nodes is not forces
                    if (unknownForced_[u1]) {
                        // Force node2 to node1-value
                        unknownValue_[u2] = unknownValue_[u1] - value;
                        unknownForced_[u2] = true;
                    } else if (unknownForced_[u2]) {
                        // Force node1 to node2+value
                        unknownValue_[u1] = unknownValue_[u2] + value;
                        unknownForced_[u1] = true;
                    } else {
                        // Force u2 to 0 and u1 to value
                        unknownValue_[u1] = value;
                        unknownValue_[u2] = 0;
                        unknownForced_[u1] = true;
                        unknownForced_[u2] = true;
                    }
                } else {
                    // As delta forces
                    deltaValue_.push_back(value);
                    deltaIndices_.push_back(std::make_tuple(u1, u2));
                }        
            }
        }
    }

    return true;
}

void Forces::clear() {
    unknownValue_.clear();
    unknownForced_.clear();
    deltaValue_.clear();
    deltaIndices_.clear();
}

void Forces::resizeUnknownForces(size_t n) {
    unknownValue_.resize(n);
    unknownForced_.resize(n, false);
}

bool Forces::formatError(Status& s) const {
    std::string txt;
    switch (lastError) {
        case Error::ConflictNode:
            s.set(Status::Force, "Conflicting forces for node '"+std::string(errorNode1->name())+"'.");
            return false;
        case Error::ConflictDelta:
            s.set(Status::Force, "Forcing delta on node pair ('"
                        +std::string(errorNode1->name())+"', '"
                        +std::string(errorNode2->name())
                        +"') conflicts previous forces."
                    );
            return false;
        default:
            s.set(Status::OK, "");
            return true;
    }
}

void Forces::dump(Circuit& circuit, std::ostream& os) const {
    auto n = unknownForced_.size();
    for(decltype(n) i=1; i<n; i++) {
        if (!unknownForced_[i]) {
            continue;
        }
        os << circuit.reprNode(i)->name() << " : " << unknownValue_[i] << "\n";
    }
    auto nd = deltaIndices_.size();
    for(decltype(nd) i=0; i<nd; i++) {
        auto [u1, u2] = deltaIndices_[i];
        os << circuit.reprNode(u1)->name() << ", " << circuit.reprNode(u2)->name()
            << " : " << deltaValue_[i] << "\n"; 
    }
}

}
