#include "ansolution.h"
#include "circuit.h"
#include "common.h"


namespace NAMESPACE {

AnnotatedSolution::AnnotatedSolution() {
}

void AnnotatedSolution::set(Circuit& circuit, Vector<double>& vec) {
    names_.clear();
    values_.clear();
    auto n = vec.size();
    // Skip ground unknown
    for(decltype(n) i=1; i<n; i++) {
        auto node = circuit.reprNode(i);
        names_.push_back(node->name());
        values_.push_back(vec[i]);
    }
}


Forces::Forces() {
}

bool Forces::setForceOnUnknown(Circuit& circuit, Node* node, double value, Status& s) {
    // Unknown
    auto u = node->unknownIndex();
    // Is it a ground node
    if (u==0) {
        s.set(Status::NotFound, "Cannot force value on ground node '"+std::string(node->name())+"'.");
        return false;
    }
    // Is it conflicting with a previous nodeset
    if (unknownForced_[u] && unknownValue_[u]!=value) {
        s.set(Status::NotFound, "Conflicting forces for node '"+std::string(node->name())+"'.");
        return false;
    }
    unknownValue_[u] = value;
    unknownForced_[u] = true;

    return true;
}

bool Forces::set(Circuit& circuit, const AnnotatedSolution& solution, bool abortOnError, Status& s) {
    // Clear forced values
    unknownValue_.clear();
    unknownForced_.clear();
    deltaValue_.clear();
    deltaIndices_.clear();

    // Number of unknowns
    auto n = circuit.unknownCount();

    // Make space for variable forces
    unknownValue_.resize(n+1);
    unknownForced_.resize(n+1, false);

    bool error = false;

    // Go through all solution components, excluding ground
    auto nSol = solution.values().size();
    for(decltype(nSol) i=1; i<nSol; i++) {
        // Node
        auto name = solution.names()[i];
        auto value = solution.values()[i];
        Node* node = circuit.findNode(name);
        if (!node) {
            // Node not found
            continue;
        }

        if (!setForceOnUnknown(circuit, node, value, s)) {
            error = true;
            if (abortOnError) {
                return false;
            }
        }
    }

    return !error;
}

void Forces::clear() {
    unknownValue_.clear();
    unknownForced_.clear();
    deltaValue_.clear();
    deltaIndices_.clear();
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
            // Have single node
            nodes.push_back(node1);
            nodeIds.push_back(id1);
            nodeValues.push_back(value);
            state = 0;
        } else if (state==4) {
            // Have node pair
            nodePairs.push_back(std::make_tuple(node1, node2)); 
            nodeIdPairs.push_back(std::make_tuple(id1, id2));
            nodePairValues.push_back(value); 

            // Check existence of extradiagonal entries
            // Skip check if an entry has been found that is not there
            if (node1 && node2 && haveAllEntries) {
                auto u1 = node1->unknownIndex();
                auto u2 = node2->unknownIndex();
                auto [dummy1, found12] = circuit.sparsityMap().find(u1, u2);
                auto [dummy2, found21] = circuit.sparsityMap().find(u2, u1);
                haveAllEntries = haveAllEntries && found12 && found21;
            }
            state = 0;
        }
        nsNdx++;
    }
    
    return std::make_tuple(true, !haveAllEntries);
}

bool Forces::set(Circuit& circuit, const PreprocessedUserForces& preprocessed, bool uicMode, bool abortOnError, Status& s) {
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
        if (!node) {
            s.set(Status::NotFound, "Cannot set forced value, node '"+std::string(preprocessed.nodeIds[i])+"' not found.");
            error = true;
            if (abortOnError) {
                return false;
            }
            // Skip this node
            continue; 
        }
        
        if (!setForceOnUnknown(circuit, node, value, s)) {
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
        std::string txt;
        if (!node1) {
            txt = "Cannot set forced value, first";
        } else if (!node2) {
            txt = "Cannot set forced value, second";
        }
        if (!node1 || !node2) {
            s.set(Status::NotFound, txt+" node in node pair ('"+
                std::string(id1)+"', "+
                std::string(id2)+"') "+
                "not found."
            );
            error = true;
            if (abortOnError) {
                return false;
            }
            // Skip this node pair
            continue;
        }
        
        // Get unknowns and value
        auto u1 = node1->unknownIndex();
        auto u2 = node2->unknownIndex();
        auto value = preprocessed.nodePairValues[i];

        // Check if both nodes are ground
        if (u1==0 && u2==0 && value!=0) {
            s.set(Status::BadArguments, "Cannot set nonzero forced delta on node pair ('"+
                std::string(id1)+"', "+
                std::string(id2)+"')."
            ); 
            error = true;
            if (abortOnError) {
                return false;
            }
            // Skip this node pair
            continue;
        }

        // Check if first node is ground, convert it to a force on an unknown
        if (u1==0) {
            // v(0,x) = value -> v(x)=-value
            if (!setForceOnUnknown(circuit, node2, -value, s)) {
                error = true;
                if (abortOnError) {
                    return false;
                }
            }
        } else if (u2==0) {
            // v(x,0) = value -> v(x)=value
            if (!setForceOnUnknown(circuit, node1, value, s)) {
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
        if (!node1 || !node2) {
            // Skip this node pair
            continue;
        }
        
        // Get unknowns and value
        auto u1 = node1->unknownIndex();
        auto u2 = node2->unknownIndex();
        auto value = preprocessed.nodePairValues[i];

        // Check if both nodes are ground
        if (u1==0 && u2==0 && value!=0) {
            // Skip this node pair
            continue;
        }

        if (u1==0) {
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
                    s.set(Status::NotFound, "Forcing delta on node pair ('"
                        +std::string(node1->name())+"', '"
                        +std::string(node2->name())
                        +"'. conflicts previous forces."
                    );
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
