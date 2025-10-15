#include <limits>
#include <cstring>
#include "osdifile.h"

namespace NAMESPACE {

typedef struct Position {
    enum PositionType { None, Nature, DisciplineFlow, DisciplinePotential };
    PositionType kind {None};
    size_t pos {0};
} Position;

struct PositionHasher {
    size_t operator()(const sim::Position& p) const noexcept {
        return std::hash<int>()(p.kind) ^ (std::hash<int>()(p.pos) << 1);
    }
};

struct PositionEqual {
    bool operator()(const sim::Position& a, const sim::Position& b) const noexcept {
        return a.kind==b.kind && a.pos==b.pos;
    }
};

void OsdiFile::processNaturesAndDisciplines() {
    size_t i, j;
    OsdiAttribute* at;

    // Collect abstols from natures
    natureAbstol.resize(naturesCount);
    std::fill(natureAbstol.begin(), natureAbstol.end(), std::numeric_limits<double>::infinity());
    for(i=0; i<naturesCount; i++) {
        auto nature = natures + i;
        at = attributes + nature->attr_start;
        for(j=0; j<nature->num_attr; j++) {
            if (!strcmp(at[j].name, "abstol")) {
                if (at[j].value_type==ATTR_TYPE_REAL) {
                    natureAbstol[i] = at[j].value.real;
                    break;
                } else if (at[j].value_type==ATTR_TYPE_INT) {
                    natureAbstol[i] = at[j].value.integer;
                    break;
                }
            }
        }
    }
    
    // Collect abstols from disciplines
    disciplineFlowAbstol.resize(disciplinesCount);
    disciplinePotentialAbstol.resize(disciplinesCount);
    std::fill(disciplineFlowAbstol.begin(), disciplineFlowAbstol.end(), std::numeric_limits<double>::infinity());
    std::fill(disciplinePotentialAbstol.begin(), disciplinePotentialAbstol.end(), std::numeric_limits<double>::infinity());
    for(i=0; i<disciplinesCount; i++) {
        auto discipline = disciplines + i;
        at = attributes + discipline->attr_start;
        for(j=0; j<discipline->num_flow_attr; j++) {
            if (!strcmp(at[j].name, "abstol")) {
                if (at[j].value_type==ATTR_TYPE_REAL) {
                    disciplineFlowAbstol[i] = at[j].value.real;
                    break;
                } else if (at[j].value_type==ATTR_TYPE_INT) {
                    disciplineFlowAbstol[i] = at[j].value.integer;
                    break;
                }
            }
        }
        at += discipline->num_flow_attr;
        for(j=0; j<discipline->num_potential_attr; j++) {
            if (!strcmp(at[j].name, "abstol")) {
                if (at[j].value_type==ATTR_TYPE_REAL) {
                    disciplinePotentialAbstol[i] = at[j].value.real;
                    break;
                } else if (at[j].value_type==ATTR_TYPE_INT) {
                    disciplinePotentialAbstol[i] = at[j].value.integer;
                    break;
                }
            }
        }
    }

    // Parent position
    auto parentPos = [this](Position pos) -> Position {
        if (pos.kind==Position::Nature) {
            // At nature
            auto nat = natures + pos.pos;
            switch (nat->parent_type) {
                case NATREF_NATURE:
                    return Position { .kind=Position::Nature, .pos=nat->parent };
                case NATREF_DISCIPLINE_FLOW: 
                    return Position { .kind=Position::DisciplineFlow, .pos=nat->parent };
                case NATREF_DISCIPLINE_POTENTIAL: 
                    return Position { .kind=Position::DisciplinePotential, .pos=nat->parent };
                default:
                    return Position { .kind=Position::None }; 
            }
        } else if (pos.kind==Position::DisciplinePotential) {
            // At discipline potential, parent is potential nature
            auto disc = disciplines + pos.pos;
            if (disc->potential!=UINT32_MAX) {
                return Position { .kind=Position::Nature, .pos=disc->potential };
            } else {
                return Position { .kind=Position::None }; 
            }
        } else if (pos.kind==Position::DisciplineFlow) {
            // At discipline potential, parent is potential nature
            auto disc = disciplines + pos.pos;
            if (disc->flow!=UINT32_MAX) {
                return Position { .kind=Position::Nature, .pos=disc->flow };
            } else {
                return Position { .kind=Position::None }; 
            }
        } else {
            // No position
            return Position { .kind=Position::None }; 
        }
    };

    // Abstol at position
    auto posToAbstol = [this](Position pos) -> std::tuple<bool, double> {
        switch (pos.kind) {
            case Position::Nature: 
                if (!std::isinf(natureAbstol[pos.pos])) {
                    return std::make_tuple(true, natureAbstol[pos.pos]);
                }
                break;
            case Position::DisciplineFlow: 
                if (!std::isinf(disciplineFlowAbstol[pos.pos])) {
                    return std::make_tuple(true, disciplineFlowAbstol[pos.pos]);
                }
                break;
            case Position::DisciplinePotential: 
                if (!std::isinf(disciplinePotentialAbstol[pos.pos])) {
                    return std::make_tuple(true, disciplinePotentialAbstol[pos.pos]);
                }
                break;
        }
        return std::make_tuple(false, 0); 
    };

    // Traverse parent chain, return abstol if found
    auto chainTraverse = [this, &parentPos, &posToAbstol](Position start) -> std::tuple<bool, double> {
        // Prepare set of visited positions
        std::unordered_set<Position, PositionHasher, PositionEqual> visited;
        
        // Start
        Position pos = start;
        visited.insert(pos);
        while (true) {
            pos = parentPos(pos);
            if (pos.kind==Position::None) {
                // End of chain, stop
                break;
            }
            auto [_, inserted] = visited.insert(pos);
            if (!inserted) {
                // In a loop, stop
                break;
            }
            auto [given, abstol] = posToAbstol(pos);
            if (given) {
                // Found abstol
                return std::make_tuple(true, abstol);
            }
        }
        return std::make_tuple(false, 0);
    };

    // Collect abstol for natures via parents
    for(i=0; i<naturesCount; i++) {
        if (std::isinf(natureAbstol[i])) {
            // Traverse parents chain
            auto [found, abstol] = chainTraverse(Position {.kind=Position::Nature, .pos=i});
            if (found) {
                natureAbstol[i] = abstol;
            }
        }
    }

    // Collect abstols for discipline flows and potentials
    for(i=0; i<disciplinesCount; i++) {
        if (std::isinf(disciplineFlowAbstol[i])) {
            // Traverse parents chain
            auto [found, abstol] = chainTraverse(Position {.kind=Position::DisciplineFlow, .pos=i});
            if (found) {
                disciplineFlowAbstol[i] = abstol;
            }
        }
        if (std::isinf(disciplinePotentialAbstol[i])) {
            // Traverse parents chain
            auto [found, abstol] = chainTraverse(Position {.kind=Position::DisciplinePotential, .pos=i});
            if (found) {
                disciplinePotentialAbstol[i] = abstol;
            }
        }
    }

    // Collect abstol for nature idts
    natureIdtAbstol.resize(naturesCount);
    std::fill(natureIdtAbstol.begin(), natureIdtAbstol.end(), std::numeric_limits<double>::infinity());
    for(i=0; i<naturesCount; i++) {
        auto nat = natures + i;
        // Get idt nature
        if (nat->idt!=UINT32_MAX) {
            natureIdtAbstol[i] = natureAbstol[nat->idt];
        } else {
            // If idt is not given, the nature itself is the idt nature
            natureIdtAbstol[i] = natureAbstol[i];
        }
    }

    // Collect abstols for discipline flow and potential ddts
    disciplineFlowIdtAbstol.resize(naturesCount);
    std::fill(disciplineFlowIdtAbstol.begin(), disciplineFlowIdtAbstol.end(), std::numeric_limits<double>::infinity());
    disciplinePotentialIdtAbstol.resize(naturesCount);
    std::fill(disciplinePotentialIdtAbstol.begin(), disciplinePotentialIdtAbstol.end(), std::numeric_limits<double>::infinity());
    for(i=0; i<disciplinesCount; i++) {
        auto disc = disciplines + i;
        if (disc->flow!=UINT32_MAX) {
            // Have flow nature, get its idt abstol
            disciplineFlowIdtAbstol[i] = natureIdtAbstol[disc->flow];
        }
        if (disc->potential!=UINT32_MAX) {
            // Have potential nature, get its idt abstol
            disciplinePotentialIdtAbstol[i] = natureIdtAbstol[disc->potential];
        }
    }
}

}
