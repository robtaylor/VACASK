#include <limits>
#include <cstring>
#include <functional>
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

    // Traverse parent chain, return abstol if found
    // Traversal includes the starting point
    auto chainTraverse = [this, &parentPos](Position start, std::function<bool(struct Position)> checker) -> Position {
        // Prepare set of visited positions
        std::unordered_set<Position, PositionHasher, PositionEqual> visited;
        
        // Start
        Position pos = start;
        visited.insert(pos);
        while (true) {
            if (checker(pos)) {
                return pos;
            }
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
        }
        return Position { kind: Position::None };
    };

    auto hasAbstol = [this](struct Position pos) -> bool {
        switch (pos.kind) {
            case Position::Nature: 
                return !std::isinf(natureAbstol[pos.pos]);
            case Position::DisciplineFlow: 
                return !std::isinf(disciplineFlowAbstol[pos.pos]);
            case Position::DisciplinePotential: 
                return !std::isinf(disciplinePotentialAbstol[pos.pos]);
        }
        return false;
    };

    // Abstol at position
    auto getAbstol = [this](Position pos) -> double {
        switch (pos.kind) {
            case Position::Nature: 
                if (!std::isinf(natureAbstol[pos.pos])) {
                    return natureAbstol[pos.pos];
                }
                break;
            case Position::DisciplineFlow: 
                if (!std::isinf(disciplineFlowAbstol[pos.pos])) {
                    return disciplineFlowAbstol[pos.pos];
                }
                break;
            case Position::DisciplinePotential: 
                if (!std::isinf(disciplinePotentialAbstol[pos.pos])) {
                    return disciplinePotentialAbstol[pos.pos];
                }
                break;
        }
        return std::numeric_limits<double>::infinity(); 
    };

    auto hasIdt = [this](Position pos) -> bool {
        switch (pos.kind) {
            case Position::Nature: 
                return natures[pos.pos].idt!=UINT32_MAX;
            case Position::DisciplineFlow: {
                    if (auto ndx = disciplines[pos.pos].flow; ndx!=UINT32_MAX) {
                        return natures[ndx].idt!=UINT32_MAX;
                    } else {
                        return false;
                    }
                }
            case Position::DisciplinePotential: {
                if (auto ndx = disciplines[pos.pos].potential; ndx!=UINT32_MAX) {
                        return natures[ndx].idt!=UINT32_MAX;
                    } else {
                        return false;
                    }
                }
        }
        return false;
    };

    auto getIdt = [this](Position pos) -> uint32_t {
        switch (pos.kind) {
            case Position::Nature:
                return natures[pos.pos].idt;
            case Position::DisciplineFlow: 
                return natures[disciplines[pos.pos].flow].idt;
            case Position::DisciplinePotential: 
                return natures[disciplines[pos.pos].potential].idt;
            default:
                return UINT32_MAX;
        }
    };

    // Collect abstol for natures via parents
    for(i=0; i<naturesCount; i++) {
        if (std::isinf(natureAbstol[i])) {
            // Traverse parents chain
            auto pos = chainTraverse(Position {.kind=Position::Nature, .pos=i}, hasAbstol);
            if (pos.kind!=Position::None) {
                natureAbstol[i] = getAbstol(pos);
            }
        }
    }

    // Collect abstols for discipline flows and potentials
    for(i=0; i<disciplinesCount; i++) {
        if (std::isinf(disciplineFlowAbstol[i])) {
            // Traverse parents chain
            auto pos = chainTraverse(Position {.kind=Position::DisciplineFlow, .pos=i}, hasAbstol);
            if (pos.kind!=Position::None) {
                disciplineFlowAbstol[i] = getAbstol(pos);
            }
        }
        if (std::isinf(disciplinePotentialAbstol[i])) {
            // Traverse parents chain
            auto pos = chainTraverse(Position {.kind=Position::DisciplinePotential, .pos=i}, hasAbstol);
            if (pos.kind!=Position::None) {
                disciplinePotentialAbstol[i] = getAbstol(pos);
            }
        }
    }

    // Collect abstol for nature idts
    natureIdtId.resize(naturesCount);
    natureIdtAbstol.resize(naturesCount);
    std::fill(natureIdtId.begin(), natureIdtId.end(), NatureRegistry::noNature);
    std::fill(natureIdtAbstol.begin(), natureIdtAbstol.end(), std::numeric_limits<double>::infinity());
    for(i=0; i<naturesCount; i++) {
        // Find idt nature
        auto pos = chainTraverse(Position {.kind=Position::Nature, .pos=i}, hasIdt);
        auto idt = getIdt(pos);
        if (idt!=UINT32_MAX) {
            // Idt nature's abstol
            natureIdtId[i] = natureId[idt];
            natureIdtAbstol[i] = natureAbstol[idt];
        } else {
            // If idt is not given, the nature itself is the idt nature
            natureIdtId[i] = natureId[i];
            natureIdtAbstol[i] = natureAbstol[i];
        }
    }

    // At this point all nature idt abstols are determined, no need to traverse chain for idt nature
    // Collect idt abstols for discipline flows and potentials
    disciplineFlowIdtId.resize(naturesCount);
    disciplineFlowIdtAbstol.resize(naturesCount);
    std::fill(disciplineFlowIdtId.begin(), disciplineFlowIdtId.end(), NatureRegistry::noNature);
    std::fill(disciplineFlowIdtAbstol.begin(), disciplineFlowIdtAbstol.end(), std::numeric_limits<double>::infinity());
    disciplinePotentialIdtId.resize(naturesCount);
    disciplinePotentialIdtAbstol.resize(naturesCount);
    std::fill(disciplinePotentialIdtId.begin(), disciplinePotentialIdtId.end(), NatureRegistry::noNature);
    std::fill(disciplinePotentialIdtAbstol.begin(), disciplinePotentialIdtAbstol.end(), std::numeric_limits<double>::infinity());
    for(i=0; i<disciplinesCount; i++) {
        auto disc = disciplines + i;
        if (disc->flow!=UINT32_MAX) {
            // Have flow nature, get its idt abstol
            disciplineFlowIdtId[i] = natureIdtId[disc->flow];
            disciplineFlowIdtAbstol[i] = natureIdtAbstol[disc->flow];
        }
        if (disc->potential!=UINT32_MAX) {
            // Have potential nature, get its idt abstol
            disciplinePotentialIdtId[i] = natureIdtId[disc->potential];
            disciplinePotentialIdtAbstol[i] = natureIdtAbstol[disc->potential];
        }
    }
}

}
