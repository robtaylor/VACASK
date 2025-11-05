#include "parseroutput.h"
#include "common.h"


namespace NAMESPACE {

std::ostream& operator<<(std::ostream& os, const PTParsedIdentifier& obj) {
    os << std::string(obj.id_);
    return os;
}

std::ostream& operator<<(std::ostream& os, const PTIdentifierList& obj) {
    for(auto it=obj.cbegin(); it!=obj.cend(); ++it) {
        if (it!=obj.cbegin()) 
            os << " ";
        if (it->name()!=Id::none) {
            os << it->name();
        } else {
            os << "<null>";
        }
    }
    return os;
}


void PTParameterValue::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << std::string(id_) << "=" << val_;
}


void PTParameterExpression::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << std::string(id_) << "=" << rpn_.str();
}


void PTParameters::add(PTParameterValue&& v) {
    values_.push_back(std::move(v));
}

void PTParameters::set(std::vector<PTParameterValue>&& v) {
    values_ = std::move(v);
}

void PTParameters::add(PTParameterExpression&& e) {
    expressions_.push_back(std::move(e));
}

void PTParameters::add(PTParameters&& p) {
    for(auto it=p.values_.begin(); it!=p.values_.end(); ++it) {
        values_.push_back(std::move(*it));
    }
    for(auto it=p.expressions_.begin(); it!=p.expressions_.end(); ++it) {
        expressions_.push_back(std::move(*it));
    }
}

bool PTParameters::verify(Status& s) const {
    std::unordered_map<Id,Loc> puniq;
    for(auto& it : values_) {
        auto [itPrev, inserted] = puniq.insert({it.name(), it.location()});
        if (!inserted) {
            s.set(Status::Redefinition, "Parameter '"+std::string(it.name())+"' redefinition.");
            s.extend(it.location());
            if (itPrev->second) {
                s.extend("Parameter was first defined here");
                s.extend(itPrev->second);
            }
            return false;
        }
    }
    for(auto& it : expressions_) {
        auto [itPrev, inserted] = puniq.insert({it.name(), it.location()});
        if (!inserted) {
            s.set(Status::Redefinition, "Parameter '"+std::string(it.name())+"' redefinition.");
            s.extend(it.location());
            if (itPrev->second) {
                s.extend("Parameter was first defined here");
                s.extend(itPrev->second);
            }
            return false;
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const PTParameters& obj) {
    for(auto it=obj.values_.cbegin(); it!=obj.values_.cend(); ++it) {
        os << (it->name()) << "=" << it->val() << " " ;
    }
    for(auto it=obj.expressions_.begin(); it!=obj.expressions_.end(); ++it) {
        os << (it->name()) << "=" << it->rpn() << " " ;
    }
    
    return os;
}

void PTModel::add(PTParameters&& par) {
    parameters_.add(std::move(par));
}

void PTModel::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "model " << (modelName) << " " << (deviceName) << " ";
    os << parameters_ << "\n";
}


void PTInstance::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << (instanceName_) << " (" << connections_ << ") ";
    os << masterName_ << " " << parameters_ << "\n";
}


void PTBlockSequence::add(Rpn&& cond, PTBlock&& block, const Loc& l) {
    entries_.push_back(std::move(std::make_tuple(l, std::move(cond), std::move(block))));
}

void PTBlockSequence::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    bool first = true;
    for(auto& entry : entries_) {
        auto& [loc, rpn, block] = entry;
        if (rpn.size()>0) {
            if (first) {
                os << pfx << "@if " << rpn << "\n";
            } else {
                os << pfx << "@elseif " << rpn << "\n";
            }
        } else {
            os << pfx << "@else\n";
        }
        block.dump(indent+2, os);
        first = false;
    }
    os << pfx << "@end\n";
}


PTBlock& PTBlock::add(PTModel&& mod) {
    models_.push_back(std::move(mod));
    return *this;
}

PTBlock& PTBlock::add(PTInstance&& inst) {
    instances_.push_back(std::move(inst));
    return *this;
}

PTBlock& PTBlock::add(PTBlockSequence&& seq) {
    if (blockSequences_==nullptr) {
        blockSequences_ = std::make_unique<std::vector<PTBlockSequence>>();
    }
    blockSequences_->push_back(std::move(seq));
    return *this;
}

void PTBlock::dump(int indent, std::ostream& os) const {
    for(auto& mod : models_) {
        mod.dump(indent, os);
    }
    for(auto& inst : instances_) {
        inst.dump(indent, os);
    }
    if (hasBlockSequences()) {
        for(auto& seq: *blockSequences_) {
            seq.dump(indent, os);
        }
    }
}


void PTSubcircuitDefinition::add(PTIdentifierList&& terms) {
    terminals_ = std::move(terms);
}

void PTSubcircuitDefinition::add(PTSubcircuitDefinition&& subDef) {
    auto* ptr = new PTSubcircuitDefinition;
    *ptr = std::move(subDef);
    subDefs_.push_back(std::unique_ptr<PTSubcircuitDefinition>(ptr));
}

// Subcircuit terminal can have the same name as a global node. 
// In that case it is simply a local node name. It does not 
// behave as a global node. 
bool PTSubcircuitDefinition::verifyTerminals(Status& s) const {
    // Create terminal set for duplicates check
    std::unordered_set<Id> tset;
    for(auto& term : terminals_) {
        auto [exIt, inserted] = tset.insert(term.name());
        if (!inserted) {
            s.set(Status::Conflicting, "Terminal '"+std::string(term.name())+"' is not unique.");
            s.extend(location());
            return false;
        }
    }
    return true; 
}

void PTSubcircuitDefinition::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');

    bool isToplevel = !modelName;

    if (!isToplevel) {
        os << pfx << "subckt " << (modelName) << " (" << terminals_ << ")\n";
    }
    
    if (parameters_.count()>0) {
        os << pfx << (isToplevel ? "" : "  ") << "parameters " << parameters_ << "\n";
    }
    if (subDefs_.size()>0) {
        os << "\n";
    }
    for(auto it=subDefs_.begin(); it!=subDefs_.end(); ++it) {
        it->get()->dump(isToplevel ? indent : indent+2, os);
    }
    if (root_.models().size()>0) {
        os << "\n";
    }
    root_.dump(isToplevel ? indent : indent+2, os);
    
    if (!isToplevel) {
        os << pfx << "ends\n";
    }
    
    os << "\n";
}

void PTLoad::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "load \"" << file_ << "\"";
    os << " " << (parameters_);
    os << "\n";
}


std::ostream& operator<<(std::ostream& os, const PTSave& s) {
    os << std::string(s.typeName()) << "(";
    if (s.id[0]) {
        os << "\"" << std::string(s.id[0]) << "\"";
    }
    if (s.id[1]) {
        os << "," << "\"" << std::string(s.id[1]) << "\"";
    }
    os << ")";
    
    return os;
}


void PTSaves::add(PTSave&& s) {
    saves_.push_back(std::move(s));
}

void PTSaves::add(std::vector<PTSave>&& s) {
    for(size_t i=0; i<s.size(); i++) {
        saves_.push_back(std::move(s[i]));
    }
}

std::ostream& operator<<(std::ostream& os, const PTSaves& s) {
    if (s.saves_.size()>0) {
        for(auto it=s.saves_.cbegin(); it!=s.saves_.cend(); ++it) {
            os << *it << " ";
        }
    }
    
    return os;
}


std::ostream& operator<<(std::ostream& os, const PTSweep& s) {
    os << "sweep " << std::string(s.name_) << " " << s.parameters_;
    
    return os;
}


void PTSweeps::add(PTSweep&& s) {
    sweeps_.push_back(std::move(s));
}

void PTSweeps::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    for(auto it=data().cbegin(); it!=data().cend(); ++it) {
        os << pfx << *it << "\n";
    }
}


PTAnalysis& PTAnalysis::add(PTParameters&& par) {
    parameters_.add(std::move(par));
    return *this;
}

PTAnalysis& PTAnalysis::add(PTSweeps&& sw) {
    sweeps_ = std::move(sw);
    return *this;
}


void PTAnalysis::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    if (sweeps_.data().size()>0) {
        sweeps_.dump(indent, os);
    }
    os << pfx << (sweeps_.data().size()>0 ? "  " : "");
    os << "analysis " << std::string(name_) << " " << std::string(typeName_) << " ";
    os << parameters_ << "\n";
}


std::ostream& operator<<(std::ostream& os, const PTEmbed& e) {
    os << "embed \"" << e.filename_ << "\" <<<FILE\n" << e.contents_ << ">>>FILE";
    return os;
}


void PTCommand::add(PTIdentifierList&& kw) {
    keywords_ = std::move(kw);
}

void PTCommand::add(PTParameters&& args) {
    args_ = std::move(args);
}

void PTCommand::add(std::vector<Rpn>&& exprs) {
    expressions_ = std::move(exprs);
}

void PTCommand::add(PTSaves&& s) {
    saves_ = std::move(s);
}

std::ostream& operator<<(std::ostream& os, const PTCommand& c) {
    os << c.name_;
    if (c.keywords_.size()>0) {
        os << " " << c.keywords_;
    }
    if (c.saves_.saves().size()==0) {
        if (c.expressions_.size()>0) {
            os << " (";
            bool first = true;
            for(auto& it : c.expressions_) {
                if (!first) {
                    os << ", ";
                }
                os << it;
                first = false;
            }
            os << ")";
        }
        os << " " << c.args_;
    } else {
        os << " " << c.saves_;
    }
    return os;
}


ParserTables& ParserTables::setTitle(const std::string t) {
    title_ = t;
    return *this;
}

const std::string& ParserTables::title() const {
    return title_;
}

ParserTables& ParserTables::addDefaultSubDef(PTSubcircuitDefinition&& def) {
    defaultSubDef_ = std::move(def);
    return *this;
}

ParserTables& ParserTables::addLoad(PTLoad&& o) {
    loads_.push_back(std::move(o));
    return *this;
}

ParserTables& ParserTables::addGround(PTParsedIdentifier parsedId) {
    groundNodes_.push_back(parsedId);
    return *this;
}

ParserTables& ParserTables::addGlobal(PTParsedIdentifier parsedId) {
    globalNodes_.push_back(parsedId);
    return *this;
}

ParserTables& ParserTables::defaultGround() {
    if (groundNodes_.size()==0) {
        groundNodes_.push_back(PTParsedIdentifier(Id("0"), Loc::bad));
    }
    return *this;
}

void ParserTables::addCommand(PTAnalysis&& a) {
    control_.push_back(std::move(a));
}

void ParserTables::addCommand(PTCommand&& c) {
    control_.push_back(std::move(c));
}

void ParserTables::addEmbed(PTEmbed&& e) {
    embed_.push_back(std::move(e));
}

bool ParserTables::verify(Status& s) const {
    // Check duplicate analyses
    std::unordered_map<Id,const PTAnalysis*> anmap;
    for(auto& it : control_) {
        if (!std::holds_alternative<PTAnalysis>(it)) {
            continue;
        }
        auto& an = std::get<PTAnalysis>(it);
        auto [itEx, inserted] = anmap.insert({an.name(), &an});
        if (!inserted) {
            s.set(Status::Redefinition, "Analysis '"+std::string(an.name())+"' redefinition.");
            s.extend(an.location());
            s.extend("Analysis was first defined here");
            s.extend(itEx->second->location());
            return false;
        }
    }

    return true;
}

void ParserTables::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');

    os << pfx << title_ << "\n\n";

    if (loads_.size()>0) {
        for(auto it=loads_.begin(); it!=loads_.end(); ++it) {
            it->dump(indent, os);
        }
        os << "\n";
    }
    
    if (groundNodes_.size()>0) {
        os << pfx << "ground";
        for(auto it=groundNodes_.begin(); it!=groundNodes_.end(); ++it) {
            os << " " << *it;
        }
        os << "\n\n";
    }
    
    if (globalNodes_.size()>0) {
        os << pfx << "global";
        for(auto it=globalNodes_.begin(); it!=globalNodes_.end(); ++it) {
            os << " " << *it;
        }
        os << "\n\n";
    }

    defaultSubDef_.dump(indent, os);
    os << "\n";

    if (embed_.size()>0) {
        for(auto& it : embed_) {
            os << pfx << it << "\n";
        }
        os << "\n";
    }

    if (control_.size()>0) {
        os << pfx << "control\n";
        for(auto& it : control_) {
            if (std::holds_alternative<PTAnalysis>(it)) {
                std::get<PTAnalysis>(it).dump(indent+2, os); 
            } else {
                os << pfx << "  " << std::get<PTCommand>(it) << "\n";
            }
        }
        os << pfx << "endc\n";
    }
}

}
