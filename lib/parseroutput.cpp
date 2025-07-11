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


void PTParameterValue::dump(int indent, std::ostream& os) {
    std::string pfx = std::string(indent, ' ');
    os << pfx << std::string(id_) << "=" << val_;
}


void PTParameterExpression::dump(int indent, std::ostream& os) {
    std::string pfx = std::string(indent, ' ');
    os << pfx << std::string(id_) << "=" << rpn_.str();
}


PTParameters::PTParameters() {
}

PTParameters::PTParameters(std::vector<PTParameterValue>&& pv, std::vector<PTParameterExpression>&& pe) 
    : values_(std::move(pv)), expressions_(std::move(pe)) {
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

PTModel::PTModel() {
}

PTModel::PTModel(const Loc& l, Id name, Id device) 
    : loc(l), modelName(name), deviceName(device) {
}

void PTModel::add(PTParameters&& par) {
    parameters_.add(std::move(par));
}

void PTModel::dump(int indent, std::ostream& os) {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "model " << (modelName) << " " << (deviceName) << " ";
    os << parameters_ << "\n";
}


PTInstance::PTInstance() {
}

PTInstance::PTInstance(const Loc& l, Id name, Id master, PTIdentifierList&& conns, PTParameters&& params) 
    : loc(l), instanceName_(name), masterName_(master), connections_(std::move(conns)), parameters_(std::move(params)) {
}

void PTInstance::dump(int indent, std::ostream& os) {
    std::string pfx = std::string(indent, ' ');
    os << pfx << (instanceName_) << " (" << connections_ << ") ";
    os << masterName_ << " " << parameters_ << "\n";
}


PTBlockSequence::PTBlockSequence() {
}

void PTBlockSequence::add(const Loc& l, Rpn&& cond, PTBlockIndex blockIndex) {
    entries_.push_back(std::move(std::make_tuple(l, std::move(cond), blockIndex)));
}


PTBlock::PTBlock() {
}

void PTBlock::add(PTModel&& mod) {
    models_.push_back(std::move(mod));
}

void PTBlock::add(PTInstance&& inst) {
    instances_.push_back(std::move(inst));
}

void PTBlock::add(PTBlockSequence&& seq) {
    blockSequences_.push_back(std::move(seq));
}


PTSubcircuitDefinition::PTSubcircuitDefinition() {
    // Add root block
    add(std::move(PTBlock()));
}

PTSubcircuitDefinition::PTSubcircuitDefinition(const Loc& l, Id name, PTIdentifierList&& terms) 
    : PTModel(l, name, "__hierarchical__"), terminals_(std::move(terms)) {
    // Add root block
    add(std::move(PTBlock()));
}

void PTSubcircuitDefinition::add(PTIdentifierList&& terms) {
    terminals_ = std::move(terms);
}

PTBlockIndex PTSubcircuitDefinition::add(PTBlock&& block) {
    auto ndx = blocks_.size();
    blocks_.push_back(std::move(block));
    return ndx;
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

void PTSubcircuitDefinition::dump(int indent, std::ostream& os) {
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
    if (block(0).models().size()>0) {
        os << "\n";
    }
    for(auto it=block(0).models().begin(); it!=block(0).models().end(); ++it) {
        it->dump(isToplevel ? indent : indent+2, os);
    }
    if (block(0).instances().size()>0) {
        os << "\n";
    }
    for(auto it=block(0).instances().begin(); it!=block(0).instances().end(); ++it) {
        it->dump(isToplevel ? indent : indent+2, os);
    }
    
    if (!isToplevel) {
        os << pfx << "ends\n";
    }
    
    os << "\n";
}

PTLoad::PTLoad() {
}

PTLoad::PTLoad(const Loc& l, const std::string& file, Id module, Id asModule)
    : loc(l), file_(file), module_(module), asModule_(asModule) {
}

void PTLoad::dump(int indent, std::ostream& os) {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "load \"" << file_ << "\"";
    if (module_) {
        os << " " << (module_);
    }
    if (asModule_) {
        os << " : " << (asModule_);
    }
    
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


PTSaves::PTSaves() {
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


PTSweep::PTSweep(const Loc& l, Id name, PTParameters&& par) : loc(l), name_(name), parameters_(std::move(par)) {
}

std::ostream& operator<<(std::ostream& os, const PTSweep& s) {
    os << "sweep " << std::string(s.name_) << " " << s.parameters_;
    
    return os;
}


PTSweeps::PTSweeps() {
}

void PTSweeps::add(PTSweep&& s) {
    sweeps_.push_back(std::move(s));
}

void PTSweeps::dump(int indent, std::ostream& os) {
    std::string pfx = std::string(indent, ' ');
    for(auto it=data().cbegin(); it!=data().cend(); ++it) {
        os << pfx << *it << "\n";
    }
}


PTAnalysis::PTAnalysis() {
}

PTAnalysis::PTAnalysis(const Loc& l, Id name, Id typeName) 
    : loc(l), name_(name), typeName_(typeName) {
}

void PTAnalysis::add(PTParameters&& par) {
    parameters_.add(std::move(par));
}

void PTAnalysis::add(PTSweeps&& sw) {
    sweeps_ = std::move(sw);
}


void PTAnalysis::dump(int indent, std::ostream& os) {
    std::string pfx = std::string(indent, ' ');
    if (sweeps_.data().size()>0) {
        sweeps_.dump(indent, os);
    }
    os << pfx << (sweeps_.data().size()>0 ? "  " : "");
    os << "analysis " << std::string(name_) << " " << std::string(typeName_) << " ";
    os << parameters_ << "\n";
}


ParserTables::ParserTables() {
}

ParserTables::~ParserTables() {
}

void ParserTables::setTitle(const std::string t) {
    title_ = t;
}

const std::string& ParserTables::title() const {
    return title_;
}

void ParserTables::addDefaultSubDef(PTSubcircuitDefinition&& def) {
    defaultSubDef_ = std::move(def);
}

void ParserTables::addLoad(PTLoad&& o) {
    loads_.push_back(std::move(o));
}

void ParserTables::addGround(PTParsedIdentifier parsedId) {
    groundNodes_.push_back(parsedId);
}

void ParserTables::addGlobal(PTParsedIdentifier parsedId) {
    globalNodes_.push_back(parsedId);
}

void ParserTables::defaultGround() {
    if (groundNodes_.size()==0) {
        groundNodes_.push_back(PTParsedIdentifier(Id("0"), Loc::bad));
    }
}

bool ParserTables::verify(Status& s) const {
    return true;
}

void ParserTables::dump(int indent, std::ostream& os) {
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
}

}
