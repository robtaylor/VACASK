#ifndef __PARSER_OUTPUT
#define __PARSER_OUTPUT

#include <iostream>
#include <string>
#include <unordered_set>
#include <memory>
#include "value.h"
#include "rpnexpr.h"
#include "filestack.h"
#include "identifier.h"
#include "acct.h"
#include "common.h"


namespace NAMESPACE {

// Parsed identifier 
class PTParsedIdentifier {
public:
    PTParsedIdentifier(const char* name, Loc location=Loc::bad) : id_(name), loc_(location) {};
    PTParsedIdentifier(Id name, Loc location=Loc::bad) : id_(name), loc_(location) {};

    // This is a simple type, allow copy constructor
    PTParsedIdentifier           (const PTParsedIdentifier&)  = default;
    PTParsedIdentifier           (      PTParsedIdentifier&&) = default;
    PTParsedIdentifier& operator=(const PTParsedIdentifier&)  = default;
    PTParsedIdentifier& operator=(      PTParsedIdentifier&&) = default;

    // Getters
    Id name() const { return id_; }; 
    Loc location() const { return loc_; };

    friend std::ostream& operator<<(std::ostream& os, const PTParsedIdentifier& obj);

private:
    Id id_;
    Loc loc_;
};


// Identifier list, when used for terminal connections Id::none identifier
// means that the corresponding terminal is not connected
typedef std::vector<PTParsedIdentifier> PTIdentifierList;

std::ostream& operator<<(std::ostream& os, const PTIdentifierList& obj); 


// A single parameter with name identifier and given Value
class PTParameterValue {
public:
    PTParameterValue(const Id ident, Value&& value, const Loc& l=Loc::bad) : loc_(l), id_(ident), val_(std::move(value)) {};

    PTParameterValue           (const PTParameterValue&)  = delete;
    PTParameterValue           (      PTParameterValue&&) = default;
    PTParameterValue& operator=(const PTParameterValue&)  = delete;
    PTParameterValue& operator=(      PTParameterValue&&) = default;

    // Getters
    Id name() const { return id_; }; 
    Loc location() const { return loc_; };
    const Value& val() const { return val_; };

    void dump(int indent, std::ostream& os) const;

private:
    Id id_;
    Value val_;
    Loc loc_;
};

// A single parameter with name identifier and given Expression
class PTParameterExpression {
public:
    PTParameterExpression(const Id ident, Rpn&& expr, const Loc& l=Loc::bad) : loc_(l), id_(ident), rpn_(std::move(expr)) {};

    PTParameterExpression           (const PTParameterExpression&)  = delete;
    PTParameterExpression           (      PTParameterExpression&&) = default;
    PTParameterExpression& operator=(const PTParameterExpression&)  = delete;
    PTParameterExpression& operator=(      PTParameterExpression&&) = default;

    // Getters
    Id name() const { return id_; }; 
    Loc location() const { return loc_; };
    const Rpn& rpn() const { return rpn_; };

    void dump(int indent, std::ostream& os) const;

private:
    Id id_;
    Rpn rpn_;
    Loc loc_;
};


// Parameters (constants and expressions) for instances, models, and subcircuit definitions
class PTParameters {
public:
    PTParameters() {};
    PTParameters(std::vector<PTParameterValue>&& pv)
        : values_(std::move(pv)) {};
    PTParameters(std::vector<PTParameterExpression>&& pe)
        : expressions_(std::move(pe)) {};
    PTParameters(std::vector<PTParameterValue>&& pv, std::vector<PTParameterExpression>&& pe)
        : values_(std::move(pv)), expressions_(std::move(pe)) {};
    
    PTParameters           (const PTParameters&)  = delete;
    PTParameters           (      PTParameters&&) = default;
    PTParameters& operator=(const PTParameters&)  = delete;
    PTParameters& operator=(      PTParameters&&) = default;

    // Getters
    inline size_t valueCount() const { return values_.size(); };
    inline size_t expressionCount() const { return expressions_.size(); };
    inline size_t count() const { return values_.size()+expressions_.size(); };
    inline const auto& values() const { return values_; };
    inline auto& values() { return values_; };
    inline const auto& expressions() const { return expressions_; };
    inline auto& expressions() { return expressions_; };

    // Fluent API
    PTParameters& add(PTParameterValue&& v) & { values_.push_back(std::move(v)); return *this; };
    PTParameters& add(PTParameterExpression&& e) & { expressions_.push_back(std::move(e)); return *this; };
    PTParameters& add(PTParameters&& p) & {
        for(auto it=p.values_.begin(); it!=p.values_.end(); ++it) {
            values_.push_back(std::move(*it));
        }
        for(auto it=p.expressions_.begin(); it!=p.expressions_.end(); ++it) {
            expressions_.push_back(std::move(*it));
        }
        return *this;
    };
    PTParameters&& add(PTParameterValue&& v) && { return std::move(this->add(std::move(v))); };
    PTParameters&& add(PTParameterExpression&& e) && { return std::move(this->add(std::move(e))); };
    PTParameters&& add(PTParameters&& p) && { return std::move(this->add(std::move(p))); };

    bool verify(Status& s=Status::ignore) const;

    // Operator that prints parameters
    friend std::ostream& operator<<(std::ostream& os, const PTParameters& obj);

private:
    std::vector<PTParameterValue> values_;
    std::vector<PTParameterExpression> expressions_;
};


using PTValueOrExpression = std::variant<const PTParameterValue*, const PTParameterExpression*>;
using PTParameterMap = std::unordered_map<Id, PTValueOrExpression>;


class PTModel {
public:
    PTModel() {};
    PTModel(Id name, Id device, const Loc& l=Loc::bad)
        : loc(l), modelName(name), deviceName(device) {};
    PTModel(Id name, Id device, PTParameters&& params, const Loc& l=Loc::bad)
        : loc(l), modelName(name), deviceName(device), parameters_(std::move(params)) {};

    PTModel           (const PTModel&)  = delete;
    PTModel           (      PTModel&&) = default;
    PTModel& operator=(const PTModel&)  = delete;
    PTModel& operator=(      PTModel&&) = default;

    // Getters
    inline const Loc& location() const { return loc; };
    inline Id name() const { return modelName; };
    inline Id device() const { return deviceName; };
    inline bool isParameterized() const { return parameters_.expressionCount()>0; };
    inline const PTParameters& parameters() const { return parameters_; };

    // Fluent API
    PTModel& add(PTParameters&& par) & { parameters_.add(std::move(par)); return *this; };
    PTModel& add(PTParameterValue&& v) & { parameters_.add(std::move(v)); return *this; };
    PTModel& add(PTParameterExpression&& e) & { parameters_.add(std::move(e)); return *this; };
    PTModel&& add(PTParameters&& par) && { return std::move(this->add(std::move(par))); };
    PTModel&& add(PTParameterValue&& v) && { return std::move(this->add(std::move(v))); };
    PTModel&& add(PTParameterExpression&& e) && { return std::move(this->add(std::move(e))); };
    
    void dump(int indent, std::ostream& os) const;

protected:
    Id modelName;
    Id deviceName;
    PTParameters parameters_;
    Loc loc;
};


class PTInstance {
public:
    PTInstance() {};
    PTInstance(Id name, Id master, PTIdentifierList&& terms, const Loc& l=Loc::bad)
        : loc(l), instanceName_(name), masterName_(master), connections_(std::move(terms)) {};
    PTInstance(Id name, Id master, PTIdentifierList&& terms, PTParameters&& params, const Loc& l=Loc::bad)
        : loc(l), instanceName_(name), masterName_(master), connections_(std::move(terms)), parameters_(std::move(params)) {};

    PTInstance           (const PTInstance&)  = delete;
    PTInstance           (      PTInstance&&) = default;
    PTInstance& operator=(const PTInstance&)  = delete;
    PTInstance& operator=(      PTInstance&&) = default;

    // Getters
    inline const Loc& location() const { return loc; };
    Id name() const { return instanceName_; };
    bool isParameterized() const { return parameters_.expressionCount()>0; };
    inline Id masterName() const { return masterName_; };
    inline const PTParameters& parameters() const { return parameters_; };
    inline const PTIdentifierList& connections() const { return connections_; };

    // Fluent API
    PTInstance& add(PTParameters&& par) & { parameters_.add(std::move(par)); return *this; };
    PTInstance& add(PTParameterValue&& v) & { parameters_.add(std::move(v)); return *this; };
    PTInstance& add(PTParameterExpression&& e) & { parameters_.add(std::move(e)); return *this; };
    PTInstance&& add(PTParameters&& par) && { return std::move(this->add(std::move(par))); };
    PTInstance&& add(PTParameterValue&& v) && { return std::move(this->add(std::move(v))); };
    PTInstance&& add(PTParameterExpression&& e) && { return std::move(this->add(std::move(e))); };

    void dump(int indent, std::ostream& os) const;

private:
    Id instanceName_;
    Id masterName_;
    PTIdentifierList connections_;
    PTParameters parameters_;
    Loc loc;
};


// Block index
typedef uint32_t PTBlockIndex; 

class PTBlockSequence;

// A single netlist block comprising models, instances, and block sequences
class PTBlock {
public:
    PTBlock() {};
    
    PTBlock           (const PTBlock&)  = delete;
    PTBlock           (      PTBlock&&) = default;
    PTBlock& operator=(const PTBlock&)  = delete;
    PTBlock& operator=(      PTBlock&&) = default;

    // Getters
    bool hasBlockSequences() const { return blockSequences_!=nullptr; };
    inline const std::vector<PTModel>& models() const { return models_; };
    inline const std::vector<PTInstance>& instances() const { return instances_; };
    inline const std::vector<PTBlockSequence>& blockSequences() const { return *blockSequences_; };
    inline std::vector<PTModel>& models() { return models_; };
    inline std::vector<PTInstance>& instances() { return instances_; };
    inline std::vector<PTBlockSequence>& blockSequences() { return *blockSequences_; };

    // Fluent API
    PTBlock& add(PTModel&& mod) & { models_.push_back(std::move(mod)); return *this; };
    PTBlock& add(PTInstance&& inst) & { instances_.push_back(std::move(inst)); return *this; };
    PTBlock& add(PTBlockSequence&& seq) & {
        if (blockSequences_==nullptr) {
            blockSequences_ = std::make_unique<std::vector<PTBlockSequence>>();
        }
        blockSequences_->push_back(std::move(seq));
        return *this;
    };
    PTBlock&& add(PTModel&& mod) && { return std::move(this->add(std::move(mod))); };
    PTBlock&& add(PTInstance&& inst) && { return std::move(this->add(std::move(inst))); };
    PTBlock&& add(PTBlockSequence&& seq) && { return std::move(this->add(std::move(seq))); };

    void dump(int indent, std::ostream& os) const;

private:
    std::vector<PTModel> models_;
    std::vector<PTInstance> instances_;
    std::unique_ptr<std::vector<PTBlockSequence>> blockSequences_;
    Loc loc;
};


typedef std::tuple<Loc, Rpn, PTBlock> PTBlockSequenceEntry;

class PTBlockSequence {
public:
    PTBlockSequence() {};
    
    PTBlockSequence           (const PTBlockSequence&)  = delete;
    PTBlockSequence           (      PTBlockSequence&&) = default;
    PTBlockSequence& operator=(const PTBlockSequence&)  = delete;
    PTBlockSequence& operator=(      PTBlockSequence&&) = default;

    // Getters
    const std::vector<PTBlockSequenceEntry>& entries() const { return entries_; };
    PTBlock& back() { return std::get<2>(entries_.back()); }; 

    // Fluent API
    PTBlockSequence& add(Rpn&& cond, PTBlock&& block, const Loc& l=Loc::bad) & {
        entries_.push_back(std::move(std::make_tuple(l, std::move(cond), std::move(block))));
        return *this;
    };
    PTBlockSequence&& add(Rpn&& cond, PTBlock&& block, const Loc& l=Loc::bad) && {
        return std::move(this->add(std::move(cond), std::move(block), l));
    };

    void dump(int indent, std::ostream& os) const;

private:
    std::vector<PTBlockSequenceEntry> entries_;
};


class PTSubcircuitDefinition : public PTModel {
public:
    PTSubcircuitDefinition() {};
    PTSubcircuitDefinition(Id name, PTIdentifierList&& terms, const Loc& l=Loc::bad)
        : PTModel(name, "__hierarchical__", l), terminals_(std::move(terms)) {};

    PTSubcircuitDefinition           (const PTSubcircuitDefinition&)  = delete;
    PTSubcircuitDefinition           (      PTSubcircuitDefinition&&) = default;
    PTSubcircuitDefinition& operator=(const PTSubcircuitDefinition&)  = delete;
    PTSubcircuitDefinition& operator=(      PTSubcircuitDefinition&&) = default;

    // Getters
    inline const PTIdentifierList& terminals() const { return terminals_; };
    inline const PTBlock& root() const { return root_; };
    inline const std::vector<std::unique_ptr<PTSubcircuitDefinition>>& subDefs() const { return subDefs_; };
    
    // Fluent API
    PTSubcircuitDefinition& add(PTSubcircuitDefinition&& subDef) & {
        auto* ptr = new PTSubcircuitDefinition;
        *ptr = std::move(subDef);
        subDefs_.push_back(std::unique_ptr<PTSubcircuitDefinition>(ptr));
        return *this;
    };
    PTSubcircuitDefinition& add(PTModel&& mod) & { root_.add(std::move(mod)); return *this; };
    PTSubcircuitDefinition& add(PTInstance&& inst) & { root_.add(std::move(inst)); return *this; };
    PTSubcircuitDefinition& add(PTBlockSequence&& seq) & { root_.add(std::move(seq)); return *this; };
    PTSubcircuitDefinition& add(PTParameters&& par) & { PTModel::add(std::move(par)); return *this; };
    PTSubcircuitDefinition& add(PTParameterValue&& v) & { PTModel::add(std::move(v)); return *this; };
    PTSubcircuitDefinition& add(PTParameterExpression&& e) & { PTModel::add(std::move(e)); return *this; };
    PTSubcircuitDefinition&& add(PTSubcircuitDefinition&& subDef) && { return std::move(this->add(std::move(subDef))); };
    PTSubcircuitDefinition&& add(PTModel&& mod) && { return std::move(this->add(std::move(mod))); };
    PTSubcircuitDefinition&& add(PTInstance&& inst) && { return std::move(this->add(std::move(inst))); };
    PTSubcircuitDefinition&& add(PTBlockSequence&& seq) && { return std::move(this->add(std::move(seq))); };
    PTSubcircuitDefinition&& add(PTParameters&& par) && { return std::move(this->add(std::move(par))); };
    PTSubcircuitDefinition&& add(PTParameterValue&& v) && { return std::move(this->add(std::move(v))); };
    PTSubcircuitDefinition&& add(PTParameterExpression&& e) && { return std::move(this->add(std::move(e))); };

    void dump(int indent, std::ostream& os) const;

    // Checker
    bool verifyTerminals(Status& s=Status::ignore) const;
    
private:
    PTIdentifierList terminals_;
    PTBlock root_;
    std::vector<std::unique_ptr<PTSubcircuitDefinition>> subDefs_;
};


class PTLoad {
public:
    PTLoad() {};
    PTLoad(const std::string& file, const Loc& l=Loc::bad)
        : loc(l), file_(file) {};
    PTLoad(const std::string& file, PTParameters&& par, const Loc& l=Loc::bad)
        : loc(l), file_(file), parameters_(std::move(par)) {};
    
    PTLoad           (const PTLoad&)  = delete;
    PTLoad           (      PTLoad&&) = default;
    PTLoad& operator=(const PTLoad&)  = delete;
    PTLoad& operator=(      PTLoad&&) = default;

    // Getters
    inline const Loc& location() const { return loc; };
    inline const std::string& file() const { return file_; };
    inline const PTParameters& parameters() const { return parameters_; };

    // Fluent API
    PTLoad& add(PTParameters&& par) & { parameters_.add(std::move(par)); return *this; };
    PTLoad& add(PTParameterValue&& v) & { parameters_.add(std::move(v)); return *this; };
    PTLoad& add(PTParameterExpression&& e) & { parameters_.add(std::move(e)); return *this; };
    PTLoad&& add(PTParameters&& par) && { return std::move(this->add(std::move(par))); };
    PTLoad&& add(PTParameterValue&& v) && { return std::move(this->add(std::move(v))); };
    PTLoad&& add(PTParameterExpression&& e) && { return std::move(this->add(std::move(e))); };
    
    void dump(int indent, std::ostream& os) const;

private:
    PTParameters parameters_;
    std::string file_;
    Loc loc;
};


class PTSave {
public:
    PTSave() {};
    PTSave(const Id tname, const Loc& l=Loc::bad) : loc(l), typeName_(tname) { id[0] = Id(); id[1] = Id(); };
    PTSave(const Id tname, const Id id1, const Loc& l=Loc::bad) : loc(l), typeName_(tname) { id[0] = id1; id[1] = Id(); };
    PTSave(const Id tname, const Id id1, const Id id2, const Loc& l=Loc::bad) : loc(l), typeName_(tname) { id[0] = id1; id[1] = id2; };

    // Allow copying - it is cheap and simplifies the API
    PTSave           (const PTSave&)  = default;
    PTSave           (      PTSave&&) = default;
    PTSave& operator=(const PTSave&)  = default;
    PTSave& operator=(      PTSave&&) = default;

    // Getters
    Id typeName() const { return typeName_; };
    Id objName() const { return id[0]; };
    Id subName() const { return id[1]; }; 
    const Loc& location() const { return loc; };

    // Operator that prints a save directive
    friend std::ostream& operator<<(std::ostream& os, const PTSave& o);

private:
    Id typeName_;
    Id id[2];
    Loc loc;
};


class PTSaves {
public:
    PTSaves() {};
    
    PTSaves           (const PTSaves&)  = delete;
    PTSaves           (      PTSaves&&) = default;
    PTSaves& operator=(const PTSaves&)  = delete;
    PTSaves& operator=(      PTSaves&&) = default;

    // Getters
    const std::vector<PTSave>& saves() const { return saves_; }; 

    // Fluent API
    PTSaves& add(PTSave&& s) & { saves_.push_back(std::move(s)); return *this; };
    PTSaves& add(std::vector<PTSave>&& s) & {
        for(size_t i=0; i<s.size(); i++) {
            saves_.push_back(std::move(s[i]));
        }
        return *this;
    };
    PTSaves&& add(PTSave&& s) && { return std::move(this->add(std::move(s))); };
    PTSaves&& add(std::vector<PTSave>&& s) && { return std::move(this->add(std::move(s))); };
    
    // Operator that prints save directives
    friend std::ostream& operator<<(std::ostream& os, const PTSaves& o);

private:
    std::vector<PTSave> saves_;
};

using PTSavesVector = std::vector<PTSave>;


class PTSweep {
public:
    PTSweep(Id name, const Loc& l=Loc::bad) : loc(l), name_(name) {};
    PTSweep(Id name, PTParameters&& par, const Loc& l=Loc::bad)
         : loc(l), name_(name), parameters_(std::move(par)) {};

    PTSweep           (const PTSweep&)  = delete;
    PTSweep           (      PTSweep&&) = default;
    PTSweep& operator=(const PTSweep&)  = delete;
    PTSweep& operator=(      PTSweep&&) = default;

    // Getters
    Id name() const { return name_; }; 
    Loc location() const { return loc; }; 
    const PTParameters& parameters() const { return parameters_; };

    // Fluent API
    PTSweep& add(PTParameters&& par) & { parameters_.add(std::move(par)); return *this; };
    PTSweep& add(PTParameterValue&& v) & { parameters_.add(std::move(v)); return *this; };
    PTSweep& add(PTParameterExpression&& e) & { parameters_.add(std::move(e)); return *this; };
    PTSweep&& add(PTParameters&& par) && { return std::move(this->add(std::move(par))); };
    PTSweep&& add(PTParameterValue&& v) && { return std::move(this->add(std::move(v))); };
    PTSweep&& add(PTParameterExpression&& e) && { return std::move(this->add(std::move(e))); };
    
    // Operator that prints a sweep directive
    friend std::ostream& operator<<(std::ostream& os, const PTSweep& o);

private:
    Id name_;
    Loc loc;
    PTParameters parameters_;
};


class PTSweeps {
public:
    PTSweeps() {};

    PTSweeps           (const PTSweeps&)  = delete;
    PTSweeps           (      PTSweeps&&) = default;
    PTSweeps& operator=(const PTSweeps&)  = delete;
    PTSweeps& operator=(      PTSweeps&&) = default;

    // Getters
    std::vector<PTSweep>& sweeps() { return sweeps_; }; 
    const std::vector<PTSweep>& sweeps() const { return sweeps_; }; 
    
    // Fluent API
    PTSweeps& add(PTSweep&& s) & { sweeps_.push_back(std::move(s)); return *this; };
    PTSweeps&& add(PTSweep&& s) && { return std::move(this->add(std::move(s))); };
    
    void dump(int indent, std::ostream& os) const;

private:
    std::vector<PTSweep> sweeps_;
};


class PTAnalysis {
public:
    PTAnalysis() {};
    PTAnalysis(Id name, Id typeName, const Loc& l=Loc::bad)
        : loc(l), name_(name), typeName_(typeName) {};

    PTAnalysis           (const PTAnalysis&)  = delete;
    PTAnalysis           (      PTAnalysis&&) = default;
    PTAnalysis& operator=(const PTAnalysis&)  = delete;
    PTAnalysis& operator=(      PTAnalysis&&) = default;

    // Getters
    inline const Loc& location() const { return loc; };
    inline Id name() const { return name_; };
    inline Id typeName() const { return typeName_; };
    const PTParameters& parameters() const { return parameters_; };
    PTParameters& parameters() { return parameters_; };
    const std::vector<PTSweep>& sweeps() const { return sweeps_; };

    // Fluent API
    PTAnalysis& add(PTSweep&& sw) & { sweeps_.push_back(std::move(sw)); return *this; };
    PTAnalysis& add(PTSweeps&& swps) & { 
        for(size_t i=0; i<swps.sweeps().size(); i++) {
            sweeps_.push_back(std::move(swps.sweeps()[i]));
        }
        return *this; 
    };
    PTAnalysis& add(PTParameters&& par) & { parameters_.add(std::move(par)); return *this; };
    PTAnalysis& add(PTParameterValue&& v) & { parameters_.add(std::move(v)); return *this; };
    PTAnalysis& add(PTParameterExpression&& e) & { parameters_.add(std::move(e)); return *this; };
    PTAnalysis&& add(PTSweep&& sw) && { return std::move(this->add(std::move(sw))); };
    PTAnalysis&& add(PTSweeps&& swps) && { return std::move(this->add(std::move(swps))); };
    PTAnalysis&& add(PTParameters&& par) && { return std::move(this->add(std::move(par))); };
    PTAnalysis&& add(PTParameterValue&& v) && { return std::move(this->add(std::move(v))); };
    PTAnalysis&& add(PTParameterExpression&& e) && { return std::move(this->add(std::move(e))); };
    
    void dump(int indent, std::ostream& os) const;

private:
    Id name_;
    Id typeName_;
    PTParameters parameters_;
    std::vector<PTSweep> sweeps_; 
    Loc loc;
};

// Embedded files
class PTEmbed {
public:
    PTEmbed() {};
    PTEmbed(std::string&& filename, std::string&& contents, const Loc& l=Loc::bad)
        : loc_(l), filename_(std::move(filename)), contents_(std::move(contents)) {};
    
    PTEmbed           (const PTEmbed&)  = delete;
    PTEmbed           (      PTEmbed&&) = default;
    PTEmbed& operator=(const PTEmbed&)  = delete;
    PTEmbed& operator=(      PTEmbed&&) = default;

    // Getters
    Loc location() const { return loc_; };
    const std::string& filename() const { return filename_; }; 
    const std::string& contents() const { return contents_; }; 

    friend std::ostream& operator<<(std::ostream& os, const PTEmbed& e);

private:
    Loc loc_;
    std::string filename_;
    std::string contents_;
};

// Control block commands
class PTCommand {
public: 
    PTCommand() {};
    PTCommand(const Loc& l, Id name) : loc_(l), name_(name) {};
    
    PTCommand           (const PTCommand&)  = delete;
    PTCommand           (      PTCommand&&) = default;
    PTCommand& operator=(const PTCommand&)  = delete;
    PTCommand& operator=(      PTCommand&&) = default;

    // Getters
    Loc location() const { return loc_; };
    Id name() const { return name_; };
    const PTIdentifierList& keywords() const { return keywords_; }; 
    const std::vector<Rpn>& expressions() const { return expressions_; }; 
    const PTParameters& args() const { return args_; }; 
    const PTSaves& saves() const { return saves_; };
    PTSaves& saves() { return saves_; };
    
    // No fluent API for now
    void set(PTIdentifierList&& kw) { keywords_ = std::move(kw); };
    void set(PTParameters&& args) { args_ = std::move(args); };
    void set(std::vector<Rpn>&& exprs) { expressions_ = std::move(exprs); };
    void set(PTSaves&& s) { saves_ = std::move(s); };

    friend std::ostream& operator<<(std::ostream& os, const PTCommand& c);

private:
    Loc loc_;
    Id name_;
    std::vector<Rpn> expressions_;
    PTIdentifierList keywords_;
    PTParameters args_;
    PTSaves saves_;
};

typedef std::variant<PTAnalysis, PTCommand> ControlEntry;
typedef std::vector<ControlEntry> PTControl;


// Top level parser tables structure
class ParserTables {
public:
    ParserTables() {};
    ParserTables(const std::string& title) : title_(title) {};
    
    ParserTables           (const ParserTables&)  = delete;
    ParserTables           (      ParserTables&&) = default;
    ParserTables& operator=(const ParserTables&)  = delete;
    ParserTables& operator=(      ParserTables&&) = default;

    // Getters
    const std::string& title() const { return title_; };
    FileStack& fileStack() { return fileStack_; };
    const std::vector<PTLoad>& loads() const { return loads_; };
    const PTSubcircuitDefinition& defaultSubDef() const { return defaultSubDef_; };
    PTSubcircuitDefinition& defaultSubDef() { return defaultSubDef_; };
    const PTIdentifierList& groundNodes() const { return groundNodes_; };
    const PTIdentifierList& globalNodes() const { return globalNodes_; };
    const PTControl& control() const { return control_; }; 
    PTControl& control() { return control_; }; 
    const std::vector<PTEmbed>& embed() const { return embed_; };
    Accounting& accounting() { return acct_; };
    
    // Fluent API
    ParserTables& setTitle(const std::string t) & { title_ = t; return *this; };
    ParserTables& setDefaultSubDef(PTSubcircuitDefinition&& def) & { defaultSubDef_ = std::move(def); return *this; };
    ParserTables& add(PTLoad&& ld) & { loads_.push_back(std::move(ld)); return *this; };
    ParserTables& addGround(PTParsedIdentifier parsedId) & { groundNodes_.push_back(parsedId); return *this; };
    ParserTables& addGlobal(PTParsedIdentifier parsedId) & { globalNodes_.push_back(parsedId); return *this; }; 
    ParserTables& defaultGround() & {
        if (groundNodes_.size()==0) {
            groundNodes_.push_back(PTParsedIdentifier(Id("0"), Loc::bad));
        }
        return *this;
    };
    ParserTables& add(PTEmbed&& e) & { embed_.push_back(std::move(e)); return *this; };
    ParserTables&& setTitle(const std::string t) && { return std::move(this->add(t)); };
    ParserTables&& setDefaultSubDef(PTSubcircuitDefinition&& def) && { return std::move(this->setDefaultSubDef(std::move(def))); };
    ParserTables&& add(PTLoad&& ld) && { return std::move(this->add(std::move(ld)));; };
    ParserTables&& addGround(PTParsedIdentifier parsedId) && { return std::move(this->addGround(parsedId)); };
    ParserTables&& addGlobal(PTParsedIdentifier parsedId) && { return std::move(this->addGlobal(parsedId)); }; 
    ParserTables&& defaultGround() && { return std::move(this->defaultGround()); };
    ParserTables&& add(PTEmbed&& e) && { return std::move(this->add(std::move(e))); };
    
    // Control block, no fluent API for now
    void addCommand(PTCommand&& c) { control_.push_back(std::move(c)); };
    void addCommand(PTAnalysis&& a) {  control_.push_back(std::move(a)); };
    
    // Post-parse checks
    bool verify(Status& s=Status::ignore) const;

    // Write embedded files
    bool writeEmbedded(int debug=0, Status& s=Status::ignore);

    void dump(int indent, std::ostream& os) const;

private:
    Accounting acct_;
    std::string title_;
    FileStack fileStack_;
    PTSubcircuitDefinition defaultSubDef_;
    PTIdentifierList globalNodes_; // Order is not important
    PTIdentifierList groundNodes_; // Order matters, first ground node is the name of the ground node, rest are just aliases
    std::vector<PTLoad> loads_;  
    std::vector<PTEmbed> embed_;
    PTControl control_;
};

// Helpers for API users

// Generic vector constructor
template<class T, class... Args>
std::vector<T> make_vector(Args&&... args) {
    std::vector<T> v;
    v.reserve(sizeof...(Args));
    (v.emplace_back(std::forward<Args>(args)), ...);
    return v;
}

// Constant parameter value vector constructor
using PV = PTParameterValue;
template<class... Args>
std::vector<PV> PVv(Args&&... args) {
    return make_vector<PV>(std::forward<Args>(args)...);
}

// Parameter expression vector constructor
using PE = PTParameterExpression;
template<class... Args>
std::vector<PE> PEv(Args&&... args) {
    return make_vector<PE>(std::forward<Args>(args)...);
}

}

#endif
