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
    PTParameterValue(const Loc& l, const Id ident, Value&& value) : loc_(l), id_(ident), val_(std::move(value)) {};

    PTParameterValue           (const PTParameterValue&)  = delete;
    PTParameterValue           (      PTParameterValue&&) = default;
    PTParameterValue& operator=(const PTParameterValue&)  = delete;
    PTParameterValue& operator=(      PTParameterValue&&) = default;

    Id name() const { return id_; }; 
    Loc location() const { return loc_; };
    const Value& val() const { return val_; };

    void dump(int indent, std::ostream& os);

private:
    Id id_;
    Value val_;
    Loc loc_;
};

// A single parameter with name identifier and given Expression
class PTParameterExpression {
public:
    PTParameterExpression(const Loc& l, const Id ident, Rpn&& expr) : loc_(l), id_(ident), rpn_(std::move(expr)) {};

    PTParameterExpression           (const PTParameterExpression&)  = delete;
    PTParameterExpression           (      PTParameterExpression&&) = default;
    PTParameterExpression& operator=(const PTParameterExpression&)  = delete;
    PTParameterExpression& operator=(      PTParameterExpression&&) = default;

    Id name() const { return id_; }; 
    Loc location() const { return loc_; };
    const Rpn& rpn() const { return rpn_; };

    void dump(int indent, std::ostream& os);

private:
    Id id_;
    Rpn rpn_;
    Loc loc_;
};


// Parameters (constants and expressions) for instances, models, and subcircuit definitions
class PTParameters {
public:
    PTParameters();
    PTParameters(std::vector<PTParameterValue>&& pv, std::vector<PTParameterExpression>&& pe);

    PTParameters           (const PTParameters&)  = delete;
    PTParameters           (      PTParameters&&) = default;
    PTParameters& operator=(const PTParameters&)  = delete;
    PTParameters& operator=(      PTParameters&&) = default;

    // Move a value into the list of parameters values
    void add(PTParameterValue&& v);

    // Move parameter values and replace current ones
    void set(std::vector<PTParameterValue>&& v);

    // Move a value into the list of parameter expressions
    void add(PTParameterExpression&& e);

    // Move a set of parameters into parameter value/expression lists
    void add(PTParameters&& p);

    // Number of parametric expressions/values
    inline size_t valueCount() const { return values_.size(); };
    inline size_t expressionCount() const { return expressions_.size(); };
    inline size_t count() const { return values_.size()+expressions_.size(); };

    // Values
    inline const auto& values() const { return values_; };
    inline auto& values() { return values_; };

    // Expressions
    inline const auto& expressions() const { return expressions_; };
    inline auto& expressions() { return expressions_; };

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
    PTModel();
    PTModel(const Loc& l, Id name, Id device);

    PTModel           (const PTModel&)  = delete;
    PTModel           (      PTModel&&) = default;
    PTModel& operator=(const PTModel&)  = delete;
    PTModel& operator=(      PTModel&&) = default;

    inline const Loc& location() const { return loc; };
    inline Id name() const { return modelName; };
    inline Id device() const { return deviceName; };
    inline bool isParameterized() const { return parameters_.expressionCount()>0; };
    inline const PTParameters& parameters() const { return parameters_; };

    void add(PTParameters&& par);

    void dump(int indent, std::ostream& os);

protected:
    Id modelName;
    Id deviceName;
    PTParameters parameters_;
    Loc loc;
};


class PTInstance {
public:
    PTInstance();
    PTInstance(const Loc& l, Id name, Id master, PTIdentifierList&& terms, PTParameters&& params); 

    PTInstance           (const PTInstance&)  = delete;
    PTInstance           (      PTInstance&&) = default;
    PTInstance& operator=(const PTInstance&)  = delete;
    PTInstance& operator=(      PTInstance&&) = default;

    inline const Loc& location() const { return loc; };

    Id name() const { return instanceName_; };
    bool isParameterized() const { return parameters_.expressionCount()>0; };
    inline Id masterName() const { return masterName_; };
    inline const PTParameters& parameters() const { return parameters_; };
    inline const PTIdentifierList& connections() const { return connections_; };

    void dump(int indent, std::ostream& os);

private:
    Id instanceName_;
    Id masterName_;
    PTIdentifierList connections_;
    PTParameters parameters_;
    Loc loc;
};


// Block index
typedef uint32_t PTBlockIndex; 

// Block sequence with entries of the form (condition, block index)
// Blocks are stored in a blocks vector within the subcircuit definition

class PTBlockSequence {
public:
    PTBlockSequence();
    
    PTBlockSequence           (const PTBlockSequence&)  = delete;
    PTBlockSequence           (      PTBlockSequence&&) = default;
    PTBlockSequence& operator=(const PTBlockSequence&)  = delete;
    PTBlockSequence& operator=(      PTBlockSequence&&) = default;

    const std::vector<std::tuple<Loc, Rpn, PTBlockIndex>>& entries() const { return entries_; };

    void add(const Loc& l, Rpn&& cond, PTBlockIndex blockIndex);

    void dump(int indent, std::ostream& os);

private:
    std::vector<std::tuple<Loc, Rpn, PTBlockIndex>> entries_;
};



// A single netlist block comprising models, instances, and block sequences
class PTBlock {
public:
    PTBlock();
    
    PTBlock           (const PTBlock&)  = delete;
    PTBlock           (      PTBlock&&) = default;
    PTBlock& operator=(const PTBlock&)  = delete;
    PTBlock& operator=(      PTBlock&&) = default;

    inline const std::vector<PTModel>& models() const { return models_; };
    inline const std::vector<PTInstance>& instances() const { return instances_; };
    inline const std::vector<PTBlockSequence>& blockSequences() const { return blockSequences_; };
    inline std::vector<PTModel>& models() { return models_; };
    inline std::vector<PTInstance>& instances() { return instances_; };
    inline std::vector<PTBlockSequence>& blockSequences() { return blockSequences_; };

    void add(PTModel&& mod);
    void add(PTInstance&& inst);
    void add(PTBlockSequence&& seq);
    
    void dump(int indent, std::ostream& os);

private:
    std::vector<PTModel> models_;
    std::vector<PTInstance> instances_;
    std::vector<PTBlockSequence> blockSequences_;
    Loc loc;
};


class PTSubcircuitDefinition : public PTModel {
public:
    PTSubcircuitDefinition();
    PTSubcircuitDefinition(const Loc& l, Id name, PTIdentifierList&& terms);

    PTSubcircuitDefinition           (const PTSubcircuitDefinition&)  = delete;
    PTSubcircuitDefinition           (      PTSubcircuitDefinition&&) = default;
    PTSubcircuitDefinition& operator=(const PTSubcircuitDefinition&)  = delete;
    PTSubcircuitDefinition& operator=(      PTSubcircuitDefinition&&) = default;

    inline const PTIdentifierList& terminals() const { return terminals_; };
    inline const std::vector<PTBlock>& blocks() const { return blocks_; };
    inline std::vector<PTBlock>& blocks() { return blocks_; };
    inline const PTBlock& block(PTBlockIndex ndx) const { return blocks_[ndx]; };
    inline PTBlock& block(PTBlockIndex ndx) { return blocks_[ndx]; };
    inline const std::vector<std::unique_ptr<PTSubcircuitDefinition>>& subDefs() const { return subDefs_; };
    
    void add(PTIdentifierList&& terms);
    PTBlockIndex add(PTBlock&& block);
    void add(PTSubcircuitDefinition&& subDef);

    using PTModel::add;

    void dump(int indent, std::ostream& os);

    bool verifyTerminals(Status& s=Status::ignore) const;
    
private:
    PTIdentifierList terminals_;
    std::vector<PTBlock> blocks_;
    std::vector<std::unique_ptr<PTSubcircuitDefinition>> subDefs_;
};


class PTLoad {
public:
    PTLoad();
    PTLoad(const Loc& l, const std::string& file, Id module=Id::none, Id asModule=Id::none);

    PTLoad           (const PTLoad&)  = delete;
    PTLoad           (      PTLoad&&) = default;
    PTLoad& operator=(const PTLoad&)  = delete;
    PTLoad& operator=(      PTLoad&&) = default;


    inline const Loc& location() const { return loc; };
    inline const std::string& file() const { return file_; };
    inline Id module() const { return module_; };
    inline Id asModule() const { return asModule_; };

    void dump(int indent, std::ostream& os);

private:
    std::string file_;
    Id module_;
    Id asModule_;
    Loc loc;
};


class PTSave {
public:
    PTSave() {};
    PTSave(const Loc& l, const Id tname, const Id id1=Id(), const Id id2=Id()) : loc(l), typeName_(tname) { id[0] = id1; id[1] = id2; };

    PTSave           (const PTSave&)  = delete;
    PTSave           (      PTSave&&) = default;
    PTSave& operator=(const PTSave&)  = delete;
    PTSave& operator=(      PTSave&&) = default;

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
    PTSaves();
    
    PTSaves           (const PTSaves&)  = delete;
    PTSaves           (      PTSaves&&) = default;
    PTSaves& operator=(const PTSaves&)  = delete;
    PTSaves& operator=(      PTSaves&&) = default;

    void add(PTSave&& s);
    void add(std::vector<PTSave>&& s);
    
    const std::vector<PTSave>& saves() const { return saves_; }; 

    // Operator that prints save directives
    friend std::ostream& operator<<(std::ostream& os, const PTSaves& o);

private:
    std::vector<PTSave> saves_;
};

using PTSavesVector = std::vector<PTSaves*>;


class PTSweep {
public:
    PTSweep(const Loc& l, Id name, PTParameters&& par);

    PTSweep           (const PTSweep&)  = delete;
    PTSweep           (      PTSweep&&) = default;
    PTSweep& operator=(const PTSweep&)  = delete;
    PTSweep& operator=(      PTSweep&&) = default;

    Id name() const { return name_; }; 
    Loc location() const { return loc; }; 
    const PTParameters& parameters() const { return parameters_; };

    // Operator that prints a sweep directive
    friend std::ostream& operator<<(std::ostream& os, const PTSweep& o);

private:
    Id name_;
    Loc loc;
    PTParameters parameters_;
};


class PTSweeps {
public:
    PTSweeps();

    PTSweeps           (const PTSweeps&)  = delete;
    PTSweeps           (      PTSweeps&&) = default;
    PTSweeps& operator=(const PTSweeps&)  = delete;
    PTSweeps& operator=(      PTSweeps&&) = default;

    void add(PTSweep&& s);
    const std::vector<PTSweep>& data() const { return sweeps_; }; 

    void dump(int indent, std::ostream& os);

private:
    std::vector<PTSweep> sweeps_;
};


class PTAnalysis {
public:
    PTAnalysis();
    PTAnalysis(const Loc& l, Id name, Id typeName);

    PTAnalysis           (const PTAnalysis&)  = delete;
    PTAnalysis           (      PTAnalysis&&) = default;
    PTAnalysis& operator=(const PTAnalysis&)  = delete;
    PTAnalysis& operator=(      PTAnalysis&&) = default;

    void add(PTParameters&& par);
    void add(PTSweeps&& sw);
    
    inline const Loc& location() const { return loc; };
    inline Id name() const { return name_; };
    inline Id typeName() const { return typeName_; };
    const PTParameters& parameters() const { return parameters_; };
    PTParameters& parameters() { return parameters_; };
    const PTSweeps& sweeps() const { return sweeps_; };

    void dump(int indent, std::ostream& os);

private:
    Id name_;
    Id typeName_;
    PTParameters parameters_;
    PTSweeps sweeps_; 
    Loc loc;
};

// Top level parser tables structure
class ParserTables {
public:
    ParserTables();
    ~ParserTables();

    ParserTables           (const ParserTables&)  = delete;
    ParserTables           (      ParserTables&&) = default;
    ParserTables& operator=(const ParserTables&)  = delete;
    ParserTables& operator=(      ParserTables&&) = default;

    void setTitle(const std::string t);

    const std::string& title() const;
    FileStack& fileStack() { return fileStack_; };
    const std::vector<PTLoad>& loads() const { return loads_; };
    const PTSubcircuitDefinition& defaultSubDef() const { return defaultSubDef_; };
    PTSubcircuitDefinition& defaultSubDef() { return defaultSubDef_; };
    void addDefaultSubDef(PTSubcircuitDefinition&& def);
    void addLoad(PTLoad&& o);
    void addGround(PTParsedIdentifier parsedId);
    void addGlobal(PTParsedIdentifier parsedId); 
    
    const PTIdentifierList& groundNodes() const { return groundNodes_; };
    const PTIdentifierList& globalNodes() const { return globalNodes_; };
    
    void defaultGround();

    // Post-parse checks
    bool verify(Status& s=Status::ignore) const;

    // Accounting
    Accounting& accounting() { return acct_; };

    void dump(int indent, std::ostream& os);

private:
    Accounting acct_;
    std::string title_;
    FileStack fileStack_;
    PTSubcircuitDefinition defaultSubDef_;
    PTIdentifierList globalNodes_; // Order is not important
    PTIdentifierList groundNodes_; // Order matters, first ground node is the name of the ground node, rest are just aliases
    std::vector<PTLoad> loads_;  
};

}

#endif
