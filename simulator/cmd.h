#ifndef __CMD_DEFINED 
#define __CMD_DEFINED

#include <unordered_set>
#include <limits>
#include <variant>
#include <unordered_map>
#include "parseroutput.h"
#include "circuit.h"
#include "status.h"
#include "common.h"

namespace NAMESPACE {

class CommandInterpreter;

typedef bool (*CommandFuncPtr)(CommandInterpreter& interpreter, PTCommand& cmd, Status& s);

template <typename T> bool evaluateExpressions(RpnEvaluator& e, const PTCommand& cmd, std::vector<T>& out, Status& s=Status::ignore);

class CommandInterpreter {
public:
    typedef struct CmdDesc {
        static const size_t many = std::numeric_limits<size_t>::max();

        size_t minKw {0};
        size_t maxKw {0};
        size_t minExpr {0};
        size_t maxExpr {0};
        bool limitArgs {false};
        std::unordered_set<Id> allowedArgs;
        CommandFuncPtr func;
    } CmdDesc;

    CommandInterpreter(ParserTables& tables, PTControl& control, Circuit& circuit);
    ~CommandInterpreter();

    CommandInterpreter           (const CommandInterpreter&)  = delete;
    CommandInterpreter           (      CommandInterpreter&&) = delete;
    CommandInterpreter& operator=(const CommandInterpreter&)  = delete;
    CommandInterpreter& operator=(      CommandInterpreter&&) = delete;

    bool postprocessingAllowed() { return runPostprocess_; };

    void setPrintProgress(bool b) { printProgress_ = b; };
    void setRunPostprocess(bool b) { runPostprocess_ = b; };

    bool printProgress() const { return printProgress_; }; 
    bool runPostprocess() const { return runPostprocess_; }; 

    void clearSaves() { commonSaves_.clear(); };
    void addSaves(PTSaves& s) { commonSaves_.insert(commonSaves_.end(), s.saves().begin(), s.saves().end()); }; 

    void clearUserOptions() { userOptions_.clear(); };
    void addUserOption(const PTParameterValue& pv);
    void addUserOption(const PTParameterExpression& pe);
    
    // If circuit is not elaboarated, perform default elaboration, otherwise elaborate only changes
    bool minimalElaboration(Status& s);
    // If circuit is not elaboarated, perform default elaboration, otherwise do nothing
    bool defaultElaboration(Status& s);
    // Elaborate circuit from given toplevel definitions
    bool elaborate(const std::vector<Id>& names, const std::string& topDefName, const std::string& topInstName, Status& s=Status::ignore);
    
    bool run(Status& s=Status::ignore);

    bool clearVariables(Status& s=Status::ignore);
    Circuit& circuit() { return circuit_; }; 
    ParserTables& tables() { return tables_; };
    RpnEvaluator& variableEvaluator() { return circuit_.variableEvaluator(); }; 

    bool addAbort(Id cmd);
    void clearAborts();
    void setAbortOnMatch(bool b);
    bool mustAbort(Id cmd);

    void dumpOptionsMap(int indent, std::ostream& os) const;
    void dumpSaves(int indent, std::ostream& os) const;

private:
    bool printProgress_;
    bool runPostprocess_;
    std::vector<PTSave> commonSaves_;
    PTParameterMap userOptions_;
    
    std::unordered_set<Id> abortCommands;
    bool abortOnMatch;

    static std::unordered_map<Id, CmdDesc> commandDescriptors;

    ParserTables& tables_;
    PTControl& control_;
    Circuit& circuit_;
};

}

#endif
