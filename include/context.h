#ifndef __CONTEXT_DEFINED
#define __CONTEXT_DEFINED

#include "status.h"
#include "identifier.h"
#include "value.h"
#include "rpnexpr.h"
#include "rpnbuiltin.h"
#include <unordered_map>
#include <vector>
#include <variant>
#include "common.h"


namespace NAMESPACE {

typedef struct Builtin {
    Rpn::Arity minArity;
    Rpn::Arity maxArity;
    bool pure;
    RpnBuiltinFunc func;
} Builtin;


class Context {
public:
    Context();
    Context(std::initializer_list<std::pair<const Id, Value>> inl);

    size_t size() const { return data.size(); }; 

    // Returns true on insertion, false on replacement
    bool insert(Id name, const Value& v);
    bool insert(Id name, Value&& v);

    // Return value: inserted, changed
    std::tuple<bool, bool> insertAndCheck(Id name, const Value& v);
    std::tuple<bool, bool> insertAndCheck(Id name, Value&& v);

    // Returns nullptr if not found
    const Value* get(Id name) const;
    Value* get(Id name);

    // Clears context
    void clear();
    
    void dump(int indent, std::ostream& os) const;

private:
    std::unordered_map<Id,Value> data;
};


class ContextStack {
public:
    ContextStack();

    ContextStack           (const ContextStack&)  = delete;
    ContextStack           (      ContextStack&&) = default;
    ContextStack& operator=(const ContextStack&)  = delete;
    ContextStack& operator=(      ContextStack&&) = default;

    void clear();

    inline int depth() const { return stack.size(); };
    inline bool isGlobal() const { return depth()==1; };
    
    // Creates a new context
    inline void enter(Context* external=nullptr, bool addToPath=false) { 
        if (external) {
            stack.emplace_back(external); 
        } else {
            stack.emplace_back(std::move(Context())); 
        }
        if (addToPath) {
            searchPath.push_back(stack.size()-1);
        }
    };

    // Returns context stack at level i
    Context& at();
    const Context& at() const;
    Context& at(size_t i);
    const Context& at(size_t i) const;

    // Exits last context and returns true on success, global context cannot be exited (returns false)
    // External context is unlinked. 
    // Internal context is deleted. 
    // If context was in search path it is removed from search path. 
    inline bool exit() { 
        if (stack.size()>1) { 
            if (searchPath.size()>0 && searchPath.back()==stack.size()-1) {
                searchPath.pop_back();
            }
            stack.pop_back(); 
            return true; 
        } else { 
            return false; 
        } 
    };

    // Exits contexts until stack size is equal to stackMarker
    inline bool exit(size_t stackMarker) { 
        if (stackMarker>stack.size()) {
            return false;
        }
        stack.resize(stackMarker);
        while (searchPath.size()>0 && searchPath.back()>stackMarker-1) {
            searchPath.pop_back();
        }
        return true;
    };

    // Adds a value to the last context, prevents insertion of values that mask constants
    bool insert(Id name, const Value& v, Status& s=Status::ignore);
    bool insert(Id name, Value&& v, Status& s=Status::ignore);

    // Adds a value to the last context, checks for change, prevents insertion of values that mask constants
    std::tuple<bool, bool> insertAndCheck(Id name, const Value& v, Status& s=Status::ignore);
    std::tuple<bool, bool> insertAndCheck(Id name, Value&& v, Status& s=Status::ignore);

    // Retrieves a value, returns nullptr if value is not found
    const Value* get(Id name) const;

    // Retrieves a builtin
    static const Builtin* getBuiltin(Id name);

    // Checks if name is a constant
    static bool isConstant(Id name) { return consts.get(name)!=nullptr; };

    void dump(int indent, std::ostream& os) const;

    typedef std::unordered_map<Id,Builtin> Builtins;
    
private:
    typedef std::variant<Context, Context*> StackEntry;
    static Builtins builtins;
    static Context consts;
    std::vector<StackEntry> stack;
    std::vector<size_t> searchPath;
};

}

#endif
