#ifndef __RPNSTACK_DEFINED
#define __RPNSTACK_DEFINED

#include "value.h"
#include <vector>
#include "common.h"


namespace NAMESPACE {

class RpnStack {
public:
    RpnStack() {};
    
    RpnStack           (const RpnStack&)  = delete;
    RpnStack           (      RpnStack&&) = default;
    RpnStack& operator=(const RpnStack&)  = delete;
    RpnStack& operator=(      RpnStack&&) = default;

    Value* get(size_t i=0) { 
        if (i<stack.size())
            return &(stack[stack.size()-1-i]); 
        else 
            return nullptr; 
    };
    const Value* get(size_t i=0) const { 
        if (i<stack.size())
            return &(stack[stack.size()-1-i]); 
        else 
            return nullptr; 
    };
    inline bool pop() { 
        if (stack.size()>0) { 
            stack.pop_back(); 
            return true; 
        } else {
            return false; 
        }
    };
    inline bool pop(Value& v) { 
        if (stack.size()>0) { 
            swap(v, stack.back());
            stack.pop_back();
            return true;
        }
        return false;
    };
    inline bool pop(size_t n) {
        if (stack.size()>=n) {
            if (n>0) {
                stack.resize(stack.size()-n);
            }
            return true;
        } else {
            return false;
        }
    };
    inline bool push(const Value& v) {
        stack.push_back(v);
        return true;
    };
    inline bool push(Value&& v) {
        stack.push_back(std::move(v));
        return true;
    };
    size_t size() const { return stack.size(); };
    void clear() { stack.clear(); };
    
private:
    std::vector<Value> stack;
};

}

#endif
