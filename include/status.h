#ifndef __STATUS_DEFINED
#define __STATUS_DEFINED

#include <string>
#include "filestack.h"
#include "sourceloc.h"
#include "common.h"


namespace NAMESPACE {

class Status {
public:
    enum Code : char {
        OK = 0, 
        NotFound, 
        Conflicting, 
        BadVersion, 
        BadArguments, 
        BadConversion, 
        CreationFailed, 
        Redefinition, 
        Syntax, 
        Unsupported, 
        Exception, 
        DivZero, 
        Range, 
        LinearSolver, 
        NonlinearSolver, 
        Save, 
        Domain, 
        Recursion,
        Force,  
        Analysis, 
        Homotopy, 
        AbortRequested, 
        FinishRequested, 
        StopRequested, 
        Process
    };

    inline Status(bool ignoreFlag=false) : ignoreFlag(ignoreFlag), error(false) {};
    
    Status           (const Status&)  = delete;
    Status           (      Status&&) = default;
    Status& operator=(const Status&)  = delete;
    Status& operator=(      Status&&);

    // Clears error code and message
    void clear();
    
    // Appends error code(s) and message(s)
    void set(Code c, const std::string& msg);
    void set(const Status& other);

    // Extends last message without setting a new error code
    void extend(const std::string& msg);
    void extend(const Loc& point);
    void prefix(const std::string& msg);

    bool ignored() { return ignoreFlag; };

    // Checks if error is set
    // operator bool() {
    //    return !error;
    //};

    // Retrieves message
    const std::string& message() const;
    
    // Singleton for ignoring status
    static Status ignore;

private:
    bool ignoreFlag = false;
    std::string message_;
    bool error;
};

}

#endif
