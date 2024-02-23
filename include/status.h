#ifndef __STATUS_DEFINED
#define __STATUS_DEFINED

#include <string>
#include "filestack.h"
#include "sourceloc.h"
#include "common.h"


namespace NAMESPACE {

// API level 0 provides only error messages
// API level 1 adds location to error message

class Status {
public:
    enum Code : char {
        OK = 0, 
        NotFound, 
        Missing, 
        Conflicting, 
        Empty, 
        BadVersion, 
        Unknown, 
        BadArguments, 
        BadConversion, 
        SizeMismatch, 
        CreationFailed, 
        Redefinition, 
        Unmatched, 
        Syntax, 
        Unsupported, 
        Evaluation, 
        NoMem, 
        Internal, 
        Exception, 
        DivZero, 
        Range, 
        NoGround, 
        BadNode, 
        BadTerminal, 
        MatrixCreate, 
        MatrixAnalysis, 
        MatrixOrdering, 
        MatrixType, 
        MatrixLU,
        MatrixSolve, 
        MatrixRcond, 
        MatrixRgrowth, 
        MatrixResidual, 
        Save, 
        Domain, 
        HierarchicalRecursion, 
        NotConverged, 
        NoAlg, 
        AnalysisFailed, 
        HomotopyFailed, 
        AbortRequested, 
        FinishRequested, 
        StopRequested, 
        NotFinite, 
        Compile, 
        Breakpoint, 
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
