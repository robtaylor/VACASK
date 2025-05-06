#ifndef __SIMULATOR_DEFINED
#define __SIMULATOR_DEFINED

#include <vector>
#include <string>
#include "status.h"
#include "identifier.h"
#include "introspection.h"
#include "options.h"
#include "an.h"
#include "common.h"
#include <iostream>


namespace NAMESPACE {

class Simulator {
public:
    static void setStreams(std::ostream& output, std::ostream& error, std::ostream& debug);
    static bool setup(        
        const std::string& moduleFilePath, 
        const std::string& includeFilePath, 
        Status& s=Status::ignore
    );
    
    static const std::string& startupPath() { return startupPath_; };
    static const std::vector<std::string>& modulePath() { return modulePath_; };
    static const std::vector<std::string>& includePath() { return includePath_; };
    
    static const int majorVersion = 0;
    static const int minorVersion = 1;

    static std::ostream& out() { return *out_; };
    static std::ostream& err() { return *err_; };
    static std::ostream& dbg() { return *dbg_; };
    static std::ostream& wrn() { return *wrn_; };
    
    static bool fileDebug() { return fileDebug_; }; 
    static void setFileDebug(bool val) { fileDebug_ = val; };

    static bool noOutput() { return noOutput_; }; 
    static void setNoOutput(bool val) { noOutput_ = val; };

private:
    template<typename T> static bool registerAnalysis(Id anType, Status& s=Status::ignore);
    static std::ostream* out_;
    static std::ostream* err_;
    static std::ostream* dbg_;
    static std::ostream* wrn_;

    static std::ostream nullStream;

    static std::vector<std::string> modulePath_;
    static std::vector<std::string> includePath_;
    static std::string startupPath_;
    static bool fileDebug_;
    static bool noOutput_;
};

}

#endif
