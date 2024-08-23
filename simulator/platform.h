#ifndef __PLATFORM_DEFINED
#define __PLATFORM_DEFINED

#include <vector>
#include <string>
#include <filesystem>
#include "common.h"

namespace NAMESPACE {

class Platform {
public:
    static bool setup(
        const std::string& openVafName="", 
        const std::vector<std::string>& openVafArgs={}
    );
    static const std::string& openVafName() { return openVafName_; };
    static const std::vector<std::string>& openVafArgs() { return openVafArgs_; };
    static const std::filesystem::path& libraryPath();
    static const std::string& pythonExecutable();
    static const std::string& pythonPath();

    static bool isTty(std::ostream& os);
    static int ttyColumns(std::ostream& os);

    // Program info
    static const std::string programName; 
    static const std::string programVersion; 
    static const std::string programCopyright; 

private:
    static const char* defaultOpenVafBinaryName();
    static std::string openVafName_;
    static std::vector<std::string> openVafArgs_;
};

}

#endif
