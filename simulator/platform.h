#ifndef __PLATFORM_DEFINED
#define __PLATFORM_DEFINED

#include <vector>
#include <string>
#include <filesystem>
#include "common.h"

namespace NAMESPACE {

class Platform {
public:
    static void setup();

    static void setOpenVaf(std::string openVaf);
    static void setOpenVafArgs(std::vector<std::string>&& openVafArgs);
    static void setPythonExecutable(std::string pythonExecutable);

    static const std::string& openVaf() { return openVaf_; };
    static const std::vector<std::string>& openVafArgs() { return openVafArgs_; };
    static const std::string& pythonExecutable() { return pythonExecutable_; };
    
    static const std::filesystem::path& libraryPath();
    static const std::string& pythonPath();
    static const std::filesystem::path& homeDir();

    static const std::string& systemConfig();
    static const std::string& userConfig();
    static const std::string& localConfig();

    static bool isTty(std::ostream& os);
    static int ttyColumns(std::ostream& os);

    // Program info
    static const std::string programName; 
    static const std::string programVersion; 
    static const std::string programCopyright; 
    static const std::string programHomepage; 

private:
    static std::string openVaf_;
    static std::vector<std::string> openVafArgs_;
    static std::string pythonExecutable_;
};

}

#endif
