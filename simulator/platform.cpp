#include "platform.h"
#include "libplatform.h"
#include "filesystem.h"
#include "processutils.h"
#include "common.h"
#include <tuple>


namespace NAMESPACE {

std::string Platform::openVafName_;
std::vector<std::string> Platform::openVafArgs_;

bool Platform::setup(
    const std::string& openVafName, 
    const std::vector<std::string>& openVafArgs
) {
    if (openVafName.size()>0) {
        openVafName_ = openVafName;
    } else {
        // Default
        openVafName_ = Platform::defaultOpenVafBinaryName();
    }
    openVafArgs_ = openVafArgs;

    return true;
}

const char* Platform::defaultOpenVafBinaryName() {
#ifdef SIMWINDOWS
    static const char binary[] = "openvaf-r.exe"; 
#else
    static const char binary[] = "openvaf-r"; 
#endif
    return binary;
}

static std::string initPythonExecutable() {
#ifdef SIMWINDOWS
    auto [found, pythonExecutable_] =  findFileInSystemPath("python.exe");
    if (found) {
        return pythonExecutable_;
    } else {
        return "";
    }
#else
    // Try python3
    auto [found, pythonExecutable_] =  findFileInSystemPath("python3");
    if (found) {
        return pythonExecutable_;
    } else {
        // Try python
        auto [found, pythonExecutable_] =  findFileInSystemPath("python");
        if (found) {
            return pythonExecutable_;
        } else {
            return "";
        }
    } 
#endif
}

const std::string& Platform::pythonExecutable() {
    static std::string pythonExecutable_ = initPythonExecutable();
    return pythonExecutable_;
}

static std::string initPythonPath() {
    auto lib = Platform::libraryPath();
    if (!lib.empty()) {
        auto pythonPath_ = lib / "python";
        return pythonPath_.string();
    }
    return "";
}

const std::string& Platform::pythonPath() {
    static std::string pythonPath_ = initPythonPath();
    return pythonPath_;
}

const std::filesystem::path& Platform::libraryPath() {
#ifdef SIMWINDOWS
    static auto libPath = std::filesystem::path(executableFile()).parent_path().parent_path() / "lib";
#else
    static auto libPath = std::filesystem::path(executableFile()).parent_path().parent_path() / "lib" / programName;
#endif
    return libPath;
}

#define str_expand(s) #s
#define to_str(s) str_expand(s)
const std::string Platform::programName = to_str(PROGRAM_NAME);
const std::string Platform::programVersion = to_str(PROGRAM_VERSION);
const std::string Platform::programCopyright = "(c)" to_str(PROGRAM_COPYRIGHT) " EDA Lab FE Uni-Lj, Arpad Buermen";
#undef str_expand
#undef to_str
} 
