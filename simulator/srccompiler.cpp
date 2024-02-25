#include <filesystem>
#include "processutils.h"
#include "simulator.h"
#include "srccompiler.h"
#include "platform.h"
#include "common.h"

namespace NAMESPACE {

// loadDirectiveCanonicalPath is the canonical path to the file with the load directive
// fileName is the name of file from the load directive
// canonicalPath is the canonical path of the found file
// outputCanonicalPath is the canonical path of the output file
// 
// if extension is .va
//   look for .osdi in the current directory
//   if found 
//     if .osdi file is older than .va file
//       compile, store in the same directory
//   else
//     compile, store in the same directory
std::tuple<bool, bool> OpenvafCompiler::compile(const std::string& loadDirectiveCanonicalPath, const std::string& fileName, const std::string& canonicalPath, std::string& outputCanonicalPath, Status& s) {
    auto extension = std::filesystem::path(fileName).extension();
    if (extension==".va" || extension==".VA") {
        // Found a .va file, see if we need to compile it
        // Look for .osdi file in the same directory
        auto pVa = std::filesystem::path(canonicalPath);
        
        // Directory of .va file
        auto pVaDir = pVa.parent_path();
        
        // va file and osdi file name
        auto vaFile = pVa.filename();
        auto osdiFile = pVa.filename();
        osdiFile.replace_extension(".osdi");
        
        decltype(pVaDir) outputPath;

        // Current directory
        auto cwd = std::filesystem::current_path();

        // Store .osdi file in the same directory as .va file
        // outputPath = pVaDir;
        // Store .osdi file in the current directory
        outputPath = cwd;
        
        // Full path to .osdi file
        outputPath /= osdiFile;
        
        // Switch to .va file directory 
        // (in case it includes stuff from the same directory)
        std::filesystem::current_path(pVaDir);

        // Do we need to compile 
        bool compile = false;
        if (!std::filesystem::exists(outputPath)) {
            compile = true;
        } else { 
            auto vaModificationTime = std::filesystem::last_write_time(vaFile);
            auto osdiModificationTime = std::filesystem::last_write_time(outputPath);
            if (vaModificationTime>osdiModificationTime) {
                compile = true;
            }
        }

        // Compile
        if (compile) {
            if (Simulator::fileDebug()) {
                Simulator::dbg() << "Compiling file '" << pVa.string() << "'.\n";
            }
            auto args = Platform::openVafArgs();
            args.push_back("-o");
            args.push_back(outputPath);
            args.push_back(vaFile);
            auto [ok, out, err] = runProcess(Platform::openVafName(), args, nullptr, true, s);
            if (!ok) {
                // Failure, error
                s.extend("Failed to compile file '"+vaFile.string()+"'.");
                s.extend(err);
                // Restore current directory
                std::filesystem::current_path(cwd);
                return std::make_tuple(false, false);
            }
        }

        // Set canonical path of .osdi file
        outputCanonicalPath = std::filesystem::canonical(outputPath).string();

        std::filesystem::current_path(cwd);

        return std::make_tuple(true, true);
    }
    // Unsupported extension, do nothing, return the original file path
    return std::make_tuple(true, false);
}

}
