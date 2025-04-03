#ifdef SIMMSC
#define _WIN32_WINNT 0x0601
#endif
#include <boost/process.hpp>
#include <boost/process/args.hpp>
#include <boost/dll.hpp>

#include <exception>
#include <filesystem>
#include "processutils.h"
#include "filesystem.h"
#include "libplatform.h"
#include "simulator.h"
#include "status.h"


namespace NAMESPACE {

bool findProgram(const std::string& prog, std::string& path) {
    static std::unordered_map<std::string, std::string> cache;

    // Try cache
    auto it = cache.find(prog);
    if (it!=cache.end()) {
        path = it->second;
        return true;
    } else {
        // First look in the directory of the executable file
        // Then look in the system path
        if (findFile((std::filesystem::path(executableFile()).parent_path() / prog).string(), path)) {
            cache.insert({prog, path});
            return true;
        } else if (auto [ok, sysPath] = findFileInSystemPath(prog); ok) {
            path = std::move(sysPath);
            cache.insert({prog, path});
            return true;
        } 
    }
    return false;
}

// TODO: replace this with something better
std::tuple<bool, std::string, std::string> runProcess(
    const std::string& prog, const std::vector<std::string>& args, 
    const std::string* pythonPath, bool collect, bool debugFiles, Status& s
) {
    
    bool ok;
    std::string path;

    if (!findProgram(prog, path)) {
        s.set(Status::NotFound, "Executable '"+prog+"' not found.");
        return std::make_tuple(false, "", "");
    }

    if (debugFiles) {
        Simulator::dbg() << "Running: '" << path;
        for(auto& arg : args) {
            Simulator::dbg() << " \"" << arg << "\"";
        }
        Simulator::dbg() << "'.\n";
    }

    auto procEnv = boost::this_process::environment();
    boost::process::environment customEnv = procEnv;
    if (pythonPath) {
        customEnv["PYTHONPATH"] += pathSeparator()+(*pythonPath);
    }

    if (collect) {
        boost::asio::io_context ios;
        std::future<std::string> fout;
        std::future<std::string> ferr;
        
        boost::process::child c(
            path, 
            boost::process::args(args), 
            customEnv, 
            boost::process::std_out > fout, 
            boost::process::std_err > ferr, 
            ios
        );

        ios.run();
        c.wait();

        bool ok = c.exit_code()==0;
        if (!ok) {
            s.set(Status::Process, "Error running '"+prog+"', exit code="+std::to_string(c.exit_code())+".");
        }

        return std::make_tuple(ok, fout.get(), ferr.get()); 
    } else {
        boost::process::child c(path, boost::process::args(args), customEnv);

        // while (c.running())
        //    do_some_stuff();

        c.wait(); //wait for the process to exit   
        
        bool ok = c.exit_code()==0;
        if (!ok) {
            s.set(Status::Process, "Error running '"+prog+"', exit code="+std::to_string(c.exit_code())+".");
        }

        return std::make_tuple(ok, "", "");
    }
}

std::string executableFile() {
    static auto str = boost::dll::program_location().string();
    return str;
}

}
