#include "parser.h"
#include "status.h"
#include "srccompiler.h"
#include "circuit.h"
#include "simulator.h"
#include "processutils.h"
#include "platform.h"
#include "cmd.h"
#include "common.h"
#include <filesystem>


using namespace sim;

char helpText[] = 
    /* Usage: programName */ "[options] [<filename>]\n"
    "\nOptions:\n"
    "  -h, --help          print help and exit\n"
    "  -dp, --dump-paths   print location of simulator's components\n"
    "  -df, --debug-files  turn on debug output for file operations\n"
    "  -dt, --dump-tables  dump parser tables\n"
    "  -se, --skip-embed   do not dump embedded files\n"
    "  -sp, --skip-postprocess\n"
    "                      do not run postprocessing steps\n"
    // "  -qw, --quiet-warnings\n"
    // "                      turn off warning messages\n"
    "  -qp, --quiet-progress\n"
    "                      turn off progress messages\n"
    ; 

#include "generator.h"

Generator<int> counter(int from, int to) {
    for(int i=from; i<to; i++) {
        co_yield i;
    }
}


int main(int argc, char**argv) {
    // Workaround for boost crash
    // setenv("LC_ALL", "C", 1);

    /*
    // Create a generator, run it until initial suspend
    Generator<int> gen;
    gen = std::move(counter(0, 10));
    gen = std::move(counter(0, 5));
    
    while (gen) {
        std::cout << gen() << "\n";
    }

    return 0;
    */

    Status status;

    // Setup platform
    //   openVafName, openVafArgs
    Platform::setup("", {});

    // Get path to simulator binary
    auto simulatorBinary = executableFile();
    
    // Default module directory and include directory
    std::string defaultModuleDirectory = (Platform::libraryPath() / "mod").string();
    std::string defaultIncludeDirectory = (Platform::libraryPath() / "inc").string(); 

    // Get SIM_MODULES_PATH and SIM_SOURCES_PATH variables, separator is ":"
    auto modCstring = std::getenv("SIM_MODULE_PATH"); 
    // Defaults when env var is not defined
    auto modStr = modCstring ? std::string(modCstring) : defaultModuleDirectory;
    auto incCstring = std::getenv("SIM_INCLUDE_PATH"); 
    // Defaults when env var is not defined
    auto incStr = incCstring ? std::string(incCstring) : defaultIncludeDirectory;

    // Setup simulator
    if (!Simulator::setup(modStr, incStr, status)) {
        Simulator::err() << status.message() << "\n";
        return 1;
    }

    Simulator::out() << 
        "This is "+Platform::programName+" "+Platform::programVersion+".\n"+Platform::programCopyright+"\n";
#ifdef SIMDEBUG
    Simulator::out() << "\n" << "Warning! This is a debug build. Simulator will be slow.\n";
#endif
#ifdef SIMPROFILE
    Simulator::out() << "\n" << "Warning! This binary is instrumeted for profiling and code coverage.\n";
#endif
    Simulator::out() << "\n";

    bool needFile = true;
    bool paths = false;
    bool dumpTables = false;
    bool fileDebug = false;
    bool dumpEmbed = true;
    bool runPostprocess = true;
    bool progress = true;

    Parser parser;

    if (argc<2) {
        // No arguments, print a hint on help
        Simulator::err() << "To get help, run as: " << Platform::programName << " -h\n"; 
        return 1;
    } else {
        // Parse arguments
        int i;
        for(i=1; i<argc; i++) {
            if (argv[i][0]!='-') {
                if (i==argc-1) {
                    // Last argument, OK
                    break;
                } else {
                    // Not last argument, error
                    Simulator::err() << "Unrecognized argument '" << argv[i] << "'.\n";
                    return 1;
                }
            }
            std::string arg = argv[i];
            if (arg=="-h" || arg=="--help") {
                Simulator::out() << "Usage: " << Platform::programName << " " << helpText;
                return 0;
            } else if (arg=="-dp" || arg=="--dump-paths") {
                paths = true;
            } else if (arg=="-df" || arg=="--debug-files") {
                fileDebug = true;
            } else if (arg=="-dt" || arg=="--dump-tables") {
                dumpTables = true;
            } else if (arg=="-se" || arg=="--skip-embed") {
                dumpEmbed = false;
            } else if (arg=="-sp" || arg=="--skip-postprocess") {
                runPostprocess = false;
            } else if (arg=="-qp" || arg=="--quiet-progress") {
                progress = false;
            } else {
                Simulator::err() << "Unrecognized argument '"+arg+"'.\n";
                return 1;
            }
        }
        
        Simulator::setFileDebug(fileDebug); 

        int fileArgIndex = i;

        // Dump paths
        if (paths) {
            Simulator::dbg() << "Simulator binary: " << simulatorBinary << "\n";
            Simulator::dbg() << "Startup directory: " << Simulator::startupPath() << "\n";
            Simulator::dbg() << "Module path:\n";
            for(auto& d : Simulator::modulePath()) {
                Simulator::dbg() << "  " << d << "\n";
            }
            Simulator::dbg() << "Include path:\n";
            for(auto& d : Simulator::includePath()) {
                Simulator::dbg() << "  " << d << "\n";
            }
            std::string openVafPath;
            if (findProgram(Platform::openVafName(), openVafPath)) {
                Simulator::dbg() << "OpenVAF compiler: " << openVafPath << "\n";
            } else {
                Simulator::dbg() << "OpenVAF compiler not found.\n";
            }
            if (Platform::pythonExecutable().size()>0) {
                Simulator::dbg() << "Python interpreter: " << Platform::pythonExecutable() << "\n";
            }
            if (Platform::pythonPath().size()>0) {
                Simulator::dbg() << "Python path addition: " << Platform::pythonPath() << "\n";
            }
            Simulator::dbg() << "\n";
        }

        if (!needFile) {
            return 0;
        }

        if (fileArgIndex>=argc) {
            Simulator::err() << "No input file specified.\n";
            return 1;
        }

        // Load file passed as argument
        ParserTables tab;
        ParserExtras extras;

        if (!parser.parseNetlistFile(argv[fileArgIndex], tab, extras, status)) {
            Simulator::err() << status.message() << "\n";
            return 1;
        }

        if (!tab.verify(status) || !extras.verify(status)) {
            Simulator::err() << status.message() << "\n";
            return 1;
        }

        if (progress) {
            Simulator::dbg() << "Simulating: " << tab.title() << "\n";
        }

        if (dumpTables) {
            Simulator::dbg() << "---- Parser tables ----\n";
            tab.dump(0, Simulator::dbg());
            extras.dump(0, Simulator::dbg()); 
            Simulator::dbg() << "---- Parser tables end ----\n\n";
        }

        // Dump embedded files
        if (dumpEmbed) {
            for(auto& e : extras.embed()) {
                // Get canonical path of file with the embed directive
                auto [fs, pos, line, offset] = e.location().data();
                auto timeRefCanonicalPath = fs->canonicalName(pos);

                // Check if the file to be dumped exists
                bool dump = false;
                if (!std::filesystem::exists(e.filename())) {
                    // Does not exist, needs dumping
                    dump = true;
                } else {
                    // File with embed directive is newer than the dumped file, needs dumping
                    auto refModificationTime = std::filesystem::last_write_time(timeRefCanonicalPath);
                    auto fileModificationTime = std::filesystem::last_write_time(e.filename());
                    if (refModificationTime>fileModificationTime) {
                        dump = true;
                    }
                }

                // TODO: multiple .sim files can create an embedded file with the same name
                // Therefore we must also check the embedded file's origin. 
                // Until this is implemented, always dump file. 
                dump = true;
                
                if (dump) {
                    if (Simulator::fileDebug()) {
                        Simulator::dbg() << "Dumping embedded file '" << e.filename() << "'.\n";
                    }
                    std::ofstream fs;
                    fs.open(e.filename(), std::ios::out);
                    fs << e.contents();
                    if (fs.fail()) {
                        status.set(Status::CreationFailed, "Failed to write file '"+e.filename()+"'.");
                        status.extend(e.location());
                        Simulator::err() << status.message() << "\n";
                        return 1;
                    }
                    fs.close();
                }
            }
        }
        
        // Create circuit
        OpenvafCompiler comp;
        Circuit cirObj(tab, comp, status);
        if (!cirObj.isValid()) {
            Simulator::err() << status.message() << "\n";
            return 1;
        }
        
        // Command interpreter
        CommandInterpreter interp(tab, extras, cirObj);
        interp.setPrintProgress(progress),
        interp.setRunPostprocess(runPostprocess);
        
        // Run iterpreter
        if (!interp.run(tab, extras, status)) {
            Simulator::err() << status.message() << "\n";
            return 1;
        }
    }
    
    return 0;
}
