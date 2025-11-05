#include "parser.h"
#include "status.h"
#include "openvafcomp.h"
#include "circuit.h"
#include "simulator.h"
#include "processutils.h"
#include "platform.h"
#include "cmd.h"
#include "config.h"
#include "common.h"
#include <filesystem>
#include <unordered_set>

#include "corehb.h"
#include "anop.h"

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
    "  --no-output         suppress output of result files\n"
    ; 

int main(int argc, char**argv) {
    // IntegratorCoeffs::test();
    // DenseMatrix<double>::test();
    // HBCore::test();
    // return 0;

    bool dumpPaths = false;
    bool dumpTables = false;
    bool fileDebug = false;
    bool dumpEmbed = true;
    bool runPostprocess = true;
    bool progress = true;
    bool noOutput = false;

    // Simulator information
    Simulator::out() << 
        "This is "+Platform::programName+" "+Platform::programVersion+".\n"+Platform::programCopyright+"\n";
    Simulator::out() << 
        Platform::programHomepage+"\n";
    #ifdef SIMDEBUG
    Simulator::out() << "\n" << "Warning! This is a debug build. Simulator will be slow.\n";
    #endif
    #ifdef SIMPROFILE
    Simulator::out() << "\n" << "Warning! This binary is instrumeted for profiling and code coverage.\n";
    #endif
    Simulator::out() << "\n";
    
    if (argc<2) {
        // No arguments, print a hint on help
        Simulator::err() << "To get help, run as: " << Platform::programName << " -h\n"; 
        return 1;
    } 
        
    // Parse arguments
    char* fileArg = nullptr;
    for(int i=1; i<argc; i++) {
        if (argv[i][0]!='-') {
            // Not an option
            if (i==argc-1) {
                // Last argument is file name, OK
                fileArg = argv[i];
                break;
            } else {
                // Not last argument, error
                Simulator::err() << "Unrecognized argument '" << argv[i] << "'.\n";
                return 1;
            }
        }
        // Must be an option
        std::string arg = argv[i];
        if (arg=="-h" || arg=="--help") {
            Simulator::out() << "Usage: " << Platform::programName << " " << helpText;
            return 0;
        } else if (arg=="-dp" || arg=="--dump-paths") {
            dumpPaths = true;
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
        } else if (arg=="--no-output") {
            noOutput = true;
        } else {
            Simulator::err() << "Unrecognized argument '"+arg+"'.\n";
            return 1;
        }
    }
    
    // FileDebug and noOutput simulator flags
    Simulator::setFileDebug(fileDebug); 
    Simulator::setNoOutput(noOutput); 

    //Status
    Status status;

    // Setup platform defaults
    Platform::setup();

    // SIM_OPENVAF environmental variable
    auto openvafCstring = std::getenv("SIM_OPENVAF"); 
    if (openvafCstring) {
        Platform::setOpenVaf(openvafCstring);
    }

    // Get path to simulator binary
    auto simulatorBinary = executableFile();

    // Default module directory and include directory
    std::string defaultModuleDirectory = (Platform::libraryPath() / "mod").string();
    std::string defaultIncludeDirectory = (Platform::libraryPath() / "inc").string(); 

    // Overrides via SIM_MODULES_PATH and SIM_SOURCES_PATH environmental variables
    auto modCstring = std::getenv("SIM_MODULE_PATH"); 
    auto incCstring = std::getenv("SIM_INCLUDE_PATH"); 
    // Defaults when env var is not defined
    auto modStr = modCstring ? std::string(modCstring) : defaultModuleDirectory;
    auto incStr = incCstring ? std::string(incCstring) : defaultIncludeDirectory;

    // Setup simulator
    if (!Simulator::setup(modStr, incStr, status)) {
        Simulator::err() << status.message() << "\n";
        return 1;
    }

    // Load file passed as argument
    ParserTables tab;
    
    // Candidate config files
    std::vector<std::string> configFiles = {
        Platform::systemConfig(), 
        Platform::userConfig(), 
        Platform::localConfig()
    };

    // Add input file to file stack
    FileStackFileIndex stackPosition;
    if (fileArg) {
        // Add input file to filestack. Search only in local directory 
        // if fileArg is not an absolute path. Stop if file not found. 
        stackPosition = tab.fileStack().addFile(fileArg);
        if (stackPosition==FileStack::badFileId) {
            Simulator::err() << std::string("File '")+fileArg+"' not found.\n";
            return 1;
        }
        
        // Get canonical name
        auto canonical = tab.fileStack().canonicalName(stackPosition);

        // Get input file directory
        auto inputFileDir = std::filesystem::path(canonical).parent_path();
        
        // Add another cadidate config file
        configFiles.push_back((inputFileDir / ".vacaskrc.toml").string());
    }

    // Read configuration files, make sure you read each one only once
    std::unordered_set<std::string> processed_configs;
    for(auto& cfg : configFiles) {
        // Was it processed already, skip if yes
        if (processed_configs.contains(cfg)) {
            continue;
        }
        processed_configs.insert(cfg);
        // Open file
        std::ifstream f(cfg);
        if (f) {
            // File found, read it
            if (fileDebug) {
                Simulator::dbg() << "Reading config file: " << cfg << "\n";
            }
            if (!readConfig(f, cfg, status)) {
                Simulator::err() << status.message() << "\n";
                return 1;
            }
        }
    }

    // Dump paths
    if (dumpPaths) {
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
        if (findProgram(Platform::openVaf(), openVafPath)) {
            Simulator::dbg() << "OpenVAF compiler: " << openVafPath << "\n";
            if (Platform::openVafArgs().size()>0) {
                Simulator::dbg() << "OpenVAF arguments: ";
                for(auto& arg : Platform::openVafArgs()) {
                    Simulator::dbg() << arg << " ";
                }
                
                Simulator::dbg() << "\n";
            }
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

    // Stop if no input file
    if (!fileArg) {
        Simulator::out() << "No input file specified.\n";
        return 0;
    }

    // Parser 
    Parser parser(tab);

    if (!parser.parseNetlistFile(stackPosition, status)) {
        Simulator::err() << status.message() << "\n";
        return 1;
    }

    if (!tab.verify(status)) {
        Simulator::err() << status.message() << "\n";
        return 1;
    }

    if (progress) {
        Simulator::dbg() << "Simulating: " << tab.title() << "\n";
    }

    if (dumpTables) {
        Simulator::dbg() << "---- Parser tables ----\n";
        tab.dump(0, Simulator::dbg());
        Simulator::dbg() << "---- Parser tables end ----\n\n";
    }

    // Dump embedded files
    if (dumpEmbed) {
        for(auto& e : tab.embed()) {
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
    OpenvafCompiler comp(Platform::openVaf(), Platform::openVafArgs());
    Circuit cirObj(tab, &comp, status);
    if (!cirObj.isValid()) {
        Simulator::err() << status.message() << "\n";
        return 1;
    }
    
    // Command interpreter
    CommandInterpreter interp(tab, tab.control(), cirObj);
    interp.setPrintProgress(progress),
    interp.setRunPostprocess(runPostprocess);
    
    // Run iterpreter
    if (!interp.run(status)) {
        Simulator::err() << status.message() << "\n";
        return 1;
    }
    
    
    return 0;
}
