#include <filesystem>
#include "simulator.h"
#include "anop.h"
#include "andcinc.h"
#include "andctf.h"
#include "anac.h"
#include "anactf.h"
#include "antran.h"
#include "annoise.h" 
#include "libplatform.h"
#include "common.h"


namespace NAMESPACE { 

template<typename T> bool Simulator::registerAnalysis(Id anType, Status& s) {
    Analysis::registerFactory(anType, T::create);
    return true;
}

class NullBuffer : public std::streambuf {
public:
    int overflow(int c) { return c; }
};

static NullBuffer nullBuff; 
std::ostream Simulator::nullStream(&nullBuff);

std::string Simulator::startupPath_;
std::ostream* Simulator::out_ = &std::cout;
std::ostream* Simulator::err_ = &std::cout;
std::ostream* Simulator::dbg_ = &std::cout;
std::ostream* Simulator::wrn_ = &std::cout;

std::vector<std::string> Simulator::modulePath_;
std::vector<std::string> Simulator::includePath_;

bool Simulator::fileDebug_ = false;

void Simulator::setStreams(std::ostream& output, std::ostream& error, std::ostream& debug) {
    Simulator::out_ = &output;
    Simulator::err_ = &error;
    Simulator::dbg_ = &debug;
}

bool Simulator::setup(
    const std::string& moduleFilePathString, 
    const std::string& includeFilePathString, 
    Status& s
) {
    std::vector<std::string> modPathVec;
    std::vector<std::string> incPathVec;

    splitString(pathSeparator(), moduleFilePathString, modPathVec);
    splitString(pathSeparator(), includeFilePathString, incPathVec);

    modulePath_ = modPathVec;
    includePath_ = incPathVec;

    startupPath_ = std::filesystem::current_path().string();

    bool ok = true;
    ok &= registerAnalysis<OperatingPoint>("op", s);
    ok &= registerAnalysis<DcIncremental>("dcinc", s);
    ok &= registerAnalysis<DcTf>("dctf", s);
    ok &= registerAnalysis<Ac>("ac", s);
    ok &= registerAnalysis<AcTf>("actf", s);
    ok &= registerAnalysis<Tran>("tran", s);
    ok &= registerAnalysis<Noise>("noise", s);
    
    return ok;
}

}
