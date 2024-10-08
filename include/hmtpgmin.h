#ifndef __HMTPGMIN_DEFINED
#define __HMTPGMIN_DEFINED

#include "homotopy.h"
#include "common.h"

namespace NAMESPACE {

// Based on Ngspice dynamic_gmin() and new_gmin()
class GminStepping : public Homotopy {
public:
    GminStepping(Circuit& circuit, AnalysisCore& core, bool gdev=false);

    std::tuple<bool, bool> run();
    
    std::string formatProgress() const;

private:
    bool gdev;
};


// Based on Ngspice spice3_gmin()
class Spice3GminStepping : public Homotopy {
public:
    Spice3GminStepping(Circuit& circuit, AnalysisCore& core);

    std::tuple<bool, bool> run();
    
    std::string formatProgress() const;
};

}

#endif

