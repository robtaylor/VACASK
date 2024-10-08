#ifndef __HMTPSRC_DEFINED
#define __HMTPSRC_DEFINED

#include "homotopy.h"
#include "common.h"

namespace NAMESPACE {

// Based on Ngspice gillespie_src()
class SourceStepping : public Homotopy {
public:
    SourceStepping(Circuit& circuit, AnalysisCore& core, const std::vector<Id>& alternativeHomotopy={});

    std::tuple<bool, bool> run();
    
    std::string formatProgress() const;

private:
    const std::vector<Id>& alternativeHomotopy;
};


// Based on Ngspice spice3_src()
class Spice3SourceStepping : public Homotopy {
public:
    Spice3SourceStepping(Circuit& circuit, AnalysisCore& core);

    std::tuple<bool, bool> run();
    
    std::string formatProgress() const;
};

}

#endif
