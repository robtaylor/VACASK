#ifndef __HMTPSRC_DEFINED
#define __HMTPSRC_DEFINED

#include "homotopy.h"
#include "common.h"

namespace NAMESPACE {

// Based on Ngspice gillespie_src()
class SourceStepping : public Homotopy {
public:
    SourceStepping(Circuit& circuit, AnalysisCore& core, const std::vector<Id>& alternativeHomotopy={});

    SourceStepping           (const SourceStepping&)  = delete;
    SourceStepping           (      SourceStepping&&) = delete;
    SourceStepping& operator=(const SourceStepping&)  = delete;
    SourceStepping& operator=(      SourceStepping&&) = delete;

    std::tuple<bool, bool> run();
    
    std::string formatProgress() const;

private:
    const std::vector<Id>& alternativeHomotopy;
};


// Based on Ngspice spice3_src()
class Spice3SourceStepping : public Homotopy {
public:
    Spice3SourceStepping(Circuit& circuit, AnalysisCore& core);

    Spice3SourceStepping           (const Spice3SourceStepping&)  = delete;
    Spice3SourceStepping           (      Spice3SourceStepping&&) = delete;
    Spice3SourceStepping& operator=(const Spice3SourceStepping&)  = delete;
    Spice3SourceStepping& operator=(      Spice3SourceStepping&&) = delete;

    std::tuple<bool, bool> run();
    
    std::string formatProgress() const;
};

}

#endif
