#ifndef __HOMOTOPY_DEFINED
#define __HOMOTOPY_DEFINED

#include <string>
#include "circuit.h"
#include "core.h"
#include "common.h"

namespace NAMESPACE {

// TODO: Bypass can cause problems with homotopy because small steps taken by homotopy
//       can be mistaken for instances converging and consequently being bypassed. 
//       This can lead to homotopy failure. 

class Homotopy {
public:
    Homotopy(Circuit& circuit, AnalysisCore& core)
        : circuit(circuit), core(core) {};

    Homotopy           (const Homotopy&)  = delete;
    Homotopy           (      Homotopy&&) = delete;
    Homotopy& operator=(const Homotopy&)  = delete;
    Homotopy& operator=(      Homotopy&&) = delete;

    // Run homotopy
    // Return value: coverged, abort
    virtual std::tuple<bool, bool> run() { return std::make_tuple(false, false); };
    
    // Format progress message
    virtual std::string formatProgress() const { return "Unknown homotopy"; };

    // Number of steps performed 
    Int stepCount() const { return itCount; };

    static Id gdev;
    static Id gshunt;
    static Id spice3Gmin;
    static Id src;
    static Id spice3Src;

protected:
    Circuit& circuit;
    AnalysisCore& core;
    Int itCount;
};

}

#endif