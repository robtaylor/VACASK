#ifndef __EVALDUMP_DEFINED
#define __EVALDUMP_DEFINED

#include <fstream>
#include <set>
#include <string>
#include "common.h"

//
// Spike instrumentation (GPU device-eval f32-accuracy validation).
//
// EvalDumpSink is a lightweight, flag-gated sink for per-instance device-eval
// I/O. The transient core creates one when VACASK_EVAL_DUMP names an output
// file, then runs a forced (no-bypass) eval-only pass at selected accepted
// timepoints. OsdiInstance::dumpEvalIO() writes, per instance, the input
// voltages and the f64 resistive residual + Jacobian computed by the OSDI
// model. Those f64 values are the ground truth the generated f32 Metal kernel
// is validated against (see docs/spikes/openvaf-gpu-codegen.md, Q2 Step 4).
//
// Not part of normal simulation: the hot NR path only ever sees a null
// EvalSetup::evalDump pointer (one branch). Gated by an env var so no options
// or parser changes are needed; remove with the rest of the spike scaffolding.
//

namespace NAMESPACE {

class EvalDumpSink {
public:
    explicit EvalDumpSink(const std::string& path);

    bool good() const { return file_.is_open(); }
    std::ofstream& os() { return file_; }

    // Records the descriptor name; returns true the first time it is seen so
    // the caller emits the static MODEL header block exactly once per model.
    bool firstSeen(const std::string& name) { return seen_.insert(name).second; }

    void setPoint(double t, long step) { time_ = t; step_ = step; }
    double time() const { return time_; }
    long step() const { return step_; }

private:
    std::ofstream file_;
    std::set<std::string> seen_;
    double time_ {0};
    long step_ {0};
};

}

#endif
