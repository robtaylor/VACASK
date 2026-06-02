#include "evaldump.h"

namespace NAMESPACE {

EvalDumpSink::EvalDumpSink(const std::string& path) : file_(path) {
    if (file_.is_open()) {
        // Self-describing header; harness keys on this version line.
        file_ << "# VACASK eval dump v1\n";
        file_ << "# Per-instance f64 device-eval I/O at accepted transient timepoints.\n";
        file_ << "# MODEL blocks (static structure) precede the INST blocks that reference them.\n";
    }
}

}
