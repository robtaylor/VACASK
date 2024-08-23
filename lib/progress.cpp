#include "progress.h"
#include "common.h"


namespace NAMESPACE {

ProgressReporter::ProgressReporter(ParameterSweeper* s, AnalysisCore* c) 
    : sweeper(s), core(c) {
    tLast = Accounting::wclk();
};

}