#ifndef __PROGRESS_DEFINED
#define __PROGRESS_DEFINED

#include "core.h"
#include "answeep.h"
#include "acct.h"
#include "common.h"

namespace NAMESPACE {

class ProgressReporter {
public:
    ProgressReporter(ParameterSweeper* s, AnalysisCore* c);

    virtual void begin() = 0;
    virtual void reportProgress(bool force=false) = 0;
    virtual void end() = 0;

protected:
    Accounting::Timepoint tLast;
    ParameterSweeper* sweeper;
    AnalysisCore* core;
};

}

#endif
