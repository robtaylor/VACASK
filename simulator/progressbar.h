#ifndef __PROGRESSBAR_DEFINED
#define __PROGRESSBAR_DEFINED

#include <string>
#include <iomanip>
#include "core.h"
#include "answeep.h"
#include "progress.h"
#include "common.h"

namespace NAMESPACE { 

class AnalysisProgress : public ProgressReporter {
public:
    AnalysisProgress(int indent, std::ostream& os, double dt=0.25);
    
    virtual bool reporter();

private:
    std::string renderBar(double norm, double val);

    std::ostream& os;
    int indent;
    int columns; 
    int barSize;
    
    std::string lastStr;
    size_t lastLen;

    Accounting::Timepoint tBegin;
};


}

#endif
