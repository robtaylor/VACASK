#ifndef __PROGRESS_DEFINED
#define __PROGRESS_DEFINED

#include <tuple>
#include <optional>
#include "acct.h"
#include "common.h"

namespace NAMESPACE {

// Objects of this class render a progress report
// Information on progress is stored here (extent, position). 
// The optional displayed value is also stored here. 
// A ProgressTracker object (sweep or analysis core) sets these
// values via the initProgress() and setProgress() methods. 
// Progress output is triggered by calling the report() method. 
class ProgressReporter {
public:
    enum class ValueFormat { Fixed, Scientific, Default };

    ProgressReporter(double dt);

    // For disabling in debug mode
    void disable() { enabled_ = false; };

    // For checking if the reporter was disabled
    bool enabled() const { return enabled_; };
    
    // Setting how a value should be displayed
    void setValueFormat(ValueFormat f, int precision=-1) { valueFormat_ = f; precision_ = precision; };
    void setValueDecoration(const std::string& prefix, const std::string& postfix) { prefix_ = prefix; postfix_ = postfix; };
    
    // Start reporting
    void begin() { tBegin = Accounting::wclk(); tLast = tBegin; };

    // Initialize progress (set extent and position)
    void initProgress(double extent, double pos, std::optional<double> value = std::nullopt) { extent_ = extent, pos_ = pos; value_ = value; };
    
    // Set progress
    void setProgress(double pos, std::optional<double> value = std::nullopt) { pos_ = pos; value_ = value; };

    // Return position
    double position() const { return pos_; }; 
    
    // Report if it is time
    bool report(bool force=false) {
        if (force || timeForReport()) {
            updateReportTime();
            return reporter();
        }
        return false;
    };

    // End of reporting
    void end() { updateReportTime(); };

    // Override in derived class to implement reporting
    virtual bool reporter() = 0;
    
    // Time between last update and begin()
    double time() { 
        return std::chrono::duration_cast<std::chrono::duration<double>>(tLast-tBegin).count(); 
    };

protected:
    Accounting::Timepoint tBegin;
    Accounting::Timepoint tLast;
    double dt;
    bool enabled_;

    bool timeForReport() { return Accounting::wclkDelta(tLast)>=dt; };
    void updateReportTime() { tLast = Accounting::wclk(); }

    ValueFormat valueFormat_;
    int precision_;
    std::string prefix_;
    std::string postfix_;

    double extent_;
    double pos_;
    std::optional<double> value_;
};


// Objects of this class track their progress. Usually this class is inherited
// by analysis cores and the multidimensional sweeper. Classes that inherit
// ProgressTracker use its initProgress() and setProgress() methods to set 
// the progress and trigger a progress report output. 
class ProgressTracker {
public:
    ProgressTracker() : progressReporter(nullptr) {};
    void install(ProgressReporter* p) { progressReporter = p; };
    // Default extent (1), position (0)
    void initProgress() { 
        if (!progressReporter) {
            return;
        }
        progressReporter->initProgress(1, 0);
        progressReporter->report();
    };
    // Set extent and position
    void initProgress(double extent, double pos, std::optional<double> value = std::nullopt) { 
        if (!progressReporter) {
            return;
        }
        progressReporter->initProgress(extent, pos, value);
        progressReporter->report(); 
    };
    // Set position
    void setProgress(double pos, std::optional<double> value = std::nullopt) {
        if (!progressReporter) {
            return;
        }
        progressReporter->setProgress(pos, value);
        progressReporter->report(); 
    }
    
protected:
    ProgressReporter* progressReporter;
};

}

#endif
