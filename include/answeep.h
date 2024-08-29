#ifndef __ANSWEEP_DEFINED
#define __ANSWEEP_DEFINED

#include "value.h"
#include "flags.h"
#include "identifier.h"
#include "options.h"
#include "parseroutput.h"
#include "status.h"
#include "progress.h"
#include "common.h"
#include <memory>


namespace NAMESPACE {

// Sweep settings
typedef struct SweepSettings  {
    Id name {Id()};          // Name of sweep, will appear as a vector 
                             // holding the swept values in the raw file
    Loc location;            // Not exposed 
    Id instance {Id()};      // Name of the instance if sweeping an instance parameter
    Id model {Id()};         // Name of the model if sweeping a model parameter
    Id parameter {Id()};     // Name of the parameter to sweep
                             // if instance/model are not given, sweeps a toplevel
                             // instance parameter
    Id option {Id()};        // Name of simulator option to sweep
    Id variable {Id()};      // Name of circuit variable to sweep
    Int component {-1};      // Not supported yet. 
    Real from {0};           // Starting point
    Real to {0};             // End point
    Real step {0};           // Step size for linear stepped sweep (when mode not given)
    Id mode {Id()};          // Sweep with given number of points 
                             // can be lin (linear), dec (points per decade), 
                             // or oct (points per octave)
    Int points {0};          // Number of sweep points when mode is given
    Value values {0};        // Sweep given values (must be a vector or a list), 
                             // used when from, to, step, mode, and points are not given
    Int continuation {1};    // Use continuation mode for speeding up the sweep 
                             // (1=enabled, 0=disabled), enabled by default
    
    SweepSettings();
    
} SweepSettings;


// TODO: sweeping of vector parameters, should work only with values= which should be a list
class ScalarSweep {
public:
    enum class SweepType {
        Stepped, Lin, Log, Value
    };

    ScalarSweep();

    // Index corresponding to point at which the sweep is now (0-based)
    Int valueIndex() const;

    // Reset index to 0
    void reset();

    // Return position (0-based)
    Int at() const;

    // Return number of values
    Int count() const;

    // Advance index by 1, returns true when sweep is exhausted
    bool advance();

    // Computes the sweep value
    bool compute(Value& v, Status& s=Status::ignore) const;

    // Format progress
    std::string progress() const;

    // Set up a sweep
    bool setupSteppedSweep(Real from_, Real to_, Real step_, Status& s=Status::ignore);
    bool setupValueSweep(const Value& values, Status& s=Status::ignore);
    bool setupLinearSweep(Real from_, Real to_, Int points, Status& s=Status::ignore); 
    bool setupLogSweep(Real from_, Real to_, Real factor_, Int pointsPerFactor, Status& s=Status::ignore); 

    // Set up a scalar sweep based on settings structure
    template<typename A> bool setup(const A& settings, Status& s=Status::ignore);

protected: 
    // Common fields
    Int at_;
    Int end;
    SweepType sweepType;

    // For range sweep (stepped, lin, dec, oct)
    Real from;
    Real to;

    // For stepped sweep
    Real step;

    // For value sweep
    const Value* vals;

    // For log sweep
    Real factor;

private:
    static Id modeLin;
    static Id modeDec;
    static Id modeOct;
};

template<typename A> bool ScalarSweep::setup(const A& settings, Status& s) {
    // Check if any sweep is pecified
    int specCount=0;
    if (settings.values.isVector()) {
        specCount++;
    }
    if (settings.mode) {
        specCount++;
    }
    if (settings.step!=0) {
        specCount++;
    }
    if (specCount>1) {
        s.set(Status::Conflicting, "Sweep needs to specify only one of the following: values, mode, step.");
        return false;
    }

    if (settings.values.isVector()) {
        return setupValueSweep(settings.values, s);
    } else if (settings.mode) {
        if (settings.mode==ScalarSweep::modeLin) {
            return setupLinearSweep(settings.from, settings.to, settings.points, s);
        } else if (settings.mode==ScalarSweep::modeDec) {
            return setupLogSweep(settings.from, settings.to, 10, settings.points, s);
        } else if (settings.mode==ScalarSweep::modeOct) {
            return setupLogSweep(settings.from, settings.to, 2, settings.points, s);
        } else {
            s.set(Status::BadArguments, "Unknown sweep mode.");
            return false;
        }
    } else if (settings.step!=0) {
        return setupSteppedSweep(settings.from, settings.to, settings.step, s);
    }
    s.set(Status::NotFound, "Sweep needs to specify values, mode, or step.");
    return false;
}


// Order method invocation:
//   bind()
//   storeState()
//   reset()
//   repeat
//     write() sweep
//     invoke analysis
//     advance()
//   write() stored state

class Circuit;

class ParameterSweeper : public ProgressTracker {
public:
    enum class WriteValues { StoredState, Sweep };
    enum class ParameterFamily { Instance=1<<0, Model=1<<1, Option=1<<2, Variable=1<<3 };
    
    ParameterSweeper(Circuit& circuit, const PTSweeps& ptSweeps);

    // Setup sweeper (evaluate expressions, fill settings structures)
    bool setup(Status& s=Status::ignore);

    // Update
    bool update(int advancedSweepIndex, Status& s=Status::ignore);

    // Number of sweeps
    int count() const { return settings.size(); };

    // Bind (lookup instances, models, and simulator options)
    bool bind(Circuit& circuit, IStruct<SimulatorOptions>& opt, Status& s=Status::ignore);

    // Store parameters corresponding to current circuit state
    bool storeState(Status& s=Status::ignore);

    // Reset
    void reset();

    // Does i-th sweep use continuation
    bool continuation(int i) { return settings[i].continuation; };

    // Advance, return value: sweep done, index of sweep that was incremented (resets do not count)
    std::tuple<bool, Int> advance();

    // Position of the innermost sweep
    Int innermostSweepPosition() const { return scalarSweeps.back().at(); }; 

    // Format progress
    std::string progress() const;

    // Write swept parameters or stored parameters to circuit
    // Options are written to opt structure
    // Return value: ok, at least one instance or model parameter changed
    std::tuple<bool, bool> write(ParameterFamily types, WriteValues what, Status& s=Status::ignore);

    // Sweep name
    Id sweepName(Int ndx) const;
    
    // Return current index
    Int valueIndex(Int ndx) const;

    // Compute current value
    bool compute(Int ndx, Value& v, Status& s=Status::ignore) const;

private:
    Circuit& circuit;
    const PTSweeps& ptSweeps;
    std::vector<SweepSettings> settings;
    std::vector<ScalarSweep> scalarSweeps;
    std::vector<ParameterFamily> parameterFamily;
    std::vector<Parameterized*> parameterizedObject;
    std::vector<ParameterIndex> parameterIndex;
    std::vector<Value> storedValues;
    Int incrementedSweepIndex;
    Circuit* circuit_;
    size_t sweepPos;
};
DEFINE_FLAG_OPERATORS(ParameterSweeper::ParameterFamily);

}

#endif
