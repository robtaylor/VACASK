#ifndef __ANSWEEP_DEFINED
#define __ANSWEEP_DEFINED

#include "value.h"
#include "flags.h"
#include "identifier.h"
#include "options.h"
#include "status.h"
#include "common.h"
#include <memory>


namespace NAMESPACE {

// TODO: sweeping of vector parameters, should work only with values= which should be a list
class ScalarSweep {
public:
    ScalarSweep();

    // Is sweep valid
    bool isValid() const;

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
    virtual bool compute(Value& v, Status& s=Status::ignore) const = 0;

    // Format progress
    std::string progress() const;

    // Factory function, uses any kind of struct that has all the neccessary components
    template<typename A> static ScalarSweep* create(const A& settings, Status& s=Status::ignore);

protected: 
    bool valid;
    Int at_;
    Int end;

private:
    static Id modeLin;
    static Id modeDec;
    static Id modeOct;
};


class SteppedScalarSweep : public ScalarSweep {
public:
    SteppedScalarSweep(Real from, Real to, Real step, Status& s=Status::ignore);

    virtual bool compute(Value& v, Status& s=Status::ignore) const;

protected:
    Real from;
    Real to;
    Real step;
};


class ValueScalarSweep : public ScalarSweep {
public:
    ValueScalarSweep(const Value& values, Status& s=Status::ignore);

    virtual bool compute(Value& v, Status& s=Status::ignore) const;

protected:
    const Value& vals;
};


class LinearScalarSweep : public ScalarSweep {
public:
    LinearScalarSweep(Real from, Real to, Int points, Status& s=Status::ignore);

    virtual bool compute(Value& v, Status& s=Status::ignore) const;

protected:
    Real from;
    Real to;
};


class LogScalarSweep : public ScalarSweep {
public:
    LogScalarSweep(Real from, Real to, Real factor, Int pointsPerFactor, Status& s=Status::ignore);

    virtual bool compute(Value& v, Status& s=Status::ignore) const;

protected:
    Real from;
    Real to;
};

template<typename A> ScalarSweep* ScalarSweep::create(const A& settings, Status& s) {
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
    if (specCount<=0) {
        s.set(Status::NotFound, "Sweep needs to specify values, mode, or step.");
        return nullptr;
    }
    if (specCount>1) {
        s.set(Status::Conflicting, "Sweep needs to specify only one of the following: values, mode, step.");
        return nullptr;
    }

    ScalarSweep* sweep = nullptr;
    if (settings.values.isVector()) {
        sweep = new ValueScalarSweep(settings.values, s);
        if (!sweep->isValid()) {
            delete sweep;
            return nullptr;
        }
    } else if (settings.mode) {
        if (settings.mode==ScalarSweep::modeLin) {
            sweep = new LinearScalarSweep(settings.from, settings.to, settings.points, s);
        } else if (settings.mode==ScalarSweep::modeDec) {
            sweep = new LogScalarSweep(settings.from, settings.to, 10, settings.points, s);
        } else if (settings.mode==ScalarSweep::modeOct) {
            sweep = new LogScalarSweep(settings.from, settings.to, 2, settings.points, s);
        } else {
            s.set(Status::BadArguments, "Unknown sweep mode.");
            return nullptr;
        }
        if (!sweep->isValid()) {
            delete sweep; 
            return nullptr;
        }
    } else if (settings.step!=0) {
        sweep = new SteppedScalarSweep(settings.from, settings.to, settings.step, s);
        if (!sweep->isValid()) {
            delete sweep;
            return nullptr;
        }
    }
    return sweep;
}


// Sweep settings
typedef struct SweepSettings  {
    Id name;
    Loc location;
    Id instance;
    Id model;
    Id parameter;
    Id option;
    Id variable;
    Int component;
    Real from;
    Real to;
    Real step;
    Id mode;
    Int points;
    Value values;
    Int continuation;
    
    SweepSettings();
    
} SweepSettings;


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

class ParameterSweeper {
public:
    enum class WriteValues { StoredState, Sweep };
    enum class ParameterFamily { Instance=1<<0, Model=1<<1, Option=1<<2, Variable=1<<3 };
    
    ParameterSweeper(const std::vector<SweepSettings>& settings, Status& s=Status::ignore);

    // Check if constructor was successfull
    bool isValid() const; 

    // Bind (lookup instances, models, and simulator options)
    bool bind(Circuit& circuit, IStruct<SimulatorOptions>& opt, Status& s=Status::ignore);

    // Store parameters corresponding to current circuit state
    bool storeState(Status& s=Status::ignore);

    // Reset
    void reset();

    // Advance, return value: sweep done, index of sweep that was incremented (resets do not count)
    std::tuple<bool, Int> advance();

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
    bool valid;
    const std::vector<SweepSettings>& settings;
    std::vector<std::unique_ptr<ScalarSweep>> scalarSweeps;
    std::vector<ParameterFamily> parameterFamily;
    std::vector<Parameterized*> parameterizedObject;
    std::vector<ParameterIndex> parameterIndex;
    std::vector<Value> storedValues;
    Int incrementedSweepIndex;
    Circuit* circuit_;
};
DEFINE_FLAG_OPERATORS(ParameterSweeper::ParameterFamily);

}

#endif
