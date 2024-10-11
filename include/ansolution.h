#ifndef __ANSOLUTION_DEFINED
#define __ANSOLUTION_DEFINED

#include "ansupport.h"
#include "node.h"
#include "value.h"
#include "common.h"


namespace NAMESPACE {

class Circuit;

class AnnotatedSolution {
public:
    AnnotatedSolution();

    AnnotatedSolution           (const AnnotatedSolution&)  = delete;
    AnnotatedSolution           (      AnnotatedSolution&&) = default;
    AnnotatedSolution& operator=(const AnnotatedSolution&)  = delete;
    AnnotatedSolution& operator=(      AnnotatedSolution&&) = default;

    
    // Actual data
    const Vector<double>& values() const { return std::get<std::vector<double>>(values_); };
    const Vector<Complex>& cxValues() const { return std::get<std::vector<Complex>>(values_); };
    void setValues(const Vector<double>& vec) { values_ = vec; };
    void setCxValues(const Vector<Complex>& vec) { values_ = vec; };
    
    // Names of unknowns
    const std::vector<Id>& names() const { return names_; };
    void setNames(Circuit& circuit);
    void clearNames() { names_.clear(); };

    // States (DC), frequencies (HB)
    const Vector<double>& auxData() const { return auxData_; };
    void setAuxData(const Vector<double>& vec) { auxData_ = vec; };
    
private:
    typedef std::variant<Vector<double>, Vector<Complex>> VectorVariant;
    // Solution vector
    // - dc: one real component per unknown, index 0 is ground (bucket)
    // - hb: nf complex components per unknown. 
    //       no bucket - index 0 is first unknown
    VectorVariant values_;
    
    // Names of unknowns for cross matching across slightly different circuits
    std::vector<Id> names_;

    // Vector of auxiliary data
    // - hb: list of frequencies including DC (first component)
    Vector<double> auxData_;
};

}

#endif
