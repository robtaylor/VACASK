#ifndef __OUTPUT_DEFINED
#define __OUTPUT_DEFINED

#include "identifier.h"
#include "value.h"
#include "ansupport.h"
#include "answeep.h"
#include <complex>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <complex>
#include "common.h"


namespace NAMESPACE {

// One save directive -> zero, one, or more output descriptors
// Conversion is taken care of by analysis. 

// One output descriptor -> one output source
// Conversion is taken care of by analysis. 

// Bound data source, contains full information for extracting a scalar value
// This must be efficient as there will be many of these
//   char (1) + ptr_int (8) + ptr (8) = 17 (24 bytes)
class OutputSource {
public:
    enum class Type : char { 
        // Real value
        IntPtr, 
        RealPtr, 
        RealArray, 
        RealVec, 
        RealVecRepo, 

        // Complex value
        ComplexPtr, 
        ComplexArray, 
        ComplexVec, 
        ComplexVecRepo, 

        // Sweep variable
        Sweep 
    };
    typedef size_t Index;
    
    OutputSource() : cPtr(&defaultBad), type(Type::ComplexPtr) {};

    OutputSource(const int* ptr) : iPtr(ptr), type(Type::IntPtr) {};

    OutputSource(const double* ptr) : rPtr(ptr), type(Type::RealPtr) {};
    OutputSource(const double* ptr, Index ndx) : rPtr(ptr), index(ndx), type(Type::RealArray) {};
    OutputSource(const Vector<double>* ptr, Index ndx) : rVec(ptr), index(ndx), type(Type::RealVec) {};
    OutputSource(const VectorRepository<double>* ptr, Index ndx) 
        : rRepoPtr(ptr), index(ndx), type(Type::RealVecRepo) {};
    
    OutputSource(const Complex* ptr) : cPtr(ptr), type(Type::ComplexPtr) {};
    OutputSource(const Complex* ptr, Index ndx) : cPtr(ptr), index(ndx), type(Type::ComplexArray) {};
    OutputSource(const Vector<Complex>* ptr, Index ndx) : cVec(ptr), index(ndx), type(Type::ComplexVec) {};
    OutputSource(const VectorRepository<Complex>* ptr, Index ndx) 
        : cRepoPtr(ptr), index(ndx), type(Type::ComplexVecRepo) {};
    
    OutputSource(const ParameterSweeper* ptr, Index ndx) : sweeperPtr(ptr), index(ndx), type(Type::Sweep) {};
    
    double getR() const {
        switch (type) {
            // Real value -> real value
            case Type::IntPtr:
                return *iPtr;
            case Type::RealPtr:
                return *rPtr;
            case Type::RealArray:
                return rPtr[index];
            case Type::RealVec:
                return (*rVec)[index];
            case Type::RealVecRepo:
                // Valid solution of NR algorithm is always the current one
                // Latest opvars are computed at this solution
                return rRepoPtr->data()[index]; 
            
            // Complex value -> real value (return real part)
            case Type::ComplexPtr:
                return cPtr->real();
            case Type::ComplexArray:
                return cPtr[index].real();
            case Type::ComplexVec:
                return (*cVec)[index].real();
            case Type::ComplexVecRepo:
                // Return current solution from repository
                return cRepoPtr->data()[index].real();
            
            // Sweep variable (real) -> real value
            case Type::Sweep:
                Value v;
                if (sweeperPtr->compute(index, v)) {
                    if (v.type()==Value::Type::Int) {
                        return v.val<Int>();
                    } else if (v.type()==Value::Type::Real) {
                        return v.val<Real>();
                    } else {
                        // Return index of value in all other cases 
                        // that cannot be converted to a real scalar 
                        return sweeperPtr->valueIndex(index);
                    }
                } else {
                    // Sweep value computation failed. Should never end up here because anaylsis will abort in this case
                    return 0.0;
                }
        }
        // Should never be reached
        return 0.0;
    };
    
    Complex getC() const {
        switch (type) {
            // Real value -> complex value
            case Type::IntPtr:
                return *iPtr;
            case Type::RealPtr:
                return *rPtr;
            case Type::RealArray:
                return rPtr[index];
            case Type::RealVec:
                return (*rVec)[index];
            case Type::RealVecRepo:
                // Valid solution of NR algorithm is always the current one
                // Latest opvars are computed at this solution
                return rRepoPtr->data()[index]; 
            
            // Complex value -> complex value
            case Type::ComplexPtr:
                return *cPtr;
            case Type::ComplexArray:
                return cPtr[index];
            case Type::ComplexVec:
                return (*cVec)[index];
            case Type::ComplexVecRepo:
                // Return current solution from repository
                return cRepoPtr->data()[index];
            
            // Sweep variable (real) -> complex value
            case Type::Sweep:
                Value v;
                if (sweeperPtr->compute(index, v)) {
                    if (v.type()==Value::Type::Int) {
                        return v.val<Int>();
                    } else if (v.type()==Value::Type::Real) {
                        return v.val<Real>();
                    } else {
                        // Return index of value in all other cases 
                        // that cannot be converted to a real scalar 
                        return sweeperPtr->valueIndex(index);
                    }
                } else {
                    // Sweep value computation failed. Should never end up here because anaylsis will abort in this case
                    return 0.0;
                }
        }
        // Should never be reached
        return 0.0;
    };

private:
    union {
        const int* iPtr;
        
        const double* rPtr;
        const Vector<double>* rVec;
        const VectorRepository<double>* rRepoPtr;

        const Complex* cPtr;
        const Vector<Complex>* cVec;
        const VectorRepository<Complex>* cRepoPtr;
        
        const ParameterSweeper* sweeperPtr;
    };
    Index index;
    Type type;
    static Complex defaultBad;
};


// Output descriptor 
// This must be efficient as there will be many of these
//   type (1) + Id (4) + Id or ndx (4)+ Id (4) = 13 (16 bytes)

// Abstract analysis stores descriptors vector and output sources vector. 
// Descriptors are created only once, just before core analysis is run for the first time. 
// Conversion to output sources follows, can happen more than once (when required). 
// This structure only holds the basic information for creating an output source. 
// Its interpretation depends on the analysis. 
// The same output descriptor can result in different output sources in two different analyses. 
// 
// Type                 Content      Coverted to output source by
// solution component       id           analysis
// opvar                    id id        analysis (via circuit)
// sweep variable           ndx          abstract analysis
// simulator internal       id           abstract analysis (via circuit)
//   (time, freq, ...)    
// tf                       id (src)     analysis
// input impedance          id (src)     analysis
// noise contribution       id id        analysis
// output noise                          analysis
// noise contribution gain  id id        analysis
// noise gain                            analysis

struct OutputDescriptor {
    using Type = uint8_t;
    using NdxType = uint32_t;
    static constexpr Type CommonMask() { return 2<<6; };
    static constexpr Type CustomMask() { return 1<<6; }; 

    OutputDescriptor(Type t, Id outputName);
    OutputDescriptor(Type t, Id outputName, size_t n);
    OutputDescriptor(Type t, Id outputName, Id id); 
    OutputDescriptor(Type t, Id outputName, Id id1, Id id2); 
    OutputDescriptor(Type t, Id outputName, Id id, size_t n); 

    Type type;
    union {
        struct {
            Id id;
            NdxType ndx;
        } idNdx;
        struct {
            Id id1;
            Id id2;
        } idId;
        Id id;
        NdxType ndx;
    };
    Id name;
};

// Common descriptor types (64 possible values)
const OutputDescriptor::Type OutdSweepvar    = OutputDescriptor::CommonMask() | 0;
const OutputDescriptor::Type OutdSimInternal = OutputDescriptor::CommonMask() | 1;
const OutputDescriptor::Type OutdSimStat     = OutputDescriptor::CommonMask() | 2;

// Descriptor types used by several analyses (64 possible values)
const OutputDescriptor::Type OutdSolComponent             = 0;
const OutputDescriptor::Type OutdPinCurrent               = 1; 
const OutputDescriptor::Type OutdOpvar                    = 2;
const OutputDescriptor::Type OutdTf                       = 3;
const OutputDescriptor::Type OutdZin                      = 4;
const OutputDescriptor::Type OutdYin                      = 5;
const OutputDescriptor::Type OutdNoiseContribInst         = 6;
const OutputDescriptor::Type OutdNoiseContribInstPartial  = 7;

const OutputDescriptor::Type OutdOutputNoise              = 60;
const OutputDescriptor::Type OutdPowerGain                = 61;
const OutputDescriptor::Type OutdFrequency                = 62;
const OutputDescriptor::Type OutdTime                     = 63;


// Common interface of all output classes
class Output {
public:
    using DescriptorList = std::vector<OutputDescriptor>;
    using SourcesList = std::vector<OutputSource>;

    Output(const std::string& baseName, DescriptorList& descriptors, SourcesList& sources);
    
    void setTitle(std::string title) { title_= title; };
    void setPlotname(std::string name) { plotname_= name; };

    // Open file
    virtual bool prologue(Status& s=Status::ignore) = 0;
    
    // Add point
    virtual bool addPoint(Status& s=Status::ignore) = 0;

    // FInalize and close file
    virtual bool epilogue(Status& s=Status::ignore) = 0;

    // Remove file
    virtual bool remove(Status& s=Status::ignore) = 0;

protected:
    std::string baseName_;
    std::ofstream outStream;
    std::string title_;
    std::string plotname_;
    DescriptorList& descrs;
    SourcesList& srcs;
};

}

#endif
