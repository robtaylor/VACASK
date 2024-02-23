#ifndef __OUTRAWFILE_DEFINED
#define __OUTRAWFILE_DEFINED

#include "value.h"
#include "identifier.h"
#include "circuit.h"
#include "output.h"
#include "flags.h"
#include <complex>
#include <vector>
#include <string>
#include <iostream>
#include "common.h"


namespace NAMESPACE {

enum class RawfileFlags {
    None = 0, 
    Binary = 1,  // ASCII otherwise
    Complex = 2, // Real otherwise
    Padded = 4,  // Unpadded otherwise
};
DEFINE_FLAG_OPERATORS(RawfileFlags);

class OutputRawfile : public Output, public FlagBase<RawfileFlags> {
public:
    OutputRawfile(const std::string& baseName, Output::DescriptorList& descriptors, Output::SourcesList& sources, Flags f=Flags::Binary|Flags::Padded);
    ~OutputRawfile();

    virtual bool prologue(Status& s=Status::ignore);
    virtual bool addPoint(Status& s=Status::ignore);
    virtual bool epilogue(Status& s=Status::ignore);
    virtual bool remove(Status& s=Status::ignore);

private:
    std::streampos pointCountPos;
    size_t count;
    std::string fileName;
};

}

#endif
