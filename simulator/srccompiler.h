#ifndef __SRCCOMPILER_DEFINED
#define __SRCCOMPILER_DEFINED

#include <string>
#include "circuit.h"
#include "common.h"

namespace NAMESPACE {

class OpenvafCompiler : public SourceCompiler {
public: 
    virtual std::tuple<bool, bool> compile(const std::string& loadDirectiveCanonicalPath, const std::string& fileName, const std::string& canonicalPath, std::string& outputCanonicalPath, Status& s=Status::ignore);
};

}

#endif
