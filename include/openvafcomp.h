#ifndef __SRCCOMPILER_DEFINED
#define __SRCCOMPILER_DEFINED

#include <string>
#include "circuit.h"
#include "common.h"

namespace NAMESPACE {

class OpenvafCompiler : public SourceCompiler {
public: 
    OpenvafCompiler(std::optional<const std::string> compiler=std::nullopt, std::optional<const std::vector<std::string>> compilerArgs=std::nullopt)
        : compiler(compiler), compilerArgs(compilerArgs) {};

    virtual std::tuple<bool, bool> compile(const std::string& loadDirectiveCanonicalPath, const std::string& fileName, const std::string& canonicalPath, std::string& outputCanonicalPath, Status& s=Status::ignore);

private:
    std::optional<const std::string> compiler;
    std::optional<const std::vector<std::string>> compilerArgs;
};

}

#endif
