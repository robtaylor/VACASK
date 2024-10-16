#ifndef __PROCESS_DEFINED
#define __PROCESS_DEFINED

#include "common.h"
#include <tuple>
#include <string>
#include <vector>
#include "status.h"


namespace NAMESPACE {

bool findProgram(const std::string& prog, std::string& path);

std::tuple<bool, std::string, std::string> runProcess(
    const std::string& prog, const std::vector<std::string>& args, const std::string* pythonPath, bool collect=true, bool debugFiles=false, Status& s=Status::ignore
);

std::string executableFile();

}

#endif

