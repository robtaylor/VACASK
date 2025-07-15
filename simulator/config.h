#ifndef __CONFIG_DEFINED
#define __CONFIG_DEFINED

#include "status.h"
#include "common.h"

namespace NAMESPACE {

bool readConfig(std::ifstream& input, const std::string& filename, Status& s=Status::ignore);

}

#endif
