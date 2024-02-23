#include "devbuiltins.h"
#include "devvisrc.h"
#include "devctlsrc.h"
#include "common.h"

namespace NAMESPACE {

void createBuiltins(std::vector<Device*>& devices) {
    devices.push_back(new BuiltinVSource("vsource"));
    devices.push_back(new BuiltinISource("isource"));
    devices.push_back(new BuiltinVccs("vccs"));
    devices.push_back(new BuiltinVcvs("vcvs"));
    devices.push_back(new BuiltinCccs("cccs"));
    devices.push_back(new BuiltinCcvs("ccvs"));
    devices.push_back(new BuiltinMutual("mutual"));
}

}
