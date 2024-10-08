#include "homotopy.h"
#include "hmtpgmin.h"
#include "hmtpsrc.h"
#include "common.h"

namespace NAMESPACE {

Id Homotopy::gdev = Id::createStatic("gdev");
Id Homotopy::gshunt = Id::createStatic("gshunt");
Id Homotopy::spice3Gmin = Id::createStatic("spice3gmin");
Id Homotopy::src = Id::createStatic("src");
Id Homotopy::spice3Src = Id::createStatic("spice3src");

}
