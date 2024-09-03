#include <vector>
#include <algorithm>
#include "core.h"
#include "corehb.h"
#include "simulator.h"
#include "common.h"

namespace NAMESPACE {

Id HbCore::truncateRaw = Id::createStatic("raw");
Id HbCore::truncateBox = Id::createStatic("box");
Id HbCore::truncateDiamond = Id::createStatic("diamond");

Id HbCore::sampleUniform = Id::createStatic("uniform");
Id HbCore::sampleRandom = Id::createStatic("random");

HbParameters::HbParameters() {
    truncate = HbCore::truncateDiamond;
    sample = HbCore::sampleRandom;
}

}
