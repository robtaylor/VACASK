#ifndef __LIBPLATFORM_DEFINED
#define __LIBPLATFORM_DEFINED

#include <vector>
#include <string>
#include <filesystem>
#include "common.h"

namespace NAMESPACE {

const std::string& pathSeparator();

size_t splitString(const std::string& separator, const std::string s, std::vector<std::string>& stringVec); 

class LibPlatform {
public:
    static const std::vector<std::string>& systemPath();
};

// std::aligned_alloc() is not available in MSVC
void* alignedAlloc(size_t alignment, size_t size);

// MSVC requires us to use special free for aligned allocations
void alignedFree(void* ptr);

// MSVC does not allow std::isinf(), isstd::isnan(), and std::isfinite() on non-floating point arguments
template<typename T> bool isInf(T x) {
    if constexpr(std::is_floating_point<T>::value) {
        return std::isinf(x);
    } else {
        return false;
    }
}

template<typename T> bool isNaN(T x) {
    if constexpr(std::is_floating_point<T>::value) {
        return std::isnan(x);
    } else {
        return false;
    }
}

template<typename T> bool isFinite(T x) {
    if constexpr(std::is_floating_point<T>::value) {
        return std::isfinite(x);
    } else {
        return false;
    }
}

std::string findPythonExecutable();
const char* defaultOpenVafBinaryName();

}

#endif
