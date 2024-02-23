#include "libplatform.h"
#include "filesystem.h"
#include "processutils.h"
#include "common.h"
#include <tuple>
#include <cstdlib>

namespace NAMESPACE {

size_t splitString(const std::string& separator, const std::string s, std::vector<std::string>& stringVec) {
    // Empty string
    if (s.size()==0) {
        return 0;
    }

    // Go through string
    auto sepSize = separator.size();
    decltype(sepSize) start = 0;
    bool exit = false;
    size_t count = 0;
    do {
        auto pos = s.find(separator, start);
        decltype(pos) len;
        if (pos == std::string::npos) {
            exit = true;
            len = s.size() - start;
        } else {
            len = pos - start;
        }
        if (len>0) {
            stringVec.push_back(s.substr(start, len));
            count++;
        }
        start += len + separator.size();
    } while (!exit);

    return count;
}

const std::string& pathSeparator() {
    static const std::string separator =
#ifdef SIMWINDOWS
        ";";
#else
        ":";
#endif
    return separator;
}

static std::vector<std::string> initSystemPath() {
    std::vector<std::string> pathStrings;
    auto sysPath = std::getenv("PATH"); 
    static auto dummy = splitString(pathSeparator(), sysPath, pathStrings);
    return pathStrings;
}

const std::vector<std::string>& LibPlatform::systemPath() {
    static std::vector<std::string> path = initSystemPath();
    return path;
}

#ifdef SIMWINDOWS
void* alignedAlloc(size_t alignment, size_t size) {
    // MSVC does not have std::aligned_alloc()
    return _aligned_malloc(size, alignment);
}

void alignedFree(void* ptr) {
    _aligned_free(ptr);
}
#else
void* alignedAlloc(size_t alignment, size_t size) {
    return std::aligned_alloc(alignment, size);
}

void alignedFree(void* ptr) {
    free(ptr);
}
#endif

}
