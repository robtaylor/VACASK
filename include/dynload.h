#ifndef __DYNLOAD_DEFINED
#define __DYNLOAD_DEFINED

#include "common.h"

namespace NAMESPACE {

// Open library, return handle (nullptr on failure)
void* openDynamicLibrary(const char* fileName);

// Close library, return true on success
bool closeDynamicLibrary(void* handle);

// Get symbol, returns nullptr if not found
void* dynamicLibrarySymbol(void* handle, const char* name);

// Last error (empty string if no error)
std::string dynamicLibraryError();

}

#endif
