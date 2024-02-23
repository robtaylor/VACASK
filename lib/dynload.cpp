#include "dynload.h"
#include "common.h"

#ifdef SIMWINDOWS
#include <windows.h> 
#else
#include <dlfcn.h>
#endif

namespace NAMESPACE {
#ifdef SIMWINDOWS
void* openDynamicLibrary(const char* fileName) {
    return LoadLibrary(fileName);
}

bool closeDynamicLibrary(void* handle) {
    if (FreeLibrary((HMODULE)handle)) {
        return true;
    } else {
        return false;
    }
}

void* dynamicLibrarySymbol(void* handle, const char* name) {
    return (void*)GetProcAddress((HMODULE)handle, name); 
}

std::string dynamicLibraryError() {
    DWORD errCode = GetLastError();
    
    LPVOID lpMsgBuf = nullptr;
    DWORD retval = FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, 
        NULL, errCode, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR)&lpMsgBuf, 0, NULL
    );

    std::string errMsg;
    if (retval == 0) { 
        errMsg = std::string("Failed to format error message with code ")+std::to_string(errCode)+".";
    } else {
        // TODO: maybe convert from wchar_t* to utf8
        errMsg = std::string(static_cast<char*>(lpMsgBuf)); 
    }
    if (lpMsgBuf) {
        LocalFree(lpMsgBuf);
    }

    return errMsg;
  }
#else
void* openDynamicLibrary(const char* fileName) {
    return dlopen(fileName, RTLD_NOW);
}

bool closeDynamicLibrary(void* handle) {
    if (dlclose(handle)==0) {
        return true;
    } else {
        return false;
    }
}

void* dynamicLibrarySymbol(void* handle, const char* name) {
    return dlsym(handle, name);
}

std::string dynamicLibraryError() {
    auto err = dlerror();
    if (err) {
        return err;
    } else {
        return "";
    }
}
#endif

}
