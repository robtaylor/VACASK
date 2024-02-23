#ifndef __OSDICALLBACK_DEFINED
#define __OSDICALLBACK_DEFINED

#include "osdi.h"
#include "common.h"


namespace NAMESPACE {

// Callbacks
typedef struct OsdiCallbackHandle {
    uint32_t kind;
    char *name;
} OsdiCallbackHandle;

void osdiLogMessage(void *handle, char *msg, uint32_t level);

typedef struct OsdiLimitFunction {
    const char *name;
    int nArgs;
    void *ptr;
} OsdiLimitFunction;

}

#endif
