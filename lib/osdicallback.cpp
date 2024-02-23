#include <stdlib.h>
#include <iostream>
#include "osdi.h"
#include "limitfunctions.h"
#include "osdicallback.h"
#include "osdifile.h"
#include "simulator.h"
#include "common.h"


namespace NAMESPACE {

// std::endl is slow because it flushes the stream, use \n

void osdiLogMessage(void *handle, char *msg, uint32_t level) {
    OsdiCallbackHandle *h = (OsdiCallbackHandle *)handle;
    std::ostream *dst = &Simulator::out();

    switch (level & LOG_LVL_MASK) {
        case LOG_LVL_DEBUG:
            dst = &Simulator::dbg();
            *dst << "OSDI(debug) " << h->name << ": ";
            break;
        case LOG_LVL_DISPLAY:
            *dst << "OSDI " << h->name << ": ";
            break;
        case LOG_LVL_INFO:
            *dst << "OSDI(info) " << h->name << ": ";
            break;
        case LOG_LVL_WARN:
            dst = &Simulator::wrn();
            *dst << "OSDI(warn) " << h->name << ": ";
            break;
        case LOG_LVL_ERR:
            dst = &Simulator::err();
            *dst << "OSDI(err) " << h->name << ": ";
            break;
        case LOG_LVL_FATAL:
            dst = &Simulator::err();
            *dst << "OSDI(fatal) " << h->name << ": ";
            break;
        default:
            dst = &Simulator::err();
            *dst << "OSDI(unknown) " << h->name << ": ";
            break;
    }

    if (level & LOG_FMT_ERR) {
        *dst << "failed to format \"" << msg << "\"\n";
        // Not allowed to free msg because it is a constant string
    } else {
        *dst << msg; 
    }
    // Free msg because we own it now, use libc free
    free(msg);
}

double osdiPnjlim(
    bool init, bool *check, double vnew, double vold, 
    double vt, double vcrit
) {
    if (init) {
        *check = true;
        return vcrit;
    }
    int icheck = 0;
    double res = DEVpnjlim(vnew, vold, vt, vcrit, &icheck);
    *check = icheck != 0;
    return res;
}

double osdiTypedpnjlim(
    bool init, bool *check, double vnew, double vold, 
    double vt, double vcrit, double type
) {
    if (init) {
        *check = true;
        return vcrit;
    }
    int icheck = 0;
    double res = DEVpnjlim(type*vnew, vold, vt, vcrit, &icheck);
    *check = icheck != 0;
    return res;
}

double osdiLimvds(bool init, bool *check, double vnew, double vold) {
    if (init) {
        *check = true;
        return 0.1;
    }
    double res = DEVlimvds(vnew, vold);
    if (res != vnew) {
        *check = true;
    }
    return res;
}

double osdiFetlim(
    bool init, bool *check, double vnew, 
    double vold, double vto
) {
    if (init) {
        *check = true;
        return vto + 0.1;
    }
    double res = DEVfetlim(vnew, vold, vto);
    if (res != vnew) {
        *check = true;
    }
    return res;
}

double osdiLimitlog(
    bool init, bool *check, double vnew, 
    double vold, double LIM_TOL
) {
    if (init) {
        *check = true;
        return 0.0;
    }
    int icheck = 0;
    double res = DEVlimitlog(vnew, vold, LIM_TOL, &icheck); 
    *check = icheck != 0;
    return res;
}

// OsdiLimitFunction *OsdiFile::limitFunctionTable
const OsdiLimitFunction OsdiFile::limitFunctionTable[] = {
    { "pnjlim", 2, (void *)osdiPnjlim }, 
    { "typedpnjlim", 3, (void *)osdiTypedpnjlim }, 
    { "limvds", 0, (void *)osdiLimvds }, 
    { "fetlim", 1, (void *)osdiFetlim }, 
    { "limitlog", 1, (void *)osdiLimitlog }, 
    { NULL, 0, NULL }
}; 

}
