#ifndef __LIMITFUNCTIONS_DEFINED
#define __LIMITFUNCTIONS_DEFINED

#include "common.h"


namespace NAMESPACE {

double DEVlimvds(double vnew, double vold);
double DEVpnjlim(double vnew, double vold, double vt, double vcrit, int *icheck);
double DEVpnjlim_old(double vnew, double vold, double vt, double vcrit, int *icheck);
double DEVfetlim(double vnew, double vold, double vto);
double DEVlimitlog(double deltemp, double deltemp_old, double LIM_TOL, int *check);

}

#endif
