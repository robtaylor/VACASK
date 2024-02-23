#include <math.h>
#include <iostream>
#include "limitfunctions.h"
#include "simulator.h"
#include "common.h"


namespace NAMESPACE {

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

double voltagelimit_option;	// For voltage limiting.

// Function limits the value to voltagelimit option value. Voltage limiting is
// performed only is voltagelimit is greater than zero.
static double voltagelimit(double value) {
	if(voltagelimit_option > 0) {
        if(voltagelimit_option < 1e100)	{
            // Soft limiting.
		    value = voltagelimit_option * tanh(value / voltagelimit_option);
	    } else if(fabs(value) > voltagelimit_option) {
            // Sharp limiting.
            if(value > 0) 
                value = 1e100;
	        else 
                value = -1e100;
        }
    }
	
    return value;
}

/* DEVlimvds(vnew,vold)
 * limit the per-iteration change of VDS
 */
double DEVlimvds(double vnew, double vold) {
    if(vold >= 3.5) {
        if(vnew > vold) {
            vnew = MIN(vnew,(3 * vold) +2);
        } else {
            if (vnew < 3.5) {
                vnew = MAX(vnew,2);
            }
        }
    } else {
        if(vnew > vold) {
            vnew = MIN(vnew,4);
        } else {
            vnew = MAX(vnew,-.5);
        }
    }
	vnew = voltagelimit(vnew);	// Limit voltage value.
    return(vnew);
}


/* New version of
 * DEVpnjlim(vnew,vold,vt,vcrit,icheck)
 * limit the per-iteration change of PN junction  voltages 
 * A fix by A. Buermen based on ngspice. Limiting for vnew>vcrit in ngspice
 * is weird, if |vnew-vold|<3vt. In that case vnew is limited to 
 * a value in a direction opposite to vnew-vold since
 * log((vnew-vold)/vt-2) is negative. 
*/
double DEVpnjlim(double vnew, double vold, double vt, double vcrit, int *icheck) {
    double arg;

    if((vnew > vcrit) && (fabs(vnew - vold) > (vt + vt))) {
        // Positive vnew above vcrit, change greater than 2*vt
        if(vold > 0) {
            arg = (vnew - vold) / vt;
            if(arg > 0) {
                vnew = vold + vt * log(1+arg);
            } else {
                vnew = vold - vt * log(1-arg);
            }
        } else {
            vnew = vt *log(vnew/vt);
        }
        *icheck = 1;
    } else if (vnew<0) {
        // Negative vnew, use ngspice limiting
        if (vold > 0) {
            arg = -1*vold-1;
        } else {
            arg = 2*vold-1;
        }
        if (vnew<arg) {
            vnew = arg;
            *icheck = 1;
        } else {
            *icheck = 0;
        }
    } else {
        // Positive vnew, but not above vcrit or the change is smaller than 2*vt
        *icheck = 0;
    }

	// Zero should not be limited since the code below produces NaN in that case.
    // Last resort limiting added by J. Puhan
    // Improved by A. Buermen
	if(fabs(vnew)>1e-3) {
        // Below 1e-3 the change in tanh() limiting is below 0.3ppm
        double tol;

        // vnew must be nonzero
		arg = voltagelimit(vnew);	// Limit voltage value.

        tol = fabs(arg/vnew);
        if (tol>1.01 || tol<0.99) {
            // There is more that 1% difference.
            vnew = arg;
            *icheck = 1;
		}
	}
    return(vnew);
}


/* DEVpnjlim(vnew,vold,vt,vcrit,icheck)
 * limit the per-iteration change of PN junction  voltages 
 */
double DEVpnjlim_old(double vnew, double vold, double vt, double vcrit, int *icheck) {
    double arg;

    if((vnew > vcrit) && (fabs(vnew - vold) > (vt + vt))) {
        if(vold > 0) {
            arg = 1 + (vnew - vold) / vt;
            if(arg > 0) {
                vnew = vold + vt * log(arg);
            } else {
                vnew = vcrit;
            }
        } else {
            vnew = vt *log(vnew/vt);
        }
        *icheck = 1;
    } else {
        *icheck = 0;
    }
	// Zero should not be limited since the code below produces NaN in that case.
	if(vnew)
	{
		arg = voltagelimit(vnew);	// Limit voltage value.
		if(arg / vnew < 0.99) *icheck = 1;	// There is more that 1% difference.
		vnew = arg;
	}
    return(vnew);
}

/*
 * DEVfetlim(vnew,vold,vto)
 *
 * limit the per-iteration change of FET voltages 
 */
double DEVfetlim(double vnew, double vold, double vto) {
    double vtsthi;
    double vtstlo;
    double vtox;
    double delv;
    double vtemp;

    vtsthi = fabs(2*(vold-vto))+2;
    // vtstlo = vtsthi/2 +2;
    vtstlo = fabs(vold-vto)+1;
    vtox = vto + 3.5;
    delv = vnew-vold;

    if (vold >= vto) {
        if(vold >= vtox) {
            if(delv <= 0) {
                /* going off */
                if(vnew >= vtox) {
                    if(-delv >vtstlo) {
                        vnew =  vold - vtstlo;
                    }
                } else {
                    vnew = MAX(vnew,vto+2);
                }
            } else {
                /* staying on */
                if(delv >= vtsthi) {
                    vnew = vold + vtsthi;
                }
            }
        } else {
            /* middle region */
            if(delv <= 0) {
                /* decreasing */
                vnew = MAX(vnew,vto-.5);
            } else {
                /* increasing */
                vnew = MIN(vnew,vto+4);
            }
        }
    } else {
        /* off */
        if(delv <= 0) {
            if(-delv >vtsthi) {
                vnew = vold - vtsthi;
            } 
        } else {
            vtemp = vto + .5;
            if(vnew <= vtemp) {
                if(delv >vtstlo) {
                    vnew = vold + vtstlo;
                }
            } else {
                vnew = vtemp;
            }
        }
    }
	vnew = voltagelimit(vnew);	// Limit voltage value.
    return(vnew);
}

// std::endl is slow because it flushes the stream, use \n

/* DEVlimitlog(deltemp, deltemp_old, LIM_TOL, check)
 * Logarithmic damping the per-iteration change of deltemp beyond LIM_TOL.
 */
double DEVlimitlog(double deltemp, double deltemp_old, double LIM_TOL, int *check) {
    static int shown = 0;
    *check = 0;
    if (!shown && (isnan (deltemp) || isnan (deltemp_old))) {
        Simulator::wrn() << "\n\nThe temperature limiting function received NaN.\n";
        Simulator::wrn() << "Please check your power dissipation and improve your heat sink Rth!\n";
        Simulator::wrn() << "    This message will be shown only once.\n\n";
        deltemp = 0.0;
        *check = 1;
        shown = 1;
    }
    /* Logarithmic damping of deltemp beyond LIM_TOL */
    if (deltemp > deltemp_old + LIM_TOL) {
        deltemp = deltemp_old + LIM_TOL + log10((deltemp-deltemp_old)/LIM_TOL);
        *check = 1;
    }
    else if (deltemp < deltemp_old - LIM_TOL) {
        deltemp = deltemp_old - LIM_TOL - log10((deltemp_old-deltemp)/LIM_TOL);
        *check = 1;
    }
    return deltemp;
}

}
