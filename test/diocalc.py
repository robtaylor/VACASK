import numpy as np
from scipy import optimize

# Solve diode voltage when diode is in series with a resistor
# and powered by a voltage source. 
# This module is used in testing. 

P_K = 1.3806503e-23
P_Q = 1.602176462e-19
P_CtoK = 273.15

def dioMod(u, Is, n, temp):
	vt = P_K*(temp+P_CtoK)/P_Q
	gd = Is/(n*vt)*np.exp(u/(n*vt))
	return Is*(np.exp(u/(n*vt))-1), gd

def f(x, Is=1e-12, n=2, temp=27, r=1000, v=0.8):
	# OpenVAF uses NIST1998 values of P_K and P_Q
	i, gd = dioMod(x, Is, n, temp)
	return r * i + x - v

def dioSolve(Is, n, temp, r, v):
	return optimize.bisect(f, 0, 2*v, args=(Is, n, temp, r, v), xtol=1e-16)
