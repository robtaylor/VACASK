import sys, os
import subprocess
from pathlib import Path
from functools import reduce
import numpy as np
from pprint import pprint

__all__ = [ "absDiff", "relDiff", "isTest" ]

def absDiff(a, b):
    return np.abs(a-b)

def relDiff(a, b, abstol):
    aabs = np.abs(a)
    babs = np.abs(b)
    ref = np.where(aabs>=babs, aabs, babs)
    ref = np.where(ref>abstol, ref, abstol)
    return np.abs(a-b)/ref

def isTest():
	return "SIM_TEST" in os.environ and os.environ["SIM_TEST"]=="yes"
