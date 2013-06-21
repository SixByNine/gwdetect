#!/usr/bin/python

import sys
from math import *


if len(sys.argv) > 1:
    f0=float(sys.argv[1])
else:
    f0=500.0

if len(sys.argv) > 2:
    f1=float(sys.argv[2])
else:
    f1=7e-16

f1=abs(f1)


sectoyear=365.25*86400.0


A=pow(sectoyear,-4) * 4 * pow(pow(f0,-1.4) * pow(f1/1e-15,1.1)  * pow(10,1.6),2)

alpha=5
fc=0.05
print alpha,A*pow(fc,-alpha),fc


