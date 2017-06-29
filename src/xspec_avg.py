#!/usr/bin/env python


import numpy as np
from matplotlib import pyplot as plt
import argparse
from sys import argv
from math import pi


A = float(argv[1])

A2=A*A


freqs=None
reals=None
imags=None
pwrA=None
pwrB=None
count=0.0

for fn in argv[2:]:
    data=np.loadtxt(fn).T
    if count==0.0:
        freqs = data[0]
        reals=data[1]
        imags=data[2]
        pwrA=data[3]
        pwrB=data[4]
    else:

        if data[3][0] > 10*(pwrA[0]/count):
            print "A", data[3][0], (pwrA[0]/count), fn

        if data[4][0] > 10*(pwrB[0]/count):
            print "B",data[4][0], (pwrB[0]/count), fn
        reals+=data[1]
        imags+=data[2]
        pwrA+=data[3]
        pwrB+=data[4]

    count+=1.0

    

reals/=count
imags/=count
pwrA/=count
pwrB/=count

print freqs

plt.loglog(freqs,pwrA,color='black')
plt.loglog(freqs,pwrB, color='green')

plt.loglog(freqs,reals,'o',color='red')
plt.loglog(freqs,-reals,'o',color='magenta')

plt.loglog(freqs,imags,'o',color='cyan')
plt.loglog(freqs,-imags,'o',color='blue')


zeta=0.5

plt.loglog(freqs,np.power(freqs,-13.0/3.0)*A2/12.0/pi/pi,'-',color='orange')
plt.loglog(freqs,np.power(freqs,-13.0/3.0)*zeta*A2/12.0/pi/pi,':',color='orange')

plt.title("%s %s"%(argv[1],argv[2]))
plt.xlabel("Frequency (yr^-1)")
plt.ylabel("PSD (yr^3)")


plt.savefig("xspec_avg.png")

plt.show()





