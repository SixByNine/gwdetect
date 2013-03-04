#!/usr/bin/python
from sys import argv
from numpy import *
from matplotlib.pyplot import *


avg=list()
W=list()
avgE2=list()
freqs=list()
first=True
n=0

A=float(argv[1])
for a in argv[2:]:
    print a
    f=open(a)
    i=0
    for line in f:
        elems=line.split()
        a=float(elems[0])
        v=float(elems[1])
        e=float(elems[2])
        w=float(elems[3])
        if first:
            freqs.append(a)
            avg.append(v*w)
            avgE2.append(e*e*w*w)
            W.append(w)
        else:
            avg[i]+=v*w
            avgE2[i]+=e*e*w*w
            W[i]+=w
        i+=1
    if first:
        print "1st"
        first=False
        avg=array(avg)
        avgE2=array(avgE2)
        freqs=array(freqs)
        W=array(W)
    n+=1.0

avg/=W
avgE2/=W*W
avgE = sqrt(avgE2)


f=open("ff.avg","w")
for a,v,e in zip(freqs,avg,avgE):
    f.write("%f %g %g\n"%(a,v,e))
f.close()


A2_guess=A*A
errorbar(freqs,avg,avgE,fmt='x',color='red')
plot(freqs,freqs*0+A2_guess,color='red')
plot(freqs,freqs*0,color='black')


show()
