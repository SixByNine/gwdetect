#!/usr/bin/python
from sys import argv
from numpy import *
from matplotlib.pyplot import *


avg=list()
W=list()
avgE2=list()
angles=list()
pairs=list()
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
        w=float(elems[7])
        if first:
            angles.append(a)
            avg.append(v*w)
            avgE2.append(e*e*w*w)
            W.append(w)
            p=elems[8]+" "+elems[9]
            pairs.append(p)
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
        angles=array(angles)
        W=array(W)
    n+=1.0

avg/=W
avgE2/=W*W
avgE = sqrt(avgE2)

aa=linspace(0.01,pi,128)
x=(1.0-cos(aa))/2.0
azeta=(3.0/2.0)*x*log(x) - x/4.0 + 0.5
x=(1.0-cos(pi*angles/180.0))/2.0
zeta=(3.0/2.0)*x*log(x) - x/4.0 + 0.5



f=open("hd.avg","w")
for a,v,e in zip(angles,avg,avgE):
    f.write("%f %g %g\n"%(a,v,e))
f.close()


A2_guess=A*A
p=0
while p < len(angles):
    if abs(avg[p]-A2_guess*zeta[p]) > 4*avgE[p]:
        print "BAD",pairs[p], angles[p],avg[p],A2_guess*zeta[p],avgE[p],abs(avg[p]-A2_guess*zeta[p])/avgE[p]
    p+=1

errorbar(angles,avg,avgE,fmt='x',color='red')
plot(180*aa/pi,azeta*A2_guess,color='red')
plot(180*aa/pi,azeta*0,color='black')


show()
