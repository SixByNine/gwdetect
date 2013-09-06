#! /usr/bin/python
import sys
from matplotlib.pyplot import *
from numpy import *


xs=list()
ys=list()
n=0

for a in sys.argv[1:]:
    with open(a) as f:
        i=0
        for line in f:
            elems=line.split()
            x=float(elems[0])
            y=float(elems[1])
            if n==0:
                xs.append(x)
                ys.append(y)
            else:
                ys[i]+=y
            i+=1

    n+=1.0


ys=array(ys)
xs=array(xs)
ys/=float(n)


f=open("av.spec","w")
for x,y in zip(xs,ys):
    f.write("%g %g\n"%(x,y))
f.close()
