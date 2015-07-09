#!/usr/bin/python
from sys import argv
from numpy import *
import argparse

parser = argparse.ArgumentParser(description='Mean, variance, other useful things.')
parser.add_argument('file')
parser.add_argument('col',type=int)
parser.add_argument('-m','--median',action="store_const", const='median')
parser.add_argument('-w','--weights',type=int,default=-1)
args=parser.parse_args()

w=list()
vals=list()
f=open(args.file)
col=int(args.col)-1
wcol=0
if args.weights > 0:
    wcol = args.weights

for line in f:
    elems=line.split()
    vals.append(float(elems[col]))
    if wcol > 0:
        s=float(elems[wcol-1])
        w.append(1.0/pow(s,2))
    else:
        w.append(1.0)


w=array(w)
w/=sum(w)
ovals=array(vals)
vals=ovals*w
m= sum(vals)
vals = (ovals - m)
rms=sum(w*w*vals*vals)/sum(w*w)


if args.median == "median":
    p=percentile(ovals,84)
    med=median(ovals)
    print med,p,p-med
else:
    print m,rms,sqrt(rms)

