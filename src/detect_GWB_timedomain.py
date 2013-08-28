#!/usr/bin/python
from sys import argv,stdout,stderr
from math import *
from numpy import *
from numpy import linalg as lalg
from scipy import linalg as slalg
from matplotlib import pyplot as plt
import random
import argparse


ALPHA=-13.0/3.0

def detect_GWB_td():
    parse = argparse.ArgumentParser()
    parse.add_argument('-A','--gwamp',type=float,default=1e-15)
    parse.add_argument('-e','--ext',default=".cm")
    parse.add_argument('-f','--angles_file',required=True)
    parse.add_argument('-p','--psrlist',required=True)
    args=parse.parse_args()

    A=args.gwamp
    A2=A*A

    print "Read pulsars"
    psrs=list()
    f=open(args.psrlist)
    for line in f:
        e = line.split()
        if e[0].startswith("J"):
            psrs.append(e[0])
    f.close()
    print "Read angles"
    pairs=list()
    angles=list()
    f=open(args.angles_file)
    for line in f:
        e = line.split()
        p1=e[0]
        p2=e[1]
        angle=float(e[2])
        if p1 in psrs and p2 in psrs:
            pairs.append((p1,p2))
            angles.append(angle)
    f.close()

    Npsr = len(psrs)
    Npair= len(pairs)
    if Npsr*(Npsr-1)/2 == Npair:
        print "Npsr = %d Npair= %d"%(Npsr,Npair)
    else:
        print "Error, Npsr = %d Npair= %d, expected Npair = %d"%(Npsr,Npair,Npsr*(Npsr-1)/2)
        exit(1)

    angles=array(angles)
    print "Generate zetas"
    x=(1.0-cos(radians(angles)))/2.0
    zetas=(3.0/2.0)*x*log(x) - x/4.0 + 0.5
    for a in range(Npair):
        for b in range (Npair):
            getZetas(a,b,pairs,zetas)
    

    print "Read timeseries"
    timeseries=dict()

    for psr in psrs:
        fn = "%s%s"%(psr,args.ext)
        print "Read: %s\r"%fn,
        stdout.flush()
        f=open(fn)
        mjd=list()
        vals=list()
        errs=list()
        for line in f:
            elems=line.split()
            if elems[0]=="_CM":
                mjd.append(float(elems[1]))
                vals.append(float(elems[2]))
                errs.append(float(elems[3]))

        mjd=array(mjd)
        vals=array(vals)
        errs=array(errs)
        timeseries[psr] = ((mjd,vals,errs))
        f.close()

    print "%d files read                       "%Npsr


    print "Compute correlations"
    for pair in range(Npair):
        psr1=pairs[pair][0]
        psr2=pairs[pair][1]
        m1,v1,e1 = timeseries[psr1]
        m2,v2,e2 = timeseries[psr2]
        




zeta_lookup=dict()
def getZetas(p1, p2, pairs,zetas):
    key="%d %d"%(p1,p2)
    if key not in zeta_lookup:
        zeta_lookup[key] = _getZetas(p1,p2,pairs,zetas)
    return zeta_lookup[key]

def _getZetas(p1, p2, pairs,zetas):
    psr_i = pairs[p1][0]
    psr_j = pairs[p1][1]
    psr_l = pairs[p2][0]
    psr_m = pairs[p2][1]

    zeta_ij = zetas[p1]
    zeta_lm = zetas[p2]
    zeta_il=1
    zeta_jm=1
    zeta_im=1
    zeta_jl=1
    idx=0
    if psr_i != psr_l:
        for pp in pairs:
            if pp[0] == psr_i and pp[1]==psr_l:
                zeta_il = zetas[idx]
                break
            elif pp[1] == psr_i and pp[0]==psr_l:
                zeta_il = zetas[idx]
                break
            idx+=1

    idx=0
    if psr_j != psr_m:
        for pp in pairs:
            if pp[0] == psr_j and pp[1]==psr_m:
                zeta_jm = zetas[idx]
                break
            elif pp[1] == psr_j and pp[0]==psr_m:
                zeta_jm = zetas[idx]
                break
            idx+=1

    idx=0
    if psr_i != psr_m:
        for pp in pairs:
            if pp[0] == psr_i and pp[1]==psr_m:
                zeta_im = zetas[idx]
                break
            elif pp[1] == psr_i and pp[0]==psr_m:
                zeta_im = zetas[idx]
                break
            idx+=1

    idx=0
    if psr_j != psr_l:
        for pp in pairs:
            if pp[0] == psr_j and pp[1]==psr_l:
                zeta_jl = zetas[idx]
                break
            elif pp[1] == psr_j and pp[0]==psr_l:
                zeta_jl = zetas[idx]
                break

            idx+=1

    return zeta_ij,zeta_lm,zeta_il,zeta_jm,zeta_im,zeta_jl


class pwr_model:
    def __init__(self):
        self.white=1e-30
#        self.GWB=1e-15
        self.red=[(1e-30,-4)]

    def value(self,freq):
        ret=freq*0+self.white
#        ret+=self.GWB*self.GWB*power(freq,ALPHA)/(12.0*pi*pi)
        for red in self.red:
            ret += red[0]/power(1.0+power(freq/red[2],2.0),red[1]/2.0)
        return ret



    
def P_gw_nrm(freq):
    K=power(0.0335955,ALPHA)
    return K*power(1.0 + power(freq/0.0335955,2),ALPHA/2.0)/(12*pi*pi)
    #return power(0.0335955+freq,ALPHA)/(12.0*pi*pi)


if __name__ == "__main__":
    detect_GWB_td()
