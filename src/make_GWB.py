#! /usr/bin/python

import toasim
import sys
import os
from math import *
from numpy import *
import random
import argparse
import ifunc_simulator
from multiprocessing import Pool
import matplotlib.pyplot as plt
from scipy import interpolate


class make_GWB:
    def __init__(self):
        pass
    def main(self,argv):
        parser = argparse.ArgumentParser(description='Detect the GWB.')
        parser.add_argument('-A','--gwamp',type=float,default=1e-15)
        parser.add_argument('-a','--angles',required=True)
        parser.add_argument('-p','--psrlist',default=None)
        parser.add_argument('--par',default=".par")
        parser.add_argument('--tim',default=".tim")
        parser.add_argument('-n','--npts',type=int,default=1024)
        parser.add_argument('-N','--nreal',type=int,default=1000)
        parser.add_argument('-U','--uncorr',default=False, action='store_true')
        

        
        args=parser.parse_args(argv)

        psrs=list()
        pairs=list()
        angles=dict()
        DO_GW_CORR = True
        if args.uncorr:
            DO_GW_CORR = False
            print "WARNING: GW correlation is DISABLED."
        
        Agwb = args.gwamp

        fixedpsrlist=False
        if args.psrlist!=None:
            fixedpsrlist=True
            with open(args.psrlist) as f:
                for line in f:
                    elems=line.split()
                    psrs.append(elems[0])
            
        
        with open(args.angles) as f:
            for line in f:
                elems = line.split()
                p1=elems[0]
                p2=elems[1]
                a=float(elems[2])
                if fixedpsrlist:
                    if p1 not in psrs or p2 not in psrs:
                        continue
                else:
                    if p1 not in psrs:
                        psrs.append(p1)
                    if p2 not in psrs:
                        psrs.append(p2)
                if p1 not in angles:
                    angles[p1]=dict()
                if p2 not in angles:
                    angles[p2]=dict()
                
                angles[p1][p2]=a
                angles[p1][p1]=0
                angles[p2][p1]=a
                angles[p2][p2]=0
                pairs.append((p1,p2))
                
        Npsr=len(psrs)
        Npair=len(pairs)
        if (Npsr*(Npsr-1))/2 != Npair:
            print "ERROR: Number of pairs is not N(N-1)/2"
            sys.exit(1)
            
        startMjd=1e99
        endMjd=0
        psrnfo=dict()
        for p in psrs:
            psrnfo[p] = psrinfo(p,"%s%s"%(p,args.tim))
            startMjd = min(psrnfo[p].start(),startMjd)
            endMjd = max(psrnfo[p].finish(),endMjd)
        
        # Extend the window by one day to ensure validity
        startMjd-=0.5
        endMjd+=0.5
        
        npts1=args.npts
        tspan=endMjd - startMjd
        nreal=args.nreal
        nptsA = npts1*nreal
        
        gridMjds = linspace(startMjd,endMjd,npts1)
        dt = gridMjds[1]-gridMjds[0]
        print "dt =",dt
        print "nreal = ",nreal
        print "npts = ",npts1
        print "tspan = ",tspan
        print "start = ",startMjd
        print "end = ",endMjd
        
        rawData = zeros((Npsr,nptsA))
        
        
        pool = Pool(4)
        procs = list()
        for p in psrs:
            procs.append(pool.apply_async(ifunc_simulator.make_red_noise,args=[npts1,nreal,-13.0/3.0,dt/365.25,Agwb,True,0.001]))
        
        i=0
        for result in procs:
            rawData[i] += result.get()
            i+=1
        if DO_GW_CORR:
            C = ifunc_simulator.make_GW_covar(angles,psrs)
            L = linalg.cholesky(C)
            timeseries = dot(L,rawData)
        else:
            timeseries=rawData
        
        
        i=0
        for psr in psrs:
            print "Write","%s%s.addGWB"%(psr,args.tim)
            with open("%s%s.addGWB"%(psr,args.tim),"w") as f:
                toasim_header = toasim.header()
                toasim_header.ntoa=len(psrnfo[psr].T)
                toasim_header.nrealisations=nreal
                toasim_header.description="Gravitational Wave background"
                toasim_header.short_desc="GWB"
                toasim_header.timfile_name="%s%s"%(psr,args.tim)
                if not DO_GW_CORR:
                    toasim_header.short_desc="GWB-uncor"
                toasim_header.rparam_len=0
                toasim_header.rparam_desc=""
                toasim_header.write(f)
                real=0
                while real < nreal:
                    print "%04d\r"%real,
                    sys.stdout.flush()
                    nn = real*npts1
                    gw = interpolate.interp1d(gridMjds,timeseries[i][nn:nn+npts1],kind='cubic')
                    offsets=list()
                    for t in psrnfo[psr].T:
                        offsets.append(gw(t))
                    offsets=array(offsets)
                    poly = poly1d(polyfit(psrnfo[psr].T,offsets,3))
                    offsets-=poly(psrnfo[psr].T)

                    toasim_corr = toasim.correction(toasim_header,offsets,0,0,0,"")
                    toasim_corr.write(f)
                    real+=1
                print "%04d"%real
            i+=1
        
        

class psrinfo:
    def __init__(self,name,timfile):
        self.name=name
        t=list()
        with open(timfile) as f:
            for line in f:
                e=line.split()
                if len(e) > 5:
                    mjd=float(e[2])
                    t.append(mjd)
        self.T = array(t)
        self.T.sort()
    def start(self):
        return self.T[0]
    def finish(self):
        return self.T[-1]

if __name__ == "__main__":
    make_GWB().main(sys.argv[1:])

















