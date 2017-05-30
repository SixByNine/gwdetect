#!/usr/bin/env python
import os
from math import *
from numpy import *
import matplotlib.pyplot as plt
from sys import argv,exit
import argparse
import matplotlib.mlab as mlab
stat_desc=["a=0","a=0.5","a=1.0","a=itr","a"]

ON="on"
OFF="off"

plotstats=zeros((10,5))

parser = argparse.ArgumentParser(description='Histogram of detection stats')
parser.add_argument('-A','--gwamp',type=float,default=1e-15)
parser.add_argument('--test',default=False, action='store_true')
parser.add_argument('-c','--oncols',nargs='*',type=int,default=[1])
parser.add_argument('-C','--offcols',nargs='*',type=int)
parser.add_argument('-t','--title')
parser.add_argument('-x','--xmin',type=float,default=None)
parser.add_argument('-X','--xmax',type=float,default=None)
parser.add_argument('--noplot',default=False, action='store_true')
parser.add_argument('-N','--nbins',type=int,default=40)
args=parser.parse_args()
A2=args.gwamp*args.gwamp

if args.offcols == None:
    args.offcols=args.oncols
stats=dict()
for onoff in [ON,OFF]:
    N=0
    for line in open(onoff+"stats"):
        N+=1
    i=0
    for line in open(onoff+"stats"):
        e=line.split()
        M=len(e)
        if onoff not in stats:
            stats[onoff] = zeros((M,N))
        for j in range(len(e)):
            if j < 4:
                stats[onoff][j][i] = float(e[j])/A2
            else:
                stats[onoff][j][i] = float(e[j])

        i+=1



print "Red/Blue"
print var(stats[ON][0])/var(stats[ON][3])
print var(stats[OFF][0])/var(stats[OFF][3])

print "Red/Green"
print var(stats[ON][0])/var(stats[ON][2])
print var(stats[OFF][0])/var(stats[OFF][2])
YMAX=0
XMAX=-1e999
XMIN=1e999

for onoff in [ON,OFF]:
    cols=args.oncols
    if onoff==OFF:
        cols=args.offcols
    for c in cols:
        XMAX = max(amax(stats[onoff][c-1]),XMAX)
        XMIN = min(amin(stats[onoff][c-1]),XMIN)



plt.figure(figsize=(12,8),dpi=100,facecolor='w',edgecolor='w')
rng=(XMIN,XMAX)


colours=['red','orange','green','blue']
for onoff in [ON,OFF]:
    if onoff==ON:
        colours=['#FF2020','orange','green','#2020FF','#20FF20']
    else:
        colours=['#B00000','orange','green','#0000B0','#00B000']

    colours.extend(colours)

    print onoff
    cols=args.oncols
    if onoff==OFF:
        cols=args.offcols
    for c in cols:
        print "%s %d %%1"%(onoff,c),
        c-=1
        A2g= mean(stats[onoff][c])
        sig3=sqrt(var(stats[OFF][c]))*3.0
        print A2g,colours[c],sqrt(var(stats[onoff][c])),sig3
        n,b,p = plt.hist(stats[onoff][c],histtype='stepfilled',alpha=0.3,color=colours[c],normed=True,bins=args.nbins,range=rng)
        sigma=sqrt(var(stats[onoff][c]))
        if onoff==ON:
            plotstats[c][0]=sigma
            plotstats[c][4]=sig3
        else:
            plotstats[c][1]=sigma
        bstep=b[1]-b[0]
        x = linspace(b[0]-bstep*5,b[-1]+bstep*5,len(n)*10)
        y = mlab.normpdf(x, A2g, sigma)
        if amax(n) > YMAX:
            YMAX=amax(n)
        #if amax(b) > XMAX:
        #    XMAX=amax(b)
        #if amin(b) < XMIN:
        #    XMIN=amin(b)

        plt.plot(x,y,'--',color=colours[c],linewidth=2,label="correlation:%s statistic:%s"%(onoff,stat_desc[c]))
        plt.vlines(A2g,0,YMAX*10,linestyle="-",color=colours[c])
        fout=open("plot.%d.%s"%(c+1,onoff),"w")
        for vx,vy in zip(x,y):
            fout.write("%g %g\n"%(vx,vy))
        fout.close()

        if onoff==OFF:
            plt.vlines(A2g+3*sigma,0,YMAX*10,linestyle=':',color=colours[c],linewidth=2)
        else:
            sig3=sqrt(var(stats[OFF][c]))*3.0
            nn=0.0
            for i in stats[onoff][c]:
                if i > sig3:
                    nn+=1.0
            mm=0.0
            ss=0.0
            dx=x[1]-x[0]
            for ix,iy in zip(x,y):
                if ix > sig3:
                    mm+=iy*dx
                ss+=iy*dx

            print "Detected count=%.1f%%, gauss=%.1f%% gauss100=%.1f%%"%(100*nn/N,100*mm,100*ss)
            print "%s %d %%2 %.4f %.4f %.4f"%(onoff,c+1, nn/N,mm,ss)
            plotstats[c][2]=nn/N
            plotstats[c][3]=mm


try:
    for cc in args.oncols:
        c=cc-1
        print cc,c,plotstats[c]
        f=open("plotstats.%d"%cc,"w")
        f.write("% 8.4f % 8.4f % 8.4f % 8.4f % 8.4f\n"%(plotstats[c][0],plotstats[c][1],plotstats[c][2],plotstats[c][3],plotstats[c][4]))
        f.close()
except Exception as e:
    print e

if args.noplot == False:

    plt.vlines([0,A2/A2],0,YMAX*2,'k',linewidth=2)
    XR=XMAX-XMIN
    xmin=XMIN-0.1*XR
    xmax=XMAX+0.1*XR
    if not args.xmin == None:
        xmin=args.xmin

    if not args.xmax == None:
        xmax=args.xmax
    plt.xlim((xmin,xmax))
    plt.ylim((0,YMAX+0.1*YMAX))


    plt.xlabel("Statistic value (A^2)")
    plt.ylabel("Probability density")
    plt.legend()
    if args.title!=None:
        plt.title(args.title)
    plt.show()



