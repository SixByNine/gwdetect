#!/usr/bin/python
import sys
from math import sqrt



after=0
band=0
bw=1e99
efac=1
equad=0
tg_eq=0
tgt=0
psd=False
for i in sys.argv:
    elems=i.split("=",2)
    k=elems[0]
    if len(elems) > 1:
        v=elems[1]
    else:
        v=""
    if k=="-psd":
        psd=True
    if k=="-equad":
        equad=float(v)
    if k=="-efac":
        efac=float(v)
    if k=='-tpsd':
        tgt=float(v)

    if k=='-tpsd-eq':
        tgt=float(v)
        tg_eq=1

    if k=='-freq':
        band=float(v)
    if k=='-bw':
        bw=float(v)
    if k=='-after':
        after=float(v)



def doit(fn,efac,equad):

    f=open(fn)
    s=0
    sw=0
    mode=0
    if fn.endswith(".ifuncDGW") or fn.endswith(".cm"):
        mode=1


    t0=9000000
    t1=0

    sss=0

    t2efac=list()
    t2equad=list()
    skipon=False
    sk=0
    for line in f:
        try:
            elems=line.split()
            if mode==0 and (len(elems) < 5 or elems[0]=="C" or elems[0]=="#"):
                if len(elems) > 0:
                    if elems[0]=="SKIP":
                        skipon=True
                    if elems[0]=="NOSKIP":
                        skipon=False
                    if skipon:
                        continue
                    if elems[0]=="EQUAD":
                        equad+=float(elems[1])
                        pass
                    if elems[0]=="EFAC":
                        efac*=float(elems[1])
                    if elems[0]=="T2EFAC":
                        t2efac.append((elems[1],elems[2],float(elems[3])))
                    if elems[0]=="T2EQUAD":
                        t2equad.append((elems[1],elems[2],float(elems[3])))
                continue
            if skipon:
                sk+=1
                continue
            if mode==0 and abs(float(elems[1])-band) > bw:
                continue
            EF=efac
            EQ=equad
            i=0
            for a in elems:
                i+=1
                for fl,v,e in t2efac:
                    if a==fl:
                        if elems[i]==v:
                            EF*=e
                for fl,v,e in t2equad:
                    if a==fl:
                        if elems[i]==v:
                            EQ+=e

            if mode==0:
                err=float(elems[3])
                t=float(elems[2])
            else:
                err=float(elems[2])*1e6
                t=float(elems[0])
            err=sqrt(err*err + EQ*EQ)
            err*=EF

            if t < after:
                continue
            
            w=1.0/(err*err*1e-12)
            

            eee=err/1e6/86400.0/365.25
            sss+=1.0/(eee*eee)

            sw+=w
            s+=1
            if t < t0:
                t0=t
            if t > t1:
                t1=t
        except:
            print "ERROR"
            print line


    tobs=(t1-t0)/365.25
    if sw==0:
        sw=1
    if sss==0:
        sss=1
    f.close()

    return 2*tobs/sss,sqrt(s/sw)



p,r = doit(sys.argv[1],efac,equad)

if tgt != 0 :
    if tg_eq==1:
        while p < tgt:
            equad += 1e-2
            p,r=doit(sys.argv[1],efac,equad)
        print "%.4g"%equad
    else:
        while abs(p-tgt)/tgt > 0.01:
            efac*=sqrt(tgt/p)
            p,r=doit(sys.argv[1],efac,equad)
        print "%.4g"%efac
else:
    if psd:
        print "%.5g"%p
    else:
        print "%.5g"%r


