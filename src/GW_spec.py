#!/usr/bin/python
import sys
import os
from numpy import *
import gtk
import matplotlib as mpl
import copy
mpl.use("GTKAgg")

import matplotlib.pyplot as plt


class gwspec():
    def __init__(this):
        this.pulsars=list()
        this.fig=plt.figure(1)
        this.psrid=0
        this.ALPHA=-13.0/3.0
        this.mode=""
        this.GWB=0
    def go(t,args):
        xmin=1e99
        xmax=0
        ymin=1e99
        ymax=0
        files=list()
        avgfiles=list()
        gw=0
        t.AVG=False
        inf=None
        for a in args:
            if a.startswith("-i"):
                inf=a[2:]
            elif a.startswith("-g"):
                gw=float(a[2:])
            elif a.startswith("-A"):
                t.AVG=True
            else:
                if t.AVG:
                    avgfiles.append(a)
                else:
                    files.append(a)
        params=dict()
        if inf!=None:
            f=open(inf)
            for line in f:
                elems=line.split()
                name=elems[0]
                pw=float(elems[1])
                pr=float(elems[2])
                a=float(elems[3])
                fc=float(elems[4])
                params[name]=(pw,pr,a,fc)
                print name,pw,pr,a,fc

        for fn in files:
            f=open(fn)
            freqs=list()
            powers=list()
            for line in f:
                e=line.split()
                freq=float(e[0])*365.25
                pwr=float(e[1])
                freqs.append(freq)
                powers.append(pwr)
            freqs=array(freqs)
            powers=array(powers)
            xmin=min(amin(freqs),xmin)
            xmax=max(amax(freqs),xmax)
            ymin=min(amin(powers),ymin)
            ymax=max(amax(powers),ymax)
            psrn=os.path.basename(fn)[:-5]
            fc=freqs[0]
            a=0
            pr=0
            pw=mean(powers[-4:])
            if psrn in params:
                pw,pr,a,fc=params[psrn]

            avgp=zeros(len(powers))
            avgpe=zeros(len(powers))
            if t.AVG:
                ps=list()
                for f in freqs:
                    ps.append(list())
                navg=0
                for ffn in avgfiles:
                    if os.path.basename(ffn)[:-5]==psrn:
                         ff=open(ffn)
                         i=0
                         for line in ff:
                            e=line.split()
                            apwr=float(e[1])
                            ps[i].append(apwr)
                            i+=1
                         ff.close()
                         navg+=1.0
                i=0
                for f in freqs:
                    ppp=array(ps[i])
                    avgp[i]=mean(ppp)
                    avgpe[i]=sqrt(var(ppp)/navg)
                    i+=1

            print psrn,pw,pr,a,fc
            t.pulsars.append(dict(N=psrn, F=freqs,P=powers,M=pwr_model(pw,GWB=gw,ALPHA=t.ALPHA,red=[[pr,a,fc]]),A=avgp,Ae=avgpe))

        t.fig.canvas.mpl_disconnect(3)
        t.fig.canvas.mpl_connect('button_press_event', t.click)
        t.fig.canvas.mpl_connect('key_press_event', t.keypress)

        t.xmin=xmin
        t.xmax=xmax
        t.ymin=ymin
        t.ymax=ymax
        t.draw()
        plt.show()

    def draw(t):
        while t.psrid >= len(t.pulsars):
            t.psrid-=len(t.pulsars)
        while t.psrid < 0:
            t.psrid+=len(t.pulsars)
        m=t.pulsars[t.psrid]['M'].value(t.pulsars[t.psrid]['F'])
        t.fig.clear()
        plt.loglog(t.pulsars[t.psrid]['F'],t.pulsars[t.psrid]['P'],color='b')
        plt.loglog(t.pulsars[t.psrid]['F'],abs(m-t.pulsars[t.psrid]['P']),'--',color='r')
        plt.loglog(t.pulsars[t.psrid]['F'],t.pulsars[t.psrid]['P']-m,color='r')
        plt.loglog(t.pulsars[t.psrid]['F'],m,color='k')
        plt.title("%s [%s] %.3f %.3f"%(t.pulsars[t.psrid]['N'],t.mode,t.gof()/float(len(m)),t.ggof()[1]))
        plt.xlim(t.xmin/2.0,t.xmax*2.0)
        plt.ylim(t.ymin/2.0,t.ymax*2.0)
        if t.AVG:
            plt.plot(t.pulsars[t.psrid]['F'],t.pulsars[t.psrid]['A'],color='gray')
            plt.errorbar(t.pulsars[t.psrid]['F'],t.pulsars[t.psrid]['A'],t.pulsars[t.psrid]['Ae'],color='gray')

        plt.draw()
        print "redraw"


    def fit(t,white=[1],GW_A=[1],r_A=[1],r_a=[1],r_fc=[1]):
        m=t.pulsars[t.psrid]['M']
        orig_m=copy.deepcopy(m)
        best_gof=t.gof()
        best_m=copy.deepcopy(m)
        for A in GW_A:
            for w in white:
                for rA in r_A:
                    for ra in r_a:
                        for rfc in r_fc:
                            m.GWB=orig_m.GWB*A
                            m.red[0][0]=orig_m.red[0][0]*rA
                            m.red[0][1]=orig_m.red[0][1]*ra
                            m.red[0][2]=orig_m.red[0][2]*rfc
                            m.white=orig_m.white*w
                            gof=t.gof()
                            if gof < best_gof:
                                best_gof=gof
                                best_m=copy.deepcopy(m)
        print best_m.GWB,best_m.white,best_m.red[0]
        t.pulsars[t.psrid]['M']=best_m
        for psr in t.pulsars:
            psr['M'].GWB=best_m.GWB

    def ggof(t):
        s=0
        k=0
        for i in range(len(t.pulsars)):
            s+=t.gof(i)
            k+=len(t.pulsars[i]['F'])
        return s,s/float(k)

    def gof(t,id=-1):
        if id < 0:
            id=t.psrid
        m=t.pulsars[id]['M'].value(t.pulsars[id]['F'])
        v=t.pulsars[id]['P']
        gof=sum(log(m/mean(v))+v/m)
        return gof

    def click(t,evt):
        #print evt
        clickX=evt.xdata
        clickY=evt.ydata
        if t.mode=='g':
            t.GWB=sqrt(clickY*pow(clickX,-t.ALPHA)*(12.0*pi*pi))
            print t.GWB
            for psr in t.pulsars:
                psr['M'].GWB=t.GWB

        elif t.mode=='w':
            t.pulsars[t.psrid]['M'].white=clickY

        elif t.mode=='r':
            t.r1=(clickX,clickY)
            t.mode='rrr'
        elif t.mode=='rrr':
            fc=t.r1[0]
            x1=log10(clickX)
            x2=log10(t.r1[0])
            y1=log10(clickY)
            y2=log10(t.r1[1])
            m=-(y1-y2)/(x1-x2)
            alpha=m
            A=clickY*power(1.0+power(clickX/fc,2.0),alpha/2.0)
            print A,alpha,fc
            t.pulsars[t.psrid]['M'].red[0]=[A,alpha,fc]
            t.mode='r'
        t.draw()

    def printit(t):
        rn=open("psr.noise","w")
        print "GWB=%.3g"%t.pulsars[0]['M'].GWB
        for psr in t.pulsars:
            m=psr['M']
            n=psr['N']
            print n,m
            rn.write("%s %.6g %.6g %.6f %.6f\n"%(n,m.white,m.red[0][0],m.red[0][1],m.red[0][2]))
        rn.close()
            


    def keypress(t,evt):
        #print evt.key
        if evt.key=='escape' or evt.key=='q':
            t.mode=""
        elif evt.key==',':
            t.psrid-=1
        elif evt.key=='.':
            t.psrid+=1
        elif evt.key=='Q':
            sys.exit(0)
        elif evt.key=='g':
            t.mode="g"
        elif evt.key=='w':
            t.mode='w'
        elif evt.key=='r':
            t.mode='r'
        elif evt.key=='p':
            t.printit()
         # set the GWB
        elif t.mode=='g':
            if evt.key=='f':
                t.fit(white=[0.5,0.8,0.9,1,1.1,1.2,2],GW_A=[0.5,0.8,0.9,1,1.1,1.2,2])
        elif t.mode=='r':
            if evt.key=='f':
                t.fit(white=[0.5,0.8,0.9,1,1.1,1.2,2],r_A=[0.5,0.8,0.9,1,1.1,1.2,2],r_a=[0.5,0.8,0.9,1,1.1,1.2,2])
            if evt.key=='d':
                t.pulsars[t.psrid]['M'].red[0]=[0,0,t.pulsars[t.psrid]['F'][0]]
        elif t.mode=='w':
            if evt.key=='f':
                t.fit(white=[0.5,0.8,0.9,1,1.1,1.2,2])


        t.draw();

class pwr_model:
    def __init__(self,white,GWB=0,red=[],ALPHA=-13.0/3.0):
        self.white=white
#        self.GWB=1e-15
        self.GWB=GWB
        self.red=red
        self.ALPHA=ALPHA

    def value(self,freq):
        ret=freq*0+self.white
        ret+=self.GWB*self.GWB*power(freq,self.ALPHA)/(12.0*pi*pi)
        for red in self.red:
            ret+= red[0]/power(1.0+power(freq/red[2],2.0),red[1]/2.0)
        return ret

    def __str__(self):
        r="Ag=%.4g Pw=%.4g"%(self.GWB,self.white)
        i=0
        for red in self.red:
            r+=" R[%d]={A=%.3g a=%.3f fc=%.3f}"%(i,red[0],red[1],red[2])
            i+=1
        return r



if __name__=="__main__":
    psrid=0
    g=gwspec()
    g.go(sys.argv[1:])




