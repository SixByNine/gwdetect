#!/usr/bin/python
import sys
import os
from math import *
from numpy import *
import random
from matplotlib import pyplot as plt

    
def make_DIPOLE_covar(angles, psrs):
    np=len(psrs)
    DPCovar=zeros((np,np))
    i=0
    for p1 in psrs:
        j=0
        for p2 in psrs:
            if p1==p2:
                zeta=1
            else:
                zeta=0.5*cos(radians(angles[p1][p2]))
            DPCovar[i][j]=zeta
            j+=1
        i+=1
    return DPCovar



    
def make_GW_covar(angles, psrs):
    np=len(psrs)
    GWCovar=zeros((np,np))
    i=0
    for p1 in psrs:
        j=0
        for p2 in psrs:
            if p1==p2:
                zeta=1
            else:
                x=(1.0-cos(radians(angles[p1][p2])))/2.0
                zeta=(3.0/2.0)*x*log(x) - x/4.0 + 0.5
            GWCovar[i][j]=zeta
            j+=1
        i+=1
    return GWCovar



def make_red_noise(N,M,alpha,dt,gwamp,gw=True,fc=0.05):
    L=M*N
    spec=zeros(L,dtype=dtype(complex))
    freq=fft.fftfreq(L,dt)
    if gw:
        amp = gwamp*gwamp / 12.0 / pi / pi # yr3
    else:
        print "TEST",fc,alpha,gwamp
        amp=gwamp/pow(fc,alpha) # For non-GW red noise
    print "A (yr3) =",amp
    print "A (fc) =",amp*pow(fc,alpha)
    amp = amp / (4*dt*L)
    amp=sqrt(amp)
    i=0
    for f in freq:
        if f > 0:
            A=amp*pow(f,alpha/2.0)
            spec[i]=A*complex(random.gauss(0,1),random.gauss(0,1))
            spec[L-i]=spec[i].conjugate()
        i+=1
    spec[0]=0
    spec[L/2]=0
    ts = fft.ifft(spec)
    ts=real(ts)*L
    ts *=365.25*86400.0
    return ts


class ifunc_simulator:
    def __init__(self,args):
        Nreal=100
        Agwb=1.2e-15
        GWcor=True
        noop=False
        if "-1000" in args:
            Nreal=1000
        if "-4000" in args:
            Nreal=4000
        if "-400" in args:
            Nreal=400
        if "-uncor" in args:
            GWcor=False
        if "-noop" in args:
            noop=True
        if "-nogw" in args:
            Agwb=0
        if "-A=2" in args:
            Agwb=2e-15

 
        print "Nreal=",Nreal
        self.xtra=list()
        f=open("xtra.psrs")
        for line in f:
            self.xtra.append(line.strip())

        print self.xtra

        f.close()
        f=open("psd.list")
        self.curr = dict()
        self.red = dict()
        self.tsamp = 60.0

        if "-14" in args:
            self.tsamp=14.0
        for line in f:
            if line.startswith("#"):
                continue
            elems=line.split()
            psr=elems[0]
            if elems[1] == "FILE":
                psdfile = open(elems[2])
                psd=list()
                lda=list()
                rte=list()
                for l2 in psdfile:
                    e2=l2.split()
                    freq=float(e2[0])*1e6
                    psd.append(float(e2[1]))
                    lda.append(3e8/freq)
                    rte.append(self.tsamp/15)
                psd=array(psd)
                lda=array(lda)
                rte=array(rte)
                self.curr[psr]=dict(psd=psd,lda=lda,rte=rte)
            elif elems[1] == "CM":
                psd=float(elems[2]) 
                self.curr[psr]=dict(psd=array([psd]),lda=array([0.0]),rte=array([1]))
            else:
                psd10=float(elems[1])
                psd20=float(elems[2])
                psd40=float(elems[3])
                self.curr[psr]=dict(psd=array([psd10,psd20,psd40]),lda=array([0.1,0.2,0.4]),rte=self.tsamp/array([15,10,15]))
            self.psd2sma(self.curr[psr])

        self.psrs=self.curr.keys()
        self.psrs.sort()
        npsr=len(self.psrs)
        print npsr
        self.curr2wbr()
        f.close()


        f=open("red.noise")
        for line in f:
            if line.startswith("#"):
                continue
            elems=line.split()
            psr=elems[0]
            alpha=float(elems[1])
            amp=float(elems[2])
            fc=float(elems[3])
            if alpha > 0:
                alpha = -alpha
            self.red[psr] = dict(amp=amp,alpha=alpha,fc=fc)
            pass
        f.close()

        f=open("toa.list","w")
        ff=open("toa.tex","w")
        f.write("#PSRJ  (us)  10cm    20cm    40cm    UWL   ;   W10cm   W13cm   W17cm   W20cm   W30cm   W40cm R:CM R:DM\n")
        ff.write("PSRJ  &  10cm  &  20cm &   40cm &   UWL   &   W10cm &  W13cm  & W17cm &  W20cm &  W30cm  & W40cm &R:CM &R:DM\\\\\n")
        for psr in self.psrs:
            c10=c20=c40=0
            for i in range(len(self.curr[psr]['sma'])):
                if self.curr[psr]['lda'][i]==0.1:
                    c10=self.curr[psr]['sma'][i]
                if self.curr[psr]['lda'][i]==0.2:
                    c20=self.curr[psr]['sma'][i]
                if self.curr[psr]['lda'][i]==0.4:
                    c40=self.curr[psr]['sma'][i]
            uwl=array(self.wbr[psr]['sma'][0:6])
            ww= sqrt(1.0/sum(1.0/(uwl*uwl)))
            #ww= sqrt(sum((uwl*uwl)))/float(len(uwl))
            uwl*=1e6
            ww*=1e6
            c10*=1e6
            c20*=1e6
            c40*=1e6
            vals=self.curr[psr]
            C_DM,C_CM = self.getCMDM(vals['sma'],vals['lda'])
            vals=self.wbr[psr]
            W_DM,W_CM = self.getCMDM(vals['sma'],vals['lda'])
            print psr,C_CM,C_DM
            flag=" "
            if W_CM < 30e-9:
                flag="*"

            f.write("%s % 7.3f % 7.3f % 7.3f % 7.3f%s; % 7.3f % 7.3f % 7.3f % 7.3f % 7.3f % 7.3f %4.2f %4.2f % 7.3f % 7.3f\n"%(psr,c10,c20,c40,ww,flag,uwl[0],uwl[1],uwl[2],uwl[3],uwl[4],uwl[5],C_CM/W_CM,C_DM/W_DM,C_CM*1e6,C_DM*1e6))
            ff.write("%s &% 7.3f &% 7.3f &% 7.3f &% 7.3f%s& % 7.3f &% 7.3f &% 7.3f &% 7.3f &% 7.3f &% 7.3f &%4.2f &%4.2f\\\\\n"%(psr,c10,c20,c40,ww,flag,uwl[0],uwl[1],uwl[2],uwl[3],uwl[4],uwl[5],C_CM/W_CM,C_DM/W_DM))

        f.close()

        if noop:
            print "noop"
            sys.exit(0)
        angles=dict()
        pairs=list()
        f=open("angles")
        for line in f:
            elems=line.split()
            p1=elems[0]
            p2=elems[1]
            a=float(elems[2])
            if p1 not in angles:
                angles[p1]=dict()
            if p2 not in angles:
                angles[p2]=dict()
            angles[p1][p2]=a
            angles[p1][p1]=0
            angles[p2][p1]=a
            angles[p2][p2]=0
        

        self.angles=angles
               
        Tcur=10
        Twbr=15
        Ncur=int(Tcur*365.25/self.tsamp)
        Nwbr=int(Twbr*365.25/self.tsamp)
        
        N=Ncur+Nwbr
        print N
        data=zeros((npsr,N*Nreal))
        print "Make GW signal ",Agwb
        i=0
        for psr in self.psrs:
            print psr
            ts = make_red_noise(N,Nreal,-13.0/3.0,self.tsamp/365.25,Agwb)
            data[i]+=ts
            i+=1

        if GWcor:
            self.GWCovar = make_GW_covar(angles,self.psrs)
            C=self.GWCovar
            L = linalg.cholesky(C)
            timeseries = dot(L,data)
        else:
            print "NO GW correlation!!!!!"
            print "NO GW correlation!!!!!"
            print "NO GW correlation!!!!!"
            timeseries=data

        print "RED noise"
        i=0
        for psr in self.psrs:
            if psr in self.red:
                print psr," RED ",self.red[psr]["alpha"],self.red[psr]["fc"],self.red[psr]["amp"]
                ts = self.make_red_noise(N,Nreal,self.red[psr]["alpha"],self.tsamp/365.25,self.red[psr]["amp"],gw=False,fc=self.red[psr]["fc"]) # Make non-GW red noise
                timeseries[i]+=ts
            i+=1


        psrts_cur=dict()
        psrts_wbr=dict()

        psrts_cur_e=dict()
        psrts_wbr_e=dict()

        x=linspace(0,N*self.tsamp,N)+53371
        print "Create timeseries"
        ipsr=0
        for psr in self.psrs:
            vals=self.curr[psr]
            C_DM,C_CM = self.getCMDM(vals['sma'],vals['lda'])
            vals=self.wbr[psr]
            W_DM,W_CM = self.getCMDM(vals['sma'],vals['lda'])
            print psr,C_DM/W_DM,C_CM/W_DM
            W_CM=max(30e-9,W_CM)
            C_CM=max(30e-9,C_CM)
            psrts_cur[psr] = list()
            psrts_wbr[psr] = list()
            psrts_cur_e[psr] = list()
            psrts_wbr_e[psr] = list()
            real=0
            while real < Nreal:
                psrts_cur_e[psr].append(zeros(N))
                psrts_wbr_e[psr].append(zeros(N))
                psrts_cur[psr].append(zeros(N))
                psrts_wbr[psr].append(zeros(N))
                psrts_cur[psr][real] += timeseries[ipsr][real*N:(real+1)*N]
                psrts_wbr[psr][real] += timeseries[ipsr][real*N:(real+1)*N]
                o=0
                while o < Ncur:
                    v=random.gauss(0,C_CM)
                    psrts_cur[psr][real][o] += v
                    psrts_wbr[psr][real][o] += v
                    psrts_cur_e[psr][real][o] = C_CM
                    psrts_wbr_e[psr][real][o] = C_CM
                    o+=1
                while o < N:
                    psrts_cur[psr][real][o] += random.gauss(0,C_CM)
                    psrts_wbr[psr][real][o] += random.gauss(0,W_CM)
                    psrts_cur_e[psr][real][o] = C_CM
                    psrts_wbr_e[psr][real][o] = W_CM
                    o+=1

                real+=1
            ipsr+=1
            
        print "write files"
        real=0
        while real < Nreal:
            try:
                os.mkdir("r.%04d"%real)
            except:
                pass
            real+=1
        ipsr=0
        for psr in self.psrs:
            print psr
            real=0
            while real < Nreal:
                y=psrts_cur[psr][real]
                w=1.0/power(psrts_cur_e[psr][real],2)
                z=polyfit(x,y,2)
                m=x*z[1]+z[2]+x*x*z[0]
                psrts_cur[psr][real]-=m
                y=psrts_wbr[psr][real]
                w=1.0/power(psrts_wbr_e[psr][real],2)
                z=polyfit(x,y,2)
                m=x*z[1]+z[2]+x*x*z[0]
                psrts_wbr[psr][real]-=m

                f_cur=open("r.%04d/%s.cur.cm"%(real,psr),"w")
                f_wbr=open("r.%04d/%s.wbr.cm"%(real,psr),"w")
                f_cur.write("DMMODEL DM 0\n")
                f_wbr.write("DMMODEL DM 0\n")
                o=0
                while o < Ncur:
                    if psr not in self.xtra:
                        f_cur.write("_CM %f %g %g\n"%(x[o],psrts_cur[psr][real][o],psrts_cur_e[psr][real][o]))
                        f_wbr.write("_CM %f %g %g\n"%(x[o],psrts_wbr[psr][real][o],psrts_wbr_e[psr][real][o]))
                    o+=1
                while o < N:
                    if psr not in self.xtra:
                        f_cur.write("_CM %f %g %g\n"%(x[o],psrts_cur[psr][real][o],psrts_cur_e[psr][real][o]))
                    f_wbr.write("_CM %f %g %g\n"%(x[o],psrts_wbr[psr][real][o],psrts_wbr_e[psr][real][o]))
                    o+=1


                f_wbr.close()
                f_cur.close()
                real+=1
            ipsr+=1



#UWL bands
#
#700 - 880  OK
#880 - 960  GSM mobile
#960 - 1260  OK
#1260 - 1300 Galileo etc
#1300 - 1560 OK
#1560 - 1590 Galileo etc
#1590 - 1920 OK
#1920 - 1980  3G mobile
#1980 - 2110 OK
#2110 - 2170 3G mobile
#2170 - 2300 OK
#2300 - 2380 NBN Wireless
#2380 - 4000 OK
#
#Total clear:  2950 = 89.4%
#Total rfi:   350 = 10.6%
#
#Note: some of these will be free some of the time, and conversely there are some bands not listed where RFI is either narrowband or only occasionally present. Also an RFI adaptive filter system could in principle remove some (or even all) of these RFI signals (while leaving the pulsar signal intact).
#
#System temperatures:
#
#Measured values (near Hydra A, from PPTA1 paper)  Tsys = Ssys / 1.6 K/Jy
#Band  Ssys   Tsys  Tsky Tsys-Tsky
#10cm  60 Jy  38 K  0 K   38 K
#20cm  36 Jy  22 K  1 K   21 K 
#50cm  50 Jy  31 K  6 K   25 K
#UWL      -        -        -       21 K
    def curr2wbr(self):

        self.wbr=dict()
        for psr in self.psrs: 
            LREF=0.2
            lda=list()
            sma=list()
            clda=self.curr[psr]['lda']
            csma=self.curr[psr]['sma']
            s10=0
            s20=0
            s40=0
            for l,s in zip(clda,csma):
                if l==0.1:
                    s10=s
                if l==0.2:
                    s20=s
                if l==0.4:
                    s40=s

            if s10==0:
                s10=1
            if s20==0:
                s20=1
                LREF=clda[0]
            if s40==0:
                s40=1
            s10*=sqrt(900)*(21./38.)
            s20*=sqrt(256)
            s40*=sqrt(64)*(28./31.)
            s13=(s10+s20)/2.0
            s17=(s10*0.5 + s20)/1.5
            s30=(s20+s40)/2.0
            print "UWL",psr,s10,s13,s17,s20,s40

            for l in clda:
                if l==LREF or (len(self.xtra) == 0 and l==0.4):
                    # 2380 - 4000
                    sig = s10 / sqrt(1620)
                    lam = 0.094
                    lda.append(lam)
                    sma.append(sig)
                    # 2170 - 2300 AND 1980 - 2110
                    sig = s13 / sqrt(260)
                    lam = 0.13
                    lda.append(lam)
                    sma.append(sig)
                    # 1590 - 1920 
                    sig = s17 / sqrt(330)
                    lam = 0.17
                    lda.append(lam)
                    sma.append(sig)
                    # 1300 - 1560
                    sig = s20 / sqrt(260)
                    lam=0.2
                    lda.append(lam)
                    sma.append(sig)
                    # 960 - 1260
                    sig = s30 / sqrt(300)
                    lam=0.3
                    lda.append(lam)
                    sma.append(sig)
                    #700 - 880
                    sig=s40/sqrt(180)
                    lam=0.4
                    lda.append(lam)
                    sma.append(sig)

            self.wbr[psr]=dict(lda=array(lda),sma=array(sma))
            print psr,len(lda),len(clda)
#            print sma[0:6],lda[0:6]
#            print csma[0:3],clda[0:3]



    def psd2sma(self,d):
        lda=d['lda']
        psd=d['psd']
        rte=d['rte']

        newlda=list()
        sma=list()
        band=0
        while band < len(lda):
            dt=self.tsamp/365.25
            N=rte[band]
            print band,dt,N
            sigma = sqrt(psd[band]*N/(2.0*dt))*365.25*86400.0
            i=0
            while i < N:
                sma.append(sigma)
                newlda.append(lda[band])
                i+=1
            band+=1
        d['lda']=array(newlda)
        d['sma']=array(sma)

    def getCMDM(self,sma, lda):
        if lda[0]==0.0:
            return 0,sma[0]
        LAMBDA_REF=0.2
        lR2=LAMBDA_REF*LAMBDA_REF
        l2=lda*lda
        l4=l2*l2
        s2=sma*sma

        DELTA = sum(1.0/s2) * sum(l4/s2) - pow(sum(l2/s2),2)


        a = lR2 * ((l2/s2)*sum(1/s2) - (1.0/s2) * sum(l2/s2))/DELTA
        b = ((1.0/s2)*sum(l4/s2) - (l2/s2) * sum(l2/s2))/DELTA

        a2=a*a
        b2=b*b

        s2_TDM = sum(a2*s2)
        s2_TCM = sum(b2*s2)

        return sqrt(s2_TDM),sqrt(s2_TCM)




if __name__ == "__main__":
    sim=ifunc_simulator(sys.argv)
