#!/usr/bin/python
from sys import argv,stdout,stderr
from math import *
from numpy import *
from numpy import linalg as lalg
from scipy import linalg as slalg
from matplotlib import pyplot as plt
import random


ALPHA=-13.0/3.0

def detect_GWB(A_guess,infile,path,args):
    A2_guess=A_guess*A_guess

    diagPlots=False
    plotSpec=False
    saveCVM=False
    loadCVM=False
    writeSpec=False
    skipem=list()
    max_channels=10
    start_channel=0
    fmax=10.0
    white_factor=1.0
    red_file=None
    factor_file=None
    ch1f=False
    DEBUG_W=False
    NIT=3
    for arg in args:
        if arg=="-p":
            diagPlots=True
        if arg.startswith("-S"):
            saveCVM=arg[2:]+".npz"
        if arg.startswith("-L"):
            loadCVM=arg[2:]+".npz"
        if arg.startswith("-k"):
            skipem.append(arg[2:])
        if arg.startswith("-F"):
            fmax=float(arg[2:])
        if arg.startswith("-n"):
            max_channels=float(arg[2:])
        if arg.startswith("-s"):
            start_channel=float(arg[2:])
        if arg.startswith("--writespec"):
            writeSpec=True
        if arg.startswith("--plotspec"):
            plotSpec=True
        if arg.startswith("-w"):
            white_factor=float(arg[2:])
        if arg.startswith("-R"):
            red_file=arg[2:]
        if arg.startswith("-f"):
            factor_file=arg[2:]+".npy"
        if arg.startswith("--ch1factor"):
            ch1f=True
        if arg.startswith("-D"):
            DEBUG_W=True
        if arg.startswith("-I"):
            NIT=int(arg[2:])


    pairs=list()
    angles=list()

    power_spectra=dict()
    xspec=list()
    Ts=list()
    white=list()
    Pn=list()
    Pgn=list()
    skip=True
    p1=False
    ii=0
    inf=open(infile)
    for line in inf:
        elems=line.split()
        if line.startswith("#"):
            skip=False
            if ii > 0:
                closeoff(p1,p2,freq,cross_r,cross_i,power_1,power_2,power_spectra,xspec,Pgn,Ts)
            else:
                if p1:
                    print "Skip",p1,p2,"No Freq"
            # now set up for the new pulsar.
            p1=elems[1]
            p2=elems[2]
            if p1 in skipem or p2 in skipem:
                print "Skip",p1,p2
                skip=True
                ii=0
                continue
            a=pi*float(elems[3])/180.0
            pw1=float(elems[4])
            pw2=float(elems[5])
            
            freq=list()
            cross_r=list()
            cross_i=list()
            power_1=list()
            power_2=list()
            ii=0
            iii=0
            continue
        if skip:
            continue

        if ii >= max_channels:
            continue
        if iii < start_channel:
            iii+=1
            continue
        fff=float(elems[0])
        if fff > fmax:
            continue
        if ii==0:
            pairs.append((p1,p2))
            angles.append(a)
            white.append((pw1,pw2))

        elems=line.split()
        freq.append(fff)
        cross_r.append(float(elems[1]))
        cross_i.append(float(elems[2]))
        power_1.append(float(elems[3]))
        power_2.append(float(elems[4]))
        ii+=1


    closeoff(p1,p2,freq,cross_r,cross_i,power_1,power_2,power_spectra,xspec,Pgn,Ts)
    angles=array(angles)
    np=len(pairs)

    # Compute zeta, the expected correlation between pairs.
    x=(1.0-cos(angles))/2.0
    zeta=(3.0/2.0)*x*log(x) - x/4.0 + 0.5


# red the red noise models.
    red=dict()
    if red_file != None:
        fle=open(red_file)
        for line in fle:
            elems=line.split()
            psr=elems[0]
            alp=float(elems[1])
            amp=float(elems[2])
            fc=float(elems[3])
            if not psr in red:
                red[psr]=list()
            print "%s A=%g a=%g fc=%g"%(psr,amp,alp,fc)
            red[psr].append((amp,alp,fc))
        fle.close()

    if writeSpec:
        fspecout=open("allspec","w")

    pwr_models=list()
# Model the power-spectra:
    p=0
    for pair in pairs:
        m1 = pwr_model()
        m1.white=white[p][0]*white_factor
        if pair[0] in red:
            m1.red=red[pair[0]]
        else:
            m1.red=[]

        m2 = pwr_model()
        m2.white=white[p][1]*white_factor
        if pair[1] in red:
            m2.red=red[pair[1]]
        else:
            m2.red=[]

        if writeSpec:
            freq=xspec[p][0]
            mp1=m1.value(freq)
            mp2=m2.value(freq)
            fle=open("%s/%s-%s.ss"%(path,pair[0],pair[1]),"w")
            f=0
            while f < len(freq):
                if writeSpec:
                    fspecout.write("%g\n"%xspec[p][1][f])
                fle.write("%g %g %g %g %g %g %g\n"%(freq[f],xspec[p][1][f],xspec[p][3][f],xspec[p][4][f],mp1[f],mp2[f],A2_guess*Pgn[p][f]))
                f+=1
            fle.close()


        pwr_models.append((m1,m2))
        p+=1

    if writeSpec:
        fspecout.close()




    DM = None

    top=0
    bottom=0
    p=0
    for psr1,psr2 in pairs:
        stderr.write("%s %s   \n"%(psr1,psr2))
        stderr.flush()
        fzero=xspec[p][0][0]
        MX=1000
        r_i = list()
        r_j = list()
        w_i = list()
        w_j = list()
        ff=open("%s.cur.cm"%psr1)
        for line in ff:
            elems=line.split()
            if len(elems) >3 and float(elems[1]) <= 58850:
                r_i.append(float(elems[2])/86400.0/365.25)
                w_i.append(1.0/pow(float(elems[3]),1))
                if len(r_i) >= MX:
                    break
        ff.close()

        ff=open("%s.cur.cm"%psr2)
        for line in ff:
            elems=line.split()
            if len(elems) >3 and float(elems[1]) <= 58850:
                r_j.append(float(elems[2])/86400.0/365.25)
                w_j.append(1.0/pow(float(elems[3]),1))
                if len(r_j) >= MX:
                    break
        ff.close()
        n=len(r_i)
        w_i=array(w_i)
        w_j=array(w_j)
        r_i=array(r_i)
        r_j=array(r_j)

        if DM==None:
            DM=zeros((n,1))
            for i in range(0,n):
                DM[i][0]=1
                #DM[i][1]=60.0*float(i)
                #DM[i][2]=pow(DM[i][1],2)
        A=DM
        ii = lalg.inv(dot(A.T,A))
        R_i = eye(n) - dot(dot(A,ii),A.T)
        R_j = eye(n) - dot(dot(A,ii),A.T)
        #R_i = eye(n)
        #R_j = eye(n)


        nn=256
        V=zeros(nn*2)
        delta = 60.0/365.25
        for r in range(0,nn+1):
            f=r/(nn*2*delta)
            V[r] = P_gw_nrm(f)
            if r > 0:
                V[2*nn-r] = V[r]

        S=real(fft.ifft(V))/delta
        SS=zeros((n,n))
        for i in range(0,n):
            for j in range(0,n):
                SS[i][j] = S[abs(i-j)]

        St = zeta[p] * dot(dot(R_i,SS),R_j.T)

        mi=pwr_models[p][0]
        mj=pwr_models[p][0]
        for r in range(0,nn+1):
            f=r/(nn*2*delta)
            V[r] = mi.value(f) + A2_guess * P_gw_nrm(f)
            if r > 0:
                V[2*nn-r] = V[r]
        P=real(fft.ifft(V))/delta


        Pi=zeros((n,n))
        for i in range(0,n):
            for j in range(0,n):
                Pi[i][j] = P[abs(i-j)]

        Pi = dot(dot(R_i,Pi),R_i.T)
        for r in range(0,nn+1):
            f=r/(nn*2*delta)
            V[r] = mj.value(f)  + A2_guess * P_gw_nrm(f)
            if r > 0:
                V[2*nn-r] = V[r]
        P=real(fft.ifft(V))/delta

        Pj=zeros((n,n))
        for i in range(0,n):
            for j in range(0,n):
                Pj[i][j] = P[abs(i-j)]
        Pj = dot(dot(R_j,Pj),R_j.T)

        #Pi+=eye(n)*1e-40
        #Pj+=eye(n)*1e-40
        Pi_inv = lalg.inv(Pi)
        Pj_inv = lalg.inv(Pj)
        #top += dot(dot(dot(dot(r_i,Pi_inv),St),Pj_inv),r_j)
        A = dot(r_i.T,Pi_inv)
        B = dot (Pj_inv,r_j)
        t = dot(dot(A,St),B)
        b = trace(dot(dot(dot(Pi_inv,St),Pj_inv),St))
        top+=t
        bottom+=b
        print t,b
        print t/b

        p+=1
    

    stderr.write("Done                  \n")
    print "Done                       "
    print "top",top
    print "bottom",bottom
    A2=top.real/bottom
    E=1.0/sqrt(bottom)
    A=0
    if A2 > 0:
        A=sqrt(A2)
    print "A2=",A2,E
    print "A =",A,sqrt(E)
    print "S =",A2/E, A2_guess/E, top.real/sqrt(bottom)
    E=0
    W=0


    return A2,E,W







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



def closeoff(p1,p2,freq,cross_r,cross_i,power_1,power_2,power_spectra,xspec,Pgn,Ts):
    # Close off the previous file.
    freq=array(freq)
    ns=len(freq)
    cross_r=array(cross_r)
    cross_i=array(cross_i)
    power_1=array(power_1)
    power_2=array(power_2)
    xspec.append((freq,cross_r,cross_i,power_1,power_2))
    Ts.append((1.0/freq[0],freq[0]))

    # get the best power spectrum for each pulsar
    if p1 in power_spectra:
        if len(power_spectra[p1][1]) < len(power_1):
            #print p1,p2,"0"
            power_spectra[p1]=[freq,power_1]
    else:
        power_spectra[p1]=[freq,power_1]

    if p2 in power_spectra:
        if len(power_spectra[p2][1]) < len(power_2):
            #print p1,p2,"1"
            power_spectra[p2]=[freq,power_2]
    else:
        power_spectra[p2]=[freq,power_2]

    # The normalised GWB signal
    Pgn.append(P_gw_nrm(freq))


def P_gw_nrm(freq):
    K=power(0.0335955,ALPHA)
    return K*power(1.0 + power(freq/0.0335955,2),ALPHA/2.0)/(12*pi*pi)
    #return power(0.0335955+freq,ALPHA)/(12.0*pi*pi)



def crossW(f,T):
    return T *power(sin(pi*f*T)/(pi*f*T),2.0)

crossCovar_lookup=dict()
def crossCovar(delta,T1,T2,bw):
    key="%g %g %g %g"%(delta,T1,T2,bw)
    if key not in crossCovar_lookup:
        crossCovar_lookup[key] = _crossCovar(delta,T1,T2,bw)
    return crossCovar_lookup[key]

def _crossCovar(delta,T1,T2,bw):
    delta/=2.0
    rnge = bw*2.5
    fs=linspace(-rnge,rnge,30)

    df=fs[1]-fs[0]
    X = min(1.0,sum(crossW(fs+delta,T1)*crossW(fs-delta,T2))*df / (0.666*sqrt(T1*T2)))
    return X


if __name__ == "__main__":

    infile=argv[1]
    path=argv[2]
    A_guess=float(argv[3])


    A2,E,W=detect_GWB(A_guess,infile,path,argv)



    #plt.plot(W)
    #plt.figure()
    #plt.show()



