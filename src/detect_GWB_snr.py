#!/usr/bin/python
from sys import argv,stdout
from math import *
from numpy import *
from numpy import linalg as lalg
from scipy import linalg as slalg
from matplotlib import pyplot as plt
import random


ALPHA=-13.0/3.0
#GW_FC=0.0335955
GW_FC=1.0/40.0
GW_K=power(GW_FC,ALPHA)

def gw_snr(A_guess,T,np,sig):
    A2_guess = A_guess*A_guess
    NF=10
    MAXPSR=np
    N=20*T
    sig*=1e-9
    white=2*T*sig*sig/float(N)/pow(86400.0*365.25,2)
    f=open("angles.r")

    angles=list()
    pairs=list()
    psrs=list()
    for line in f:
        e=line.split()
        p1=e[0]
        p2=e[1]
        if not p1 in psrs:
            if len(psrs) < MAXPSR:
                psrs.append(p1)
            else:
                continue
        if not p2 in psrs:
            if len(psrs) < MAXPSR:
                psrs.append(p2)
            else:
                continue
        theta=radians(float(e[2]))
        angles.append(theta)
        pairs.append((p1,p2))

    angles=array(angles)
    npsr=len(psrs)
    np=len(angles)
    x=(1.0-cos(angles))/2.0
    zeta=(3.0/2.0)*x*log(x) - x/4.0 + 0.5

    model=pwr_model()
    model.white=white


    GWcovar=zeros((np,np))
    # compute the GW covariances
    ij=0
    while ij < np:
        psr1=pairs[ij][0]
        psr2=pairs[ij][1]
        stdout.flush()
        lm=0
        while lm <= ij:
            zeta_ij,zeta_lm,zeta_il,zeta_jm,zeta_im,zeta_jl = getZetas(ij,lm,pairs,zeta)
            C=0.5 * pow(A_guess,4) * (zeta_il*zeta_jm + zeta_im*zeta_jl)/(zeta_ij*zeta_lm)
            GWcovar[ij][lm]=C
            GWcovar[lm][ij]=C
            lm+=1
        ij+=1


    C = zeros((np*(NF-1),np*(NF-1)))
    fzero=1.0/T

    a=0
    f=1
    while f < NF:
        freq=f*fzero
        p1=0
        N = model.value(freq)
        while p1 < np:
            i=pairs[p1][0]
            j=pairs[p1][1]
            p2=0
            while p2 <= p1:
                zeta_ij,zeta_pq,zeta_ip,zeta_jq,zeta_iq,zeta_jp = getZetas(p1,p2,pairs,zeta)
                p=pairs[p2][0]
                q=pairs[p2][1]
                a = (f-1)*np + p1
                b = (f-1)*np + p2
                C[a][b] += GWcovar[p1][p2]
                if j==q:
                    C[a][b] += A2_guess * zeta_iq*N/(zeta_ij*zeta_pq*P_gw_nrm(freq))
                if i==q:
                    C[a][b] += A2_guess * zeta_jp*N/(zeta_ij*zeta_pq*P_gw_nrm(freq))
                if i==p:
                    C[a][b] += A2_guess * zeta_jq*N/(zeta_ij*zeta_pq*P_gw_nrm(freq))
                if j==q:
                    C[a][b] += A2_guess * zeta_ip*N/(zeta_ij*zeta_pq*P_gw_nrm(freq))
                if p1==p2:
                    C[a][b] += 0.5 * N*N/pow(zeta_ij*P_gw_nrm(freq),2)
                else:
                    C[b][a] = C [a][b]
                
                p2+=1
            p1+=1
        f+=1
 

    d= diagonal(C)
    sm=0
    p=0
    while p < np:
        f=1
        while f < NF:
            freq=f*fzero
            P = model.value(freq) + A2_guess*P_gw_nrm(freq)
            Vijk = 0.5 * P*P/pow(zeta[p]*P_gw_nrm(freq),2)
            sm+=1.0 / Vijk
            f+=1
        p+=1

    var=1.0/sm
    sn=A2_guess/sqrt(var)

    

    M=zeros(np*(NF-1))+1
    eq2=slalg.solve(C,M,sym_pos=True,check_finite=False,lower=True)
    eq1=1.0/dot(eq2,M)
    var2 = 1.0/sum(eq2)


    sn2=A2_guess/sqrt(var2)


    rms = sqrt(white * N / (2*T))*86400.0*365.25*1e9
    print "%d %d %d %.1f %.2f %.2f %.2g %.2g"%(T,npsr,np,rms, sn,sn2, sqrt(var),var)





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
        self.red=[]

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
    return GW_K*power(1.0 + power(freq/GW_FC,2),ALPHA/2.0)/(12*pi*pi)
    #return power(freq,ALPHA)/(12.0*pi*pi)


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

    A_guess=float(argv[1])
    T=float(argv[2])
    np=float(argv[3])
    w=float(argv[4])


    gw_snr(A_guess,T,np,w)



    #plt.plot(W)
    #plt.figure()
    #plt.show()



