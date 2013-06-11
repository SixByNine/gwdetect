#!/usr/bin/python
from sys import argv,stdout
from math import *
from numpy import *
from numpy import linalg as lalg
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
    diag=False
    skipem=list()
    max_channels=10
    start_channel=0
    fmax=10.0
    white_factor=1.0
    red_file=None
    factor_file=None
    for arg in args:
        if arg=="-p":
            diagPlots=True
        if arg.startswith("-S"):
            saveCVM=arg[2:]+".npy"
            saveCVMi=arg[2:]+"_i.npy"
        if arg.startswith("-L"):
            loadCVM=arg[2:]+".npy"
            loadCVMi=arg[2:]+"_i.npy"
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
        if arg.startswith("-D"):
            diag=True


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
            fc=float(elems[2])
            amp=float(elems[3])
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

        if saveCVM or writeSpec:
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
    print "Compute estimates of A2_ijk"
    covar_idx=list()
    covar_ridx=list()
    A2ijk=list()


    p=0
    i=0
    while p < np:
        f=0
        covar_ridx.append(list())
        mx=min(len(xspec[p][0]),max_channels)
        while f < mx:
            covar_idx.append((p,f))
            covar_ridx[p].append(i)
            i+=1
            f+=1
        p+=1

    N=len(covar_idx)

    factors=zeros(N)+1.0
    if factor_file != None:
        print "Loaded gain factors"
        factors=load(factor_file)


    p=0
    i=0
    while p < np:
        f=0
        mx=min(len(xspec[p][0]),max_channels)
        while f < mx:
            A2ijk.append(xspec[p][1][f]/zeta[p]/Pgn[p][f]/factors[i]) 
            i+=1
            f+=1
        p+=1

    A2ijk=array(A2ijk)

    print "N =",N

    if loadCVM != False:
        print "Load precomputed Covariance matrix"
        covar=load(loadCVM)
    else:
        covar=zeros((N,N))

        print "Compute GW covariances"

        
        GWcovar=zeros((np,np))
        # compute the GW covariances
        ij=0
        while ij < np:
            psr1=pairs[ij][0]
            psr2=pairs[ij][1]
            print "\r% 4d % 12s % 12s"%(ij,psr1,psr2),
            stdout.flush()
            lm=0
            while lm <= ij:
                if diag and ij!=lm:
                    lm+=1
                    continue
                #psr3=pairs[lm][0]
                #psr4=pairs[lm][1]
                zeta_ij,zeta_lm,zeta_il,zeta_jm,zeta_im,zeta_jl = getZetas(ij,lm,pairs,zeta)
                C=0.5 * pow(A_guess,4) * (zeta_il*zeta_jm + zeta_im*zeta_jl)/(zeta_ij*zeta_lm)
                #print "%s %s %s %s %.4f %.4f %.4f %.4f %.4f %.4f %g"%(psr1,psr2,psr3,psr4,zeta_ij,zeta_lm,zeta_il,zeta_jm,zeta_im,zeta_jl,C)
                GWcovar[ij][lm]=C
                GWcovar[lm][ij]=C
                lm+=1
            ij+=1

        print "Done"
        print "Compute covariance matrix",N
    # Create the covariance matrix for the noise.
        a=0
        while a < N:
            ij=covar_idx[a][0]
            psr1=pairs[ij][0]
            psr2=pairs[ij][1]
            print "\r% 4d % 12s % 12s"%(a,psr1,psr2),
            stdout.flush()
            b=0
            while b <= a:
                lm=covar_idx[b][0]
                psr3=pairs[lm][0]
                psr4=pairs[lm][1]
                f1 = xspec[ij][0][covar_idx[a][1]]
                f2 = xspec[lm][0][covar_idx[b][1]]
                bw=max(Ts[ij][1],Ts[lm][1])
                delta=abs(f1-f2)
                if delta > 2*bw:
                    b+=1
                    continue
#                if not (psr1==psr3 or psr1==psr4 or psr2==psr3 or psr2==psr4):
#                    b+=1
#                    continue

                factor=factors[ij]*factors[lm]
                X=crossCovar(delta,Ts[ij][0],Ts[lm][0],bw)

                #X=0
                #if covar_idx[a][1]==covar_idx[b][1]:
                #    X=1.0
                if X < 0.1:
                    b+=1
                    continue

#                print X,psr1,psr2,psr3,psr4,f1,f2,Ts[ij][1],Ts[lm][1]
                
                # The noise is only relevent for the same pulsar pair.
                if ij==lm:
                    psr1=pairs[ij][0]
                    psr2=pairs[ij][1]
                    fa=(f1+f2)/2.0
                    #fa=min(f1,f2)
                    P1 = pwr_models[ij][0].value(fa)
                    P2 = pwr_models[ij][1].value(fa)

                    C = 0.5 * (P1*P2)/pow(zeta[ij]*P_gw_nrm(fa),2)/factor
                    C += 0.5 * (P1*A2_guess)/pow(zeta[ij]*factor,2)/P_gw_nrm(fa)/factor
                    C += 0.5 * (P2*A2_guess)/pow(zeta[ij]*factor,2)/P_gw_nrm(fa)/factor
                    covar[a][b] += X*C
                    if a!=b:
                        covar[b][a] += X*C

                if diag and ij!=lm:
                    b+=1
                    continue
                C=GWcovar[ij][lm]/factor
                covar[a][b] += X*C
                if a!=b:
                    covar[b][a] += X*C
                b+=1
            a+=1

        print "Done"



        if diag:
            print "DIAG GW"
                        

    if saveCVM != False:
        print "Save Covariance matrix"
        save(saveCVM,covar)


    M=zeros(N)+1
    if loadCVM != False:
        print "Load precomputed inverse of CVM"
        cinv=load(loadCVMi)
    else:
        print "Inverting matrix"
        stdout.flush()
        cinv = lalg.inv(covar)
        print "Done"

    if saveCVM != False:
        print "Save inverse matrix"
        save(saveCVMi,cinv)


    eq1=1.0/dot(dot(M,cinv),M)

    eq2=dot(M,cinv)

    W=eq1*eq2
    print "Weights sum =",sum(W)
    if False:
        fff=open("W","w")
        i=0
        while i < len(W):
            ij=covar_idx[i][0]
            psr1=pairs[ij][0]
            psr2=pairs[ij][1]
            f=covar_idx[i][1]
            fff.write("%d %s %s %g %g %g\n"%(i,psr1,psr2,f,W[i],A2ijk[i]))
            i+=1

        fff.close()

    A2=dot(W,A2ijk)

    varA2 = sum(W*W*power(A2-A2ijk,2))/pow(sum(W),2)
    E = sqrt(varA2)

    MX=A2-A2ijk
    chisq = dot(MX,dot(cinv,MX))
    print "CHISQ(A2)",chisq,chisq/float(len(A2ijk))

    MX=A2_guess-A2ijk
    chisq = dot(MX,dot(cinv,MX))
    print "CHISQ(A2_guess)",chisq,chisq/float(len(A2ijk))

    MX=A2+E-A2ijk
    chisq = dot(MX,dot(cinv,MX))
    print "CHISQ(A2+E)",chisq,chisq/float(len(A2ijk))

    MX=A2-E-A2ijk
    chisq = dot(MX,dot(cinv,MX))
    print "CHISQ(A2-E)",chisq,chisq/float(len(A2ijk))
    
    F=zeros(N)+A2_guess
    ryan_parameter = dot(F,dot(cinv,A2ijk))
    print "RYAN",ryan_parameter

    
    pair_W = list()
#    pair_W = zeros((np,max_channels))
    

    all_freq = list()
    w=0
    p=0
    while p < np:
        mx=min(max_channels,len(xspec[p][1]))
        pair_W.append(zeros(mx))
        f=0
        while f < mx:
            fff=floor(xspec[p][0][f]*100)
            if not fff in all_freq:
                all_freq.append(fff)
            pair_W[p][f] = W[w]
            w+=1
            f+=1
        p+=1

    all_freq.sort()
    nf=len(all_freq)
    ff=zeros(nf)
    fC=zeros(nf)
    fE=zeros(nf)


    p=0
    while p < np:
        mx=min(len(xspec[p][0]),max_channels)
        f=0
        while f < mx:
            fff=floor(xspec[p][0][f]*100)
            fidx= all_freq.index(fff)
            i=covar_ridx[p][f]
            fC[fidx]+=pair_W[p][f]
            fE[fidx]+=pair_W[p][f]*pair_W[p][f]*covar[i][i]
            ff[fidx]+=pair_W[p][f]*xspec[p][1][f]/Pgn[p][f]/zeta[p]
            f+=1
        p+=1

    ff/=fC
    fE /= fC*fC
    fE = sqrt(fE)

    fle=open("%s/ff.plot"%path,"w")
    f=0
    while f < nf:
        fle.write("%g\t%g\t%g\t%g\n"%(float(all_freq[f])/100.0,ff[f],fE[f],fC[f]))
        f+=1

    fle.close()


    pp=zeros(np)
    pC=zeros(np)
    pE=zeros(np)
    p_A=zeros(3)
    p_Aa=zeros(3)
    p_AC=zeros(3)
    p_AE=zeros(3)

    p=0
    while p < np:
        mx=min(len(xspec[p][0]),max_channels)
        a=0
        if angles[p] > pi*49.0/180.0:
            a=1
        if angles[p] > pi*122.0/180.0:
            a=2

        pp[p]=sum(pair_W[p]*xspec[p][1][0:mx]/Pgn[p][0:mx]/zeta[p])
        p_A[a]+=sum(pair_W[p]*xspec[p][1][0:mx]/Pgn[p][0:mx])
        p_Aa[a]+=sum(pair_W[p])*angles[p]
        f=0
        while f < mx:
            i=covar_ridx[p][f]
            pE[p]+=pair_W[p][f]*pair_W[p][f]*covar[i][i]
            p_AE[a]+=pair_W[p][f]*pair_W[p][f]*covar[i][i]*zeta[p]*zeta[p]
            f+=1
        pC[p]=sum(pair_W[p])
        p_AC[a]+=sum(pair_W[p])
        p+=1

    pp/=pC

#    pE = 1.0/(sqrt(abs(pC/eq1)))
    pE /= pC*pC
    pE = sqrt(pE)

    p_A/=p_AC
    p_Aa/=p_AC
    p_AE/=p_AC*p_AC
    p_AE=sqrt(p_AE)

    f=open("%s/hd.plot"%path,"w")
    p=0
    while p < np:
        f.write("%f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%s %s\n"%(180.0*angles[p]/pi,pp[p]*zeta[p],abs(pE[p]*zeta[p]),pp[p],pE[p],A2_guess*zeta[p],A2*zeta[p],pC[p],pairs[p][0],pairs[p][1]))
        p+=1

    f.close()

    print "A2=%g %g"%(A2,E)
    if A2 > 0 :
        print "A=%g %g"%(sqrt(A2),E/(2*sqrt(A2)))
        print "sigma(A2)=%f"%(A2/E)
        print "sigma(A)=%f"%(2*A2/E)



    if diagPlots:
        aa=linspace(0,pi,128)
        x=(1.0-cos(aa))/2.0
        azeta=(3.0/2.0)*x*log(x) - x/4.0 + 0.5

        plt.figure()
        plt.errorbar(180.0*angles/pi,pp*zeta,pE*zeta,fmt='x',color='#ffaaaa')
        plt.plot(180.0*aa/pi,azeta*A2,color='blue')
        plt.plot(180.0*aa/pi,azeta*(A2-E),'--',color='blue')
        plt.plot(180.0*aa/pi,azeta*(A2+E),'--',color='blue')
        plt.plot(180.0*aa/pi,azeta*A2_guess,':',color='grey')
        plt.plot(180.0*aa/pi,azeta*0,color='black')
        plt.errorbar(180.0*p_Aa/pi,p_A,p_AE,fmt='x',color='green')
        v=2*max(A2_guess,abs(A2))
        plt.ylim(-v,v)
        plt.xlim(0,180)
        plt.figure()
        plt.plot(W)
        plt.show()

    return A2,E,W







def getZetas(p1, p2, pairs,zetas):
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
    for pp in pairs:
        if pp[0] == psr_i and pp[1]==psr_l:
            zeta_il = zetas[idx]
        if pp[1] == psr_i and pp[0]==psr_l:
            zeta_il = zetas[idx]

        if pp[0] == psr_j and pp[1]==psr_m:
            zeta_jm = zetas[idx]
        if pp[1] == psr_j and pp[0]==psr_m:
            zeta_jm = zetas[idx]

        if pp[0] == psr_i and pp[1]==psr_m:
            zeta_im = zetas[idx]
        if pp[1] == psr_i and pp[0]==psr_m:
            zeta_im = zetas[idx]

        if pp[0] == psr_j and pp[1]==psr_l:
            zeta_jl = zetas[idx]
        if pp[1] == psr_j and pp[0]==psr_l:
            zeta_jl = zetas[idx]

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
            ret+= red[0]/power(1.0+power(freq/red[2],2.0),red[1]/2.0)
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
    return power(freq,ALPHA)/(12.0*pi*pi)


def crossW(f,T):
    return T *power(sin(pi*f*T)/(pi*f*T),2.0)

def crossCovar(delta,T1,T2,bw):
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



