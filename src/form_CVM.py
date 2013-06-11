#!/usr/bin/python
from sys import argv,stdout
from math import *
from numpy import *
from numpy import linalg as lalg
from matplotlib import pyplot as plt
import random

ALPHA=-13.0/3.0

def formCVM(args):
    diagPlots=False
    saveCVM=False
    loadCVM=False
    writeSpec=False
    skipem=list()
    files=list()
    max_channels=10
    start_channel=0
    fmax=10.0
    white_factor=1.0
    red_file=None
    factor_file=None
    A=float(args[0])
    A2=A*A
    for arg in args[1:]:
        if arg=="-p":
            diagPlots=True
        elif arg.startswith("-L"):
            loadCVM=arg[2:]+".npy"
            loadCVMi=arg[2:]+"_i.npy"
        elif arg.startswith("-S"):
            saveCVM=arg[2:]+".npy"
            saveCVMi=arg[2:]+"_i.npy"
        elif arg.startswith("-k"):
            skipem.append(arg[2:])
        elif arg.startswith("-F"):
            fmax=float(arg[2:])
        elif arg.startswith("-n"):
            max_channels=float(arg[2:])
        elif arg.startswith("-s"):
            start_channel=float(arg[2:])
        elif arg.startswith("--writespec"):
            writeSpec=True
        elif arg.startswith("-f"):
            factor_file=arg[2:]+".npy"
        else:
            files.append(arg)
    

    spectra=list()

    J=5
    for fle in files:
        print "Loading % 30s\r"%fle,
        spectra.append(readSpec(start_channel,max_channels,fmax,skipem,fle))

    print "Loaded all data"
    nr = len(spectra)
    np = len(spectra[0])
    # Go pairwise.
    N=0
    p=0
    while p < np:
        N+=len(spectra[0][p][0])
        p+=1

    factors=zeros(N)+1.0
    if factor_file != None:
        factors=load(factor_file)

    print "N =",N
    M=zeros((N,N))
    sumspec=zeros(N)
    r=0
    while r < nr:
        print "Compute means % 4d\r"%r,
        i=0
        p=0
        while p < np:
            nf=len(spectra[r][p][J])
            f=0
            while f < nf:
                spectra[r][p][J][f]/=factors[i]
                sumspec[i]+=spectra[r][p][J][f]
                f+=1
                i+=1
            p+=1
        r+=1
    meanspec=sumspec/float(nr)
    if loadCVM != False:
        print "Load precomputed Covariance matrix"
        M=load(loadCVM)
    else:


        print ""
        r=0
        while r < nr:
            print "Compute CVM % 4d\r"%r,
            spec=zeros(N)
            i=0
            p=0
            while p < np:
                nf=len(spectra[r][p][J])
                f=0
                while f < nf:
                    spec[i]=spectra[r][p][J][f]
                    f+=1
                    i+=1
                p+=1
            spec-=meanspec
            M+=outer(spec,spec)
            r+=1

        print ""

        M/=float(nr)


    V=diagonal(M)
    C=zeros((N,N))


    
    print ""
    print "Invert"
    cinv=lalg.inv(M)
    print "Done"

    if saveCVM != False:
        print "Save Covariance matrix"
        save(saveCVM,M)
    if saveCVM != False:
        print "Inverting Matrix"
        print "Save inverse matrix"
        save(saveCVMi,cinv)

    M=zeros(N)+1
    eq1=1.0/dot(dot(M,cinv),M)
    eq2=dot(M,cinv)
    W=eq1*eq2
    newfactors=zeros(N)
    """    p=0
    i=0
    while p < np:
        f=0
        nf=len(spectra[0][p][0])
        while f < nf:
            if f==0:
                newfactors[i]=0.66
            f+=1
            i+=1
        p+=1"""

    WW=A2/meanspec
    WW/=sum(WW)
    newfactors=factors*sum(W*meanspec/A2)

    if factor_file==None:
        save("factors",newfactors)
    print sqrt(sum(W*meanspec)),W*meanspec,mean(factors),mean(newfactors)


    sumspec=zeros(N*2)
    r=0
    while r < nr:
        i=0
        p=0
        while p < np:
            nf=len(spectra[r][p][J])
            f=0
            while f < nf:
                sumspec[i]+=spectra[r][p][2][f]
                i+=1
                f+=1
            f=0
            while f < nf:
                sumspec[i]+=spectra[r][p][3][f]
                i+=1
                f+=1
            p+=1
        r+=1
    meanpwr=sumspec/float(nr)
    save("meanpwr",meanpwr)


    if diagPlots:
        print ""
        print "Form correlation matrix"
        i=0
        while i < N:
            print "Compute C % 4d\r"%i,
            j=0
            while j < N:
                C[i][j]=M[i][j]/sqrt(V[i]*V[j])
                j+=1
            i+=1


        print mean(meanspec)
        print sum(W*meanspec)
        plt.plot(meanspec)
        plt.title("%.4g -- %.4g"%((sum(W*meanspec),sqrt(sum(W*meanspec)))))
        plt.figure()


        plt.imshow(abs(C))
        i=0
        p=0
        while p < np:
            nf=len(spectra[0][p][J])
            i+=nf
            #plt.axvline(i,0,N,color='grey')
            #plt.axhline(i,0,N,color='grey')
            p+=1


        plt.xlim(0,N)
        plt.ylim(0,N)
    #    plt.figure()
    #    plt.plot(newfactors)
    #    plt.plot(factors)
    #    plt.errorbar(range(len(factors)),newfactors,nf_e)
        plt.show()


def readSpec(start_channel,max_channels,fmax,skipem,infile):

    pairs=list()
    angles=list()

    power_spectra=dict()
    xspec=list()
    Ts=list()
    white=list()
    skip=True
    p1=False
    ii=0
    inf=open(infile)
    for line in inf:
        elems=line.split()
        if line.startswith("#"):
            skip=False
            if ii > 0:
                closeoff(p1,p2,freq,cross_r,cross_i,power_1,power_2,power_spectra,xspec,Ts,zeta)
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
            x=(1.0-cos(a))/2.0
            zeta=(3.0/2.0)*x*log(x) - x/4.0 + 0.5
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


    closeoff(p1,p2,freq,cross_r,cross_i,power_1,power_2,power_spectra,xspec,Ts,zeta)
    angles=array(angles)
    np=len(pairs)
    
    return xspec





def closeoff(p1,p2,freq,cross_r,cross_i,power_1,power_2,power_spectra,xspec,Ts,zeta):
    # Close off the previous file.
    freq=array(freq)
    ns=len(freq)
    cross_r=array(cross_r)
    cross_i=array(cross_i)
    power_1=array(power_1)
    power_2=array(power_2)
#    factor=cross_r*0+1.0
#    factor[0]=0.66
#    factor[1]=0.9
    A2ijk = (cross_r/P_gw_nrm(freq)/zeta)
#    A2ijk = power_1/factor

    xspec.append((freq,cross_r,cross_i,power_1,power_2,A2ijk,zeta))
    Ts.append((1.0/freq[0],freq[0]))

def P_gw_nrm(freq):
    return power(freq,ALPHA)/(12.0*pi*pi)



if __name__ == "__main__":

    formCVM(argv[1:])


