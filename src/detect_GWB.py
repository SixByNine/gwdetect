#!/usr/bin/env python
import gc
import argparse
from sys import stdout
import os
from math import pi,  floor
from numpy import zeros, array, power,dot, linspace, log,cos,sin,diag,diagonal,eye,argmin,argmax, abs, amin,amax, transpose,mean,var,sqrt
import numpy
from scipy import sparse
#from numpy import linalg as lalg
from scipy import linalg as lalg
from matplotlib import pyplot as plt

class SpecAddr:
    def __init__(self,name):
        self.Nadd=0
        self.name=name

class GravWavBkgrdDetect:
    def __init__(self,A_guess):
        self.ALPHA=-13.0/3.0
        self.GW_FC=0.005
        self.GW_K=power(self.GW_FC,self.ALPHA)
        self.zeta_lookup=dict()
        self.crossCovar_lookup=dict()
        self.A2_guess=A_guess*A_guess
        self.A_guess=A_guess
        self.XCOL=1
        self.factors=zeros(1000)+1
        self.max_channels=10
        self.little_as=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        self.skipem=list()
        self.fmax = 1000.0
        self.cinvs=None
        self.start_channel=0
        self.skipped_pairs = list()
        self.QUICK=False
        self.Nadd=0
        self.corn=None
        self.wf=dict()
        self.yrskip=0.1
    
    def add_spec(self,addr):
        if addr.Nadd==0:
            # Initialise the arrays
            N=0
            ip=0
            for p in self.pairs:
                N+=len(self.xspec[ip][0])*2
            addr.av2 = zeros(N)
            addr.avg_spec = zeros(N)
            addr.avg2_spec = zeros(N)
            addr.spec_copy = list()
            addr.A2copy = list()
            addr.model_spec = zeros(N)
            addr.m2 = zeros(N)
            ip=0
            i=0
            for p in self.pairs:
                for m in self.pwr_models[ip]:
                    for freq in self.xspec[ip][0]:
                        G=self.A2_guess * self.P_gw_nrm(freq)
                        P= m.value(freq)
                        addr.model_spec[i] += P+G
                        i+=1
                ip+=1
            ip=0
            i=0
            for p in self.pairs:
                for freq in self.xspec[ip][0]:
                    m1 = self.pwr_models[ip][0]
                    m2 = self.pwr_models[ip][1]
                    G1=self.A2_guess * self.P_gw_nrm(freq)
                    P1= m1.value(freq)
                    G2=self.A2_guess * self.P_gw_nrm(freq)
                    P2= m2.value(freq)
                    addr.m2[i] += sqrt((P1+G1)*(P2+G2))
                    i+=1
                ip+=1
        ip=0
        i=0
        for p in self.pairs:
            j=0
            for freq in self.xspec[ip][0]:
                addr.av2[i] += sqrt(self.xspec[ip][3][j]*self.xspec[ip][4][j])
                i+=1
                j+=1
            ip+=1
        addr.spec_copy.append(zeros(len(addr.avg_spec)))
        addr.A2copy.append(zeros(len(self.A2ijk)))
        ip=0
        i=0
        for p in self.pairs:
            for idx in [3,4]:
                for v in self.xspec[ip][idx]:
                    addr.avg_spec[i] += v
                    addr.avg2_spec[i] += v*v
                    addr.spec_copy[-1][i] = v
                    i+=1
            ip+=1
        i=0
        for v in self.A2ijk:
            addr.A2copy[-1][i]=v
            i+=1
        addr.Nadd+=1
    
    def write_avspec(self,addr,outpath="."):
        f=open("%s/white.%s"%(outpath,addr.name),"w")
        N=float(addr.Nadd)
        ip=0
        i=0
        wf=dict()
        for p1,p2 in self.pairs:
            for p in [p1,p2]:
                if p not in wf:
                    wf[p]=list()
            for m,psr in zip(self.pwr_models[ip],[p1,p2]):
                nf=len(self.xspec[ip][0])
                for ii in range(nf):
                    if ii == nf-1:
                        v=addr.avg_spec[i]/N
                        m=addr.model_spec[i]
                        wf[psr].append(v/m)
                    i+=1

            ip+=1

        for psr in wf:
            f.write("%s %f\n"%(psr,mean(wf[psr])))
            wf[psr] = mean(wf[psr])
        f.close()

        new_model_spec = zeros(len(addr.avg_spec))
        fff = zeros(len(addr.avg_spec))

        ppp=list()
        ip=0
        i=0
        for p1,p2 in self.pairs:
            for m,psr in zip(self.pwr_models[ip],[p1,p2]):
                for freq in self.xspec[ip][0]:
                    newm = pwr_model()
                    newm.white=m.white*wf[psr]
                    newm.red=m.red
                    G=self.A2_guess * self.P_gw_nrm(freq)
                    P= newm.value(freq)
                    new_model_spec[i] += P+G
                    fff[i] = freq
                    ppp.append("%s-%s"%(p1,p2))
                    i+=1
            ip+=1

        f=open("%s/avspec.%s"%(outpath,addr.name),"w")
        N=float(addr.Nadd)
        specs=zeros((addr.Nadd,len(addr.avg_spec)))
        for i in range(addr.Nadd):
            for j in range(len(addr.avg_spec)):
                specs[i][j] = addr.spec_copy[i][j]

        specs= specs.T

        for i in range(len(addr.avg_spec)):
            v=addr.avg_spec[i]
            err=sqrt(addr.avg2_spec[i]/N - v*v/N/N)
            m=addr.model_spec[i]
            new_m=new_model_spec[i]
            f.write("%g %g %g %g %s %g %g %g\n"%(v/N,m,new_m,fff[i],ppp[i],mean(specs[i]),sqrt(var(specs[i])/N),err/sqrt(N)))
        f.close()

        f=open("%s/a2avspec.%s"%(outpath,addr.name),"w")
        ms=mean(addr.A2copy,0)
        vs=var(addr.A2copy,0)
        for i in range(len(ms)):
            p=int(i/self.max_channels)
            f.write("%g %g %g %g %s %s %g %g %g\n"%(ms[i],vs[i],sqrt(vs[i]),self.A2_guess/ms[i],self.pairs[p][0],self.pairs[p][1],self.xspec[p][0][i%self.max_channels],addr.av2[i]/N,addr.m2[i]))
        f.close()


        f.close()
    
    def getA2(self,ai):
        M=zeros(self.N)+1.0
        
        self.Wijk = self.Ww[ai]
        
        if len(self.Wijk) != len(self.A2ijk):
            print ""
            print "ai=%d len(wijk)=%d len(a2ijk)=%d"%(ai,len(self.Wijk),len(self.A2ijk))
            print "ERROR"
        
        self.A2 = dot(self.Wijk,self.A2ijk)
        self.varA2 = sum(self.Wijk*self.Wijk*power(self.A2ijk,2))/pow(sum(self.Wijk),2)
        
        if self.QUICK:
            self.chisq=1
            self.chisq2=2
        else:
            self.R = self.A2ijk - self.A2

            C=self.cinvs[ai]
            
            self.chisq = dot(dot(self.R.T,C),self.R) / (self.N-1)
            
            R2 = self.A2ijk - self.A2_guess
            
            self.chisq2 = dot(dot(R2.T,self.cinvs[ai]),R2) / (self.N-1)
        
        return self.A2,self.varA2, self.chisq, self.chisq2
    
    def clear_big_matricies(self):
        self.uinvs=None
        self.cinvs=None
        return
    
    def get_whitened_residuals(self,ai):
        white_res = dot(self.uinvs[ai].T,self.R)
        return white_res
    
    
    def form_Uinvs(self):
        self.uinvs=list()
        for a in self.little_as:
            n=1-a
            covar = self.covar_k \
                  + self.covar_aa * a * a  \
                  + self.covar_a  * a      \
                  + self.covar_n  * n      \
                  + self.covar_nn * n * n  \
                  + self.covar_an * a * n
            covar_full = covar + covar.T - diag(covar.diagonal())
    
            U = lalg.cholesky(covar_full)
            Uinv = lalg.inv(U)
            self.uinvs.append(Uinv)
        return
    
    def clear_covar(self):
        self.covar_k=None
        self.covar_aa=None
        self.covar_a=None
        self.covar_n=None
        self.covar_nn=None
        self.covar_an=None
        gc.collect()

    def form_Cinvs(self,a2varfix=None):
        self.cinvs=list()
        self.Ww = list()
        self.theory_var = list()
        self.eq1 = list()
        self.covarDiag=list()
        for a in self.little_as:
            print a
            gc.collect()
            n=1-a
            covar = self.covar_k \
                  + self.covar_aa * a * a  \
                  + self.covar_a  * a      \
                  + self.covar_n  * n      \
                  + self.covar_nn * n * n  \
                  + self.covar_an * a * n
            if a2varfix !=None:
                covar *= diag(a2varfix)
                  
            self.covarDiag.append(covar.diagonal())
            M=zeros(self.N)+1
            if self.QUICK:
                eq2= lalg.solve(covar,M,sym_pos=True,check_finite=False,lower=True,overwrite_a=True,overwrite_b=True)
            else:
                Cinv = lalg.solve(covar,eye(self.N),sym_pos=True,check_finite=False,lower=True,overwrite_a=True, overwrite_b=True)
                eq2=dot(Cinv,M)
                self.cinvs.append(Cinv)
            ### TODO: Check that this is the same as inverting after filling in the upper triangle!!! TODO

            
            M=zeros(self.N)+1
            eq1=1.0/dot(M,eq2)
            self.Ww.append(eq1*eq2)
            self.eq1.append(eq1)
            self.theory_var.append(sqrt(1.0/sum(eq2)))

            # Encourage garbage collection
            Cinv=None
            covar=None
    
    def get_A2ijks(self):
        if self.corn != None:
            self.correct_xspec()
        A2ijk=list()
        p=0
        i=0
        while p < self.np:
            f=0
            mx=min(len(self.xspec[p][0]),self.max_channels)
            while f < mx:
                A2ijk.append(self.xspec[p][self.XCOL][f]/self.zeta[p]/self.Pgn[p][f]/self.factors[f]) 
                i+=1
                f+=1
                
            p+=1
    

        self.A2ijk=array(A2ijk)
        
        
    def correct_xspec(self):
        if self.xspec_corrected:
            return
        i=0
        ip=0
        for p in self.pairs:
            nf=len(self.xspec[ip][0])
            for ii in range(nf):
                F=sqrt(self.corn[i]*self.corn[i+nf])
                self.xspec[ip][1][ii] *= F
                self.xspec[ip][2][ii] *= F
                self.xspec[ip][3][ii] *= self.corn[i]
                self.xspec[ip][4][ii] *= self.corn[i+nf]
                i+=1
            i+=nf
            ip+=1
        self.xspec_corrected=True
        
    def model_noise(self,red_file=None,white_factor=1.0):
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
            
        pwr_models=list()
        p=0
        for pair in self.pairs:
            wf1=1.
            wf2=1.
            if pair[0] in self.wf:
                wf1=self.wf[pair[0]]
            m1 = pwr_model()
            m1.white=self.white[p][0]*white_factor*wf1
            if pair[0] in red:
                m1.red=red[pair[0]]
            else:
                m1.red=[]
    
            if pair[1] in self.wf:
                wf2=self.wf[pair[1]]
            m2 = pwr_model()
            m2.white=self.white[p][1]*white_factor*wf2
            if pair[1] in red:
                m2.red=red[pair[1]]
            else:
                m2.red=[]
    
            pwr_models.append((m1,m2))
            p+=1    
        self.pwr_models=pwr_models
   
    def compute_covaraince_matrix(self):
        covar_idx=list()
        covar_ridx=list()
        
        np=self.np
        p=0
        i=0
        while p < np:
            f=0
            covar_ridx.append(list())
            mx=min(len(self.xspec[p][0]),self.max_channels)
            while f < mx:
                covar_idx.append((p,f))
                covar_ridx[p].append(i)
                i+=1
                f+=1
            p+=1
        self.covar_idx=covar_idx
        N=len(covar_idx)
        self.N=N
        pairs=self.pairs
        zeta=self.zeta
        
        
        covar_a = zeros((N, N))
        covar_aa = zeros((N, N))
        covar_k = zeros((N, N))
        covar_nn = zeros((N, N))
        covar_an = zeros((N, N))
        covar_n = zeros((N, N))
        print "Compute GW covariances"
        GWcovar = zeros((np, np))
    # compute the GW covariances
        ij = 0
        while ij < np:
            psr1 = pairs[ij][0]
            psr2 = pairs[ij][1]
            print "\r% 4d % 12s % 12s" % (ij, psr1, psr2),
            stdout.flush()
            lm = 0
            while lm <= ij: #psr3=pairs[lm][0]
        #psr4=pairs[lm][1]
                zeta_ij, zeta_lm, zeta_il, zeta_jm, zeta_im, zeta_jl = self.getZetas(ij, lm, pairs, zeta)
                #C = 0.5 * pow(self.A_guess, 4) * (zeta_il * zeta_jm + zeta_im * zeta_jl) / (zeta_ij * zeta_lm)
                C = 0.5 * pow(self.A_guess, 4) * (zeta_il * zeta_jm + zeta_im * zeta_jl) / (zeta_ij * zeta_lm)

                #print "%s %s %s %s %.4f %.4f %.4f %.4f %.4f %.4f %g"%(psr1,psr2,psr3,psr4,zeta_ij,zeta_lm,zeta_il,zeta_jm,zeta_im,zeta_jl,C)
                GWcovar[ij][lm] = C
                GWcovar[lm][ij] = C
                lm += 1
            
            ij += 1
        
        print "Compute covariance matrix", N
        # Create the covariance matrix for the noise.
        a = 0
        while a < N:
            ij = covar_idx[a][0]
            psr1 = pairs[ij][0]
            psr2 = pairs[ij][1] 
            print "\r% 4d % 12s % 12s"%(a,psr1,psr2),
            stdout.flush()
            b = 0
            while b <= a:
                lm = covar_idx[b][0]
                psr3 = pairs[lm][0]
                psr4 = pairs[lm][1]
                if1 = covar_idx[a][1]
                if2 = covar_idx[b][1]
                f1 = self.xspec[ij][0][if1]
                f2 = self.xspec[lm][0][if2]
                bw = max(self.Ts[ij][1], self.Ts[lm][1])
                delta = abs(f1 - f2)
                X = 0
                if delta < 2 * bw:
                    X = self.crossCovar(delta, self.Ts[ij][0], self.Ts[lm][0], bw)
                if X < 0.01:
                    b += 1
                    continue
                
                # noise only corrlates if at least one pulsar appears twice
                if psr2 == psr3 or psr1 == psr4 or psr1 == psr3 or psr2 == psr4:
                    fa = (f1 + f2) / 2.0
                    zeta_ij, zeta_lm, zeta_il, zeta_jm, zeta_im, zeta_jl = self.getZetas(ij, lm, pairs, zeta)
                #fa=min(f1,f2)
                    P1 = self.pwr_models[ij][0].value(fa)
                    P2 = self.pwr_models[ij][1].value(fa)
                    C = 0
                    D = 0
                    if psr2 == psr3:
                        C += 0.5 * (P2 * self.A2_guess) * zeta_im / zeta_ij / zeta_lm / self.P_gw_nrm(fa)
                        D += 0.5 * (self.A2_guess * self.A2_guess) * zeta_im / zeta_ij / zeta_lm
                    if psr1 == psr4:
                        C += 0.5 * (P1 * self.A2_guess) * zeta_jl / zeta_ij / zeta_lm / self.P_gw_nrm(fa)
                        D += 0.5 * (self.A2_guess * self.A2_guess) * zeta_jl / zeta_ij / zeta_lm
                    if psr1 == psr3:
                        C += 0.5 * (P1 * self.A2_guess) * zeta_jm / zeta_ij / zeta_lm / self.P_gw_nrm(fa)
                        D += 0.5 * (self.A2_guess * self.A2_guess) * zeta_jm / zeta_ij / zeta_lm
                    if psr2 == psr4:
                        C += 0.5 * (P2 * self.A2_guess) * zeta_il / zeta_ij / zeta_lm / self.P_gw_nrm(fa)
                        D += 0.5 * (self.A2_guess * self.A2_guess) * zeta_il / zeta_ij / zeta_lm
                    
                    covar_an[a][b] += X * D
                    covar_a[a][b] += X * C
                    if ij == lm:
                        G = 0.5 * (P1 * P2) / pow(zeta[ij] * self.P_gw_nrm(fa), 2)
                        F = 2 * 0.5 * (P1 * self.A2_guess) / pow(zeta[ij],2) / self.P_gw_nrm(fa)
                        E = 0.5 * (self.A2_guess * self.A2_guess) / pow(zeta[ij],2)
                        covar_k[a][b] += X * G
                        covar_n[a][b] += X * F
                        covar_nn[a][b] += X * E
                        
                    
                        
                C = GWcovar[ij][lm]
                covar_aa[a][b] += X * C
                b += 1
            
            a += 1
        
        print "Done"
        

        
        self.covar_aa = covar_aa
        self.covar_a = covar_a
        self.covar_nn = covar_nn
        self.covar_n = covar_n
        self.covar_an = covar_an
        self.covar_k = covar_k
        self.covar_ridx = covar_ridx
        
    def getZetas(self,p1, p2, pairs,zetas):
        key="%d %d"%(p1,p2)
        if key not in self.zeta_lookup:
            self.zeta_lookup[key] = self._getZetas(p1,p2,pairs,zetas)
        return self.zeta_lookup[key]
    
    def _getZetas(self,p1, p2, pairs,zetas):
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
    
    


    def P_gw_nrm(self,freq):
        return self.GW_K*power(1.0 + power(freq/self.GW_FC,2),self.ALPHA/2.0)/(12*pi*pi)
        #return power(freq,ALPHA)/(12.0*pi*pi)
    
    
    def crossW(self,f,T):
        return T *power(sin(pi*f*T)/(pi*f*T),2.0)
    
    def crossCovar(self,delta,T1,T2,bw):
        key="%g %g %g %g"%(delta,T1,T2,bw)
        if key not in self.crossCovar_lookup:
            self.crossCovar_lookup[key] = self._crossCovar(delta,T1,T2,bw)
        return self.crossCovar_lookup[key]
    
    def _crossCovar(self,delta,T1,T2,bw):
        delta/=2.0
        rnge = bw*2.5
        fs=linspace(-rnge,rnge,30)
    
        df=fs[1]-fs[0]
        X = min(1.0,sum(self.crossW(fs+delta,T1)*self.crossW(fs-delta,T2))*df / (0.666*sqrt(T1*T2)))
        return X
    

    def read_xspec(self,infile):
        self.xspec_corrected=False
        self.pairs=None
        self.xspec=None
        self.Pgn=None
        self.angles=None
        self.zeta=None
        gc.collect()

        pairs = list()
        angles = list()
        power_spectra = dict()
        xspec = list()
        Ts = list()
        white = list()
        Pn = list()
        Pgn = list()
        skip = True
        p1 = False
        p2=freq=cross_r=cross_i=power_1=power_2=None
        ii = 0
        root=os.path.dirname(infile)+"/"
        inf = open(infile)
        for aline in inf:
            data=[aline]
            if aline.startswith("I"):
                elems = aline.split()
                FF=open(root+elems[1])
                data=FF.readlines()
                FF.close()
            for line in data:
                elems=line.split()
                if line.startswith("#"):
                    skip = False
                    if ii > 0:
                        self.closeoff(p1, p2, freq, cross_r, cross_i, power_1, power_2, power_spectra, xspec, Pgn, Ts)
                    elif p1:
                        pr="%s-%s"%(p1,p2)
                        if pr not in self.skipped_pairs:
                            print "Skip", pr
                            self.skipped_pairs.append(pr)
                    # now set up for the new pulsar.
                    p1 = elems[1]
                    p2 = elems[2]
                    if p1 in self.skipem or p2 in self.skipem:
                        skip = True
                        ii = 0
                        continue
                    a = pi * float(elems[3]) / 180.0
                    pw1 = float(elems[4])
                    pw2 = float(elems[5])
                    freq = list()
                    cross_r = list()
                    cross_i = list()
                    power_1 = list()
                    power_2 = list()
                    ii = 0
                    iii = 0
                    continue
                if skip:
                    continue
                if ii >= self.max_channels:
                    continue
                if iii < self.start_channel:
                    iii += 1
                    continue
                fff = float(elems[0])
                if fff > self.fmax or abs(fff-1.0) < self.yrskip:
                    continue
                if ii == 0:
                    pairs.append((p1, p2))
                    angles.append(a)
                    white.append((pw1, pw2))
                elems = line.split()
                freq.append(fff)
                cross_r.append(float(elems[1]))
                cross_i.append(float(elems[2]))
                power_1.append(float(elems[3]))
                power_2.append(float(elems[4]))
                ii += 1
            
        self.closeoff(p1, p2, freq, cross_r, cross_i, power_1, power_2, power_spectra, xspec, Pgn, Ts)
        angles = array(angles)
        self.np = len(pairs)
        
        self.pairs=pairs
        self.xspec=xspec
        self.Pgn=Pgn
        self.angles=angles
        self.Ts=Ts
        self.white=white
        x=(1.0-cos(angles))/2.0
        self.zeta=(3.0/2.0)*x*log(x) - x/4.0 + 0.5
        return

    def write_plotfiles(self, outpath, ia):
        np=self.np
        zeta=self.zeta
        pairs=self.pairs
        W = self.Ww[ia]
        pair_W = list()        
        all_freq = list()
        w=0
        p=0
        while p < np:
            mx=min(self.max_channels,len(self.xspec[p][self.XCOL]))
            pair_W.append(zeros(mx))
            f=0
            while f < mx:
                fff=floor(self.xspec[p][0][f]*100)
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
            mx=min(len(self.xspec[p][0]),self.max_channels)
            f=0
            while f < mx:
                fff=floor(self.xspec[p][0][f]*100)
                fidx= all_freq.index(fff)
                i=self.covar_ridx[p][f]
                fC[fidx]+=pair_W[p][f]
                fE[fidx]+=pair_W[p][f]*pair_W[p][f]*self.covarDiag[ia][i]
                ff[fidx]+=pair_W[p][f]*self.xspec[p][self.XCOL][f]/self.Pgn[p][f]/zeta[p]
                f+=1
            p+=1
    
        ff/=fC
        fE /= fC*fC
        fE = sqrt(fE)
    
        fle=open("%s/ff.plot"%outpath,"w")
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
            mx=min(len(self.xspec[p][0]),self.max_channels)
            a=0
            if self.angles[p] > pi*49.0/180.0:
                a=1
            if self.angles[p] > pi*122.0/180.0:
                a=2
    
            pp[p]=sum(pair_W[p]*self.xspec[p][self.XCOL][0:mx]/self.Pgn[p][0:mx]/zeta[p])
            p_A[a]+=sum(pair_W[p]*self.xspec[p][self.XCOL][0:mx]/self.Pgn[p][0:mx])
            p_Aa[a]+=sum(pair_W[p])*self.angles[p]
            f=0
            while f < mx:
                i=self.covar_ridx[p][f]
                pE[p]+=pair_W[p][f]*pair_W[p][f]*self.covarDiag[ia][i]
                p_AE[a]+=pair_W[p][f]*pair_W[p][f]*self.covarDiag[ia][i]*zeta[p]*zeta[p]
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
    
        f=open("%s/hd.plot"%outpath,"w")
        p=0
        while p < np:
            f.write("%f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%s %s\n"%(180.0*self.angles[p]/pi,pp[p]*zeta[p],abs(pE[p]*zeta[p]),pp[p],pE[p],self.A2_guess*zeta[p],self.A2*zeta[p],pC[p],pairs[p][0],pairs[p][1]))
            p+=1
    
        f.close()
 
        f=open("%s/avhd.plot"%outpath,"w")
        p=0
        while p < 3:
            f.write("%f\t%g\t%g\n"%(180.0*p_Aa[p]/pi,p_A[p],p_AE[p]))
            p+=1
    
        f.close()




    
    def closeoff(self,p1,p2,freq,cross_r,cross_i,power_1,power_2,power_spectra,xspec,Pgn,Ts):
        # Close off the previous file.
        freq=array(freq)
        cross_r=array(cross_r)
        cross_i=array(cross_i)
        power_1=array(power_1)
        power_2=array(power_2)
        xspec.append([freq,cross_r,cross_i,power_1,power_2])
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
        Pgn.append(self.P_gw_nrm(freq))




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


        #return self.GW_K*power(1.0 + power(freq/self.GW_FC,2),self.ALPHA/2.0)/(12*pi*pi)


def get_stats(detector):
    Na=len(detector.little_as)
    chisq=zeros(Na)
    chisq2 = zeros(Na)
    A2=zeros(Na)
    varA2=zeros(Na)
    white_rms=zeros(Na)
    
    
    zeroA2 = -1
    oneA2 = -1
    halfA2 = -1
    ai=0
    for a in detector.little_as:
        tA2,tvarA2,tchisq,tchisq2 = detector.getA2(ai)
        #white_res = detector.get_whitened_residuals(ai)
        chisq[ai] = tchisq
        chisq2[ai] = tchisq2
        A2[ai] = tA2
        varA2[ai] = tvarA2
        #white_rms[ai]=sqrt(var(white_res))
        if a == 0.0:
            zeroA2 = tA2
        if a == 1.0:
            oneA2 = tA2
        if a == 0.5:
            halfA2 = tA2

        ai+=1
    best_chi_a = detector.little_as[argmin(chisq)]
    best_chi = amin(chisq)
    best_chi_A2 = A2[argmin(chisq)]
    best_chi2_a = detector.little_as[argmin(chisq2)]
    
    adaptive1_a = zeroA2 / detector.A2_guess
    adaptive1_ai = (abs(detector.little_as-adaptive1_a)).argmin()
    adaptive1_A2 = A2[adaptive1_ai]
    
    
    ret= [zeroA2, halfA2, oneA2, adaptive1_A2, adaptive1_ai, best_chi_a, best_chi,best_chi2_a, best_chi_A2]
    ret.extend(chisq)
    ret.extend(chisq2)
    return ret

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Detect the GWB.')
    parser.add_argument('-A','--gwamp',type=float,default=1e-15)
    parser.add_argument('-R','--redfile',default=None)
    parser.add_argument('-f','--filename',default="GW.sum")
    parser.add_argument('-w','--white_factor',type=float,default=1.0)
    parser.add_argument('-n','--maxchans',type=int,default=10)
    parser.add_argument('-Z','--zerodir',default=None)
    parser.add_argument('-k','--skipfile',default=None)
    parser.add_argument('-d','--outdir',default=".")
    parser.add_argument('-Q','--quick',default=False, action='store_true')
    parser.add_argument('-P','--plotfiles',default=False, action='store_true')
    parser.add_argument('--test',default=False, action='store_true')
    parser.add_argument('-a','--write_avspec',default=False, action='store_true')
    parser.add_argument('-c','--cfile',default=None)
    parser.add_argument('-C','--a2fix',default=None)
    parser.add_argument('-W','--whitefile',default=None)
    parser.add_argument('-y','--yrskip',type=float,default=-1)
    parser.add_argument('--start_chan',type=int,default=0)

    parser.add_argument('dirs',nargs="+")

    args=parser.parse_args()

    try:
        os.mkdir(args.outdir)
    except OSError as err:
        pass
    try:
        test=open("%s/onstats"%args.outdir,"w")
        test.close()
        test=open("%s/offstats"%args.outdir,"w")
        test.close()
    except OSError as err:
        print err
        exit(1)



    detector = GravWavBkgrdDetect(args.gwamp)
    detector.yrskip = args.yrskip
    detector.start_channel=args.start_chan
    if args.whitefile != None:
        print "WF"
        wf=dict()
        with open(args.whitefile) as f:
            for line in f:
                elems=line.split()
                wf[elems[0]]=float(elems[1])
        detector.wf=wf
        
        
    
    if args.cfile != None:
        print "CORN"
        iii=1
        if args.whitefile!=None:
            iii=2
        corn=list()
        with open(args.cfile) as f:
            for line in f:
                elems=line.split()
                corn.append(float(elems[iii])/float(elems[0]))
        detector.corn=array(corn)
        print "mean corn=",mean(detector.corn),detector.corn[0],detector.corn[10]


    a2varfix=None
    if args.a2fix != None:
        print "FIX A2"
        a2fix=list()
        a2varfix=list()
        with open(args.a2fix) as f:
            for line in f:
                elems=line.split()
                a2v = float(elems[0])
                a2e = float(elems[2])/sqrt(1000.0)
                if abs(a2v-detector.A2_guess)/a2e<2:
                    a2fix.append(1.0)
                    a2varfix.append(1.0)
                    continue
                factor = detector.A2_guess/a2v
                if 1.0/factor < 0.1:
                    a2fix.append(0)
                    a2varfix.append(100.0)
                    continue
                a2fix.append(factor)
                rv=float(elems[7])
                rm=float(elems[8])
                rfactor = rm/rv
                if rfactor > factor:
                    a2varfix.append(rfactor/factor)
                else:
                    a2varfix.append(1.0)

            a2fix=array(a2fix)
            a2varfix=array(a2varfix)

            ff=open("tt","w")
            for a,b in zip(a2fix,a2varfix):
                ff.write("%g %g\n"%(a,b))
            ff.close()
    detector.QUICK=args.quick
    if detector.QUICK:
        print "QUICK MODE - NO chisq computaiton possible"

    if args.skipfile != None:
        f = open (args.skipfile)
        for line in f:
            e=line.split()
            detector.skipem.append(e[0])

    detector.max_channels = args.maxchans
    if args.test:
        TESTY=list()
        print "*** TEST MODE ***"
    
    if args.write_avspec:
        onaddr=SpecAddr("on")
        offaddr=SpecAddr("off")
    ii=0
    onstats=None
    offstats=None
    #detector.little_as=[0.0]
    for d in args.dirs:
        infile=d+"/"+args.filename
        print "Processing '%s'                [%d]\r"%(infile,ii),
        stdout.flush()
        detector.read_xspec(infile)

        if detector.cinvs==None:
            print ""
            npsr=1
            np=0
            while np < detector.np:
                npsr+=1
                np=npsr*(npsr-1)/2
            if np != detector.np:
                print "Error: number of pairs is not valid"
                print "npsr=%d npair=%d expected=%d"%(npsr,np,detector.np)
                exit(1)
            else:
                print "Npsr = %d, Npair = %d"%(npsr,np)

            Ncov = np*detector.max_channels
            Cmem1 = Ncov*Ncov*8.0/pow(2,20) # 4 byte floats
            print "Mem per matrix %.1f MB"%Cmem1
            Nmatrix = len(detector.little_as)+6
            print "Tot Mem = %.1f MB"%(Cmem1*Nmatrix)

            detector.model_noise(args.redfile, args.white_factor)
            print "Compute covar"
            detector.compute_covaraince_matrix()
            
           
            print "Compute cinvs"
            detector.form_Cinvs(a2varfix=a2varfix)
            
            #print "Compute Uinvs"
            #detector.form_Uinvs()
            detector.clear_covar()
            
            f = open("%s/theory.var"%args.outdir,"w")
            ai=0
            for a in detector.little_as:
                print "%0.2f  %.3g"%(a,detector.theory_var[ai])
                f.write("%0.2f  %.3g\n"%(a,detector.theory_var[ai]))
                ai+=1
            f.close()


        ##### BEGIN TEST
        if args.test:
            print detector.theory_var[0], detector.theory_var[10] ,"TT"
            print detector.A2_guess/detector.theory_var[0], detector.A2_guess/detector.theory_var[10] ,"SN"
            exit(1)
       
        
        
        ##### END TEST
        detector.get_A2ijks()
        if args.a2fix != None:
            detector.A2ijk*=a2fix
        if args.write_avspec:
            detector.add_spec(onaddr)

        ss = get_stats(detector)
        if args.plotfiles:
            ai=ss[4]
            detector.write_plotfiles(d,ai)

        if onstats==None:
            onstats = zeros((len(ss),len(args.dirs)))
        si=0
        for s in ss:
            onstats[si][ii] = s
            si+=1
                
        if args.zerodir != None:
            infile="%s/%s/%s"%(args.zerodir,d,args.filename)
            print "Processing '%s/%s' [ZERO]            \r"%(d,args.filename),
            stdout.flush()
            detector.read_xspec(infile)
            detector.get_A2ijks()
            if args.a2fix != None:
                detector.A2ijk*=a2fix

            if args.write_avspec:
                detector.add_spec(offaddr)
            ss = get_stats(detector)

            if offstats==None:
                offstats = zeros((len(ss),len(args.dirs)))
            si=0
            for s in ss:
                offstats[si][ii] = s
                si+=1
        ii+=1
            
    if onstats != None:
        f = open("%s/onstats"%args.outdir,"w")
        for itr in transpose(onstats):
            for s in itr:
                f.write("%g "%s)
            f.write("\n")

        f.close()
        
    
    if offstats != None:
        f = open("%s/offstats"%args.outdir,"w")
        for itr in transpose(offstats):
            for s in itr:
                f.write("%g "%s)
            f.write("\n")

        f.close()
        
    if args.write_avspec:
        detector.write_avspec(onaddr,args.outdir)
        detector.write_avspec(offaddr,args.outdir)
        
    print "0", detector.little_as[0], mean(offstats[0]),sqrt(var(offstats[0])), detector.theory_var[0], sqrt(var(offstats[0]))/detector.theory_var[0]
    
    print "1", detector.little_as[0], mean(onstats[0]),sqrt(var(onstats[0])), detector.theory_var[0], sqrt(var(onstats[0]))/detector.theory_var[0]

    print "1", detector.little_as[5],mean(onstats[1]),sqrt(var(onstats[1])), detector.theory_var[5], sqrt(var(onstats[1]))/detector.theory_var[5]
    
    print "1", detector.little_as[10],mean(onstats[2]),sqrt(var(onstats[2])), detector.theory_var[10], sqrt(var(onstats[2]))/detector.theory_var[10]


    print mean(onstats[4]), mean(offstats[4])
