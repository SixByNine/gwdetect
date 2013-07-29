#!/usr/bin/python
import argparse
from sys import stdout
from math import pi, sqrt, floor
from numpy import zeros, array, power,dot, linspace, log,cos,sin,diag,diagonal,eye
#from numpy import linalg as lalg
from scipy import linalg as lalg
from matplotlib import pyplot as plt

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
        self.path="."
        self.cinvs=None
        self.start_channel=0
        

    
    def getA2(self,ai):
        M=zeros(self.N)+1.0
        CM = dot(self.cinvs[ai],M)
        MCM = dot(M,CM)
        self.Wijk = (1.0/MCM) * CM
        
        
        self.A2 = dot(self.Wijk,self.A2ijk)
        self.varA2 = sum(self.Wijk*self.Wijk*power(self.A2ijk,2))/pow(sum(self.Wijk),2)
        self.R = self.A2ijk - self.A2
        self.chisq = dot(dot(self.R,self.cinvs[ai]),self.R) / (self.N-1)
        return self.A2,self.varA2, self.chisq
    
    def form_Cinvs(self):
        self.cinvs=list()
        self.Ww = list()
        for a in self.little_as:
            covar=a*a*self.covar_GW4 + (1-a*a)*self.covar_GN + self.covar_PN + a*self.covar_GW2
            Cinv = lalg.solve(covar,eye(self.N),sym_pos=True,check_finite=False,lower=True)
            ### TODO: Check that this is the same as inverting after filling in the upper triangle!!! TODO
            #self.cinvs.append(lalg.inv(covar))
            self.cinvs.append(Cinv)
            M=zeros(self.N)+1
            
            eq2=lalg.solve(covar,M,sym_pos=True,check_finite=False,lower=True)
            eq1=1.0/dot(eq2,M)
            self.Ww.append(eq1*eq2)
    
    def get_Aijks(self):
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
            m1 = pwr_model()
            m1.white=self.white[p][0]*white_factor
            if pair[0] in red:
                m1.red=red[pair[0]]
            else:
                m1.red=[]
    
            m2 = pwr_model()
            m2.white=self.white[p][1]*white_factor
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
        
        
        covar_GW2 = zeros((N, N))
        covar_GW4 = zeros((N, N))
        covar_PN = zeros((N, N))
        covar_GN = zeros((N, N))
        print "Compute GW covariances"
        GWcovar = zeros((np, np))
        noGWvar = zeros(np)
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
                C = 0.5 * pow(self.A_guess, 4) * (zeta_il * zeta_jm + zeta_im * zeta_jl) / (zeta_ij * zeta_lm) #print "%s %s %s %s %.4f %.4f %.4f %.4f %.4f %.4f %g"%(psr1,psr2,psr3,psr4,zeta_ij,zeta_lm,zeta_il,zeta_jm,zeta_im,zeta_jl,C)
                GWcovar[ij][lm] = C
                GWcovar[lm][ij] = C
                if ij == lm:
                    noGWvar[ij] = 0.5 * pow(self.A_guess, 4) / (zeta_ij * zeta_lm)
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
                    if psr2 == psr3:
                        C += 0.5 * (P2 * self.A2_guess) * zeta_im / zeta_ij / zeta_lm / self.P_gw_nrm(fa)
                    if psr1 == psr4:
                        C += 0.5 * (P1 * self.A2_guess) * zeta_jl / zeta_ij / zeta_lm / self.P_gw_nrm(fa)
                    if psr1 == psr3:
                        C += 0.5 * (P1 * self.A2_guess) * zeta_jm / zeta_ij / zeta_lm / self.P_gw_nrm(fa)
                    if psr2 == psr4:
                        C += 0.5 * (P2 * self.A2_guess) * zeta_il / zeta_ij / zeta_lm / self.P_gw_nrm(fa)


                    if ij == lm:
                        C += 0.5 * (P1 * P2) / pow(zeta[ij] * self.P_gw_nrm(fa), 2)
                        covar_PN[a][b] += X * C
                    else:
                        covar_GW2[a][b] += X * C
                C = GWcovar[ij][lm]
                covar_GW4[a][b] += X * C
                if ij == lm:
                    C = noGWvar[ij]
                    covar_GN[a][b] += X * C
                b += 1
            
            a += 1
        
        print "Done"
        self.covar_GN = covar_GN
        self.covar_GW2 = covar_GW2
        self.covar_GW4 = covar_GW4
        self.covar_PN = covar_PN
    
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
        inf = open(infile)
        for line in inf:
            elems = line.split()
            if line.startswith("#"):
                skip = False
                if ii > 0:
                    self.closeoff(p1, p2, freq, cross_r, cross_i, power_1, power_2, power_spectra, xspec, Pgn, Ts)
                elif p1:
                    print "Skip", p1, p2, "No Freq"
                # now set up for the new pulsar.
                p1 = elems[1]
                p2 = elems[2]
                if p1 in self.skipem or p2 in self.skipem:
                    print "Skip", p1, p2
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
            if fff > self.fmax:
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

    def write_plotfiles(self, covar, A2, a, W):
        covar=a*a*self.covar_GW4 + pow(1-a,2)*self.covar_GN + self.covar_PN + a*self.covar_GW2

        np=self.np
        zeta=self.zeta
        pairs=self.pairs
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
                fE[fidx]+=pair_W[p][f]*pair_W[p][f]*covar[i][i]
                ff[fidx]+=pair_W[p][f]*self.xspec[p][self.XCOL][f]/self.Pgn[p][f]/zeta[p]
                f+=1
            p+=1
    
        ff/=fC
        fE /= fC*fC
        fE = sqrt(fE)
    
        fle=open("%s/ff.plot"%self.path,"w")
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
    
        f=open("%s/hd.plot"%self.path,"w")
        p=0
        while p < np:
            f.write("%f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%s %s\n"%(180.0*self.angles[p]/pi,pp[p]*zeta[p],abs(pE[p]*zeta[p]),pp[p],pE[p],self.A2_guess*zeta[p],A2*zeta[p],pC[p],pairs[p][0],pairs[p][1]))
            p+=1
    
        f.close()



    
    def closeoff(self,p1,p2,freq,cross_r,cross_i,power_1,power_2,power_spectra,xspec,Pgn,Ts):
        # Close off the previous file.
        freq=array(freq)
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
        Pgn.append(self.P_gw_nrm(freq))




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




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Detect the GWB.')
    parser.add_argument('-A','--gwamp',type=float,default=1e-15)
    parser.add_argument('-R','--redfile',default=None)
    parser.add_argument('-w','--white_factor',default=1.0)
    parser.add_argument('-n','--maxchans',type=int,default=10)

    parser.add_argument('dirs',nargs="+")

    args=parser.parse_args()

    detector = GravWavBkgrdDetect(args.gwamp)
    detector.max_channels = args.maxchans
    #detector.little_as=[0.0]
    for d in args.dirs:
        infile=d+"/GW.sum"
        print "Processing '%s'"%(infile)
        detector.read_xspec(infile)

        if detector.cinvs==None:
            detector.model_noise(args.redfile, args.white_factor)
            print "Compute covar"
            detector.compute_covaraince_matrix()
            print "Compute cinvs"
            detector.form_Cinvs()


        detector.get_Aijks()
        ai=0
        for a in detector.little_as:
            A2,varA2,chisq = detector.getA2(ai)
            
            print a, A2, varA2, chisq
            ai+=1









