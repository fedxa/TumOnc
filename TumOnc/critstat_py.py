from __future__ import print_function

import numpy as np
import math as ma

class critstat:
    def __init__(self, p, pprime=0.):
        self.C = len(p)
        self.p = p
        self.lp = np.log(p)
        self.l1p = np.log(1-p)
        self.suml1p = np.sum(self.l1p)
        # Calculate all the partial sums preventively for speed
        self.l1p_psum = np.zeros((self.C+1, self.C+1))
        for i in range(self.C+1):
            for j in range(i,self.C+1):
                self.l1p_psum[i,j] = self.l1p[i:j].sum()
        # if pprime==0.:
        #     self.pprime = np.mean(p)
        # else:
        #     self.pprime = pprime
        # self.lpprime = np.log(self.pprime)
        # self.l1pprime = np.log(1-self.pprime)

    # def lpprime2dim(self,d1,d2):
    #     if d1<0 and d2<0:
    #         return self.l1pprime*self.C
    #     elif d2<0:
    #         return self.lpprime+self.l1pprime*(self.C-1)
    #     elif d1>d2:
    #         return self.lpprime*2+self.l1pprime*(self.C-d2-2)
    #     else:
    #         print("lpprime2dim: Impossible argument")
    #         return -np.inf
        
    def lp2dim(self,d1,d2):
        """deg=[d1,d2] with d1>d2>...
        of the _least_ probable hits"""
        if d1<0 and d2<0:
            return self.suml1p
        elif d2<0:
            return self.suml1p-self.l1p[d1]+self.lp[d1]
        elif d1>d2:
            return ( self.lp[d2]+self.l1p_psum[d2+1,d1]
                     +self.lp[d1]+self.l1p_psum[d1+1,self.C] )
        else:
            print("lp2dim: Impossible argument")
            return -np.inf

    # def lpprime3dim(self,d1,d2,d3):
    #     if d1<0 and d2<0 and d3<0:
    #         return self.l1pprime*self.C
    #     elif d2<0 and d3<0:
    #         return self.lpprime+self.l1pprime*(self.C-1)
    #     elif d3<0:
    #         return self.lpprime*2+self.l1pprime*(self.C-2)
    #     elif d1>d2 and d2>d3:
    #         return self.lpprime*3+self.l1pprime*(self.C-d3-3)
    #     else:
    #         print("lpprime3dim: Impossible argument")
    #         return -np.inf

    def lp3dim(self,d1,d2,d3):
        """deg=[d1,d2,d3] with d1>d2>...
        of the _least_ probable hits"""
        if d1<0 and d2<0 and d3<0:
            return self.suml1p
        elif d2<0 and d3<0:
            return self.suml1p-self.l1p[d1]+self.lp[d1]
        elif d3<0:
            return self.suml1p-self.l1p[d1]+self.lp[d1]-self.l1p[d2]+self.lp[d2]
        elif d1>d2 and d2>d3:
            return ( self.lp[d3]+self.l1p_psum[d3+1,d2]
                     +self.lp[d2]+self.l1p_psum[d2+1,d1]
                     +self.lp[d1]+self.l1p_psum[d1+1,self.C] )
        else:
            print("lp3dim: Impossible argument")
            return -np.inf

    def lgamma2dim(self,d1,d2):
        """log(gamma) which is small for probable event under current hypothesys"""
#        return self.lpprime2dim(d1,d2)-self.lp2dim(d1,d2)
        return -self.lp2dim(d1,d2)

    def lgamma3dim(self,d1,d2,d3):
#        return self.lpprime3dim(d1,d2,d3)-self.lp3dim(d1,d2,d3)
        return -self.lp3dim(d1,d2,d3)

    def plgamma2dim(self,lgamma):
        """Probability that hypothesis H0 is true"""
        p = 0.
        for i in range(-1,self.C):
            for j in range(i,self.C):
                if i<0 or j>i:
                    if self.lgamma2dim(j,i)<lgamma:
                        p += ma.exp(self.lp2dim(j,i))
        return p

    def negplgamma2dim(self,lgamma):
        """1-probability of true H0 -- p-value"""
        p = 0.
        for i in range(-1,self.C):
            for j in range(i,self.C):
                if i<0 or j>i:
                    if self.lgamma2dim(j,i)>=lgamma:
                        p += ma.exp(self.lp2dim(j,i))
        return p

    def plgamma3dim(self,lgamma):
        p = 0.
        for i in range(-1,self.C):
            for j in range(i,self.C):
                for k in range(j,self.C):
                    if j<0 or (i<0 and k>j) or (k>j and j>i):
                        if self.lgamma3dim(k,j,i)<lgamma:
                            p += ma.exp(self.lp3dim(k,j,i))
        return p

    def negplgamma3dim(self,lgamma):
        p = 0.
        for i in range(-1,self.C):
            for j in range(i,self.C):
                for k in range(j,self.C):
                    if j<0 or (i<0 and k>j) or (k>j and j>i):
                        if self.lgamma3dim(k,j,i)>=lgamma:
                            p += ma.exp(self.lp3dim(k,j,i))
        return p

    def pvalue2dim(self,d1,d2):
        return self.negplgamma2dim(self.lgamma2dim(d1,d2))

    def pvalue3dim(self,d1,d2,d3):
        return self.negplgamma3dim(self.lgamma3dim(d1,d2,d3))

    def pvalue3(self,l):
        sl=np.sort(l)[::-1]
        if len(sl)>=3:
            return self.pvalue3dim(sl[0],sl[1],sl[2])
        elif len(sl)==2:
            return self.pvalue3dim(sl[0],sl[1],-1)
        elif len(sl)==1:
            return self.pvalue3dim(sl[0],-1,-1)
        else:
            return self.pvalue3dim(-1,-1,-1)

    def pvalue2(self,l):
        sl=np.sort(l)[::-1]
        if len(sl)>=2:
            return self.pvalue2dim(sl[0],sl[1])
        elif len(sl)==1:
            return self.pvalue2dim(sl[0],-1)
        else:
            return self.pvalue2dim(-1,-1)
