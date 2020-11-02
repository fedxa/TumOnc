#-*- Mode: python -*-
import numpy as np
import math as ma
cimport numpy as np
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

from libc.math cimport exp

cdef class critstat:
    cdef int C
    cdef np.ndarray f
    cdef np.ndarray p
    cdef np.ndarray lp
    cdef np.ndarray l1p
    cdef DTYPE_t suml1p
    cdef np.ndarray l1p_psum
    def __init__(self, p, f, pprime=0.):
        """WARNING: p array should be sorted in descending order!"""
        self.C = len(f)
        self.f = np.array(f)
        self.p = np.array(p)*self.f
        self.lp = np.log(self.p)
        self.l1p = np.log(1-self.p)
        self.suml1p = np.sum(self.l1p)
        self.l1p_psum = np.zeros((self.C+1,self.C+1), dtype=DTYPE)
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
        
    cdef DTYPE_t lp2dim(self, int d1, int d2):
        """deg=[d1,d2] with d1>d2>...
        of the _least_ probable hits"""
        if d1<0 and d2<0:
            return self.suml1p
        elif d2<0:
            return self.suml1p-self.l1p[d1]+self.lp[d1]
        elif d1>d2:
            return ( self.lp[d2]+self.l1p[d2+1:d1].sum()
                     +self.lp[d1]+self.l1p[d1+1:].sum() )
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

    cdef DTYPE_t lp3dim(self,int d1,int d2,int d3):
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

    cdef DTYPE_t lp4dim(self,int d1,int d2,int d3, int d4):
        """deg=[d1,d2,d3] with d1>d2>...
        of the _least_ probable hits"""
        if d1<0 and d2<0 and d3<0 and d4<0:
            return self.suml1p
        elif d2<0 and d3<0 and d4<0:
            return self.suml1p-self.l1p[d1]+self.lp[d1]
        elif d3<0 and d4<0:
            return self.suml1p-self.l1p[d1]+self.lp[d1]-self.l1p[d2]+self.lp[d2]
        elif d4<0:
            return (self.suml1p-self.l1p[d1]+self.lp[d1]
                    -self.l1p[d2]+self.lp[d2]
                    -self.l1p[d3]+self.lp[d3])
        elif d1>d2 and d2>d3 and d3>d4:
            return ( self.lp[d4]+self.l1p_psum[d4+1,d3]
                     +self.lp[d3]+self.l1p_psum[d3+1,d2]
                     +self.lp[d2]+self.l1p_psum[d2+1,d1]
                     +self.lp[d1]+self.l1p_psum[d1+1,self.C] )
        else:
            print("lp4dim: Impossible argument")
            return -np.inf

    cdef DTYPE_t lgamma2dim(self, int d1, int d2):
        """log(gamma) which is small for probable event under current hypothesys"""
#        return self.lpprime2dim(d1,d2)-self.lp2dim(d1,d2)
        return -self.lp2dim(d1,d2)

    cdef DTYPE_t lgamma3dim(self,int d1,int d2,int d3):
#        return self.lpprime3dim(d1,d2,d3)-self.lp3dim(d1,d2,d3)
        return -self.lp3dim(d1,d2,d3)

    cdef DTYPE_t lgamma4dim(self,int d1,int d2,int d3, int d4):
#        return self.lpprime3dim(d1,d2,d3)-self.lp3dim(d1,d2,d3)
        return -self.lp4dim(d1, d2, d3, d4)

    def plgamma2dim(self,lgamma):
        """Probability that hypothesis H0 is true"""
        p = 0.
        for i in range(-1,self.C):
            for j in range(i,self.C):
                if i<0 or j>i:
                    if self.lgamma2dim(j,i)<lgamma:
                        p += exp(self.lp2dim(j,i))
        return p

    cdef DTYPE_t negplgamma2dim(self, DTYPE_t lgamma):
        """1-probability of true H0 -- p-value"""
        cdef DTYPE_t p = 0.
        cdef int i, j
        for i in range(-1,self.C):
            for j in range(i,self.C):
                if i<0 or j>i:
                    if self.lgamma2dim(j,i)>=lgamma:
                        p += exp(self.lp2dim(j,i))
        return p

    def plgamma3dim(self,lgamma):
        p = 0.
        for i in range(-1,self.C):
            for j in range(i,self.C):
                for k in range(j,self.C):
                    if j<0 or (i<0 and k>j) or (k>j and j>i):
                        if self.lgamma3dim(k,j,i)<lgamma:
                            p += exp(self.lp3dim(k,j,i))
        return p

    cdef DTYPE_t negplgamma3dim(self,DTYPE_t lgamma):
        cdef DTYPE_t p = 0.
        cdef int i, j, k
        for i in range(-1,self.C):
            for j in range(i,self.C):
                for k in range(j,self.C):
                    if j<0 or (i<0 and k>j) or (k>j and j>i):
                        if self.lgamma3dim(k,j,i)>=lgamma:
                            p += exp(self.lp3dim(k,j,i))
        return p

    cdef DTYPE_t negplgamma4dim(self,DTYPE_t lgamma):
        cdef DTYPE_t p = 0.
        cdef int i, j, k, l
        for i in range(-1,self.C):
            for j in range(i,self.C):
                for k in range(j,self.C):
                    for l in range(k,self.C):
                        if (k<0
                            or (j<0 and l>k)
                            or (i<0 and l>k and k>j)
                            or (l>k and k>j and j>i)):
                            if self.lgamma4dim(l, k, j, i)>=lgamma:
                                p += exp(self.lp4dim(l, k, j, i))
        return p

    cpdef DTYPE_t pvalue2dim(self, int d1, int d2):
        return self.negplgamma2dim(self.lgamma2dim(d1,d2))

    cpdef DTYPE_t pvalue3dim(self, int d1, int d2, int d3):
        return self.negplgamma3dim(self.lgamma3dim(d1,d2,d3))

    cpdef DTYPE_t pvalue4dim(self, int d1, int d2, int d3, int d4):
        return self.negplgamma4dim(self.lgamma4dim(d1, d2, d3, d4))

    def pvalue3(self,l):
        cdef np.ndarray sl = np.sort(l)[::-1]
        assert sl.dtype == np.int64
        if len(sl)>=3:
            return self.pvalue3dim(sl[0],sl[1],sl[2])
        elif len(sl)==2:
            return self.pvalue3dim(sl[0],sl[1],-1)
        elif len(sl)==1:
            return self.pvalue3dim(sl[0],-1,-1)
        else:
            return self.pvalue3dim(-1,-1,-1)

    def pvalue4(self,l):
        cdef np.ndarray sl = np.sort(l)[::-1]
        assert sl.dtype == np.int64
        if len(sl)>=4:
            return self.pvalue4dim(sl[0],sl[1],sl[2],sl[3])
        elif len(sl)==3:
            return self.pvalue4dim(sl[0],sl[1],sl[2],-1)
        elif len(sl)==2:
            return self.pvalue4dim(sl[0],sl[1],-1,-1)
        elif len(sl)==1:
            return self.pvalue4dim(sl[0],-1,-1,-1)
        else:
            return self.pvalue4dim(-1,-1,-1,-1)

    def pvalue2(self,l):
        sl=np.sort(l)[::-1]
        assert sl.dtype == np.int64
        if len(sl)>=2:
            return self.pvalue2dim(sl[0],sl[1])
        elif len(sl)==1:
            return self.pvalue2dim(sl[0],-1)
        else:
            return self.pvalue2dim(-1,-1)

    def guessp(self, l):
        ppind = np.ones_like(self.f, dtype=bool)
        for i in l:
            ppind[i] = False
        sf = np.sum(self.f[ppind])
        sf2 = np.sum(self.f[ppind]*self.f[ppind])
        return (ma.sqrt(sf*sf+4*len(l)*sf2)-sf)/(2*sf2)
