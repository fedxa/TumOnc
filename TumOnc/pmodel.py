from __future__ import print_function
from __future__ import division

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import math as ma
import pandas as pd
import scipy.stats as st
import scipy.misc as sm
import scipy.special as ss
# import scipy.misc as sm


class pmodel:
    def fit(self):
        print("Use proper subclass!")
        return 0

    def searchab_step(self, a):
        fi = self.dlnP3(a)
        fij = self.d2lnP3(a)
        return -np.linalg.solve(fij, fi)

    def searchab(self, a0, err=1e-6):
        a = a0
        da = self.searchab_step(a)
        a += da
        while np.linalg.norm(da)/np.linalg.norm(a) > err*err:
            print(a, da)
            da = self.searchab_step(a)
            a += da
        return a

    # The same, but using scipy.optimize
    def searchab_opt(self, a0, err=1e-6):
        res = opt.minimize(lambda a: -self.lnP3(a),
                           a0, method='Newton-CG',
                           jac=lambda a: -self.dlnP3(a),
                           hess=lambda a: -self.d2lnP3(a),
                           options={'xtol': err, 'disp': True})
        return res.x

     
class pmodel2(pmodel):
    def fit(self, ns, N, V):
        self.ns = ns
        self.N = N
        self.V = V
        self.a = [0, 0, 0]
        self.f3 = self.f3f()
        p3 = np.outer(self.f3, np.ones_like(self.V[0]))
        self.pmatrix = p3

    def f3f(self):
        C = np.size(self.ns, 3)
        return (np.sum(self.ns, axis=(0, 2, 3))
                / np.sum(self.N, axis=1)
                / C)


class pmodel3(pmodel):
    def fit(self, ns, N, V):
        self.ns = ns
        self.N = N
        self.V = V
        a0 = np.zeros(len(V))
        self.a = self.searchab(a0)
        self.f3 = self.f3f(self.a)
        self.expcv = np.exp(sum(self.a[i]*V[i] for i in range(len(self.a))))
        self.pmatrix = np.outer(self.f3, self.expcv)

    def f3f(self, a):
        C = np.size(self.ns, 3)
        expcv = np.exp(sum(a[i]*self.V[i] for i in range(len(a))))
        return (np.sum(self.ns, axis=(0, 2, 3))
                / np.sum(self.N*expcv, axis=1)
                / C)

    def df3di(self, i, a):
        C = np.size(self.ns, 3)
        expcv = np.exp(sum(a[i]*self.V[i] for i in range(len(a))))
        return -(np.sum(self.ns, axis=(0, 2, 3))
                 / np.sum(self.N*expcv, axis=1)**2
                 / C)*np.sum(self.N*self.V[i]*expcv, axis=1)

    def lnP3(self, a):
        C = np.size(self.ns, 3)
        expcv = np.exp(sum(a[i]*self.V[i] for i in range(len(a))))
        f3ab = self.f3f(a)
        return (np.sum(np.sum(self.ns, axis=(0, 2, 3))*np.log(f3ab))
                + np.sum(np.sum(self.ns, axis=(0, 1, 3))
                         * sum(a[i]*self.V[i] for i in range(len(a))))
                - np.sum(C*f3ab*np.sum(self.N*expcv, axis=1)))

    def dlnP3(self,a):
        d = np.empty_like(a)
        C = np.size(self.ns,3)
        expcv = np.exp( sum(a[i]*self.V[i] for i in range(len(a))) )
        f3ab = self.f3f(a)
        for i in range(len(a)):
            d[i] = ( np.sum( np.sum(self.ns,axis=(0,3))*self.V[i] )
                     - np.sum( C*f3ab*np.sum( self.N*self.V[i]*expcv,axis=1) ) )
        return d
    def d2lnP3(self,a):
        d = np.empty([len(a),len(a)])
        C = np.size(self.ns,3)
        expcv = np.exp( sum(a[i]*self.V[i] for i in range(len(a))) )
        f3ab = self.f3f(a)
        for i in range(len(a)):
            for j in range(len(a)):
                d[i,j] = ( - np.sum( C*f3ab*np.sum( self.N*self.V[i]*self.V[j]*expcv,axis=1) )
                           - np.sum( C*self.df3di(j,a)*
                                     np.sum( self.N*self.V[i]*expcv,axis=1 ) ) )
        return d


class pmodel3prime(pmodel3):
    def fit(self, ns, N, V):
        pmodel3.__init__(self, ns, N, V)
        expcv = np.exp( sum(self.a[i]*V[i] for i in range(len(self.a))) )
        p3 = np.outer( self.f3, expcv )
        self.pmatrix = p3


class pmodel4(pmodel):
    def __init__(self, avg_coeff=10):
        self.avg_coeff = avg_coeff

    def fit(self, ns, N, V):
        self.ns = ns
        self.N = N
        self.V = V
        self.C = np.size(ns, 3)
        pguess = ns.sum(axis=(0, 2, 3))/N.sum(axis=1)/self.C
        Navg = N.mean(axis=1)*self.C*self.avg_coeff
        nguess = pguess*Navg
        # p3 = np.zeros_like(N)
        # for c in range(np.size(N, 0)):
        #     for g in range(np.size(N, 1)):
        #         p3[c, g] = ((ns.sum(axis=(0,3))[c, g]+nguess[c])
        #                     / (N[c, g]+Navg[c]))
        p3 = ((ns.sum(axis=(0, 3)).transpose()+nguess) /
              (N.transpose()*self.C+Navg)).transpose()
        self.pmatrix = p3


class pmodel4bis(pmodel):
    def __init__(self, avg_coeff=10):
        self.avg_coeff = avg_coeff

    def fit(self, ns, N, V):
        self.ns = ns
        self.N = N
        self.V = V
        self.C = np.size(ns, 3)
        pguess = ns.sum(axis=(0, 1, 3))/N.sum(axis=0)/self.C
        Navg = N.mean(axis=0)*self.C*self.avg_coeff
        nguess = pguess*Navg
        # p3 = np.zeros_like(N)
        # for c in range(np.size(N, 0)):
        #     for g in range(np.size(N, 1)):
        #         p3[c, g] = ((ns.sum(axis=(0,3))[c, g]+nguess[c])
        #                     / (N[c, g]+Navg[c]))
        p3 = ((ns.sum(axis=(0, 3))+nguess) /
              (N*self.C+Navg))
        self.pmatrix = p3


def Hbetabin(n1, N1, n2, N2):
    a = n2+1
    b = N2-n2+1
    return sm.comb(N1, n1)*ss.beta(n1+a, N1-n1+b)/ss.beta(a, b)


class pmodel_bagels:
    def __init__(self, nBmax=50, Qmin=0.05):
        self.nBmax = nBmax
        self.Qmin = Qmin

    def fit(self, ns, N, V):
        self.ns = ns
        self.N = N
        self.V = V
        self.C = np.size(ns, 3)
        self.ppd = self.pairwise_distance(V)
        self.nsg = ns.sum(axis=(0, 1, 3))
        self.Ng = N.sum(axis=0)*self.C
        self.xg = np.empty_like(self.nsg)
        self.Xg = np.empty_like(self.Ng)
        self.Bsize = np.empty(len(self.Ng), dtype=int)
        for i in range(len(self.Ng)):
            print("computing gene %d" % (i), flush=True)
            self.compute_bagel(i)
        self.fc = (self.ns.sum(axis=(0, 2, 3)) /
                   self.ns.sum(axis=(0, 2, 3)).mean())
        self.fp = (self.ns.sum(axis=(0, 1, 2)) /
                   self.ns.sum(axis=(0, 1, 2)).mean())
        self.pmatrix = np.outer(self.fc, self.xg/self.Xg)
        del(self.ppd)

    def Qleft(self, i, g):
        return sum(Hbetabin(n, self.Ng[i], self.nsg[g], self.Ng[g])
                   for n in range(int(self.nsg[i])))

    def Qig(self, i, g):
        qleft = self.Qleft(i, g)
        return 2*min(qleft, 1-qleft)

    def get_bagel(self, i):
        gind = np.argsort(self.ppd[i, :])
        j = 1
        while j < self.nBmax and self.Qig(i, gind[j]) > self.Qmin:
            j += 1
        return gind[0:j]

    def compute_bagel(self, i):
        bg = self.get_bagel(i)
        self.xg[i] = np.sum(self.nsg[bg])
        self.Xg[i] = np.sum(self.Ng[bg])
        self.Bsize[i] = len(bg)

    # This one is the slowest!
    def pairwise_distance(self, V):
        lc = len(V)
        lg = len(V[0])
        ppd = np.empty((lg, lg))
        for i in range(lg):
            print("computing gene distance %d" % (i), flush=True)
            ppd[i, i] = 0
            for j in range(i+1, lg):
                ppd[i, j] = ma.sqrt(sum(ma.pow(V[c][i]-V[c][j], 2)
                                        for c in range(lc)))
                ppd[j, i] = ppd[i, j]
        return ppd


def plot_pmodel_test(pmodel, minnmut=100):
    print("1", sep='')
    bd = pd.DataFrame({'n': np.sum(pmodel.ns, axis=(0, 3)).flatten(),
                       'N': pmodel.N.flatten()*pmodel.C,
                       'p': pmodel.pmodel.pmatrix.flatten()})
    print("2", sep='')
    bd.sort(columns='p', inplace=True)
    print("3", sep='')
    group_stat2(bd, minnmut)
    print("4", sep='')
    bdg = bd.groupby('gid').agg({'N': np.sum,
                                 'n': np.sum,
                                 'p': np.mean})
    print("5", sep='')
    bdg['pp'] = [np.sum(f.p*f.N)/np.sum(f.N)
                 for n, f in bd.groupby('gid')]
    print("6")
    print(st.chisquare(bdg.n, bdg.pp*bdg.N))
    plt.ylim(ymin=1e-6, ymax=5e-3)
    plt.semilogy(bdg.n/bdg.N)
    plt.semilogy(bdg.pp)
    plt.semilogy(bdg.pp+np.sqrt(bdg.pp/bdg.N))
    plt.semilogy(bdg.pp-np.sqrt(bdg.pp/bdg.N))
    plt.show()
    return bdg


# TODO -- verification over genes and over categories


def group_stat(t, deltap=0.1):
    gid = 0
    t['gid'] = 0
    sumnp = 0
    for rid, r in t.iterrows():
        t.at[rid, 'gid'] = gid
        sumnp += r.p*r.N
        if sumnp > 0 and 1/ma.sqrt(sumnp) <= deltap:
            gid += 1
            sumnp = 0


def group_stat2(t, nmin=100):
    gid = 0
    t['gid'] = 0
    sumn = 0
    for rid, r in t.iterrows():
        t.at[rid, 'gid'] = gid
        sumn += r.n
        if sumn > nmin:
            gid += 1
            sumn = 0
