import itertools
import numpy as np
import matplotlib.pyplot as plt
import bisect
import copy
import random
import collections


def gen_uniformly_events(beta, lam):
    n = 0
    d = np.exp(-beta*lam)
    p = d
    zeta = np.random.rand()
    while zeta > p:
        n += 1
        d *= beta*lam/n
        p += d
    tauk = beta*np.random.rand(n)
    tauk.sort()
    return tauk, n


class cut():
    def __init__(self, beta, L, Hx):
        self.beta = beta
        self.L = L
        cuttau = collections.namedtuple('cuttau', ['t', 's', 'null'])
        self.C = []
        for ix in range(self.L):
            tauk, n = gen_uniformly_events(self.beta, Hx/2.0)
            if n > 0:
                Ctau = [cuttau(tauk[itau], random.choice([-1, 1]), False)
                        for itau in range(n)]
            else:
                Ctau = [cuttau(0, random.choice([-1, 1]), True)]
            self.C.append(Ctau)

    def __repr__(self,):
        s = ''
        for ix in range(self.L):
            for itau in range(len(self.C[ix])):
                s += str(self.C[ix][itau].t) + ' ' + \
                    str(self.C[ix][itau].s) + '\n'
        return s

    def remove_cut(self,):
        C_romoved = []
        for Cx in self.C:
            Cx_romoved = np.zeros((len(Cx), 2))
            n = 0
            for i in range(len(Cx)):
                if Cx[i, 1] * Cx[i-1, 1] < 0:
                    Cx_romoved[n] = Cx[i]
                    n += 1
            C_romoved.append(Cx_romoved[:n])
        self.C = C_romoved

    def __cut(Hx, L):
        cut = []
        for ix in range(L):
            tauk, n = gen_uniformly_events(Hx/2.0)
            cut.append(tauk)
        return cut


class Solver():
    def __init__(self, beta, J, Hx, L, ntot, measure_interval):
        self.beta = beta
        self.J = J
        self.Hx = Hx
        self.L = L
        self.ntot = ntot
        self.measure_interval = measure_interval
        self.C = []
        self.Ccut = []
        self.flip = []
        self.cluster = []
        self.cut = []
        self.bond = []

    def gen_bond(self,):
        bond = []
        for ix in range(self.L):
            tauk, n = self.gen_uniformly_events(self.J/2.0)
            bond.append(tauk)
        bond = [bondx.tolist() for bondx in bond]
        self.bond = bond

    def get_itau(self, Cx, tau): return bisect.bisect_left(Cx[:, 0], tau) - 1

    def get_spin(self, Cx, tau):
        itau = self.get_itau(Cx, tau)
        return Cx[itau, 1]

    def gen_cutspin(self,):
        cutspin = []
        for ix in range(len(self.cut)):
            cutspinx = []
            for itau in range(len(self.cut[ix])):
                cutspinx.append(
                    [self.cut[ix][itau], self.get_spin(self.C[ix], self.cut[ix][itau])])
            cutspin.append(cutspinx)
        return cutspin

    def gen_Ccut(self,):
        cutspin = self.gen_cutspin()
        Ccut = []
        for ix in range(len(self.C)):
            cutspinx = self.C[ix].tolist() + cutspin[ix]
            Ccut.append(np.array(sorted(cutspinx, key=lambda x: x[0])))
        self.Ccut = Ccut

    def remove_bond(self,):
        bond_removed = copy.deepcopy(self.bond)
        for ix in range(len(self.bond)):
            for itau in range(len(self.bond[ix])):
                tau = self.bond[ix][itau]
                if self.get_spin(self.Ccut[ix], tau) * self.get_spin(self.Ccut[(ix + 1) % self.L], tau) < 0:
                    bond_removed[ix][itau] = np.nan
        bond_removed = [[tau for tau in bond_removed[ix] if not np.isnan(tau)]
                        for ix in range(len(bond_removed))]
        self.bond = bond_removed

    def gen_flip(self, ):
        flip = []
        for ix in range(self.L):
            ntau = self.Ccut[ix].shape[0]
            flip.append(np.random.choice([-1, 1], ntau))
        self.flip = flip

    def gen_cluster(self, ):
        cluster = []
        for ix in range(self.L):
            ntau = self.Ccut[ix].shape[0]
            cluster.append([[ix, itau] for itau in range(ntau)])
        self.cluster = cluster

    def get_cluster_number(self, index_x, index_tau):
        ix = index_x
        itau = index_tau
        while ix != self.cluster[ix][itau][0] or itau != self.cluster[ix][itau][1]:
            ix, itau = self.cluster[ix][itau]
        return ix, itau

    def connect_bond_ixtau(self, ix, tau):
        itau1 = self.get_itau(self.Ccut[ix], tau)
        itau2 = self.get_itau(self.Ccut[(ix + 1) % self.L], tau)
        ix1, itau1 = self.get_cluster_number(ix, itau1)
        ix2, itau2 = self.get_cluster_number((ix+1) % self.L, itau2)
        if ix1 < ix2:
            self.cluster[ix2][itau2] = [ix1, itau1]
        elif ix1 > ix2:
            self.cluster[ix1][itau1] = [ix2, itau2]
        else:
            if itau1 < itau2:
                self.cluster[ix2][itau2] = [ix1, itau1]
            else:
                self.cluster[ix1][itau1] = [ix2, itau2]

    def connect_bond(self,):
        for ix in range(self.L):
            for tau in self.bond[ix]:
                self.connect_bond_ixtau(ix, tau)

    def cluster_flip(self,):
        for ix in range(self.L):
            for itau in range(len(self.Ccut[ix])):
                self.Ccut[ix][itau][1] = self.flip[self.cluster[ix]
                                                   [itau][0]][self.cluster[ix][itau][1]]

    def get_Ccutp(self,):
        Ccutp = []
        for Ccutx in self.Ccut:
            Ccutp.append(np.concatenate(
                [[[0, Ccutx[-1, 1]]], Ccutx, [[self.beta, Ccutx[-1, 1]]]]))
        return Ccutp

    def plot(self,):
        Ccutp = self.get_Ccutp()
        colors = {-1: 'w', 1: 'k'}
        fig, ax = plt.subplots()
        for ix, Cpx in enumerate(Ccutp):
            for itau in range(Cpx.shape[0]-1):
                ax.fill([ix, ix+1, ix+1, ix], [Cpx[itau, 0], Cpx[itau, 0], Cpx[itau+1, 0],
                        Cpx[itau+1, 0]], color=colors[int(Cpx[itau, 1])])
