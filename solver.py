import itertools
import numpy as np
import matplotlib.pyplot as plt
from numba import njit, jit, float32, int32
from numba.experimental import jitclass

# spec = [
#    ('beta', float32),
#    ('J', float32),
#    ('Hx', float32),
#    ('L', int32),
#    ('M', int32),
#    ('ntot', int32),
#    ('measure_interval', int32),
#    ('spins', int32[:, :]),
# ]
#
#
# @jitclass(spec)


class Solver():
    def __init__(self, beta, J, Hx, L, ntot, measure_interval):
        self.beta = beta
        self.J = J
        self.Hx = Hx
        self.L = L
        self.ntot = ntot
        self.measure_interval = measure_interval

    def gen_cut(self, lam):
        n = 0
        d = np.exp(-self.beta*lam)
        p = d
        zeta = np.random.rand()
        while zeta > p:
            n += 1
            d *= self.beta*lam/n
            p += d
        tauk = self.beta*np.random.rand(n)
        tauk.sort()
        return tauk, n

    def gen_conf(self,):
        Cxtau = []
        for ix in range(self.L):
            tauk, n = self.gen_cut(1.0)
            Ctau = np.zeros((n, 2))
            for itau in range(n):
                Ctau[itau, 0] = tauk[itau]
                Ctau[itau, 1] = np.random.choice([-1, 1])
            Cxtau.append(Ctau)
        return Cxtau, n

    def remove_cut(self, C):
        C_romoved = np.zeros((len(C), 2))
        n = 0
        for i in range(len(C)):
            if C[i, 1] * C[i-1, 1] < 0:
                C_romoved[n] = C[i]
                n += 1
        return C_romoved[:n]
