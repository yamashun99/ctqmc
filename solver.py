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
        self.spins = np.random.choice([-1, 1], size=(self.L, self.M))

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
        tauk, n = self.gen_cut(1.0, self.beta)
        C = np.zeros((self.L, n, 2))
        for i in range(n):
            C[i, 0] = tauk[i]
            C[i, 1] = np.random.choice([-1, 1])
        return C, n

    def remove_cut(self, C):
        C_romoved = np.zeros((len(C), 2))
        n = 0
        for i in range(len(C)):
            if C[i, 1] * C[i-1, 1] < 0:
                C_romoved[n] = C[i]
                n += 1
        return C_romoved[:n]
