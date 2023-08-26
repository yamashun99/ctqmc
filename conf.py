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


class Cuttau():
    def __init__(self, beta, Hx):
        self.beta = beta
        self.cut = collections.namedtuple(
            'cut', ['t', 's', 'null'])
        self.Hx = Hx
        tauk, n = gen_uniformly_events(self.beta, Hx/2.0)
        self.cuttau = []
        if n > 1:
            self.cuttau = [self.cut(tauk[itau], random.choice([-1, 1]), False)
                           for itau in range(n)]
        else:
            self.cuttau = [self.cut(0, random.choice([-1, 1]), True)]

    def __len__(self):
        return len(self.cuttau)

    def __getitem__(self, position):
        return self.cuttau[position]

    def __repr__(self,):
        s = ''
        for itau in range(len(self.cuttau)):
            s += str(self.cuttau[itau].t) + ' ' + \
                str(self.cuttau[itau].s) + '\n'
        s += '--\n'
        return s

    def remove_cut(self,):
        removed_cut = []
        for itau in range(len(self.cuttau)):
            print(itau)
            if self.cuttau[itau].s*self.cuttau[itau-1].s < 0:
                removed_cut.append(self.cuttau[itau])
        if len(removed_cut) == 0:
            removed_cut.append(
                self.cut(0, self.cut[0].s, True))
        self.cut = removed_cut

    def add_cut(self,):
        for ix in range(self.L):
            tauk, n = gen_uniformly_events(self.beta, self.Hx/2.0)
            if n > 0:
                confitaus = [bisect.bisect_left()]
                self.cut[ix].append(
                    self.cuttau(ix, tauk[itau], random.choice([-1, 1]), False))
                cutt = [self.cuttau(ix, tauk[itau], random.choice([-1, 1]), False)
                        for itau in range(n)]
            self.cut[ix] += cutt
