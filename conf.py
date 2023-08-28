import numpy as np
import matplotlib.pyplot as plt
import bisect
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
    tauk = np.sort(np.random.rand(n)*beta)
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

    def __str__(self,):
        s = ''
        for itau in range(len(self.cuttau)):
            s += str(self.cuttau[itau].t) + ' ' + \
                str(self.cuttau[itau].s) + '\n'
        return s
    
    def __setitem__(self, position, value):
        self.cuttau[position] = value
    
    def get_itau(self, tau):
        return bisect.bisect_left(self.cuttau, tau, key=lambda x: x.t) - 1

    def remove_cut(self,):
        removed_cut = []
        for itau in range(len(self.cuttau)):
            if self.cuttau[itau].s*self.cuttau[itau-1].s < 0:
                removed_cut.append(self.cuttau[itau])
        if len(removed_cut) == 0:
            removed_cut.append(
                self.cut(0, self.cuttau[0].s, True))
        self.cuttau = removed_cut
    
    def get_spin(self, tau):
        itau = self.get_itau(tau)
        return self.cuttau[itau].s

    def add_cut(self,):
        tauk, n = gen_uniformly_events(self.beta, self.Hx/2.0)
        for tau in tauk:
            spin = self.get_spin(tau)
            bisect.insort(self.cuttau, self.cut(tau, spin, False), key=lambda x: x.t)
        if self.cuttau[0].null == True and len(self.cuttau) > 1:
            self.cuttau.pop(0)

    def measure_sz(self, ):
        sz = 0
        if len(self) > 1:
            for itau in range(len(self)):
                if itau == 0:
                    sz += self[itau - 1].s*(self[itau +1].t)
                elif itau == len(self) - 1:
                    sz += self[itau].s*(self.beta - self[itau].t)
                else:
                    sz += self[itau].s*(self[itau+1].t - self[itau].t)
        else:
            sz += self[0].s*self.beta
        return sz/self.beta

class Cuttaux():
    def __init__(self, beta, Hx, L):
        self.beta = beta
        self.Hx = Hx
        self.L = L
        self.cuttaux = [Cuttau(beta, Hx) for i in range(L)]

    def __getitem__(self, position):
        return self.cuttaux[position]

    def __repr__(self,):
        s = ''
        for ct in self.cuttaux:
            s += str(ct) + '--\n'
        return s
    
    def remove_cut(self, ):
        for ct in self.cuttaux:
            ct.remove_cut()

    def add_cut(self, ):
        for ct in self.cuttaux:
            ct.add_cut()
    
    def measure_sz(self, ):
        sz = 0
        for ct in self.cuttaux:
            sz += ct.measure_sz()
        return sz/self.L
    

class Bondtaux():
    def __init__(self, cuttaux, J):
        self.beta = cuttaux.beta
        self.L = cuttaux.L
        self.J = J
        self.bond = []
        for ix in range(self.L):
            tauk, n = gen_uniformly_events(self.beta, self.J/2.0)
            self.bond.append(tauk)
        self.remove_bond(cuttaux)
            
    def __getitem__(self, position):
        return self.bond[position]

    def __repr__(self,):
        s = ''
        for bt in self.bond:
            for b in bt:
                s += str(b) + '\n'
            s += '--\n'
        return s
    
    def __len__(self,):
        return len(self.bond)
    
    def remove_bond(self, cuttaux):
        removed_bond = []
        for ix in range(len(self.bond)):
            rbt = []
            for tau in self.bond[ix]:
                if cuttaux[ix].get_spin(tau) * cuttaux[(ix + 1) % self.L].get_spin(tau) > 0:
                    rbt.append(tau)
            removed_bond.append(rbt)
        self.bond = removed_bond