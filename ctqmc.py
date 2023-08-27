import itertools
import numpy as np
import matplotlib.pyplot as plt
import bisect
import copy
import random
import collections
import conf
import random


class Updater():
    def __init__(self, beta, J, Hx, L):
        self.beta = beta
        self.J = J
        self.Hx = Hx
        self.L = L
        self.cuttaux = conf.Cuttaux(beta, Hx, L)
        self.cuttaux.remove_cut()
        self.bond = conf.Bondtaux(self.cuttaux, self.J)
        self.cluster = self._gen_cluster(self.cuttaux)


    def _gen_flip(self, cuttaux):
        flip = []
        for ix in range(self.L):
            ntau = len(cuttaux[ix])
            flip.append([random.choice([-1, 1]) for itau in range(ntau)])
        return flip

    def _gen_cluster(self, cuttaux):
        cluster = []
        for ix in range(self.L):
            ntau = len(cuttaux[ix])
            cluster.append([[ix, itau] for itau in range(ntau)])
        return cluster

    def _get_cluster_number(self, index_x, index_tau):
        ix = index_x
        itau = index_tau
        while ix != self.cluster[ix][itau][0] or itau != self.cluster[ix][itau][1]:
            ix, itau = self.cluster[ix][itau]
        return ix, itau

    def _connect_cluster(self, cuttaux, ix, tau):
        #connected_cluster = copy.deepcopy(cluster)
        itau1 = cuttaux[ix].get_itau(tau)
        itau2 = cuttaux[(ix+1)%self.L].get_itau(tau)
        ix1, itau1 = self._get_cluster_number(ix, itau1)
        ix2, itau2 = self._get_cluster_number((ix+1) % self.L, itau2)
        if ix1 < ix2:
            self.cluster[ix2][itau2] = [ix1, itau1]
        elif ix1 > ix2:
            self.cluster[ix1][itau1] = [ix2, itau2]
        else:
            if itau1 < itau2:
                self.cluster[ix2][itau2] = [ix1, itau1]
            else:
                self.cluster[ix1][itau1] = [ix2, itau2]

    def _connect_clusterbybond(self, cuttaux,  bond):
        #connected_cluster = copy.deepcopy(cluster)
        for ix in range(self.L):
            for tau in bond[ix]:
                self._connect_cluster(cuttaux, ix, tau)

    def _cluster_flip(self, cluster, flip):
        for ix in range(self.L):
            for itau in range(len(self.cuttaux[ix])):
                self.cuttaux[ix][itau] = self.cuttaux[ix][itau]._replace(s = 
                    flip[cluster[ix][itau][0]][cluster[ix][itau][1]])
    
    def remove_cut(self):
        self.cuttaux.remove_cut()
    
    def add_cut(self):
        self.cuttaux.add_cut()
    
    def update_bond(self): 
        self.bond = conf.Bondtaux(self.cuttaux, self.J)

    def update(self):
        self.cluster = self._gen_cluster(self.cuttaux)
        flip = self._gen_flip(self.cuttaux)
        self._connect_clusterbybond(self.cuttaux, self.bond)
        self._cluster_flip(self.cluster, flip)
    
    def measure_sz(self, ):
        return self.cuttaux.measure_sz()