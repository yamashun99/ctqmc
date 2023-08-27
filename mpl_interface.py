import matplotlib as mpl
import matplotlib.pyplot as plt
import copy


def spinplot(self, cuttaux):
    colors = {-1: 'w', 1: 'k'}
    for ix, cuttau in enumerate(cuttaux):
        ltau = len(cuttau)
        xaxis = [ix, ix+1, ix+1, ix]
        self.fill(xaxis, [0.0, 0.0,  cuttau[0].t,
                cuttau[0].t], color=colors[int(cuttau[-1].s)])
        for itau in range(ltau-1):
            self.fill(xaxis, [cuttau[itau].t, cuttau[itau].t,  cuttau[itau+1].t,
                    cuttau[itau+1].t], color=colors[int(cuttau[itau].s)])
        self.fill(xaxis, [cuttau[ltau - 1].t,
                cuttau[ltau - 1].t, cuttaux.beta, cuttaux.beta], color=colors[int(cuttau[ltau - 1].s)])

def bondplot(self, bondtaux):
    for ix, bondtau in enumerate(bondtaux):
        for tau in bondtau:
            self.plot([ix+0.5, ix+1.5], [tau, tau], 'C0')

def cutplot(self, cuttaux):
    for ix, cuttau in enumerate(cuttaux):
        if len(cuttau) > 1:
            for cut in cuttau:
                self.plot([ix, ix+1], [cut.t, cut.t], 'C1')

mpl.axes.Axes.spinplot = spinplot
mpl.axes.Axes.bondplot = bondplot
mpl.axes.Axes.cutplot = cutplot