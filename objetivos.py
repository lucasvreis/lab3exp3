from auxiliares import *
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.odr import ODR, RealData
from scipy.optimize import curve_fit
from numpy import sqrt, diag, array, linspace
from collections import namedtuple

Ajustes = {}
Variáveis = {}

def Limpa(dados):
    Ajustes.clear()
    Variáveis.clear()


class Memório:
    def __init__(self, funcs: dict):
        self.funcs = funcs

    def armazena(self, dados):
        arg = self.funcs['arg']
        Variáveis.update({k: v(arg) for k, v in self.funcs.items() if k != 'arg'})


class AjusteLinear(namedtuple('AjusteLinear', 'tipo, x, ux, y, uy, ref')):

    def roda(self, dados):
        x, ux = Variáveis[self.x], Variáveis[self.ux]
        y, uy = Variáveis[self.y], Variáveis[self.uy]
        if self.tipo == 'odr':
            data = RealData(x, y, ux, uy)
            odr = ODR(data, linear, beta0=[1, 0],
                      ndigit=16, maxit=100)
            DadosAjuste = odr.run()
            p = DadosAjuste.beta
            u = sqrt(diag(DadosAjuste.cov_beta))
        elif self.tipo == 'cfit':
            p, cov = curve_fit(flinear, x, y, sigma=uy,
                               absolute_sigma=True, method='trf')
            u = sqrt(diag(cov))
        elif self.tipo == 'cfitse':
            p, cov = curve_fit(flinear, x, y, method='trf')
            u = sqrt(diag(cov))
        print(p, u)
        Ajustes[self.ref] = ({
            'a': p[0], 'ua': u[0],
            'b': p[1], 'ub': u[1]
        },
            x, ux, y, uy
        )


class GráficoAjuste(namedtuple('GráficoAjuste', 'tipo, ref, título, lx, ly, legendas')):
    def plota(self, dados):
        if self.tipo == 'ajuste':
            plt.figure()
            x, ux, y, uy = Ajustes[self.ref][1:]
            folha = Ajustes[self.ref][0]
            p = [folha['a'], folha['b']]
            u = [folha['ua'], folha['ub']]
            # plt.plot(x, y, ':')
            if ux[0] > 0:
                plt.errorbar(x, y, xerr=ux, yerr=uy, fmt='o', elinewidth=0, capsize=3, ecolor='#EE5747', ms=2, label='média das medições')
            else:
                plt.plot(x, y, 'o', ms=4, label='média das medições')
            uinterval = [[(p[0] + u[0]) * x[0] + p[1] + u[1], (p[0] + u[0]) * x[-1] + p[1] + u[1]],
                         [(p[0] - u[0]) * x[0] + p[1] - u[1], (p[0] - u[0]) * x[-1] + p[1] - u[1]]]
            plt.fill_between([x[0], x[-1]], uinterval[1], uinterval[0], alpha=0.3, edgecolor='#7BC814', facecolor='#8BC84F',
                             linewidth=0, linestyle=':', antialiased=True, label='incerteza do ajuste')
            plt.plot([x[0], x[-1]], [p[0] * x[0] + p[1], p[0] * x[-1] + p[1]], '-k', label='ajuste')

            handles, labels = plt.axes().get_legend_handles_labels()
            η = (Variáveis[self.ref + '-η']*1e3, Variáveis[self.ref + '-uη']*1e3)

            handles.append(Patch(color='none', label=f'a = {sigdig(p[0], u[0])} 1/m.s'))
            handles.append(Patch(color='none', label=f'b = {sigdig(p[1], u[1])} m/s'))
            handles.append(Patch(color='none', label=f'η = {sigdig(*η)} mPa.s'))  # !!!!
            plt.legend(handles=handles)
            plt.xlabel(self.lx)
            plt.ylabel(self.ly)
            plt.title(self.título)
            plt.savefig(f'ajuste-{self.ref}.png', dpi=400, bbox_inches='tight')
            plt.show()
        elif self.tipo == 'majuste':
            for k, ref in enumerate(self.ref):
                x, ux, y, uy = Ajustes[ref][1:]
                folha = Ajustes[ref][0]
                p = [folha['a'], folha['b']]
                u = [folha['ua'], folha['ub']]
                plt.plot(x, y, ':', color=f'C{k}')
                plt.errorbar(x, y, xerr=ux, yerr=uy, fmt=f'oC{k}', elinewidth=0, capsize=2, ecolor=f'C{k}', ms=2, label='médias de ' + ref)
                uinterval = [[(p[0] + u[0]) * x[0] + p[1] + u[1], (p[0] + u[0]) * x[-1] + p[1] + u[1]],
                             [(p[0] - u[0]) * x[0] + p[1] - u[1], (p[0] - u[0]) * x[-1] + p[1] - u[1]]]
                plt.fill_between([x[0], x[-1]], uinterval[1], uinterval[0], alpha=0.1, edgecolor='#7BC814', facecolor=f'C{k}',
                                 linewidth=0, linestyle=':', antialiased=True)
                plt.plot([x[0], x[-1]], [p[0] * x[0] + p[1], p[0] * x[-1] + p[1]], '-', color=f'C{k}', label='ajuste para ' + ref)
            handles, labels = plt.axes().get_legend_handles_labels()
            plt.legend(handles=handles)
            plt.xlabel(self.lx)
            plt.ylabel(self.ly)
            plt.title(self.título)
            plt.savefig(f'ajuste-{self.ref[0]}-{self.ref[-1]}.png', dpi=400)
            plt.show()

class GráficoComparação(namedtuple('GráficoComparação', 'xs, ys, refs, título, f, consts, lx, ly, legendas')):
    def plota(self, dados):
        for k, ref in enumerate(self.refs):
            xs = Variáveis[ref + self.xs]
            ys = Variáveis[ref + self.ys]
            τ = (Variáveis[ref + '-τ'], Variáveis[ref + '-uτ'])
            Ip = (Variáveis[ref + '-Ip'], Variáveis[ref + '-uIp'])
            locais = {'τ': τ[0], 'uτ': τ[1], 'Ip': Ip[0], 'uIp': Ip[1], **self.consts}
            xg = linspace(0.016, 0.12, 100)
            f = sp.lambdify('Δm', self.f.subs(locais), modules=['numpy'])
            # yg = vetórmula([{'Δm': x} for x in xg], self.f, locais)
            yg = f(xg)
            plt.plot(xg, yg, '-', color=f'C{k}', label='curva para ' + ref, zorder=0)
            plt.scatter(xs, ys, 8, c=f'C{k}', marker='o', label='média de ' + ref, zorder=1)
        plt.legend()
        plt.xlabel(self.lx)
        plt.ylabel(self.ly)
        plt.title(self.título)
        plt.savefig(f'comparação.png', dpi=400)
        plt.show()

