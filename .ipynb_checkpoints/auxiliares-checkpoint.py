import sympy as sp
import math
from sympy.parsing.sympy_parser import parse_expr, standard_transformations
from scipy.odr import Model
from numpy import mean, vectorize
from numpy.linalg import norm
from scipy.stats import sem

# Funções auxiliares

def média(dados, árvore, nível, nome, restrições=None, exclusões=()):
    for subconjunto in macaco(árvore, corte=nível, restrições=restrições, exclusões=exclusões):
        lista = []
        for galho in macaco(árvore, list(subconjunto), restrições=restrições, exclusões=exclusões):
            lista.append(cata_galho(dados, galho)[1][nome])
        yield mean(lista)

def desvio(dados, árvore, nível, nome, restrições=None, exclusões=()):
    for subconjunto in macaco(árvore, corte=nível, restrições=restrições, exclusões=exclusões):
        lista = []
        for galho in macaco(árvore, list(subconjunto), restrições=restrições, exclusões=exclusões):
            lista.append(cata_galho(dados, galho)[1][nome])
        yield sem(lista)

def utri(a):
    return a/(2*math.sqrt(6))

def uret(a):
    return a/(2*math.sqrt(3))

def combinada(*args):
    return norm(args)

def funfórmula(pares, f, subs=None):
    if subs is None:
        subs = {}
    subs.update(pares)
    return float(f.evalf(16, subs=subs))

def vetórmula(pares, f, subs=None):
    return vectorize(lambda p: funfórmula(p, f, subs))(pares)

def parse(s):
    return parse_expr(s, transformations=standard_transformations)

def sigdig(v: float, u:float):
    if u == 0:
        return f'({v} ± 0)'
    sdig = -math.floor(math.log10(abs(u)))
    if sdig > 0:
        return f'({round(v, sdig)} ± {round(u, sdig)})'
    else:
        return f'({int(round(v, sdig))} ± {int(round(u, sdig))})'

def macaco(árvore, galho=None, corte=None, restrições=None, exclusões=()):
    macaquices = []
    if galho is None:
        galho = []
    if corte is None:
        corte = len(árvore)
    if len(galho) < corte:
        for g in árvore[len(galho)][0]:
            a = (len(galho), g)
            r = restrições
            e = exclusões
            if (r is None or a[0] not in list(zip(*r))[0] or a in r) and (a not in e):
                macaquices += macaco(árvore, galho + [g], corte, restrições, exclusões)
        return macaquices
    else:
        return [(*galho,)]

def cata_galho(dados, endereço):
    galho = dados
    for k in endereço:
        galho = galho[k]
    return galho


def incerteza(f: sp.Expr) -> sp.Expr:
    uf = 0
    udict = {}
    for s in f.free_symbols:
        udict[s] = sp.Symbol('u' + str(s))
        uf += udict[s] ** 2 * sp.simplify(f.diff(s) ** 2)
    uf = sp.sqrt(uf)
    return uf


def incertezarápido(f: sp.Expr) -> sp.Expr:
    uf = 0
    udict = {}
    for s in f.free_symbols:
        udict[s] = sp.Symbol('u' + str(s))
        uf += udict[s] ** 2 * (f.diff(s) ** 2)
    uf = sp.sqrt(uf)
    return uf

def float_input(prompt=""):
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("Sinto muito, mas eu só gosto de floats...")

def int_input(prompt=""):
    while True:
        try:
            return int(input(prompt))
        except ValueError:
            print("Sinto muito, mas eu só gosto de ints...")

def flinear(x: float, a: float, b: float) -> float:
    return a * x + b

linear = Model(lambda beta, x: beta[0]*x + beta[1])
