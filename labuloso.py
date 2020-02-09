#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:08:28 2019

@author: Lucas Viana
"""
from objetivos import *
from auxiliares import *
from math import sqrt, pi
import pickle
import csv

Dados = {}  # Árvore de casos

# Incerteza dos instrumentos
uBalança = uret(1e-3)  # kg
uTermômetro = combinada(utri(0.1), uret(0.3))
uVoltímetro = uret(0.01)

# Diâmetros das esferas em mm
D = [2, 2.99, 3.49, 3.95]
r = list(d * 1e-3 / 2 for d in D)

# Dicionário das variáveis mensuradas por casos
dvar = dict.fromkeys(['r', 't'], 0.)

consts = {  # S.I. Constantes e suas incertezas
    'c': 1, # 1 caloria
    'uc': 0 
}

# Fórmulas
f = {
    'K': parse('(1 + 2.4 * (r/R))*(1 + 3.3 * (r/H))'),
    'vlim': parse('2*(ρ - ρl)*g*r**2/(9*η*K)'),
    'x': parse('r**2/(1 + 2.4 * (r/R))*(1 + 3.3 * (r/H))'),
    'y': parse('d/t'),
    'η': parse('2*(ρ - ρl)*g / (9*a)')
}
u = {n: incerteza(f) for n, f in f.items()}
from sympy.printing import pprint
# pprint(u['x'])
# print(u['x'].evalf(8, dict({'r': r[0]}, **consts)))
# print(f['x'].evalf(8, dict({'r': r[0]}, **consts)))

árvore = [
    (list(range(len(r))), lambda galho, e: f"Diâmetro = {round(D[e], 4)} mm"),
    (list(range(9)), lambda galho, e: f"Intervalo #{e}"),
    (list(range(5)), lambda galho, e: f"Lançamento #{e}, t = {galho[e][1]['t']} s")
]



objetivos = [
    Limpa,
    Memório({
        'x': lambda arg: vetórmula(
            [{'r':raio} for raio in r for i in range(1)], f['x'], consts
        ),
        'ux': lambda arg: vetórmula(
            [{'r':raio} for raio in r for i in range(1)], u['x'], consts
        ),
        'y': lambda arg: vetórmula(
            [{'t':t} for t in média(Dados, árvore, 1, 't', exclusões=((1,8),(1,7),(1,6),(1,5)))],
            f['y'],
            consts
        ),
        'uy': lambda arg: vetórmula(
            [{'t': m, 'ut': combinada(d, uCronômetro)}
                for m, d in zip(média(Dados, árvore, 1, 't', exclusões=((1,8),(1,7),(1,6),(1,5))),
                                desvio(Dados, árvore, 1, 't', exclusões=((1,8),(1,7),(1,6),(1,5))))],
            u['y'],
            consts
        ),
        'uxl': lambda arg: [0]*4,
        'uyl': lambda arg: [0]*4,
        # pessoa + 't': lambda arg: list(média(Dados, árvore, 2, 't', (1, arg))),
        # pessoa + 'Δm': lambda arg: list([round(massasMa[i] - massasMe[i], 4) for i in range(6)]),
        'arg': ""
    }).armazena,
    AjusteLinear(
        tipo='odr',
        x='x',
        ux='ux',
        y='y',
        uy='uy',
        ref="err"
    ).roda,
    AjusteLinear(
        tipo='cfitse',
        x='x',
        ux='uxl',
        y='y',
        uy='uyl',
        ref="nerr"
    ).roda,
    Memório({
        'err-η': lambda arg: f['η'].evalf(8, subs=dict(Ajustes[arg][0], **consts)),
        'err-uη': lambda arg: u['η'].evalf(8, subs=dict(Ajustes[arg][0], **consts)),
        'arg': "err"
    }).armazena,
    Memório({
        'nerr-η': lambda arg: f['η'].evalf(8, subs=dict(Ajustes[arg][0], **consts)),
        'nerr-uη': lambda arg: u['η'].evalf(8, subs=dict(Ajustes[arg][0], **consts)),
        'arg': "nerr"
    }).armazena,
    GráficoAjuste(
        tipo='ajuste',
        ref="err",
        título=f"Ajuste Linear da equação modelo, com erros (cronômetro)",
        ly='Vlim - variável dependente',
        lx='r²/K - variável independente',
        legendas=None
    ).plota,
    GráficoAjuste(
        tipo='ajuste',
        ref="nerr",
        título=f"Ajuste Linear da equação modelo, sem erros (cronômetro)",
        ly='Vlim - variável dependente',
        lx='r²/K - variável independente',
        legendas=None
    ).plota
]

def edita_galho(endereço):
    global Dados
    while True:
        galho = cata_galho(Dados, endereço)
        if type(galho) is not dict:
            edita_folha(endereço, galho)
            return edita_galho(endereço[:-1])
        print(f"\n ~ Editando o galho de endereço: /{'/'.join(str(e) for e in endereço)}")
        print("Opções:")
        mapa = {n:k for n, k in enumerate(galho.keys())}
        for e in mapa.items():
            print(f"{e[0]}: {árvore[len(endereço)][1](galho, e[1])}")
        entrada = input("\n'int': entrar | 'r': recursivo | 'c': voltar | 's': salvar | 'l': ler")
        if entrada.isdigit():
            entrada = int(entrada)
            if entrada in mapa:
                return edita_galho(endereço + [mapa[entrada]])
            else:
                print("Índice inválido")
        elif entrada == 'r':
            for e in macaco(árvore, endereço):
                edita_folha(e, cata_galho(Dados, e))
        elif entrada == 'rd':
            for e in macaco(árvore, endereço):
                folha = cata_galho(Dados, e)
                if folha[1] == dvar:
                    edita_folha(e, folha)
        elif entrada == 'c':
            return edita_galho(endereço[:-1])
        elif entrada == 's':
            try:
                with open(input("Nome do arquivo"), 'wb') as arquivo:
                   pickle.dump(Dados, arquivo)
                print("Sucesso")
            except:
                print("Erro")
        elif entrada == 'l':
            try:
                with open(input("Nome do arquivo"), 'rb') as arquivo:
                    Dados = pickle.load(arquivo)
                print("Sucesso")
            except:
                print("Erro")
        elif entrada == 'o':
            for o in objetivos:
                o(Dados)
        elif entrada == 'e':
            return


def edita_folha(endereço, folha):
    print(f"\n ~ Editando a folha de endereço: /{'/'.join(str(e) for e in endereço)}")
    raio = round(r[endereço[0]], 4)
    print(f"Diâmetro: {D[endereço[0]]}")
    print(f"Nesse caso, {', '.join(f'{e[0]} = {e[1]}' for e in folha[1].items())}")
    entrada = float_input("\nTempo de queda, em segundos")
    t = entrada
    folha[1].update({
        'r': raio,
        't': t
    })
    return


def planta_árvore():
    for endereço in macaco(árvore):
        aqui = Dados
        for i, galho in enumerate(endereço):
            if galho in aqui:
                aqui = aqui[galho]
            else:
                if i == len(endereço) - 1:
                    aqui[galho] = ('f', dvar.copy())
                else:
                    aqui[galho] = {}
                    aqui = aqui[galho]


planta_árvore()

# with open('dados.csv', newline='') as csvfile:
#     reader = csv.reader(csvfile, delimiter=',', quotechar='|')
#     data = list([list([c for c in r]) for r in reader])
#     for d in range(4):
#         for n in range(9):
#             for l in range(5):
#                 Dados[d][n][l][1]['r'] = r[d]
#                 Dados[d][n][l][1]['t'] = float(data[15 + 12*d + l][2 + 3*n])



edita_galho([])

