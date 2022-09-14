#!/usr/bin/env python
# coding: utf-8

"""
Run RBA
"""

import efmtool
import numpy as np
import pandas as pd
import plotext as plt
from general_functions import get_const_parameters, make_matrix, bisection_search_mu

k1, k2, k3, k5, c, kel, n1, n2, n3, n4, n5, Mtot, mwG, mwAA, mwNT, cell_mass, avogadro = get_const_parameters()

prot_fractions = np.arange(0.01, 1, 0.01)  # protein fractions to test
mus = np.arange(2,4,0.001) # growth rates to test

last_mus = []
prog = 0
plt.title("RBA")
plt.xlabel("Protein fraction in ribosome")
plt.ylabel("Growth rate [1/h]")
plt.theme("clear")
plt.plot_size(70,25)

for frac in prot_fractions:
    nR = Mtot*frac/mwAA
    mR = Mtot*(1-frac)/mwNT

    matrices = {}
    for mu in mus:
        S, rxns, row_names, revs = make_matrix(mwG, mwAA, mwNT, n1, n2, n3, n4, n5, nR, mR, mu, k1, k2, k3, k5, c, kel)
        matrices[str(round(mu,4))] = S    
    
    last_mu = bisection_search_mu(matrices, revs, rxns, row_names)
    last_mus.append(last_mu)
    prog += 1
    
    plt.clt()
    plt.ylim(2,4)
    plt.xlim(0,1)
    plt.scatter(prot_fractions[:prog], last_mus, marker="heart", color = 125)
    plt.show()  

res = pd.DataFrame({"x": prot_fractions, "mu": last_mus})
res.to_csv("../data/00_RBA.csv")

