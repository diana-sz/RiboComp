#!/usr/bin/env python
# coding: utf-8

"""
Run RBA with a limit on RNAP flux (vmax)
"""

import efmtool
import numpy as np
import pandas as pd
import plotext as plt
from general_functions import make_nice_results, get_const_parameters, make_matrix_v4max, bisection_search_mu

k1, k2, k3, k5, c, kel, n1, n2, n3, n4, n5, Mtot, mwG, mwAA, mwNT, cell_mass, avogadro = get_const_parameters()

prot_fractions = np.arange(0.01, 1, 0.01)  # protein fractions to test
mus = np.arange(0.01,4,0.001)  # growth rates to test
vmaxes = [1000, 0.0005, 0.0002, 0.00002, 0.000002]  # vmaxes to test

all_mus = {}
for vmax in vmaxes:
    plt.title(f"RBA+vmax {vmax}")
    plt.xlabel("Protein fraction in ribosome")
    plt.ylabel("Growth rate [1/h]")
    plt.theme("clear")
    plt.plot_size(70,25)

    prog = 0
    last_mus = []
    for frac in prot_fractions:
        nR = Mtot*frac/mwAA
        mR = Mtot*(1-frac)/mwNT

        matrices = {}
        for mu in mus:
            S, rxns, row_names, revs = make_matrix_v4max(mwG, mwAA, mwNT, n1, n2, n3, n4, n5, nR, mR, mu, k1, k2, k3, k5, c, kel, vmax, 0)
            matrices[str(round(mu,4))] = S    

        last_mu = bisection_search_mu(matrices, revs, rxns, row_names)
        last_mus.append(last_mu)
        prog += 1

        plt.clt()
        plt.ylim(0,4)
        plt.xlim(0,1)
        
        try:
            for v in all_mus:
                plt.scatter(prot_fractions, all_mus[v], marker="heart", color = 125)
        except KeyError:
            pass
        plt.scatter(prot_fractions[:prog], last_mus, marker="heart", color = 141)
        plt.show()

    all_mus[vmax] = last_mus

all_results = pd.DataFrame(columns=["x", "mu", "vmax"])
for vmax in all_mus:
    res = pd.DataFrame({"x": prot_fractions, "mu": all_mus[vmax]})
    res["vmax"] = vmax
    all_results = pd.concat([all_results, res])

all_results.to_csv("../data/09_vmax.csv")





