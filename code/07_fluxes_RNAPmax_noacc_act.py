#!/usr/bin/env python
# coding: utf-8

"""
Run RBA with:
* fixed total RNAP (c_pol_max)
* fixed ribosome composition (0.36% protein)
* no accumulation of R and RNA
* activities of R and RNAP included
* mu increased step by step 
* save RNAP fluxes
"""

import efmtool
import numpy as np
import pandas as pd
from general_functions import make_nice_results, get_const_parameters, make_matrix_cpolmax, egv_dict_to_df, run_efmtool

k1, k2, k3, k5, c, kel, n1, n2, n3, n4, n5, Mtot, mwG, mwAA, mwNT, cell_mass, avogadro = get_const_parameters()

# activities from Bremer
kel = kel*0.8
c = c*0.3

c_pol_max = 0.00004 # fitted so the trade-off is at 0.36 protein fraction
mus = np.arange(0.01,3.5,0.001)  # growth rates to test
frac = 0.36
nR = Mtot*frac/mwAA
mR = Mtot*(1-frac)/mwNT

nice_egvs = {}
for mu in mus:
    S, rxns, row_names, revs  = make_matrix_cpolmax(mwG, mwAA, mwNT, n1, n2, n3, n4, n5, nR, mR, mu, k1, k2, k3, k5, c, kel, c_pol_max)
    S = np.delete(S, [11, 18], 1)  # remove slack columns of rRNA and R
    rxns.pop(11)
    rxns.pop(18)
    revs = [0] * S.shape[1]

    egvs = run_efmtool(S, revs, rxns, row_names)

    # save data
    if egvs.shape[1] != 0:
        nice_egvs[str(round(mu,4))] = make_nice_results(egvs, rxns, False)

all_results = egv_dict_to_df(nice_egvs, rxns, frac)
all_results.to_csv("../data/07_fluxes_RNAPmax_noacc_act.csv")

