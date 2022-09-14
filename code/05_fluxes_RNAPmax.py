#!/usr/bin/env python
# coding: utf-8

"""
Run RBA with:
* fixed total RNAP (c_pol_max)
* fixed ribosome composition (0.36% protein)
* mu increased step by step 
* save RNAP fluxes
"""

import efmtool
import numpy as np
import pandas as pd
from general_functions import make_nice_results, get_const_parameters, make_matrix_cpolmax, egv_dict_to_df, run_efmtool

k1, k2, k3, k5, c, kel, n1, n2, n3, n4, n5, Mtot, mwG, mwAA, mwNT, cell_mass, avogadro = get_const_parameters()
c_pol_max = 0.0000126  # fitted so the trade-off is at 0.36 protein fraction
mus = np.arange(0.01,3.5,0.001)  # growth rates to test
frac = 0.36
nR = Mtot*frac/mwAA
mR = Mtot*(1-frac)/mwNT

nice_egvs = {}
for mu in mus:
    S, rxns, row_names, revs  = make_matrix_cpolmax(mwG, mwAA, mwNT, n1, n2, n3, n4, n5, nR, mR, mu, k1, k2, k3, k5, c, kel, c_pol_max)

    egvs = run_efmtool(S, revs, rxns, row_names)

    # save data
    if egvs.shape[1] != 0:
        nice_egvs[str(round(mu,4))] = make_nice_results(egvs, rxns, False)

all_results = egv_dict_to_df(nice_egvs, rxns, frac)
all_results.to_csv("../data/05_fluxes_RNAPmax.csv")


# bremer data -- convert to mmol/g/h
mus = [np.log(2)*mu for mu in [0.6, 1.0, 1.5, 2.0, 2.5]]  # growth rates (doublings/h)
rna_syn = [3,9.9,29,66.4,132.5]  # stable RNA synthesis rate (10^5 nt/cell/min)
masses = [mass*10**-15 for mass in [150,260,430,640,870]]  # dry masses
bremer_fluxes = [rate*10**5*60/mR/avogadro for rate in rna_syn]  # to mmol/cell/h
bremer_fluxes = [bremer_fluxes[i]/masses[i] for i in range(len(masses))]  # to mmol/g/h
pd.DataFrame({"mu": mus, "fluxes": bremer_fluxes}).to_csv("../data/fluxes_bremer.csv")

