#!/usr/bin/env python
# coding: utf-8

"""
Run RBA with:
    * fixed total RNAP (c_pol_max)
    * fixed ribosome composition (0.36% protein)
    * mu increased step by step
    * save RNAP fluxes
Three parameters sets:
    * noact: max. activities of R and RNAP NOT included (activities = False)
    * noacc_noact: max. activities of R and RNAP NOT included (activities = False)
                   no accumulation of rRNA or R
    * noacc_act: max. activities of R and RNAP included (activities = True)
                 no accumulation of rRNA or R
"""

import numpy as np
import pandas as pd
from general_functions import make_matrix, run_efmtool, make_nice_results, egv_dict_to_df, get_parameters

# Parameters
FILENAME = "fluxes_RNAPmax"
MATRIX_TYPE = "RNAPmax"
ORGANISM = "ecoli"
FRAC = 0.36

# Variable parameters
parameters = {"noact": {"activities": False,
                        "rnap_max": 0.0000126},
              "noacc_noact": {"activities": False,
                        "rnap_max": 0.0000126},
              "noacc_act": {"activities": True,
                      "rnap_max": 0.0000776}}
mus = np.arange(0.001,3.5,0.001)  # growth rates to test

# Simulations
for name, par in parameters.items():
    nice_egvs = {}
    for mu in mus:
        S, rxns, row_names, revs = make_matrix(MATRIX_TYPE, FRAC, mu, ORGANISM,
                                               activities=par["activities"],
                                               rnap_max=par["rnap_max"])
        # remove slack vars to prevent accumulation of rRNA/ribosome
        if "noacc" in name:
            S = np.delete(S, [11, 18], 1)  # remove slack columns of rRNA and R
            rxns.pop(11)
            rxns.pop(18)
            revs = [0] * S.shape[1]
        egvs = run_efmtool(S, revs, rxns, row_names)

        if egvs.shape[1] != 0:
            nice_egvs[str(round(mu,4))] = make_nice_results(egvs, rxns, False)

    all_results = egv_dict_to_df(nice_egvs, rxns, FRAC)
    all_results.to_csv(f"../data/{FILENAME}_{name}.csv")

# Bremer 2003 data -- convert to mmol/g/h
avogadro = 6.022e20
*_, mwR, mwG, mwAA, mwNT = get_parameters(ORGANISM)
nrrna = mwR*(1-FRAC)/mwNT
mus = [np.log(2)*mu for mu in [0.6, 1.0, 1.5, 2.0, 2.5]]  # growth rates (doublings/h)
rna_syn = [3,9.9,29,66.4,132.5]  # stable RNA synthesis rate (10^5 nt/cell/min)
masses = [mass*10**-15 for mass in [150,260,430,640,870]]  # dry masses
bremer_fluxes = [rate*10**5*60/nrrna/avogadro for rate in rna_syn]  # to mmol/cell/h
bremer_fluxes = [bremer_fluxes[i]/masses[i] for i in range(len(masses))]  # to mmol/g/h
pd.DataFrame({"mu": mus, "fluxes": bremer_fluxes}).to_csv("../data/fluxes_bremer.csv")
