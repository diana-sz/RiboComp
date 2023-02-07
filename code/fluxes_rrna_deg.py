#!/usr/bin/env python
# coding: utf-8

"""
Run RBA  with R degradation rate:
    * fixed ribosome composition (0.36% protein)
    * mu increased step by step
    * save RNAP fluxes
Three parameters sets:
    * three types of matrices with different degradation function: R_deg, R_deg2, R_deg_hill
    * activities = False ("noact") / True ("glc" - activities in minimal glucose medium used)
"""

import numpy as np
import pandas as pd
from general import Model

# Parameters
MEDIUM = 2
FRAC = 0.36

# fitted degradation rates (without / with activities)
deg_rates = {"R_deg2": (53.4, 16.9),
             "R_deg": (61.0, 18.15),
             "R_deg_hill": (157.2, 56.7)
            }

# activities
parameters = {"noact": {"medium": MEDIUM},
              "glc_noacc": {"parameter_set": "activities",
                            "medium": MEDIUM},
              "glc": {"parameter_set": "activities",
                      "medium": MEDIUM}}

growth_rates = np.arange(0.001,2.5,0.001)  # growth rates to test

# Simulations
for matrix_type, rates in deg_rates.items():
    all_results = pd.DataFrame()

    for name, par in parameters.items():
        par["matrix_type"] = matrix_type
        par["kdeg_max"] = rates[par.get("parameter_set") == "activities"]
        par["frac"] = FRAC
        model = Model(par)

        for mu in growth_rates:
            model.update_mu(mu)

            if "noacc" in name:
                model.matrix_c = np.delete(model.matrix_c, [13, 22, 23], 1)
                model.columns_c.remove("S0")
                model.columns_c.remove("S1")
                model.columns_c.remove("S2")

            model.run_efmtool(drop_slack = False)
            res = model.nice_fluxes

            res["growth_rate"] = mu
            res["prot_fraction"] = FRAC
            res["EGVs"] = res.index
            res["name"] = f"{matrix_type}_{name}"
            all_results = pd.concat([all_results, res])

    all_results.to_csv(f"../data/fluxes_{matrix_type}.csv")


# Bremer 2003 data -- convert to mmol/g/h
avogadro = 6.022e20
nrrna = model.mol_masses["R"]*(1-FRAC)/model.mol_masses["NT"]
mus = [np.log(2)*mu for mu in [0.6, 1.0, 1.5, 2.0, 2.5]]  # growth rates (doublings/h)
rna_syn = [3,9.9,29,66.4,132.5]  # stable RNA synthesis rate (10^5 nt/cell/min)
masses = [mass*10**-15 for mass in [150,260,430,640,870]]  # dry masses
bremer_fluxes = [rate*10**5*60/nrrna/avogadro for rate in rna_syn]  # to mmol/cell/h
bremer_fluxes = [bremer_fluxes[i]/masses[i] for i in range(len(masses))]  # to mmol/g/h
pd.DataFrame({"mu": mus, "fluxes": bremer_fluxes}).to_csv("../data/fluxes_bremer.csv")
