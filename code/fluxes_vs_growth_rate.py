#!/usr/bin/env python
# coding: utf-8

"""
Author: Diana Szeliova
Run RBA with or without rRNA degradation rate:
    * fixed ribosome composition (0.36% protein)
    * mu increased step by step
    * save RNAP fluxes
"""

import numpy as np
import pandas as pd
from model_functions import Model

# simulation done only for some parameters to save time
parameters = pd.read_csv("../data/parameters.csv")
parameter_subset = parameters.query("(name == 'glc') & (matrix_type != 'Kostinski')")

# duplicate the conditions and label them 'noacc'
# R, rRNA accumulation / excess rRNA degradation will be blocked for these
noacc = parameter_subset.assign(name=parameter_subset.name + "_noacc")
parameter_subset = pd.concat([parameter_subset, noacc])

growth_rates = np.geomspace(0.0001, 3, 500)

all_results = pd.DataFrame()
for index, row in parameter_subset.iterrows():
    name = row["name"]
    par = dict(row)

    print(f"Running {par['matrix_type']}_{par['parameter_set']}_{name}")

    model = Model(par)

    for mu in growth_rates:
        model.update_mu(mu)

        # remove rRNA & R accumulation / excess rRNA degradation
        if ("noacc" in name) and ("extended" in par['matrix_type']):
            model.matrix_c = np.delete(model.matrix_c, [13, 21, 22], 1)
            model.columns_c.remove("S0")
            model.columns_c.remove("S8")
            model.columns_c.remove("S9")

        # remove rRNA & R accumulation
        elif ("noacc" in name) and ("base" in par['matrix_type']):
            model.matrix_c = np.delete(model.matrix_c, [11, 18], 1)
            model.columns_c.remove("S0")
            model.columns_c.remove("S7")

        egvs = model.run_efmtool()
        res = model.make_nice_results(egvs, drop_slack = False)

        res["growth_rate"] = mu
        res["prot_fraction"] = model.frac
        res["EGVs"] = res.index
        res["name"] = f"{par['matrix_type']}_{par['parameter_set']}_{name}"
        all_results = pd.concat([all_results, res])

all_results.to_csv("../data/fluxes_x0.36.csv")


# Bremer 2008 data, converted to mmol/g/h
AVOGADRO = 6.022e20
nrrna = model.parameters["R"]*(1-model.frac)/model.parameters["NT"]
growth_rates = [np.log(2)*mu for mu in [0.6, 1.0, 1.5, 2.0, 2.5, 3]]  # doublings/h -> 1/h
rna_syn = [3.5, 11, 29, 65, 113, 161]  # stable RNA synthesis rate (10^5 nt/cell/min)
dry_masses = [mass*10**-15 for mass in [226, 374, 555, 774, 921, 1023]]
fluxes = [rate*10**5*60/nrrna/AVOGADRO for rate in rna_syn]  # to mmol/cell/h
fluxes = [fluxes[i]/dry_masses[i] for i in range(len(dry_masses))]  # to mmol/g/h

# Gausing 1977 data - correction for fraction of degraded rRNA
correction_factors = [1.29, 1.14, 1.11, 1.11, 1.11, 1.11]
corrected_fluxes = [flux*correction for flux,correction in zip(fluxes, correction_factors)]

df = pd.DataFrame({"mu": growth_rates,
                   "fluxes": fluxes,
                   "fluxes_corrected": corrected_fluxes})
df.to_csv("../data/fluxes_bremer.csv", index = False)
