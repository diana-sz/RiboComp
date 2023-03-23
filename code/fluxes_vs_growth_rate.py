#!/usr/bin/env python
# coding: utf-8

"""
Run RBA with rRNA degradation rate:
    * fixed ribosome composition (0.36% protein)
    * mu increased step by step
    * save RNAP fluxes
"""

import numpy as np
import pandas as pd
from general import Model

# simulation done only for some parameters to save time
parameters = pd.read_csv("../data/parameters_deg.csv")
parameter_subset = parameters.query("name == 'glc'")

# for one condition R, rRNA accumulation / excess rRNA degradation will be blocked
noacc = parameter_subset.assign(name=parameter_subset.name + "_noacc")
parameter_subset = pd.concat([parameter_subset, noacc])

growth_rates = np.geomspace(0.0001, 3, 500) #np.arange(0.0001,3,0.005)

all_results = pd.DataFrame()
for index, row in parameter_subset.iterrows():
    name = row["name"]
    par = dict(row)
    _ = par.pop("name")

    model = Model(par)

    for mu in growth_rates:
        model.update_mu(mu)

        if "noacc" in name:
            model.matrix_c = np.delete(model.matrix_c, [13, 21, 22], 1)
            model.columns_c.remove("S0")
            model.columns_c.remove("S8")
            model.columns_c.remove("S9")

        egvs = model.run_efmtool()
        res = model.make_nice_results(egvs, drop_slack = False)

        res["growth_rate"] = mu
        res["prot_fraction"] = model.frac
        res["EGVs"] = res.index
        res["name"] = f"{par['matrix_type']}_{name}"
        all_results = pd.concat([all_results, res])

all_results.to_csv(f"../data/fluxes_x0.36.csv")


# Bremer 2003 data, converted to mmol/g/h
avogadro = 6.022e20
nrrna = model.mol_masses["R"]*(1-model.frac)/model.mol_masses["NT"]
mus = [np.log(2)*mu for mu in [0.6, 1.0, 1.5, 2.0, 2.5]]  # growth rates (doublings/h)
rna_syn = [3,9.9,29,66.4,132.5]  # stable RNA synthesis rate (10^5 nt/cell/min)
masses = [mass*10**-15 for mass in [150,260,430,640,870]]  # dry masses
bremer_fluxes = [rate*10**5*60/nrrna/avogadro for rate in rna_syn]  # to mmol/cell/h
bremer_fluxes = [bremer_fluxes[i]/masses[i] for i in range(len(masses))]  # to mmol/g/h

# .... data - correction for fraction of degraded rRNA
fraction_degraded = [1.3, 1.13, 1.1, 1.1, 1.1]
corrected_fluxes = [i*f for i,f in zip(bremer_fluxes, fraction_degraded)]

df = pd.DataFrame({"mu": mus, 
                   "fluxes": bremer_fluxes,
                   "fluxes_corrected": corrected_fluxes})
df.to_csv("../data/fluxes_bremer.csv", index = False)
