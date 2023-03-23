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
import matplotlib.pyplot as plt

experimental_fluxes = pd.read_csv("../data/fluxes_bremer.csv")

# simulation done only for some parameters to save time
parameters = pd.read_csv("../data/parameters_deg.csv")
parameter_subset = parameters.query("name == 'glc'")

# R, rRNA accumulation / excess rRNA degradation will be blocked
parameter_subset = parameter_subset.assign(name=parameter_subset.name + "_noacc")

deg_dict = {"deg": np.arange(1, 5, 1),
            "deg_hill-2": np.arange(1, 5, 1),
            "deg_hill-6": np.arange(14, 24, 2)}
growth_rates = np.arange(0.25, 2.5, 0.25)

for index, row in parameter_subset.iterrows():
    name = row["name"]
    par = dict(row)
    _ = par.pop("name")
    degs = deg_dict[par["matrix_type"]]

    for deg in degs:
        all_results = pd.DataFrame()
        par["kdeg_max"] = deg

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

            if not res.empty:
                res["growth_rate"] = mu
                all_results = pd.concat([all_results, res])

        plt.plot(all_results["growth_rate"], all_results["vRNAP"])
    plt.plot(experimental_fluxes["mu"], experimental_fluxes["fluxes_corrected"])
    plt.legend(list(degs)+["experimental"])
    plt.title(par["matrix_type"])
    plt.show()
    plt.close()



