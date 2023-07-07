#!/usr/bin/env python
# coding: utf-8

"""
Author: Diana Szeliova

Run RBA with or without R degradation rate
Parameters from "parameters.csv"
Save optimal growth rate, fluxes, concentrations and ribosome allocations
"""

import numpy as np
import pandas as pd
from general import Simulation

prot_fractions = np.arange(0.005, 1, 0.005)  # protein fractions in ribosome
# add numbers close to 0 and 1 (otherwise no solution, the matrix would have to be changed)
prot_fractions = np.concatenate(([0.00001], prot_fractions, [0.99999]))
growth_rates = np.arange(0.001, 4, 0.001)
parameters = pd.read_csv("../data/parameters_v2.csv")

results = {
    "growth_rates": pd.DataFrame(),
    "allocations": pd.DataFrame(),
    "mass_fractions": pd.DataFrame(),
    "fluxes": pd.DataFrame()
}

for index, row in parameters.iterrows():
    name = row["name"]
    par = dict(row)
    _ = par.pop("name")
    par["medium"] = int(par["medium"])

    sim = Simulation(par,
                     growth_rates = growth_rates,
                     prot_fractions = prot_fractions)
    sim.test_xrps(True)

    temp_results = {
        "growth_rates": sim.max_growth_rates, # 1/h
        "allocations": sim.allocations, # ribosome allocations
        "mass_fractions": sim.mass_fractions, # g/g
        "fluxes": sim.fluxes # mmol/gh
    }

    # add a column with an identifier and concatenate data frames
    for result_type in results:
        temp_results[result_type]["name"] = f"{par['matrix_type']}_{par['parameter_set']}_{name}"
        results[result_type] = pd.concat([results[result_type], temp_results[result_type]])

for result_type, df in results.items():
    df.to_csv(f"../data/RBA_{result_type}.csv")
