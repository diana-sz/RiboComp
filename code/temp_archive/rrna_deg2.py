#!/usr/bin/env python
# coding: utf-8

"""
Run RBA with R degradation rate
    * parameter set "activities2" (transcription rate corrected for activity and allocation)
    * 6 media ("glc", "gly", "succ", "LB", "glcAA", "glyAA")
value KDEG optimized with "fit_deg_rate.py"

Save optimal growth rate, fluxes,
concentrations and ribosome allocations
"""

import numpy as np
import pandas as pd
from general import Simulation

# Parameters
media = {
    "glc": 2,
    "gly": 1,
    "succ": 0,
    "LB": 5,
    "glcAA": 4,
    "glyAA": 3
}

parameters = {
    "matrix_type": "R_deg_hill",
    "parameter_set": "activities2",
    "kdeg_max": 41.1
}

prot_fractions = np.arange(0.005, 1, 0.005)  # protein fractions to test
# add numbers close to 0 and 1 (otherwise no solution, the matrix would have to be changed)
prot_fractions = np.concatenate(([0.00001], prot_fractions, [0.99999])) 
growth_rates = np.arange(0.1, 4, 0.001)

# Simulations
results = {
    "growth_rates": pd.DataFrame(),
    "allocations": pd.DataFrame(),
    "mass_fractions": pd.DataFrame(),
    "fluxes": pd.DataFrame()
}

for name, medium in media.items():
    parameters["medium"] = medium

    sim = Simulation(parameters,
                     growth_rates = growth_rates,
                     prot_fractions = prot_fractions)
    sim.test_xrps()

    temp_results = {
        "growth_rates": sim.max_growth_rates, # 1/h
        "allocations": sim.allocations, # ribosome allocations
        "mass_fractions": sim.mass_fractions, # g/g
        "fluxes": sim.fluxes}  # mmol/gh

    # add a column with an identifier and concatenate data frames
    for result_type in results:
        temp_results[result_type]["name"] = f"{parameters['matrix_type']}2_{name}"
        results[result_type] = pd.concat([results[result_type], temp_results[result_type]])

for result_type, df in results.items():
    df.to_csv(f"../data/{parameters['matrix_type']}2_{result_type}.csv")
