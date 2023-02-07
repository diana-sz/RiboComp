#!/usr/bin/env python
# coding: utf-8

"""
Run RBA with R degradation rate
15 parameters sets (see function "make_matrix()" in "general_functions.py"):
    * three types of matrices with different degradation function: R_deg, R_deg2, R_deg_hill
    * activities = True / False
    * archaea ("arch") / ecoli ("noact") (for activities=False)
    * 3 media ("glc", "gly", "succ") (for activities=True)
values in "deg_rates" optimized with "fit_deg_rate.py"

Save optimal growth rate and in case of R_deg_hill also fluxes,
concentrations and phis (ribosome allocations)
"""

import numpy as np
import pandas as pd
from general import Simulation

# fitted degradation rates (without / with activities)
deg_rates = {"R_deg2": (53.4, 16.9),
             "R_deg": (61.0, 18.15),
            "R_deg_hill": (157.2, 56.7)
            }

# Variable parameters
parameters = {"noact": {},
              "arch": {"parameter_set": "archaea"},
              "glc": {"parameter_set": "activities",
                      "medium": 2},
              "gly": {"parameter_set": "activities",
                      "medium": 1},
              "succ": {"parameter_set": "activities",
                       "medium": 0},
              "LB": {"parameter_set": "activities",
                       "medium": 5},
              "glcAA": {"parameter_set": "activities",
                       "medium": 4},
              "glyAA": {"parameter_set": "activities",
                       "medium": 3}}

prot_fractions = np.arange(0.005, 1, 0.005)  # protein fractions to test
prot_fractions = np.concatenate(([0.00001], prot_fractions, [0.99999]))
growth_rates = np.arange(0.1, 4, 0.001) # growth rates to test

# Simulations
for matrix_type, rates in deg_rates.items():
    results = {"growth_rates": pd.DataFrame(),
               "allocations": pd.DataFrame(),
               "mass_fractions": pd.DataFrame(),
               "fluxes": pd.DataFrame()}

    for name, par in parameters.items():
        par["matrix_type"] = matrix_type
        par["kdeg_max"] = rates[par.get("parameter_set") == "activities"]

        sim = Simulation(par, growth_rates = growth_rates, prot_fractions = prot_fractions)
        sim.test_xrps()

        # get results
        temp_results = {"growth_rates": sim.max_growth_rates, # growth rates (1/h)
                        "allocations": sim.allocations, # ribosome allocations
                        "mass_fractions": sim.mass_fractions, # concentrations (g/g)
                        "fluxes": sim.fluxes}  # fluxes (mmol/gh)

        # add a column with an identifier and concatenate data frames
        for what in results:
            temp_results[what]["name"] = f"{matrix_type}_{name}"
            results[what] = pd.concat([results[what], temp_results[what]])

    for what, df in results.items():
        df.to_csv(f"../data/{matrix_type}_{what}.csv")
