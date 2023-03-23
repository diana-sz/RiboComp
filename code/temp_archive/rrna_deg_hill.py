#!/usr/bin/env python
# coding: utf-8

"""
Run RBA with R degradation rate
    * three types of matrices with different degradation function: R_deg, R_deg2, R_deg_hill
    * for each matrix, 8 parameter sets
    * archaea ("arch") / ecoli ("noact") (for activities=False)
    * 6 media ("glc", "gly", "succ", "LB", "glcAA", "glyAA") (for activities=True)
values in "deg_rates" optimized with "fit_deg_rate.py"

Save optimal growth rate and in case of R_deg_hill also fluxes,
concentrations and ribosome allocations
"""

import numpy as np
import pandas as pd
from general import Simulation

prot_fractions = np.arange(0.005, 1, 0.05)  # protein fractions in ribosome
# add numbers close to 0 and 1 (otherwise no solution, the matrix would have to be changed)
prot_fractions = np.concatenate(([0.00001], prot_fractions, [0.99999])) 
growth_rates = np.arange(0.1, 4, 0.01)

parameters = pd.read_csv("../data/parameters_deg.csv")

for index, row in parameters.iterrows():
    name = row.name
    par = dict(row)
    _ = par.pop("name")

# # fitted degradation rates
# matrix_types = {"R_deg_hill-6": {
#                     "no_activities": 619.4
#                     "activities": 112.5 },
#                 "R_deg_hill-2": {
#                     "no_activities": 210.4,
#                     "activities": 36.0},
#                 "R_deg": {
#                     "no_acttivities": 125.3
#                     "activities": 21.0}}

# for matrix_type, deg_rates in matrix_types.items():
#     parameters = {
#         "ecoli": {
#             "kdeg_max": deg_rates["no_activities"]},
#         "arch": {
#             "parameter_set": "archaea", 
#             "kdeg_max": deg_rates["no_activities"]},
#         "arch2": {
#             "kdeg_max": deg_rates["no_activities"]*2}, # ecoli parameters but higher kdeg
#         "mito": {
#             "kdeg_max": deg_rates["no_activities"]},
#         "glc": {
#             "parameter_set": "activities",
#             "medium": 2,
#             "kdeg_max": deg_rates["activities"]},
#         "gly": {
#             "parameter_set": "activities",
#             "medium": 1,
#             "kdeg_max": deg_rates["activities"]},
#         "succ": {
#             "parameter_set": "activities",
#             "medium": 0,
#             "kdeg_max": deg_rates["activities"]},
#         "LB": {
#             "parameter_set": "activities",
#             "medium": 5,
#             "kdeg_max": deg_rates["activities"]},
#         "glcAA": {
#             "parameter_set": "activities",
#             "medium": 4,
#             "kdeg_max": deg_rates["activities"]},
#         "glyAA": {
#             "parameter_set": "activities",
#             "medium": 3,
#             "kdeg_max": deg_rates["activities"]},
#     }

    results = {
        "growth_rates": pd.DataFrame(),
        "allocations": pd.DataFrame(),
        "mass_fractions": pd.DataFrame(),
        "fluxes": pd.DataFrame()
    }

    for name, par in parameters.items():
        par["matrix_type"] = matrix_type
        
        # mitochondria has a specific matrix
        if name == "mito":
            par["matrix_type"] = matrix_type + "_mito"

        sim = Simulation(par, 
                         growth_rates = growth_rates, 
                         prot_fractions = prot_fractions)
        sim.test_xrps()

        temp_results = {
            "growth_rates": sim.max_growth_rates, # 1/h
            "allocations": sim.allocations, # ribosome allocations
            "mass_fractions": sim.mass_fractions, # g/g
            "fluxes": sim.fluxes # mmol/gh
        }

        # add a column with an identifier and concatenate data frames
        for result_type in results:
            temp_results[result_type]["name"] = f"{matrix_type}_{name}"
            results[result_type] = pd.concat([results[result_type], temp_results[result_type]])

    for result_type, df in results.items():
        df.to_csv(f"../data/{matrix_type}_{result_type}.csv")
