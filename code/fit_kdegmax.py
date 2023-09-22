#!/usr/bin/env python
# coding: utf-8

"""
Author: Diana Szeliova

Find optimal RNA degradation rate (kdeg_max) such that optimal ribosome
composition (xrP) = 36% protein for different degradation functions
always done for 'glc' condition (glucose minimal medium)
"""

import numpy as np
import pandas as pd
from model_functions import Simulation

parameters = pd.read_csv("../data/parameters.csv")
prot_fractions = np.arange(0.35, 0.37, 0.005)  # protein fractions to test
growth_rates = np.arange(1,3.5,0.001)  # growth rates to test

deg_rates_dic = {
    "activities": np.arange(19, 20, 0.1),
    "hill-2_activities": np.arange(32, 33, 0.1),
    "hill-6_activities": np.arange(100.3, 101, 0.1),
    "hill-6_activities2": np.arange(37, 42, 0.1),
    "hill-6_activities_var_kel": np.arange(112, 113, 0.1)
}

for index, row in parameters.iterrows():
    par = dict(row)

    # skip models without degradation
    if par["matrix_type"] in ["Kostinski", "base"]:
        continue
    # skip everything that is not glucose minimal medium
    if row["name"] != "glc":
        continue

    for k_deg in deg_rates_dic[par["parameter_set"]]:
        par["kdeg_max"] = k_deg
        sim = Simulation(par,
                         growth_rates = growth_rates,
                         prot_fractions = prot_fractions)
        sim.test_xrps(plot = False)

        mu_opt = sim.max_growth_rates.max()["growth_rate"]
        mu_idx = sim.max_growth_rates.idxmax()["growth_rate"]
        xrp_opt = sim.max_growth_rates.loc[mu_idx]["prot_fraction"]

        if abs(xrp_opt - 0.36) <= 0.001:
            print(f"{par['parameter_set']}: optimal deg. rate {k_deg:.2f}")
            break
    else:
        print("Optimal deg. rate not found")
