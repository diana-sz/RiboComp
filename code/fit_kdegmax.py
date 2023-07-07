#!/usr/bin/env python
# coding: utf-8

"""
Author: Diana Szeliova

Find optimal RNA degradation rate such that optimal ribosome
composition (xrP) = 36% protein
5 iterations:
* 2 matrices
* without and with activities of RNAP and R
* for the matrix R_deg_hill also activites2 (=RNAP allocation included)
see "general_functions.py", function "make_matrix()")
"""

import numpy as np
from general import Simulation

# Parameters
prot_fractions = np.arange(0.35, 0.37, 0.005)  # protein fractions to test
growth_rates = np.arange(1,3.5,0.001)  # growth rates to test

deg_rates_dic = {
    "activities": np.arange(20, 22, 0.1),
    "hill-2_activities": np.arange(35.5, 37, 0.1),
    "hill-6_activities": np.arange(112, 113, 0.1),
    "hill-6_activities2": np.arange(41, 42, 0.1),
}

for parameter_set, deg_rates in deg_rates_dic.items():
    for k_deg in deg_rates:
        par = {"matrix_type": "extended",
               "kdeg_max": k_deg,
               "parameter_set": parameter_set}
        sim = Simulation(par,
                         growth_rates = growth_rates,
                         prot_fractions = prot_fractions)
        sim.test_xrps(plot = False)

        mu_opt = sim.max_growth_rates.max()["growth_rate"]
        mu_idx = sim.max_growth_rates.idxmax()["growth_rate"]
        xrp_opt = sim.max_growth_rates.loc[mu_idx]["prot_fraction"]

        if abs(xrp_opt - 0.36) <= 0.001:
            to_print = f"(parameter set={parameter_set}): {k_deg:.2f}"
            print(f"Optimal deg. rate {to_print}")
            break
    else:
        print("Optimal deg. rate not found")
