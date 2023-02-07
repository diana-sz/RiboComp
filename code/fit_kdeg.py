#!/usr/bin/env python
# coding: utf-8

"""
Find optimal RNA degradation rate such that optimal ribosome
composition (xrP) = 36% protein
3x2 iterations:
* without and with activities of RNAP and R
* three different types of degradation (matrices R_deg, R_deg2 and R_deg_hill;
see "general_functions.py", function "make_matrix()")
"""

import numpy as np
from general import Simulation

# Parameters
prot_fractions = np.arange(0.35, 0.37, 0.005)  # protein fractions to test
growth_rates = np.arange(1,3.5,0.001)  # growth rates to test
deg_rates_dic = {"R_deg": (np.arange(61, 62, 0.1), np.arange(18, 19, 0.05)),
                 "R_deg2": (np.arange(53, 54, 0.1), np.arange(16, 17, 0.1)),
                 "R_deg_hill": (np.arange(157, 158, 0.1), np.arange(56, 57, 0.1))}

for matrix_type, deg_rates_tuple in deg_rates_dic.items():
    for activities in [True, False]:
        deg_rates = deg_rates_tuple[activities]

        for k_deg in deg_rates:
            par = {"matrix_type": matrix_type,
                   "kdeg_max": k_deg}
            if activities:
                par["parameter_set"] = "activities"

            sim = Simulation(par,
                             growth_rates = growth_rates,
                             prot_fractions = prot_fractions)
            sim.test_xrps(plot = False)
            mu_opt = sim.max_growth_rates.max()["growth_rate"]
            mu_idx = sim.max_growth_rates.idxmax()["growth_rate"]
            xrp_opt = sim.max_growth_rates.loc[mu_idx]["prot_fraction"]

            print(f"Deg. rate: {k_deg:.2f}, xrP: {xrp_opt:.3f}, growth rate: {mu_opt:.5f}")

            if abs(xrp_opt - 0.36) <= 0.001:
                to_print = f"activities={activities}, matrix={matrix_type}: {k_deg:.2f}"
                print(f"Optimal deg. rate ({to_print})")
                break

        else:
            print("Optimal deg. rate not found")
