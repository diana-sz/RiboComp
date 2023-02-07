#!/usr/bin/env python
# coding: utf-8

"""
Find optimal RNAPmax such that optimal ribosome composition (xrP) = 36% protein
2 iterations: without and with activities of RNAP and R
"""

import numpy as np
from general import Simulation

# Parameters
MATRIX_TYPE = "RNAPmax"

rnap_maxes = np.arange(0.2e-05, 9e-05, 1e-7)
prot_fractions = np.arange(0.35, 0.37, 0.005)  # protein fractions to test
growth_rates = np.arange(1,4,0.001)  # growth rates to test

for activities in [False, True]:
    # bisection search initiation
    left = 0  # the starting index of the list we have to search in
    right = len(rnap_maxes)-1  # the last index of the list we have to search in
    mid = (right + left)//2  # floored division
    xrp_opt = -1
    iterations = 0

    while abs(xrp_opt - 0.36) >= 0.001:
        iterations += 1
        new_rnap = rnap_maxes[mid]

        sim = Simulation(growth_rates = growth_rates,
                         prot_fractions = prot_fractions,
                         matrix_type = MATRIX_TYPE,
                         rnapmax = new_rnap,
                         parameter_set = ["default", "activities"][activities])
        sim.test_xrps(plot = False)
        res = sim.max_growth_rates
        mu_opt = res.max()["growth_rate"]
        mu_idx = res.idxmax()["growth_rate"]
        xrp_opt = res.loc[mu_idx]["prot_fraction"]

        print(f"RNAP: {new_rnap:.7f}, xrP: {xrp_opt:.3f}, growth rate: {mu_opt:.5f}")


        # RNAP_max too high
        if xrp_opt < 0.36:
            right = mid - 1
        # RNAP_max too low
        else:
            left = mid + 1
        mid = (right + left)//2

        if iterations > len(rnap_maxes):
            print("Optimal RNAPmax not found")
            break

    CELL_MASS = 8.7e-13  # cell mass at max. growth rate [g/cell]
    n_rnap = round(new_rnap*6.022e20*CELL_MASS)  # converted to molecules per cell
    print(f"RNAPmax (activities={activities}): {new_rnap:.7f} ({n_rnap} RNAP/cell), mu: {mu_opt:.5f}")
