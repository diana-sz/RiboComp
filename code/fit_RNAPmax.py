#!/usr/bin/env python
# coding: utf-8

"""
Find optimal RNAPmax such that optimal ribosome composition (xrP) = 36% protein
2 iterations: without and with activities of RNAP and R
"""

import numpy as np
from general_functions import run_simulations

# Parameters
FILENAME = "fit_RNAP_max"
MATRIX_TYPE = "RNAPmax"
ORGANISM = "ecoli"

rnap_maxes = np.arange(0.2e-05, 9e-05, 1e-7)
prot_fractions = np.arange(0.34, 0.38, 0.005)  # protein fractions to test
mus = np.arange(1,4,0.001)  # growth rates to test

for activities in [False, True]:
    # bisection search initiation
    left = 0  # The starting index of the list we have to search in
    right = len(rnap_maxes)-1  # the last index of the list we have to search in
    mid = (right + left)//2  # // means floored division
    xrp_opt = -1
    iterations = 0

    while abs(xrp_opt - 0.36) >= 0.001:
        iterations += 1
        new_rnap = rnap_maxes[mid]
        res, egvs, conc, phis = run_simulations(prot_fractions,
                                                mus,
                                                MATRIX_TYPE,
                                                ORGANISM,
                                                activities = activities,
                                                rnap_max = new_rnap,
                                                plot = False)
        mu_opt = res.max()["mu"]
        mu_idx = res.idxmax()["mu"]
        xrp_opt = res.loc[mu_idx]["x"]

        print(new_rnap, xrp_opt, mu_opt)

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
    print(f"Optimal RNAPmax (activities={activities}): {new_rnap:.7f} ({n_rnap} molecules/cell), max. mu: {mu_opt}")
