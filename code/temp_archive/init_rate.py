#!/usr/bin/env python
# coding: utf-8

"""
Run RBA with a limit on RNAP initiation rate instead of elongation
"""

import numpy as np
import pandas as pd
from general import Simulation

# Parameters
MATRIX_TYPE = "RBA_init_rate"

kins = [6600, 5, 0.5, 0.05]  # limit on transcription inititation rate
prot_fractions = np.arange(0.005, 1, 0.005)  # protein fractions to test
prot_fractions = np.concatenate(([0.00001], prot_fractions, [0.99999]))
growth_rates = np.arange(0.01,4.5,0.001)

all_results = pd.DataFrame()
for kin_max in kins:
    sim = Simulation(growth_rates = growth_rates, 
                     prot_fractions = prot_fractions,
                     matrix_type = MATRIX_TYPE,
                     kin_max = kin_max)
    sim.test_xrps()

    res_df = pd.DataFrame(sim.max_growth_rates)
    res_df["kin"] = kin_max
    res_df["name"] = MATRIX_TYPE
    all_results = pd.concat([all_results, res_df])

all_results.to_csv(f"../data/{MATRIX_TYPE}_growth_rates.csv")
