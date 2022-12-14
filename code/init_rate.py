#!/usr/bin/env python
# coding: utf-8

"""
Run RBA with a limit on RNAP flux (vmax)
"""

import numpy as np
import pandas as pd
from general_functions import run_simulations

# Constant parameters
MATRIX_TYPE = "init_rate"
ORGANISM = "ecoli"

# Variable parameters
kins = [6600, 5, 0.5, 0.05]  # limit on transcriptio inititation rate
prot_fractions = np.arange(0.005, 1, 0.005)  # protein fractions to test
prot_fractions = np.concatenate(([0.00001], prot_fractions, [0.99999]))
mus = np.arange(0.01,4.5,0.001)  # growth rates to test

all_results = pd.DataFrame(columns=["x", "mu", "kin"])
for kin in kins:
    res, egvs, conc, phis = run_simulations(prot_fractions,
                                            mus,
                                            MATRIX_TYPE,
                                            ORGANISM,
                                            kin=kin)
    res_df = pd.DataFrame(res)
    res_df["kin"] = kin
    all_results = pd.concat([all_results, res_df])

all_results.to_csv(f"../data/{MATRIX_TYPE}_mus.csv")
