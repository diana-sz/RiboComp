#!/usr/bin/env python
# coding: utf-8

"""
RBA with fixed allocation of resources like in Kostinski & Reuveni 2020
* iterated over 6 different media
"""

import numpy as np
import pandas as pd
from general import Simulation

# Parameters
MATRIX_TYPE = "Kostinski"

prot_fractions = np.arange(0.01,1,0.005)
growth_rates = np.arange(0.01,2.5,0.001)

parameters = {"glc": {"medium": 2},
              "gly": {"medium": 1},
              "succ": {"medium": 0},
              "LB": {"medium": 5},
              "glcAA": {"medium": 4},
              "glyAA": {"medium": 3}}

growth_rates_res = pd.DataFrame()

for name, par in parameters.items():
    print(f"Medium {name}")

    sim = Simulation(growth_rates = growth_rates,
                     prot_fractions = prot_fractions,
                     matrix_type = MATRIX_TYPE,
                     medium = par["medium"],
                     parameter_set = MATRIX_TYPE)
    sim.test_xrps(plot=False)

    temp_results = sim.max_growth_rates
    temp_results["name"] = f"{MATRIX_TYPE}_{name}"
    growth_rates_res = pd.concat([growth_rates_res, temp_results])
    
    #all_last_mus[medium] = sim.max_growth_rates["growth_rate"]

#res = pd.DataFrame(all_last_mus)
#res.set_index(prot_fractions, inplace=True)
growth_rates_res.to_csv(f"../data/{MATRIX_TYPE}_reproduction_growth_rates.csv")
