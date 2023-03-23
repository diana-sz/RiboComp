#!/usr/bin/env python
# coding: utf-8

"""
Run RBA with a limit on total RNAP concentration
Three parameters sets:
    * noact: max. activities of R and RNAP NOT included (activities = False)
    * act: max. activities of R and RNAP included (activities = True)
    * arch: archaeal parameters used, activities = False
"""

import numpy as np
import pandas as pd
from general import Simulation

# Constant parameters
MATRIX_TYPE = "RNAPmax"
RNAPMAX_NOACT = 6.8e-06
RNAPMAX_ACT = 4.46e-05

# Variable parameters
parameters = {"noact": {"parameter_set": "default",
                        "rnapmax": RNAPMAX_NOACT},
              "arch": {"parameter_set": "archaea",
                      "rnapmax": RNAPMAX_NOACT},
              "glc": {"parameter_set": "activities",
                      "rnapmax": RNAPMAX_ACT,
                      "medium": 2},
              "gly": {"parameter_set": "activities",
                      "rnapmax": RNAPMAX_ACT,
                      "medium": 1},
              "succ": {"parameter_set": "activities",
                       "rnapmax": RNAPMAX_ACT,
                       "medium": 0},
              "LB": {"parameter_set": "activities",
                     "rnapmax": RNAPMAX_ACT,
                     "medium": 5},
              "glcAA": {"parameter_set": "activities",
                        "rnapmax": RNAPMAX_ACT,
                       "medium": 4},
              "glyAA": {"parameter_set": "activities",
                        "rnapmax": RNAPMAX_ACT,
                       "medium": 3}}

prot_fractions = np.arange(0.005, 1, 0.005)  # protein fractions to test
prot_fractions = np.concatenate(([0.00001], prot_fractions, [0.99999]))
growth_rates = np.arange(0.5,4,0.001) # growth rates to test

# Simulations
results = {"growth_rates": pd.DataFrame(),
           #"allocations": pd.DataFrame(),
           #"mass_fractions": pd.DataFrame(),
           #"fluxes": pd.DataFrame()
          }

for name, par in parameters.items():
    par["matrix_type"] = MATRIX_TYPE
    sim = Simulation(par, growth_rates = growth_rates, prot_fractions = prot_fractions)
    sim.test_xrps()

    # get results
    temp_results = {"growth_rates": sim.max_growth_rates,
                    #"allocations": sim.allocations,
                    #"mass_fractions": sim.mass_fractions,
                    #"fluxes": sim.fluxes
                   }

    # add a column with an identifier and concatenate data frames
    for what in results:
        temp_results[what]["name"] = f"{MATRIX_TYPE}_{name}"
        results[what] = pd.concat([results[what], temp_results[what]])

for what, df in results.items():
    df.to_csv(f"../data/{MATRIX_TYPE}_{what}.csv")
