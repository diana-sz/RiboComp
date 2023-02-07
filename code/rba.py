#!/usr/bin/env python
# coding: utf-8

"""
Run RBA
2 parameters sets:
    * standard: realistic E. coli parameters (rna_expensive = False)
    * reverse: parameters that make RNA moRe expensive (rna_expensive = True)
"""

import numpy as np
import pandas as pd
from general import Simulation

MATRIX_TYPE = "RBA"

# Variable parameters
parameters = {"standard": {},  # default pars
              "reverse": {"parameter_set": "rna_expensive"}}

prot_fractions = np.arange(0.005, 1, 0.005)  # protein fractions to test
prot_fractions = np.concatenate(([0.00001], prot_fractions, [0.99999]))
growth_rates = np.arange(0.1,4,0.001) # growth rates to test

# Simulations
results = {"growth_rates": pd.DataFrame(),  
           "allocations": pd.DataFrame(),
           "mass_fractions": pd.DataFrame(),
           "fluxes": pd.DataFrame()}

for name, par in parameters.items():
    sim = Simulation(par, growth_rates = growth_rates, prot_fractions = prot_fractions)
    sim.test_xrps()
    
    # get results
    temp_results = {"growth_rates": sim.max_growth_rates, # growth rates (1/h)
                    "allocations": sim.allocations, # ribosome allocations
                    "mass_fractions": sim.mass_fractions, # mass fractions
                    "fluxes": sim.fluxes}  # fluxes (mmol/gh)

    # add a column with an identifier and concatenate data frames
    for what in results:
        temp_results[what]["name"] = f"{MATRIX_TYPE}_{name}"
        results[what] = pd.concat([results[what], temp_results[what]])

for what, df in results.items():
    df.to_csv(f"../data/{MATRIX_TYPE}_{what}.csv")
        