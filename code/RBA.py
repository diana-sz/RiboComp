#!/usr/bin/env python
# coding: utf-8

"""
Run RBA
2 parameters sets:
    * standard: realistic E. coli parameters (rna_expensive = False)
    * reverse: parameters that make RNA moRe expensive (rna_expensive = True)
"""

import numpy as np
from general_functions import run_simulations

# Parameters
MATRIX_TYPE = "RBA"
ORGANISM = "ecoli"

# Variable parameters
parameters = {"standard": False, "reverse": True}
prot_fractions = np.arange(0.005, 1, 0.005)  # protein fractions to test
prot_fractions = np.concatenate(([0.00001], prot_fractions, [0.99999]))
mus = np.arange(2,5,0.001) # growth rates to test

# Simulations
for name, par in parameters.items():
    res, egvs, conc, phis = run_simulations(prot_fractions,
                                            mus,
                                            MATRIX_TYPE,
                                            ORGANISM,
                                            rna_expensive=par)

    res.to_csv(f"../data/{MATRIX_TYPE}_{name}_mus.csv")
    egvs.to_csv(f"../data/{MATRIX_TYPE}_{name}_fluxes.csv")
    conc.to_csv(f"../data/{MATRIX_TYPE}_{name}_concentrations.csv")
    phis.to_csv(f"../data/{MATRIX_TYPE}_{name}_phis.csv")
