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
from general_functions import run_simulations

# Constant parameters
MATRIX_TYPE = "RNAPmax"

# Variable parameters
parameters = {"noact": {"activities": False,
                        "rnap_max": 0.0000126,
                        "organism": "ecoli"},
              "act": {"activities": True,
                      "rnap_max": 0.0000776,
                      "organism": "ecoli"},
             "arch": {"activities": False,
                      "rnap_max": 0.0000126,
                      "organism": "archaea"}}

prot_fractions = np.arange(0.005, 1, 0.005)  # protein fractions to test
prot_fractions = np.concatenate(([0.00001], prot_fractions, [0.99999]))
mus = np.arange(0.5,4,0.001) # growth rates to test

# Simulations
for name, par in parameters.items():
    res, egvs, conc, phis = run_simulations(prot_fractions,
                                            mus,
                                            MATRIX_TYPE,
                                            organism=par["organism"],
                                            activities=par["activities"],
                                            rnap_max=par["rnap_max"])

    res.to_csv(f"../data/{MATRIX_TYPE}_{name}_mus.csv")
    egvs.to_csv(f"../data/{MATRIX_TYPE}_{name}_fluxes.csv")
    conc.to_csv(f"../data/{MATRIX_TYPE}_{name}_concentrations.csv")
    phis.to_csv(f"../data/{MATRIX_TYPE}_{name}_phis.csv")
