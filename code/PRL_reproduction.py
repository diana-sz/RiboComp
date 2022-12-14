#!/usr/bin/env python
# coding: utf-8

"""
RBA with fixed allocation of resources like in Kostinski & Reuveni 2020
* iterated over 6 different media
* phi_4 -- allocation of ribosome to RNAP production
* phi_rP -- allocation of ribosome to rP production
* phi_P -- allocation of RNAP to rRNA production
* f_R -- fraction of active ribosomes
* f_4 -- fraction of active RNAP
"""

import numpy as np
import pandas as pd
from general_functions import get_parameters, bisection_search_mu

# Parameters
FILENAME = "PRL_reproduction"
ORGANISM = "ecoli"
ACTIVITIES = True
RNA_EXPENSIVE = False

prot_fractions = np.arange(0.01,1,0.005)  # protein fractions to test
mus = np.arange(0.1,2.2,0.001)

# row and column names for the matrix
row_names = ["G", "AA", "NT", "rRNA", "rP",
             "C1", "C2", "C3", "C4", "CA", "CR", "phi1","phi2", "DM"]
rxns = ["vIG", "vEAA", "vENT", "vRNAP", "vAF",
        "wIG", "wEAA", "wENT", "wRNAP", "wAF", "wrP",
        "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "C"]

all_last_mus = {}
for case in range(6):
    print(f"Case {case}")

    k1, k2, k3, k5, c, kel, nig, neaa, nent, nrnap, naf, mwR, mwG, mwAA, mwNT = get_parameters(ORGANISM,
                                                                                               ACTIVITIES,
                                                                                               medium = case)
    # ribosome parameters - kel has to be split into several rows
    # so we cannot use the kel from get_parameters
    phi_4 = [0.93, 1.14, 1.35, 1.5, 1.61, 1.66][case]/100  # phi^RNAP_R, data in %
    phi_rP = [7.8, 9.4, 11.8, 15.3, 19.2, 23.1][case]/100  # phi^rP_R, data in %
    phi_o = 1-phi_4-phi_rP # phi_other

    # RNAP allocation fractions
    phi_P = [0.18, 0.28, 0.42, 0.52, 0.60, 0.65][case]  # phi^rRNA_RNAP, data as fraction
    c = c*phi_P # transcription rate (nt/h)

    last_mus = []
    for frac in prot_fractions:
        nrp = mwR*frac/mwAA
        nrrna = mwR*(1-frac)/mwNT
        naa = mwAA/mwG
        nnt = (mwNT-mwAA)/mwG

        matrices = {}
        for mu in mus:
            S = np.array([[   1, -naa, -nnt,   0,         0,     0,     0,     0,     0,     0,     0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                          [   0,  1,  -1,      0,         0,   -nig,  -neaa,   -nent,   -nrnap,   -naf,   -nrp,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                          [   0,  0,   1,    -nrrna,         0,     0,     0,     0,     0,     0,     0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                          [   0,  0,   0,      1,        -1,     0,     0,     0,     0,     0,     0, -1, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                          [   0,  0,   0,      0,        -1,     0,     0,     0,     0,     0,     1,  0,-1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                          [ -mu,  0,   0,      0,         0,    k1,     0,     0,     0,     0,     0,  0, 0, -1,  0,  0,  0,  0,  0,  0,  0,  0],
                          [   0,-mu,   0,      0,         0,     0,    k2,     0,     0,     0,     0,  0, 0,  0, -1,  0,  0,  0,  0,  0,  0,  0],
                          [   0,  0, -mu,      0,         0,     0,     0,    k3,     0,     0,     0,  0, 0,  0,  0, -1,  0,  0,  0,  0,  0,  0],
                          [   0,  0,   0, -mu*nrrna,         0,     0,     0,     0,   c,     0,     0,  0, 0,  0,  0,  0, -1,  0,  0,  0,  0,  0],
                          [   0,  0,   0,      0,       -mu,     0,     0,     0,     0,    k5,     0,  0, 0,  0,  0,  0,  0, -1,  0,  0,  0,  0],
                          [   0,  0,   0,      0, phi_o*kel,-mu*nig,-mu*neaa,-mu*nent,     0,-mu*naf,     0,  0, 0,  0,  0,  0,  0,  0, -1,  0,  0,  0],
                          [   0,  0,   0,      0, phi_4*kel,      0,    0,     0,-mu*nrnap,     0,     0,  0, 0,  0,  0,  0,  0,  0,  0, -1,  0,  0],
                          [   0,  0,   0,      0,phi_rP*kel,      0,    0,     0,     0,     0,-mu*nrp,  0, 0,  0,  0,  0,  0,  0,  0,  0, -1,  0],
                          [-mwG,  0,   0,      0,         0,      0,    0,     0,     0,     0,     0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0, mu]])
            matrices[str(round(mu,4))] = S

        # reversibilities - all irreversible
        revs = [0] * S.shape[1]

        last_mus.append(bisection_search_mu(matrices, revs, rxns, row_names))
    all_last_mus[case] = last_mus

res = pd.DataFrame(all_last_mus, index = prot_fractions)
res.to_csv(f"../data/{FILENAME}.csv")
