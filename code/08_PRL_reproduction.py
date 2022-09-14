#!/usr/bin/env python
# coding: utf-8

"""
RBA with fixed allocation of resources like in Kostinski & Reuveni 2020
* cases -- 6 different media
* phi_4 -- allocation of ribosome to RNAP production
* phi_rP -- allocation of ribosome to rP production
* phi_P -- allocation of RNAP to rRNA production
* f_R -- fraction of active ribosomes
* f_4 -- fraction of active RNAP
"""

import efmtool
import numpy as np
import pandas as pd
from general_functions import make_nice_results, get_const_parameters, bisection_search_mu

k1, k2, k3, k5, c, kel, n1, n2, n3, n4, n5, Mtot, mwG, mwAA, mwNT, cell_mass, avogadro = get_const_parameters()

# row and column names for the matrix
row_names = ["G", "AA", "NT", "rRNA", "rP", "C1", "C2", "C3", "C4", "CA", "CR", "phi1","phi2", "DM"]
rxns = ["vIG", "vEAA", "vENT", "vRNAP", "vrA", "wIG", "wEAA", "wENT", "wRNAP", "wrAS", "wRP", 
             "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "C"]
prot_fractions = np.arange(0.01,1,0.01)  # protein fractions to test
cases = [0,1,2,3,4,5]
mus = np.arange(0.1,2.2,0.001)

all_last_mus = {}
for case in cases:
    print(f"Case {case}")
    
    # ribosome parameters
    phi_4 = [0.93, 1.14, 1.35, 1.5, 1.61, 1.66][case]/100  # phi^RNAP_R, data in %
    phi_rP = [7.8, 9.4, 11.8, 15.3, 19.2, 23.1][case]/100  # phi^rP_R, data in %
    phi_o = 1-phi_4-phi_rP # phi_other
    f_R = 0.85  # fraction active R
    kel = [12, 16.83, 21, 20.17, 21, 22.25][case]*3600*f_R  # translation rate (AA/h)
    
    # RNAP parameters
    f_4 = [13.2, 14.4, 15.0, 18.8, 24.2, 31.0][case]/100  # fraction active P, data in %
    phi_P = [0.18, 0.28, 0.42, 0.52, 0.60, 0.65][case]  # phi^rRNA_RNAP, data as fraction
    c2 = c*f_4*phi_P # transcription rate (nt/h)

    last_mus = []
    for frac in prot_fractions:
        nR = Mtot*frac/mwAA
        mR = Mtot*(1-frac)/mwNT
        nAA = mwAA/mwG
        nNT = (mwNT-mwAA)/mwG

        matrices = {}
        for mu in mus:
            S = np.array([[   1, -nAA, -nNT,   0,         0,     0,     0,     0,     0,     0,     0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                          [   0,  1,  -1,      0,         0,   -n1,   -n2,   -n3,   -n4,   -n5,   -nR,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                          [   0,  0,   1,    -mR,         0,     0,     0,     0,     0,     0,     0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                          [   0,  0,   0,      1,        -1,     0,     0,     0,     0,     0,     0, -1, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                          [   0,  0,   0,      0,        -1,     0,     0,     0,     0,     0,     1,  0,-1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                          [ -mu,  0,   0,      0,         0,    k1,     0,     0,     0,     0,     0,  0, 0, -1,  0,  0,  0,  0,  0,  0,  0,  0],
                          [   0,-mu,   0,      0,         0,     0,    k2,     0,     0,     0,     0,  0, 0,  0, -1,  0,  0,  0,  0,  0,  0,  0],
                          [   0,  0, -mu,      0,         0,     0,     0,    k3,     0,     0,     0,  0, 0,  0,  0, -1,  0,  0,  0,  0,  0,  0],
                          [   0,  0,   0, -mu*mR,         0,     0,     0,     0,    c2,     0,     0,  0, 0,  0,  0,  0, -1,  0,  0,  0,  0,  0],
                          [   0,  0,   0,      0,       -mu,     0,     0,     0,     0,    k5,     0,  0, 0,  0,  0,  0,  0, -1,  0,  0,  0,  0],
                          [   0,  0,   0,      0, phi_o*kel,-mu*n1,-mu*n2,-mu*n3,     0,-mu*n5,     0,  0, 0,  0,  0,  0,  0,  0, -1,  0,  0,  0],
                          [   0,  0,   0,      0, phi_4*kel,      0,    0,     0,-mu*n4,     0,     0,  0, 0,  0,  0,  0,  0,  0,  0, -1,  0,  0],
                          [   0,  0,   0,      0,phi_rP*kel,      0,    0,     0,     0,     0,-mu*nR,  0, 0,  0,  0,  0,  0,  0,  0,  0, -1,  0],
                          [-mwG,  0,   0,      0,         0,      0,    0,     0,     0,     0,     0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0, mu]])
            matrices[str(round(mu,4))] = S
 
        # reversibilities - all irreversible
        revs = [0] * S.shape[1]

        last_mu = bisection_search_mu(matrices, revs, rxns, row_names)
        last_mus.append(last_mu)
    all_last_mus[case] = last_mus

res = pd.DataFrame(all_last_mus, index = prot_fractions)
res.to_csv("../data/08_PRL_reproduction.csv")
