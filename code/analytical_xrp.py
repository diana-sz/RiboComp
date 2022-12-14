#!/usr/bin/env python
# coding: utf-8

"""
Calculate optimal protein content in ribosome analytically
"""

import numpy as np
import pandas as pd
from general_functions import get_parameters

def calculate_xrp(kel, c, activities):
    """
    Calculate optimal protein content in ribosome analytically
    """

    k1, k2, k3, k5, _, _, nig, neaa, nent, nrnap, naf, mwR, mwG, mwAA, mwNT = get_parameters("ecoli",
                                                                                             False,
                                                                                             False)
    beta = kel*((nig/k1*mwAA/mwG) + (neaa/k2 + nent/k3 + nrnap/c)*mwAA/mwNT)*mwAA/mwR
    gamma = kel*(nig/k1*mwAA/mwG + neaa/k2)*mwAA/mwR
    
    cpolmax = 0.0000126
    if activities:
        cpolmax = 0.0000776 #0.000079 #

    a = (kel/c)*(1/(cpolmax*mwNT))*(mwAA/mwR)
    zeta = (a-beta)/(1+gamma)
    xrp = 1-a/((1+zeta)*zeta)
    
    return xrp

kels_r = np.arange(1, 90000, 200)
kels_rnap = np.arange(1, 90000, 200)
comb_array = np.array(np.meshgrid(kels_r, kels_rnap)).T.reshape(-1,2)

xrps = []
for row in comb_array:
    xrp_opt = calculate_xrp(row[0], row[1], True)
    xrps.append(xrp_opt)
res = pd.DataFrame({"kel_r": comb_array[:,0],
                    "kel_rnap": comb_array[:,1],
                    "xrp": xrps})
res.to_csv("../data/analytical_xrp_act.csv")
