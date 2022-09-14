import math
import numpy as np
import pandas as pd
import efmtool

def make_nice_results(egms, rxns, drop_slack = True):
    """Convert output of efmtool to a data frame
    normalize by the C (constraint) column

    Parameters: 
        egms: array with the output from efmtool
        rxns: list of reactions
        if drop_slack = True, drop columns with slack variables and constraints (named S* or C)
    
    Returns: 
        data frame with column and row names
    """
    
    # solve "ValueError: Big-endian buffer not supported on little-endian compiler"
    egms = egms.byteswap().newbyteorder()
    
    # remove zero flux EGVs
    zeros = (np.count_nonzero(egms, axis=0)-1).astype(bool)  # subtract one because zero EGVs have 1 for a slack variable
    egms = egms.T[zeros]

    res = pd.DataFrame(egms, index = ["EGM%s" % (i+1) for i in range(egms.shape[0])], columns = rxns)

    # normalize values by the column C
    if "C" in res.columns:
        res = res.div(res.C, axis = 0)
        
    if drop_slack:
        cols = [c for c in res.columns if c.startswith("S")]
        if "C" in res.columns:
            cols = cols+["C"]
        res = res[res.columns.drop(cols)]
    return res


# def mu_RBA(alpha, beta, gamma, x):
#     mu = 2/(alpha + x + math.sqrt((alpha+1)**2 - (beta-1)**2 + (beta-x)**2 + 4*gamma))
#     return mu


# def mu_RNApol(delta, x):
#     mu = delta*(1+math.sqrt(1+(2/(delta*(1-x)))))
#     return mu


def egv_dict_to_df(egv_dict, rxns, frac):
    """
    Convert dictionary of data frames to one big data frame

    Parameters:
        egv_dict: dictionary of data frames with EGVs at different mus 
        (EGV data frames are output of make_nice_results())
        rxns: list of reaction names which will be used as column names
        frac: protein fraction of ribosome, will be added as a new column
        
    Returns:
        data frame with all results
    
    """
    all_results = pd.DataFrame(columns=rxns)
    for mu in egv_dict:
        res = egv_dict[mu]
        res["mu"] = mu
        res["x"] = str(frac)
        res["EGVs"] = res.index
        all_results = pd.concat([all_results, res])
    all_results = all_results.set_index([all_results["mu"]+all_results.index])
    return all_results



def get_concentrations(all_results, original_matrices, mets_all, Mtot, mwG, n1, n2, n3, n4, n5):
    """Calculate metabolite concentrations"""
    all_concentrations = pd.DataFrame(columns = mets_all, index = all_results.index)
    for egv in all_results.index:
        frac = all_results.loc[egv, "x"]
        mat = original_matrices[frac]
        nR = float(frac)*Mtot
        mR = (Mtot-nR)/2 # NT have twice the mass of AA
        mw = [mwG, mwG, mwG*2, mR*mwG*2, nR*mwG, nR*mwG+mR*mwG*2, 
              n1*mwG, n2*mwG, n3*mwG, n4*mwG, n5*mwG]
        mu = all_results.loc[egv, "mu"]
        all_concentrations.loc[egv] = mat.multiply(all_results.loc[egv][:11]).sum(axis=1)/float(mu)
        all_concentrations.loc[egv] = all_concentrations.loc[egv]*mw
    all_concentrations["mu"] = all_results.mu
    all_concentrations["EGVs"] = all_results.EGVs
    all_concentrations["x"] = all_results.x
    
    return all_concentrations



def get_const_parameters():
    """get constant parameters for stoich. matrix with constraints"""
    to_h = 3600  # conversion factor from seconds to hours

    # kinetic parameters
    k1 = 180*to_h  # glucose transport rate, 1/h
    k2 = 10*to_h  # kcat for AA synthesis, 1/h
    k3 = 10*to_h  # kcat for NT synthesis, 1/h
    k5 = 1/120*to_h  # kcat for ribisome assembly, 1/h
    c = 85*to_h # transcription elongation, nt/s 
    kel = 21*to_h
    
    # stoichiometry
    n1 = 646  # AA in glucose transporter
    avg_protein = 325  # AA in a typical protein, 108986
    n2 = avg_protein*15  # estimated 15 enzymes, 9 steps in glycolysis, 5 steps from Pyr to Glu, 6 to Gln, other AAs more or less
    n3 = avg_protein*15  # estimated 15 enzymes, 6 PPP (glc -> prpp) + 5 pyrimidine/13 purine
    n4 = 3498  # AA in RNAP without sigma factor
    n5 = avg_protein*12  # https://dx.doi.org/10.1128%2FMMBR.00013-07 https://doi.org/10.3389/fmicb.2019.02982
    
    # molecular masses (g/mmol)
    Mtot = 2300  # ribosome
    mwG = 0.18  # glucose
    mwAA = 0.109  # avergage AA
    mwNT= 0.3243  # average NT
    
    avogadro = 6.022*10**20  # molecules/mmol
    cell_mass = 8.7*10**(-13)  # g/cell

    return((k1, k2, k3, k5, c, kel, n1, n2, n3, n4, n5, Mtot, mwG, mwAA, mwNT, cell_mass, avogadro))



def get_free_pol_parameters():
    """get parameters for free polymerase calculation"""
    avogadro = 6.022*10**20  # molecules/mmol
    cell_mass = 8.7*10**(-13)  # g/cell
    vmax = 110  # initiations/min/promoter
    vmax = vmax*60  # 1/h
    pro = 36  # promoters/cell
    pro = (pro/avogadro)/cell_mass  # mmol/g
    Km = 132 #1504  # RNAP molecules/cell
    Km = (Km/avogadro)/cell_mass  # RNAP Michaelis constant, mmol/g
    Kns = 4.7*10**(-3)  # non-speicifc RNAP dissociation constant, mmol/g
    n_sites = 17480000/cell_mass/avogadro  # nonspecific binding sites on DNA, mmol/g
    
    return((vmax, pro, Km, Kns, n_sites, cell_mass, avogadro))



def make_matrix_v4max(mwG, mwAA, mwNT, n1, n2, n3, n4, n5, nR, mR, mu, k1, k2, k3, k5, c, kel, v4max, c_pol_max):
    """create RBA matrix with additional vmax constraint on v4
    Parameters: 
    * mwG, mwAA, mwNT - molecular masses of glucose (or other substrate), AA and NT
    * n1, n2, n3, n4, n5 - stoichiometric coefficients of transporter (n1), AA/NT synthesis enzymes (n2/n3),
    RNAP (n4) and ribosome assembly factors (n5)
    * nR, mR - stoich coefficients of AA and NT in ribosome
    * mu - growth rate [1/h]
    * k1 - kcat of substrate importer [1/h]
    * k2/k3/k5 - kcat of AA synthesis (k2) / NT synthesis (k3) / ribosome assembly (k5) [1/h]
    * c - transcription rate [nt/h]
    * kel - translation rate [aa/h]
    * v4max -  max RNA synthesis rate [1/h]
    * c_pol_max - max. RNA polymerase [mmol/g]

    Returns:
    * numpy matrix
    """
    
    row_names = ["G", "AA", "NT", "rRNA", "rP", "C1", "C2", "C3", "C4", "CA", "CR", "v4max", "c_pol_max", "DM"]
    rxns = ["vIG", "vEAA", "vENT", "vRNAP", "vrA", "wIG", "wEAA", "wENT", "wRNAP", "wrAS", "wRP", 
             "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "C"]
    nAA = mwAA/mwG
    nNT = (mwNT-mwAA)/mwG  
    
    S = np.array([[  1, -nAA, -nNT,   0,      0,   0,   0,   0,   0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,   0],
                  [  0,  1,  -1,      0,      0, -n1, -n2, -n3, -n4, -n5, -nR,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,   0],
                  [  0,  0,   1,    -mR,      0,   0,   0,   0,   0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,   0],
                  [  0,  0,   0,      1,     -1,   0,   0,   0,   0,   0,   0, -1, 0,  0,  0,  0,  0,  0,  0,  0,  0,   0],
                  [  0,  0,   0,      0,     -1,   0,   0,   0,   0,   0,   1,  0,-1,  0,  0,  0,  0,  0,  0,  0,  0,   0],
                  [-mu,  0,   0,      0,      0,  k1,   0,   0,   0,   0,   0,  0, 0, -1,  0,  0,  0,  0,  0,  0,  0,   0],
                  [  0,-mu,   0,      0,      0,   0,  k2,   0,   0,   0,   0,  0, 0,  0, -1,  0,  0,  0,  0,  0,  0,   0],
                  [  0,  0, -mu,      0,      0,   0,   0,  k3,   0,   0,   0,  0, 0,  0,  0, -1,  0,  0,  0,  0,  0,   0],
                  [  0,  0,   0, -mu*mR,      0,   0,   0,   0,   c,   0,   0,  0, 0,  0,  0,  0, -1,  0,  0,  0,  0,   0],
                  [  0,  0,   0,      0,    -mu,   0,   0,   0,   0,  k5,   0,  0, 0,  0,  0,  0,  0, -1,  0,  0,  0,   0],
                  [  0,  0,   0,      0, kel/mu, -n1, -n2, -n3, -n4, -n5, -nR,  0, 0,  0,  0,  0,  0,  0, -1,  0,  0,   0],
                  [  0,  0,   0,    -mu,      0,   0,   0,   0, v4max, 0,   0,  0, 0,  0,  0,  0, 0,  0,  0,  -1, 0,   0],
                  [  0,  0,   0,      0,      0,   0,   0,   0, -1/mu, 0,   0,  0, 0,  0,  0,  0,  0,  0,  0,   0, -1, c_pol_max],
                  [-mwG,  0,   0,      0,      0,   0,   0,   0,   0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  mu]])
    
    revs = [0] * S.shape[1]

    return((S, rxns, row_names, revs))


def make_matrix_cpolmax(mwG, mwAA, mwNT, n1, n2, n3, n4, n5, nR, mR, mu, k1, k2, k3, k5, c, kel, c_pol_max):
    """create RBA matrix with additional vmax constraint on v4
    Parameters: 
    * mwG, mwAA, mwNT - molecular masses of glucose (or other substrate), AA and NT
    * n1, n2, n3, n4, n5 - stoichiometric coefficients of transporter (n1), AA/NT synthesis enzymes (n2/n3),
    RNAP (n4) and ribosome assembly factors (n5)
    * nR, mR - stoich coefficients of AA and NT in ribosome
    * mu - growth rate [1/h]
    * k1 - kcat of substrate importer [1/h]
    * k2/k3/k5 - kcat of AA synthesis (k2) / NT synthesis (k3) / ribosome assembly (k5) [1/h]
    * c - transcription rate [nt/h]
    * kel - translation rate [aa/h]
    * c_pol_max - max. RNA polymerase [mmol/g]

    Returns:
    * numpy matrix
    """
    
    row_names = ["G", "AA", "NT", "rRNA", "rP", "C1", "C2", "C3", "C4", "CA", "CR", "c_pol_max", "DM"]
    rxns = ["vIG", "vEAA", "vENT", "vRNAP", "vrA", "wIG", "wEAA", "wENT", "wRNAP", "wrAS", "wRP", 
             "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "C"]
    nAA = mwAA/mwG
    nNT = (mwNT-mwAA)/mwG  
    
    S = np.array([[  1, -nAA, -nNT,   0,      0,   0,   0,   0,   0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0,   0],
                  [  0,  1,  -1,      0,      0, -n1, -n2, -n3, -n4, -n5, -nR,  0, 0,  0,  0,  0,  0,  0,  0,  0,   0],
                  [  0,  0,   1,    -mR,      0,   0,   0,   0,   0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0,   0],
                  [  0,  0,   0,      1,     -1,   0,   0,   0,   0,   0,   0, -1, 0,  0,  0,  0,  0,  0,  0,  0,   0],
                  [  0,  0,   0,      0,     -1,   0,   0,   0,   0,   0,   1,  0,-1,  0,  0,  0,  0,  0,  0,  0,   0],
                  [-mu,  0,   0,      0,      0,  k1,   0,   0,   0,   0,   0,  0, 0, -1,  0,  0,  0,  0,  0,  0,   0],
                  [  0,-mu,   0,      0,      0,   0,  k2,   0,   0,   0,   0,  0, 0,  0, -1,  0,  0,  0,  0,  0,   0],
                  [  0,  0, -mu,      0,      0,   0,   0,  k3,   0,   0,   0,  0, 0,  0,  0, -1,  0,  0,  0,  0,   0],
                  [  0,  0,   0, -mu*mR,      0,   0,   0,   0,   c,   0,   0,  0, 0,  0,  0,  0, -1,  0,  0,  0,   0],
                  [  0,  0,   0,      0,    -mu,   0,   0,   0,   0,  k5,   0,  0, 0,  0,  0,  0,  0, -1,  0,  0,   0],
                  [  0,  0,   0,      0, kel/mu, -n1, -n2, -n3, -n4, -n5, -nR,  0, 0,  0,  0,  0,  0,  0, -1,  0,   0],
                  [  0,  0,   0,      0,      0,   0,   0,   0, -1/mu, 0,   0,  0, 0,  0,  0,  0,  0,  0,  0, -1, c_pol_max],
                  [-mwG, 0,   0,      0,      0,   0,   0,   0,   0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  mu]])
    
    revs = [0] * S.shape[1]

    return((S, rxns, row_names, revs))



def make_matrix(mwG, mwAA, mwNT, n1, n2, n3, n4, n5, nR, mR, mu, k1, k2, k3, k5, c, kel):
    """create RBA matrix
    Parameters: 
    * mwG, mwAA, mwNT - molecular masses of glucose (or other substrate), AA and NT
    * n1, n2, n3, n4, n5 - stoichiometric coefficients of transporter (n1), AA/NT synthesis enzymes (n2/n3),
    RNAP (n4) and ribosome assembly factors (n5)
    * nR, mR - stoich coefficients of AA and NT in ribosome
    * mu - growth rate [1/h]
    * k1 - kcat of substrate importer [1/h]
    * k2/k3/k5 - kcat of AA synthesis (k2) / NT synthesis (k3) / ribosome assembly (k5) [1/h]
    * c - transcription rate [nt/h]
    * kel - translation rate [aa/h]

    Returns:
    * numpy matrix
    """
    
    row_names = ["G", "AA", "NT", "rRNA", "rP", "C1", "C2", "C3", "C4", "CA", "CR", "DM"]
    rxns = ["vIG", "vEAA", "vENT", "vRNAP", "vrA", "wIG", "wEAA", "wENT", "wRNAP", "wrAS", "wRP", 
             "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "C"]
    
    nAA = mwAA/mwG
    nNT = (mwNT-mwAA)/mwG
    S = np.array([[  1, -nAA, -nNT,   0,      0,   0,   0,   0,   0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  1,  -1,      0,      0, -n1, -n2, -n3, -n4, -n5, -nR,  0, 0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  0,   1,    -mR,      0,   0,   0,   0,   0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  0,   0,      1,     -1,   0,   0,   0,   0,   0,   0, -1, 0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  0,   0,      0,     -1,   0,   0,   0,   0,   0,   1,  0,-1,  0,  0,  0,  0,  0,  0,  0],
                  [-mu,  0,   0,      0,      0,  k1,   0,   0,   0,   0,   0,  0, 0, -1,  0,  0,  0,  0,  0,  0],
                  [  0,-mu,   0,      0,      0,   0,  k2,   0,   0,   0,   0,  0, 0,  0, -1,  0,  0,  0,  0,  0],
                  [  0,  0, -mu,      0,      0,   0,   0,  k3,   0,   0,   0,  0, 0,  0,  0, -1,  0,  0,  0,  0],
                  [  0,  0,   0, -mu*mR,      0,   0,   0,   0,   c,   0,   0,  0, 0,  0,  0,  0, -1,  0,  0,  0],
                  [  0,  0,   0,      0,    -mu,   0,   0,   0,   0,  k5,   0,  0, 0,  0,  0,  0,  0, -1,  0,  0],
                  [  0,  0,   0,      0, kel/mu, -n1, -n2, -n3, -n4, -n5, -nR,  0, 0,  0,  0,  0,  0,  0, -1,  0],
                  [-mwG,  0,   0,      0,      0,   0,   0,   0,   0,   0,   0,  0, 0,  0,  0,  0,  0,  0, 0, mu]])
    
    revs = [0] * S.shape[1]

    return((S, rxns, row_names, revs))


def make_matrix_2rib(mwG, mwAA, mwNT, n1, n2, n3, n4, n5, nR, mR, mu, k1, k2, k3, k5, c, kel, x):
    """create RBA matrix
    Parameters: 
    * mwG, mwAA, mwNT - molecular masses of glucose (or other substrate), AA and NT
    * n1, n2, n3, n4, n5 - stoichiometric coefficients of transporter (n1), AA/NT synthesis enzymes (n2/n3),
    RNAP (n4) and ribosome assembly factors (n5)
    * nR, mR - stoich coefficients of AA and NT in ribosome
    * mu - growth rate [1/h]
    * k1 - kcat of substrate importer [1/h]
    * k2/k3/k5 - kcat of AA synthesis (k2) / NT synthesis (k3) / ribosome assembly (k5) [1/h]
    * c - transcription rate [nt/h]
    * kel - translation rate [aa/h]
    * x - fraction of protein in ribosome

    Returns:
    * tuple of numpy array with stoich. matrix, list of reactions, list of row names
    """
    
    nAA = mwAA/mwG
    nNT = (mwNT-mwAA)/mwG
    row_names = ["G", "AA", "NT", "rRNA", "rP", "C1", "C2", "C3", "C4", "CArRNA", "CArP", "CR", "Cfrac", 
                 "DM"]
    rxns = ["vIG", "vEAA", "vENT", "vRNAP", "vrArRNA", "vrArP", "wIG", "wEAA", "wENT", "wRNAP", "wrAS", "wrAS2", "wRP", 
             "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "C"]
    
    S = np.array([[  1, -nAA, -nNT,   0,      0,    0,   0,   0,   0,   0,   0,   0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  1,  -1,      0,      0,    0, -n1, -n2, -n3, -n4, -n5, -n5, -nR, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  0,   1,    -mR,      0,    0,   0,   0,   0,   0,   0,   0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  0,   0,      1,     -1,    0,   0,   0,   0,   0,   0,   0,  0, -1, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  0,   0,      0,      0,   -1,   0,   0,   0,   0,   0,   0,  1,  0,-1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [-mu,  0,   0,      0,      0,    0,  k1,   0,   0,   0,   0,   0,  0,  0, 0, -1,  0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,-mu,   0,      0,      0,    0,   0,  k2,   0,   0,   0,   0,  0,  0, 0,  0, -1,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  0, -mu,      0,      0,    0,   0,   0,  k3,   0,   0,   0,  0,  0, 0,  0,  0, -1,  0,  0,  0,  0,  0,  0],
                  [  0,  0,   0, -mu*mR,      0,    0,   0,   0,   0,   c,   0,   0,  0,  0, 0,  0,  0,  0, -1,  0,  0,  0,  0,  0],
                  [  0,  0,   0,      0,    -mu,    0,   0,   0,   0,   0,  k5,   0,  0,  0, 0,  0,  0,  0,  0, -1,  0,  0,  0,  0],
                  [  0,  0,   0,      0,      0,  -mu,   0,   0,   0,   0,   0,  k5,  0,  0, 0,  0,  0,  0,  0,  0, -1,  0,  0,  0],
                  [  0,  0,   0,      0, kel/mu,kel/mu,-n1, -n2, -n3, -n4, -n5, -n5, -nR, 0, 0,  0,  0,  0,  0,  0,  0, -1,  0,  0],
                  [  0,  0,   0,      0, -x/mu,(1-x)/mu, 0,   0,   0,   0,  0,   0,    0, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [-mwG,  0,   0,      0,      0,   0,   0,   0,   0,   0,   0,   0,   0, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0, mu]])
    
    return((S, rxns, row_names))



def make_matrix_2rib_vmax(mwG, mwAA, mwNT, n1, n2, n3, n4, n5, nR, mR, mu, k1, k2, k3, k5, c, kel, x, v4max, c_pol_max):
    """create RBA matrix
    Parameters: 
    * mwG, mwAA, mwNT - molecular masses of glucose (or other substrate), AA and NT
    * n1, n2, n3, n4, n5 - stoichiometric coefficients of transporter (n1), AA/NT synthesis enzymes (n2/n3),
    RNAP (n4) and ribosome assembly factors (n5)
    * nR, mR - stoich coefficients of AA and NT in ribosome
    * mu - growth rate [1/h]
    * k1 - kcat of substrate importer [1/h]
    * k2/k3/k5 - kcat of AA synthesis (k2) / NT synthesis (k3) / ribosome assembly (k5) [1/h]
    * c - transcription rate [nt/h]
    * kel - translation rate [aa/h]
    * x - fraction of protein in ribosome
    * v4max -  max RNA synthesis rate [1/h]
    * c_pol_max - max. RNA polymerase [mmol/g]

    Returns:
    * tuple of numpy array with stoich. matrix, list of reactions, list of row names
    """
    
    nAA = mwAA/mwG
    nNT = (mwNT-mwAA)/mwG
    row_names = ["G", "AA", "NT", "rRNA", "rP", "C1", "C2", "C3", "C4", "CArRNA", "CArP", "CR", "vmax", "cpolmax", "Cfrac", 
                 "DM"]
    rxns = ["vIG", "vEAA", "vENT", "vRNAP", "vrArRNA", "vrArP", "wIG", "wEAA", "wENT", "wRNAP", "wrAS", "wrAS2", "wRP", 
             "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14", "C"]
    
    S = np.array([[  1, -nAA, -nNT,   0,      0,    0,   0,   0,   0,   0,   0,   0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  1,  -1,      0,      0,    0, -n1, -n2, -n3, -n4, -n5, -n5, -nR, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  0,   1,    -mR,      0,    0,   0,   0,   0,   0,   0,   0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  0,   0,      1,     -1,    0,   0,   0,   0,   0,   0,   0,  0, -1, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  0,   0,      0,      0,   -1,   0,   0,   0,   0,   0,   0,  1,  0,-1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [-mu,  0,   0,      0,      0,    0,  k1,   0,   0,   0,   0,   0,  0,  0, 0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,-mu,   0,      0,      0,    0,   0,  k2,   0,   0,   0,   0,  0,  0, 0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  0, -mu,      0,      0,    0,   0,   0,  k3,   0,   0,   0,  0,  0, 0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0],
                  [  0,  0,   0, -mu*mR,      0,    0,   0,   0,   0,   c,   0,   0,  0,  0, 0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0],
                  [  0,  0,   0,      0,    -mu,    0,   0,   0,   0,   0,  k5,   0,  0,  0, 0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0],
                  [  0,  0,   0,      0,      0,  -mu,   0,   0,   0,   0,   0,  k5,  0,  0, 0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0],
                  [  0,  0,   0,      0, kel/mu,kel/mu,-n1, -n2, -n3, -n4, -n5, -n5, -nR, 0, 0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0],
                  [  0,  0,   0,    -mu,      0,    0,   0,   0,   0,v4max, 0,   0,    0, 0, 0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0],
                  [  0,  0,   0,      0,      0,   0,    0,   0,   0,-1/mu, 0,   0,    0, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0, -1, c_pol_max],
                  [  0,  0,   0,      0, -x/mu,(1-x)/mu, 0,   0,   0,   0,  0,   0,    0, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                  [-mwG,  0,   0,      0,      0,   0,   0,   0,   0,   0,   0,   0,   0, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, mu]])
    
    return((S, rxns, row_names))


def run_efmtool(S, revs, rxns, row_names):
    """run efmtool with options that disable weird warnings
    Parameters:
        S: stoichiometric matrix
        revs: list of reversibilities
        rxns: list of reactions
        row_names: list of metabolites/constraints
        """

    options = efmtool.get_default_options()
    options["arithmetic"] = "fractional"
    options['level'] = 'WARNING'
    egvs = efmtool.calculate_efms(stoichiometry = S, reversibilities = revs, 
                              reaction_names = rxns, metabolite_names = row_names, 
                              options=options, jvm_options=["--illegal-access=deny"])
    return egvs
        
        

def bisection_search_mu(one_frac, revs, rxns, row_names, get_fluxes=False):
    """quickly find the biggest possible growth rate with bisection search
    Parameters: 
        * one_frac: dictionary of matrices for each tested growth rate
        * revs: list of reversibilities
        * rxns: list of reactions
        * row_names: list of row names
    Returns: 
        If get_fluxes=False:
        * optimal mu as float
        If get_fluxes=True
        * tuple of optimal mu as a float and calculated egvs as an array
    """
    
    mu_list = list(one_frac.keys())
    left = 0  # Determines the starting index of the list we have to search in
    right = len(mu_list)-1  # Determines the last index of the list we have to search in
    mid = (right + left)//2  # // means floored division
    new_mu = mu_list[mid]
    last_mu = -1
    iterations = 1

    while new_mu != last_mu: # check if we found the last mu
        # calculate efvs
        egvs = run_efmtool(one_frac[new_mu], revs, rxns, row_names)

        # if there are no egvs/only zero egvs in non-slack columns, growth rate is too big (no solution)
        if (egvs.shape[1] == 0) or (np.count_nonzero(egvs[:, :11], axis=None) == 0): #np.count_nonzero(egvs, axis=None)/egvs.shape[1] < 3:
            right = mid - 1
            
        # otherwise growth rate is too small (or just right)
        else:
            left = mid + 1
            last_mu = new_mu

        mid = (right + left)//2
        new_mu = mu_list[mid]

        iterations += 1
        
        if iterations > len(mu_list):
            print("Optimal growth rate not found")
            return(0)
            break
        
    if get_fluxes:
        return((float(last_mu)), egvs)
    else:
        return(float(last_mu))
