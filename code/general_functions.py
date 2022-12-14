"""
Functions used by other scripts
"""

import numpy as np
import pandas as pd
import efmtool
import plotext as plt

# Default values for several functions
RNA_EXPENSIVE = False
ACTIVITIES = False
KIN = None
RNAP_MAX = None


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
    # subtract one because zero EGVs have 1 for a slack variable
    nonzeros = (np.count_nonzero(egms, axis=0)-1).astype(bool)
    egms = egms.T[nonzeros]

    res = pd.DataFrame(egms,
                       index = [f"EGM{(i+1)}" for i in range(egms.shape[0])],
                       columns = rxns)

    # normalize values by the column C
    if "C" in res.columns:
        res = res.div(res.C, axis = 0)

    if drop_slack:
        # remove columns with slack variables
        cols = [c for c in res.columns if c.startswith("S")]
        if "C" in res.columns:
            cols = cols+["C"]
        res = res[res.columns.drop(cols)]
    return res


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
    for mu, res in egv_dict.items():
        #res = egv_dict[mu]
        res["mu"] = mu
        res["x"] = str(frac)
        res["EGVs"] = res.index
        all_results = pd.concat([all_results, res])
    all_results = all_results.set_index([all_results["mu"]+all_results.index])
    return all_results


def get_parameters(organism, activities=ACTIVITIES, rna_expensive=RNA_EXPENSIVE, medium = 2):
    """
    Get constant parameters for stoich. matrix with constraints
        * By default E. coli parameters in glucose minimal medium are returned

    Parameters:
        * organism: one of ["ecoli", "archea"]
          If organism="archaea", some parameters (c, kel, mwR) are adjusted for Archaea
        * activities: True or False
          If activities=True, kel and c elongation rates are multiplied
          by their active fractions (values from Kostinski 2020)
        * rna_expensive: True or False. If rna_expensive=True, parameters are
                         artificially adjusted to make RNA
          more expensive than proteins
        * medium: one of [0,1,2,3,4,5], kel and c parameters are medium specific
                  0: succinate minimal
                  1: glycerol minimal
                  2: glucose minimal [default]
                  3: glycerol + AA
                  4: glucose + AA
                  5: LB

    Returns:
        * k1 - kcat of substrate importer [1/h]
        * k2/k3/k5 - kcat of AA synthesis (k2) / NT synthesis (k3) / R assembly (k5) [1/h]
        * c - transcription rate [nt/h]
        * kel - translation rate [aa/h]
        * nig, neaa, nent, nrnap, naf - stoichiometric coefficients of transporter (nig),
                                        AA/NT synthesis enzymes (neaa/nent),
                                        RNAP (nrnap) and ribosome assembly factors (naf)
        * mwR, mwG, mwAA, mwNT - molecular masses of R, G, AA and NT [g/mmol]
    """
    to_h = 3600  # conversion factor from seconds to hours

    # kinetic parameters
    k1 = 180*to_h  # glucose transport rate, 1/h
    k2 = 10*to_h  # kcat of AA synthesis, 1/h
    k3 = 10*to_h  # kcat of NT synthesis, 1/h
    k5 = 1/120*to_h  # kcat of ribisome assembly, 1/h
    c = 85*to_h  # transcription elongation rate, nt/h
    kels = np.array([12, 16.83, 21, 20.17, 21, 22.25])*to_h
    kel = kels[medium]  # translation elongation rate, aa/h

    # stoichiometry
    nig = 646  # AA in glucose transporter
    avg_protein = 325  # AA in a typical protein, BNID 108986
    neaa = avg_protein*15  # estimated 15 enzymes, 9 steps in glycolysis,
                           # 5 steps from Pyr to Glu, 6 to Gln, other AAs variable
    nent = avg_protein*15  # estimated 15 enzymes, 6 PPP (glc -> prpp)
                           # + 5 pyrimidine/13 purine
    nrnap = 3498  # AA in RNAP without sigma factor
    naf = avg_protein*12  # https://doi.org/10.3389/fmicb.2019.02982

    # molecular masses (g/mmol)
    mwR = 2300  # ribosome
    mwG = 0.18  # glucose
    mwAA = 0.109  # avergage AA
    mwNT= 0.3243  # average NT

    # organism specific parameters
    if organism == "archaea":
        c = 25*to_h
        kel = 25/3*to_h
        mwR = 3040

    # parameters that make RNA more expensive than proteins
    if rna_expensive:
        kel = kel*3
        nrnap = nrnap*15
        c = c/10

    # R and RNAP elongation rates multiplied by their active fractions (Bremer 2003)
    if activities:
        kel = kel*0.85
        f_RNAP = np.array([13.2, 14.4, 15.0, 18.8, 24.2, 31.0])/100  # fraction active, in %
        c = c*f_RNAP[medium]

    return((k1, k2, k3, k5, c, kel, nig, neaa, nent, nrnap, naf, mwR, mwG, mwAA, mwNT))


def make_matrix(matrix_type, frac, mu, organism,
                activities=ACTIVITIES, rna_expensive=RNA_EXPENSIVE,
                rnap_max=RNAP_MAX, kin=KIN, medium=2):
    """
    Parameters:
        * matrix_type: one of ["no_constraints", "RBA", "rnap_max", "init_rate"]
        * frac: protein mass fraction in ribosome
        * mu: growth rate
        * organism: one of ["ecoli", "archaea"] -- it will get different parameters
          from "get_parameters()", default is "ecoli"
        * activities: True or False -- it will get different parameters
          from "get_parameters()", default is False
        * rnap_max: max. RNAP concentration [mmol/g]
        * kin: max. transcription initiation rate [1/h]
        * medium: one of [0,1,2,3,4,5], specifies kel and c in "get_parameters()"
    Returns:
        * S: numpy matrix -- stoichiometric matrix
        * rxns: list of reactions
        * row_names: list of rows
        * revs: list of reversibilities
    """

    k1, k2, k3, k5, c, kel, nig, neaa, nent, nrnap, naf, mwR, mwG, mwAA, mwNT = get_parameters(organism,
                                                                                               activities,
                                                                                               rna_expensive,
                                                                                               medium)
    # calculate stoichiometric coefficients
    naa = mwAA/mwG
    nnt = (mwNT-mwAA)/mwG
    nrp = mwR*frac/mwAA
    nrrna = mwR*(1-frac)/mwNT
    rxns = ["vIG", "vEAA", "vENT", "vRNAP", "vAF",
            "wIG", "wEAA", "wENT", "wRNAP", "wAF", "wrP"]
    row_names = ["G", "AA", "NT", "rRNA", "rP",
                 "IG", "EAA", "ENT", "RNAP", "AF", "R"]

    if matrix_type == "no_constraints":
        S = np.array([[1, -naa, -nnt,      0,  0,    0,     0,     0,      0,    0,    0],
                      [0,    1,   -1,      0,  0, -nig, -neaa, -nent, -nrnap, -naf, -nrp],
                      [0,    0,    1, -nrrna,  0,    0,     0,     0,      0,    0,    0],
                      [0,    0,    0,      1, -1,    0,     0,     0,      0,    0,    0],
                      [0,    0,    0,      0, -1,    0,     0,     0,      0,    0,    1],
                      [0,    0,    0,      0,  0,    1,     0,     0,      0,    0,    0],
                      [0,    0,    0,      0,  0,    0,     1,     0,      0,    0,    0],
                      [0,    0,    0,      0,  0,    0,     0,     1,      0,    0,    0],
                      [0,    0,    0,      0,  0,    0,     0,     0,      1,    0,    0],
                      [0,    0,    0,      0,  0,    0,     0,     0,      0,    1,    0],
                      [0,    0,    0,      0,  1,    0,     0,     0,      0,    0,    0]])

    if matrix_type == "RBA":
        row_names = row_names[:5] + ["C1", "C2", "C3", "C4", "CA", "CR", "DM"]
        rxns = rxns + ["S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "C"]

        S = np.array([[  1, -naa, -nnt,     0,      0,    0,     0,     0,     0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0],
                      [  0,  1,    -1,      0,      0, -nig, -neaa, -nent,-nrnap,-naf, -nrp,  0, 0,  0,  0,  0,  0,  0,  0,  0],
                      [  0,  0,     1,    -nrrna,      0,    0,     0,     0,     0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0],
                      [  0,  0,     0,      1,     -1,    0,     0,     0,     0,   0,   0, -1, 0,  0,  0,  0,  0,  0,  0,  0],
                      [  0,  0,     0,      0,     -1,    0,     0,     0,     0,   0,   1,  0,-1,  0,  0,  0,  0,  0,  0,  0],
                      [-mu,  0,     0,      0,      0,   k1,     0,     0,     0,   0,   0,  0, 0, -1,  0,  0,  0,  0,  0,  0],
                      [  0,-mu,     0,      0,      0,    0,    k2,     0,     0,   0,   0,  0, 0,  0, -1,  0,  0,  0,  0,  0],
                      [  0,  0,   -mu,      0,      0,    0,     0,    k3,     0,   0,   0,  0, 0,  0,  0, -1,  0,  0,  0,  0],
                      [  0,  0,     0, -mu*nrrna,      0,    0,     0,     0,     c,   0,   0,  0, 0,  0,  0,  0, -1,  0,  0,  0],
                      [  0,  0,     0,      0,    -mu,    0,     0,     0,     0,  k5,   0,  0, 0,  0,  0,  0,  0, -1,  0,  0],
                      [  0,  0,     0,      0, kel/mu, -nig, -neaa, -nent,-nrnap, -naf,-nrp,  0, 0,  0,  0,  0,  0,  0, -1,  0],
                      [-mwG,  0,    0,      0,      0,    0,     0,     0,     0,   0,   0,  0, 0,  0,  0,  0,  0,  0, 0, mu]])

    if matrix_type == "RNAPmax":
        row_names = row_names[:5] + ["C1", "C2", "C3", "C4", "CA", "CR", "c_pol_max", "DM"]
        rxns = rxns + ["S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "C"]

        S = np.array([[  1, -naa, -nnt,   0,      0,   0,      0,     0,     0,     0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0,   0],
                      [  0,  1,  -1,      0,      0, -nig, -neaa, -nent, -nrnap, -naf, -nrp, 0, 0,  0,  0,  0,  0,  0,  0,  0,   0],
                      [  0,  0,   1,    -nrrna,      0,   0,      0,     0,     0,     0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0,   0],
                      [  0,  0,   0,      1,     -1,   0,      0,     0,     0,     0,   0, -1, 0,  0,  0,  0,  0,  0,  0,  0,   0],
                      [  0,  0,   0,      0,     -1,   0,      0,     0,     0,     0,   1,  0,-1,  0,  0,  0,  0,  0,  0,  0,   0],
                      [-mu,  0,   0,      0,      0,  k1,      0,     0,     0,     0,   0,  0, 0, -1,  0,  0,  0,  0,  0,  0,   0],
                      [  0,-mu,   0,      0,      0,   0,     k2,     0,     0,     0,   0,  0, 0,  0, -1,  0,  0,  0,  0,  0,   0],
                      [  0,  0, -mu,      0,      0,   0,      0,    k3,     0,     0,   0,  0, 0,  0,  0, -1,  0,  0,  0,  0,   0],
                      [  0,  0,   0, -mu*nrrna,      0,   0,      0,     0,     c,     0,   0,  0, 0,  0,  0,  0, -1,  0,  0,  0,   0],
                      [  0,  0,   0,      0,    -mu,   0,      0,     0,     0,    k5,   0,  0, 0,  0,  0,  0,  0, -1,  0,  0,   0],
                      [  0,  0,   0,      0, kel/mu, -nig, -neaa, -nent,-nrnap,  -naf, -nrp,  0, 0,  0,  0,  0,  0,  0, -1,  0,   0],
                      [  0,  0,   0,      0,      0,   0,      0,     0, -1/mu,     0,   0,  0, 0,  0,  0,  0,  0,  0,  0, -1, rnap_max],
                      [-mwG, 0,   0,      0,      0,   0,      0,     0,     0,     0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  mu]])

    if matrix_type == "init_rate":
        row_names = row_names[:5] + ["C1", "C2", "C3", "C4", "CA", "CR", "DM"]
        rxns = rxns + ["S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "C"]

        S = np.array([[  1, -naa, -nnt,   0,      0,   0,   0,   0,    0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,   0],
                      [  0,  1,  -1,      0,      0, -nig, -neaa, -nent, -nrnap, -naf, -nrp, 0, 0,  0,  0,  0,  0,  0,  0,   0],
                      [  0,  0,   1,    -nrrna,      0,   0,   0,   0,    0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,   0],
                      [  0,  0,   0,      1,     -1,   0,   0,   0,    0,   0,   0, -1, 0,  0,  0,  0,  0,  0,  0,   0],
                      [  0,  0,   0,      0,     -1,   0,   0,   0,    0,   0,   1,  0,-1,  0,  0,  0,  0,  0,  0,   0],
                      [-mu,  0,   0,      0,      0,  k1,   0,   0,    0,   0,   0,  0, 0, -1,  0,  0,  0,  0,  0,   0],
                      [  0,-mu,   0,      0,      0,   0,  k2,   0,    0,   0,   0,  0, 0,  0, -1,  0,  0,  0,  0,   0],
                      [  0,  0, -mu,      0,      0,   0,   0,  k3,    0,   0,   0,  0, 0,  0,  0, -1,  0,  0,  0,   0],
                      [  0,  0,   0,    -mu,      0,   0,   0,   0, kin,   0,   0,  0, 0,  0,  0,  0, -1,  0,  0,   0],
                      [  0,  0,   0,      0,    -mu,   0,   0,   0,    0,  k5,   0,  0, 0,  0,  0,  0,  0, -1,  0,   0],
                      [  0,  0,   0,      0, kel/mu, -nig, -neaa, -nent,-nrnap, -naf,-nrp,  0, 0,  0,  0,  0,  0,  0, -1,  0],
                      [-mwG, 0,   0,      0,      0,   0,   0,   0,    0,   0,   0,  0, 0,  0,  0,  0,  0,  0,  0,  mu]])

    revs = [0] * S.shape[1]

    return((S, rxns, row_names, revs))


def run_efmtool(S, revs, rxns, row_names):
    """run efmtool with options that disable weird warnings
    Parameters:
        S: stoichiometric matrix
        revs: list of reversibilities
        rxns: list of reactions
        row_names: list of metabolites/constraints
    Returns:
        egvs: array of calculated EGMs/EGVs
    """

    options = efmtool.get_default_options()
    options["arithmetic"] = "fractional"
    options['level'] = 'WARNING'
    egvs = efmtool.calculate_efms(stoichiometry = S,
                                  reversibilities = revs,
                                  reaction_names = rxns,
                                  metabolite_names = row_names,
                                  options = options,
                                  jvm_options = ["--illegal-access=deny"])
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
    left = 0  # The starting index of the list we have to search in
    right = len(mu_list)-1  # the last index of the list we have to search in
    mid = (right + left)//2  # // means floored division
    new_mu = mu_list[mid]
    last_mu = -1
    iterations = 1

    while new_mu != last_mu: # check if we already found the last mu
        egvs = run_efmtool(one_frac[new_mu], revs, rxns, row_names)

        # no solution (no egvs/zero egvs in non-slack columns) -- growth rate too big
        if (egvs.shape[1] == 0) or (np.count_nonzero(egvs[:, :11], axis=None) == 0):
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
            return 0

    if get_fluxes:
        # recalculate egvs for last_mu (the last feasible solution)
        egvs = run_efmtool(one_frac[last_mu], revs, rxns, row_names)
        return((float(last_mu)), egvs)

    return float(last_mu)


def initiate_progress_plot(title):
    """
    Make an empty plot in the terminal
    it will be updated with the function "fill_progress_plot()"
    """
    plt.title(title)
    plt.xlabel("Protein fraction in ribosome")
    plt.ylabel("Growth rate [1/h]")
    plt.theme("clear")
    plt.plot_size(70,25)


def fill_progress_plot(xdata, ydata, ymin=0, ymax=4):
    """
    Plot data into the plot initiated with "initiate_progress_plot()"
    """
    plt.clt()
    plt.ylim(ymin,ymax)
    plt.xlim(0,1)
    plt.scatter(xdata, ydata, marker="heart", color = 125)
    plt.show()


def calculate_concentrations(all_results, original_matrices, mets, rxns, organism,
                             activities=ACTIVITIES, rna_expensive=RNA_EXPENSIVE):
    """
    Calculates metabolite mass fractions in g/g (x = (Nv/mu)*mw).
    Parameters:
        * all_results: output of the function egv_dict_to_df
                       or a data frame with columns: rxns, mu, x
                       indices must be unique
        * original_matrices: a dictionary of matrices for each protein fraction
                             without constraints and with all catalysts
                             made with "make_matrix(matrix_type="no_constraints")"
        * mets: list of all metabolites
        * rxns: list of reactions
        * organism: one of ["ecoli", "archea"], changes parameters got from "get_parameters()"
        * activities: True or False, changes parameters got from "get_parameters()"
        * rna_expensive: True or False,  changes parameters got from "get_parameters()"
    Returns:
        * pandas DataFrame with mass fractions of the metabolites
    """

    *_, nig, neaa, nent, nrnap, naf, mwR, mwG, mwAA, mwNT = get_parameters(organism,
                                                                           activities,
                                                                           rna_expensive)

    all_concentrations = pd.DataFrame(columns = mets, index = all_results.index)
    for egv in all_results.index:
        frac = all_results.loc[egv, "x"]
        mat = pd.DataFrame(original_matrices[frac], index=mets, columns=rxns)
        nR = mwR*float(frac)/mwAA
        mR = mwR*(1-float(frac))/mwNT

        # molecular masses
        mw = [mwG, mwAA, mwNT, mR*mwNT, nR*mwAA, nig*mwAA, neaa*mwAA,
              nent*mwAA, nrnap*mwAA, naf*mwAA, mwR]

        mu = all_results.loc[egv, "mu"]
        fluxes = all_results.loc[egv][:11]
        all_concentrations.loc[egv] = mat.multiply(fluxes).sum(axis=1)/float(mu)
        all_concentrations.loc[egv] = all_concentrations.loc[egv]*mw
    all_concentrations["mu"] = all_results.mu
    all_concentrations["EGVs"] = all_results.EGVs
    all_concentrations["x"] = all_results.x

    return all_concentrations


def calculate_allocations(all_egvs,
                          organism,
                          activities=ACTIVITIES,
                          rna_expensive=RNA_EXPENSIVE):
    """
    Calculates protein allocations (phis).
    Parameters:
        * all_egvs: pandas DataFrame with calculated EGVs
                    rows: EGVs
                    columns: reactions, "mu" (growth rate),
                    "x" (rP mass fraction in R), "EGVs" (EGV number)
        * organism: one of ["ecoli", "archea"], changes parameters from "get_parameters()"
        * activities: True or False, changes parameters from "get_parameters()"
        * rna_expensive: True or False, changes parameters from "get_parameters()"

    Returns:
        * pandas DataFrame with ribosome allocations
    """

    *_, kel, nig, neaa, nent, nrnap, naf, mwR, mwG, mwAA, mwNT = get_parameters(organism,
                                                                                activities,
                                                                                rna_expensive)

    # multiply fluxes by stoichiometries
    protein_columns = ["wIG", "wEAA", "wENT", "wRNAP", "wAF"]
    protein_stoichiometries = [nig, neaa, nent, nrnap, naf]
    prot_fluxes = all_egvs[protein_columns].multiply(protein_stoichiometries)

    # rP stoichiometry calculated from 'x' (rP mass fraction in ribosome)
    prot_fluxes["wrP"] = all_egvs["wrP"]*mwR*all_egvs["x"].astype(float)/mwAA

    # multiply with growth rates and divide with assembly flux and elongation rate
    phis = prot_fluxes.apply(lambda x: x*all_egvs["mu"]/(all_egvs["vAF"]*kel))

    # add column with protein fractions
    phis["x"] = all_egvs["x"]

    return phis


def run_simulations(prot_fractions, mus, matrix_type, organism,
                    activities=ACTIVITIES, rna_expensive=RNA_EXPENSIVE,
                    rnap_max=RNAP_MAX, kin=KIN,
                    medium=2, plot=True):
    """
    Iterate over protein fractions (prot_fractions)
    For each fraction calculate max. growth rate and fluxes
    (function "bisection_search_mu()").

    Parameters:
    * matrix_type: type of matrix to create with "make_matrix",
                    one of ["no_constraints", "RBA", "RNAPmax", "init_rate"]
    * organism: one of ["ecoli", "archea"], changes parameters from "get_parameters()"
    * activities: True or False, changes parameters from "get_parameters()"
    * rna_expensive: True or False, changes parameters from "get_parameters()"
    * rnap_max: max. concentration of RNAP in mmol/g/h
    * kin: max. transcription initiation rate in 1/h
    * prot_fractions: rP mass fractions in ribosome to test
    * mus: growth rates to test
    * medium: one of [0,1,2,3,4,5], kel and c parameters are medium specific, see "get_parameters()"
    * plot: if True, progress plot is shown in the terminal

    Returns:
    res: pandas DataFrame with columns "x" (rP mass fraction in R)
         and "mu" (max. growth rate)
    all_egvs: pandas DataFrame with fluxes
              rows - EGMs
              columns - reactions, "mu" (growth rate),
                    "x" (rP mass fraction in R), "EGVs" (EGV number)
    conc: pandas DataFrame with metabolite mass fractions
          caluclated with "calculate_concentrations()"
    phis: pandas DataFrame with ribosome allocations
          calculated with "calculate_allocations()"
    """

    original_matrices = {}
    last_mus = []
    all_egvs = pd.DataFrame()

    if plot:
        initiate_progress_plot(f"{matrix_type}, {organism}, activities: {activities}, RNAPmax: {rnap_max}")

    for prog, frac in enumerate(prot_fractions):
        # matrices without constraints
        S, rxns1, mets, _ = make_matrix("no_constraints", frac, 0,
                                         organism, activities, rna_expensive,
                                         rnap_max, kin, medium)
        original_matrices[str(round(frac,5))] = S

        # matrices with constraints
        matrices = {}
        for mu in mus:
            S, rxns, row_names, revs = make_matrix(matrix_type, frac, mu,
                                                   organism, activities, rna_expensive,
                                                   rnap_max, kin, medium)
            matrices[str(round(mu,4))] = S

        # calculate max. growth rate and EGVs
        last_mu, egvs = bisection_search_mu(matrices, revs, rxns, row_names, True)
        last_mus.append(last_mu)

        # save egvs
        egv_df = make_nice_results(egvs, rxns)
        egv_df["x"] = str(round(frac,5))  # str because frac are used as keys in dicts
                                          # (floats make problems)
        egv_df["mu"] = last_mu
        egv_df["EGVs"] = egv_df.index
        egv_df = egv_df.set_index([egv_df["x"]+egv_df.index])  # make unique index
        all_egvs = pd.concat([all_egvs, egv_df])

        # plot progress in the terminal
        if plot:
            fill_progress_plot(prot_fractions[:prog], last_mus, 0.5, 4.5)

    res = pd.DataFrame({"x": prot_fractions, "mu": last_mus})

    # calculate concentrations and ribosome allocations
    conc = calculate_concentrations(all_egvs, original_matrices, mets, rxns1,
                                    organism, activities, rna_expensive)
    phis = calculate_allocations(all_egvs, organism, activities, rna_expensive)

    return(res, all_egvs, conc, phis)
