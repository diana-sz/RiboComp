"""
General classes and functions used by other scripts
* Model definitions
* Simulations with efmtool
* Calculations of mass fractions and ribosome allocations

Author: Diana Szeliova, 8.2.2023
"""

import efmtool
import numpy as np
import pandas as pd
import plotext as plt

TO_H = 3600  # factor to convert from seconds to hours

def slack(n_columns, pos = -100):
    """
    make a 1D array with zeros and optionally with a -1 at a specified position

    Parameters
    ----------
    pos: optional, specifies position where -1 is inserted

    Returns
    -------
    numpy array
    """
    zeros = np.zeros(n_columns, dtype = int)
    if pos in range(len(zeros)):
        zeros[pos] = -1
    return list(zeros)


def initiate_progress_plot(title):
    """
    Make an empty plot in the terminal
    it will be updated with the function "fill_progress_plot()"

    Parameters
    ----------
    title: title of the plot
    """
    plt.title(title)
    plt.xlabel("Protein fraction in ribosome")
    plt.ylabel("Growth rate [1/h]")
    plt.theme("clear")
    plt.plot_size(70,25)


def fill_progress_plot(xdata, ydata, ymin=0, ymax=4):
    """
    Plot data into the plot initiated with "initiate_progress_plot()"

    Parameters
    ----------
    xdata, ydata: x and y coordinates
    ymin, ymax: min and max limits of y-axis
    """
    plt.clt()
    plt.ylim(ymin,ymax)
    plt.xlim(0,1)
    plt.scatter(xdata, ydata, marker="heart", color = 125)
    plt.show()


class Model:
    """
        Class to store models and results of a single run of efmtool

        Attributes
        ----------
        parameter_set: one of ["default", "activities", "archaea", "rna_expensive"]
            - changes some parameters from defaults
            - see function "set_specific_parameters()"
        matrix_type: see funcions "make_matrix_uc", "make_matrix_c"
            - "RBA": standard RBA (default)
            - "RBA_init": RBA with RNAP capacity constraint on initiation
                           instead of elongation
            - "RNAPmax": RBA + constraint on max. RNAP concentration
            - "R_deg", "R_deg2", "R_deg_hill": matrices with ribosome degradation reaction
                                               three different functions for calculation
                                               of degradation rate
            - "PRL" - RBA with fixed ribosome and RNAP allocations
        medium: influences the values of keaa, mol_masses["C"], c (if parameter_set = "activities")
            - 0: succinate minimal medium
            - 1: glycerol minimal medium
            - 2: glucose minimal medium [default]
            - 3: glycerol + amino acids
            - 4: glucose + amino acids
            - 5: LB
        frac: protein mass fraction in the ribosome
        mu: growth rate
        rnapmax: maximum concentration of RNAP (mmol/g/h)
                 used only for matrix_type: RNAPmax
        mol_masses: a dictionary of molecular masses (g/mmol) of ribosome (R),
                    average amino acid (AA), average nucleotide (NT), carbon source (C)
        stoichiometries: a dictionary of stoichiometric coefficients
            - nig, neaa, nent, nrnap, naf, nrnase - transporter (nig),
              AA/NT synthesis enzymes (neaa/nent), RNAP (nrnap),
              ribosome assembly factors (naf) and RNase R (nrnase)
            - naa, nnt - stoich. coefficients for synthesis of AA, NT,
                         depend on the molecular mass of carbon source (mol_masses["C"])
            - nrp, nrrna - stoich. coefficients for synthesis of rP and rRNA, depend on "frac"
        kinetics: a dictionary of kinetic parameters
            - kig: kcat of substrate importer [1/h]
            - kent: kcat of NT synthesis [1/h]
            - kaf: kcat of R assembly [1/h]
            - keaa - kcat of AA synthesis (keaa), depends on the value of "medium"
            - ktr: transcription rate [nt/h]
            - krnap: transcription rate [1/h]
            - kel: translation rate [aa/h]
            - kexo: exonuclease rate (RNase R) [nt/h]
            - kdeg_max: max. degradation rate [1/h]
            - kin: transcription initiation rate [1/h]
        matrix_uc: np.array - matrix without constraints
        matrix_c: np.array - matrix with constraints
        rows_uc, rows_c: row names for unconstrained / constrained matrices
        columns_uc, columns_c: column names for unconstrained / constrained matrices
    """

    def __init__(self, *model_pars, **kwargs):
        """
        1. set default parameters - E. coli in glucose minimal medium
        2. make matrices without ("make_matrix_uc") and with constraints
           ("make_matrix_c")
        """

        # default simulation settings
        self.parameter_set = "default"
        self.matrix_type = "RBA"
        self.medium = 2
        self.frac = 0.36
        self.growth_rate = 1
        self.rnapmax = 1000
        self.kin = 1000
        self.kdeg_max = 0

        # molecular masses
        self.mol_masses = {"R": 2300,
                           "AA": 0.109,
                           "NT": 0.3243}
        # stoichiometries
        avg_protein = 325
        self.stoichiometries = {"nig": 646,
                                "neaa": avg_protein*15,
                                "nent": avg_protein*15,
                                "nrnap": 3498,
                                "naf": avg_protein*12,
                                "nrnase": 813}
        # kinetic parameters
        self.kinetics = {"kig": 180*TO_H,
                         "kent": 10*TO_H,
                         "kaf": 1/120*TO_H,
                         "ktr": 85*TO_H,
                         "kel": 21*TO_H,
                         "kexo": 88*TO_H}

        # if model_pars, kwargs are provided, update basic settings
        for dictionary in model_pars:
            for key, value in dictionary.items():
                setattr(self, key, value)
        for key, value in kwargs.items():
            setattr(self, key, value)

        # set parameters for specific medium/parameter set
        self.set_specific_parameters(*model_pars, **kwargs)

        # calculated stoichiometries - must be placed after all other parameters are set
        self.stoichiometries["naa"] = self.mol_masses["AA"]/self.mol_masses["C"]
        self.stoichiometries["nnt"] = (self.mol_masses["NT"]-self.mol_masses["AA"])/self.mol_masses["C"]
        self.stoichiometries["nrp"] = self.mol_masses["R"]*self.frac/self.mol_masses["AA"]
        self.stoichiometries["nrrna"] = self.mol_masses["R"]*(1-self.frac)/self.mol_masses["NT"]

        self.rows_c = []
        self.columns_c = []

        # create matrices
        self.make_matrix_uc()
        self.make_matrix_c()


    def set_specific_parameters(self, *model_pars, **kwargs):
        """change a specific set of parameters according to values
        of "Model.medium", "Model.parameter_set"
        """
        # make a list of custom parameters and check
        # if the specific parameters overwrite any of them
        custom_attr = []
        for dictionary in model_pars:
            for key in dictionary:
                custom_attr.append(key)
        for key in kwargs:
            custom_attr.append(key)

        # medium specific parameters - always set
        # kcat of AA synthesis, 1/h
        self.kinetics["keaa"] = [2, 3.5, 5, 7, 8.5, 10.5][self.medium]*TO_H

        # molecular mass of the carbon source
        # 1,3: glycerol; 2,4,5: glucose; 0: succinate
        carbon_sources = {5: 0.18, 4: 0.18, 3: 0.092,
                          2: 0.18, 1: 0.092, 0: 0.118}
        self.mol_masses["C"] = carbon_sources[self.medium]
        overwritten = [any(i in custom_attr for i in ["kinetics", "mol_masses"])]

        if self.parameter_set == "archaea":
            self.kinetics["ktr"] = 25*TO_H
            self.kinetics["kel"] = 25/3*TO_H
            self.stoichiometries["nrnap"] = 3338
            self.mol_masses["R"] = 3040
            overwritten.append(any(i in custom_attr for i in ["kinetics",
                                                              "stoichiometries",
                                                              "mol_masses"]))
        if self.parameter_set == "rna_expensive":
            self.stoichiometries["nrnap"] = 3498*15
            self.kinetics["kel"] = 21*TO_H*3
            self.kinetics["ktr"] = 85/10*TO_H
            overwritten.append(any(i in custom_attr for i in ["kinetics",
                                                              "stoichiometries"]))

        if self.parameter_set == "activities":
            self.kinetics["kel"] = self.kinetics["kel"]*0.85
            f_rnap = np.array([13.2, 14.4, 15.0, 18.8, 24.2, 31.0])/100  # % active RNAP
            self.kinetics["ktr"] = self.kinetics["ktr"]*f_rnap[self.medium]
            overwritten.append(any(i in custom_attr for i in ["kinetics"]))

        if self.parameter_set == "Kostinski":
            # parameters from Kostinski & Reuveni
            kel = [12, 16.83, 21, 20.17, 21, 22.25][self.medium]*3600
            self.kinetics["kel"] = kel*0.85  # % active ribosomes

            # % active polymerase
            f_rnap = [13.2, 14.4, 15.0, 18.8, 24.2, 31.0][self.medium]/100

            # RNAP allocation fractions
            phi_rnap = [0.18, 0.28, 0.42, 0.52, 0.60, 0.65][self.medium]

            self.kinetics["ktr"] = self.kinetics["ktr"]*f_rnap*phi_rnap
            overwritten.append(any(i in custom_attr for i in ["kinetics"]))

        if any(overwritten):
            print(f"Custom parameters may be overwritten with parameter set: {self.parameter_set}")


    def caluclate_degradation_rate(self):
        """
        Calculate ribosome degradation rate.
        The rate decreases with protein fraction in the ribosome (frac)
        either linearly, with a square of frac or according to a Hill function

        Returns
        -------
        float of the degradation rate
        """
        # linear increase
        if self.matrix_type == "R_deg":
            kdeg = self.kdeg_max*(1-self.frac)

        # squared
        if self.matrix_type == "R_deg2":
            kdeg = self.kdeg_max*(1-self.frac)**2

        # hill function
        if self.matrix_type == "R_deg_hill":
            hill = 6
            kdeg = self.kdeg_max*(1-(self.frac**hill/(self.frac**hill+0.2**hill)))

        return kdeg


    def __str__(self):
        return f"""
        Parameter set: {self.parameter_set}
        Matrix type: {self.matrix_type}
        Medium: {self.medium}
        """


    def make_matrix_uc(self):
        """Make matrix without constraints"""
        naa, nnt, = self.stoichiometries["naa"], self.stoichiometries["nnt"]
        nig, neaa, = self.stoichiometries["nig"], self.stoichiometries["neaa"]
        nent, nrnap = self.stoichiometries["nent"], self.stoichiometries["nrnap"]
        naf, nrp = self.stoichiometries["naf"], self.stoichiometries["nrp"]
        nrrna, nrnase = self.stoichiometries["nrrna"], self.stoichiometries["nrnase"]

        self.matrix_uc = np.array(
             [[1,-naa,-nnt,     0,  0,   0,    0,    0,     0,   0,   0],
              [0,   1,  -1,     0,  0,-nig,-neaa,-nent,-nrnap,-naf,-nrp],
              [0,   0,   1,-nrrna,  0,   0,    0,    0,     0,   0,   0],
              [0,   0,   0,     1, -1,   0,    0,    0,     0,   0,   0],
              [0,   0,   0,     0, -1,   0,    0,    0,     0,   0,   1],
              [0,   0,   0,     0,  0,   1,    0,    0,     0,   0,   0],
              [0,   0,   0,     0,  0,   0,    1,    0,     0,   0,   0],
              [0,   0,   0,     0,  0,   0,    0,    1,     0,   0,   0],
              [0,   0,   0,     0,  0,   0,    0,    0,     1,   0,   0],
              [0,   0,   0,     0,  0,   0,    0,    0,     0,   1,   0],
              [0,   0,   0,     0,  1,   0,    0,    0,     0,   0,   0]]
        )

        self.columns_uc = ["vIG", "vEAA", "vENT", "vRNAP", "vAF",
                           "wIG", "wEAA", "wENT", "wRNAP", "wAF", "wrP"]
        self.rows_uc = ["G", "AA", "NT", "rRNA", "rP",
                        "IG", "EAA", "ENT", "RNAP", "AF", "R"]

        # matrix with ribosome degradation rate
        if "deg" in self.matrix_type:
            self.matrix_uc = np.array(
             [[1,-naa,-nnt,     0,  0,    0,   0,    0,    0,     0,   0,     0,   0],
              [0,   1,  -1,     0,  0,    0,-nig,-neaa,-nent,-nrnap,-naf,-nrnase,-nrp],
              [0,   0,   1,-nrrna,  0,nrrna,   0,    0,    0,     0,   0,     0,   0],
              [0,   0,   0,     1, -1,    0,   0,    0,    0,     0,   0,     0,   0],
              [0,   0,   0,     0, -1,    1,   0,    0,    0,     0,   0,     0,   1],
              [0,   0,   0,     0,  1,   -1,   0,    0,    0,     0,   0,     0,   0],
              [0,   0,   0,     0,  0,    0,   1,    0,    0,     0,   0,     0,   0],
              [0,   0,   0,     0,  0,    0,   0,    1,    0,     0,   0,     0,   0],
              [0,   0,   0,     0,  0,    0,   0,    0,    1,     0,   0,     0,   0],
              [0,   0,   0,     0,  0,    0,   0,    0,    0,     1,   0,     0,   0],
              [0,   0,   0,     0,  0,    0,   0,    0,    0,     0,   1,     0,   0],
              [0,   0,   0,     0,  0,    0,   0,    0,    0,     0,   0,     1,   0]]
            )
            self.columns_uc = ["vIG", "vEAA", "vENT", "vRNAP", "vAF", "vdeg",
                               "wIG", "wEAA", "wENT", "wRNAP", "wAF", "wRNase", "wrP"]
            self.rows_uc = ["G", "AA", "NT", "rRNA", "rP", "R",
                            "IG", "EAA", "ENT", "RNAP", "AF", "RNase"]


    def make_matrix_c(self):
        """
        Create stoichiometric matrix with constraints (matrix_c), list of rows (rows_c),
        list of columns (columns_c) based on the value of "matrix_type" attribute
        """
        mu = self.growth_rate
        nig, neaa, = self.stoichiometries["nig"], self.stoichiometries["neaa"]
        nent, nrnap = self.stoichiometries["nent"], self.stoichiometries["nrnap"]
        naf, nrp = self.stoichiometries["naf"], self.stoichiometries["nrp"]
        nrrna, nrnase = self.stoichiometries["nrrna"], self.stoichiometries["nrnase"]
        kig, keaa = self.kinetics["kig"], self.kinetics["keaa"]
        kent, kaf, kel = self.kinetics["kent"], self.kinetics["kaf"], self.kinetics["kel"]
        mwc = self.mol_masses["C"]
        krnap = self.kinetics["ktr"]/nrrna  # transcripts/h
        if "_init_rate" in self.matrix_type:
            krnap = self.kin  # initiation rate limiting

        if self.matrix_type in ["RBA", "RBA_init_rate", "RNAPmax"]:
            # take the metabolite rows from the unconstrained stocihiometric matrix
            met_matrix = self.matrix_uc[0:5,:]

            # add slack columns
            ncol = 9
            if self.matrix_type == "RNAPmax":
                ncol = 10  # RNAPmax matrix needs one extra slack column
            slacks = np.array([slack(ncol), slack(ncol), slack(ncol),
                               slack(ncol, 0), slack(ncol, 1)])
            met_matrix = np.concatenate((met_matrix, slacks), axis = 1)

            # rba constraints
            constraints = np.array(
              [[-mu,  0,  0,  0,     0, kig,    0,    0,     0,   0,   0]+slack(ncol,2),
               [  0,-mu,  0,  0,     0,   0, keaa,    0,     0,   0,   0]+slack(ncol,3),
               [  0,  0,-mu,  0,     0,   0,    0, kent,     0,   0,   0]+slack(ncol,4),
               [  0,  0,  0,-mu,     0,   0,    0,    0, krnap,   0,   0]+slack(ncol,5),
               [  0,  0,  0,  0,   -mu,   0,    0,    0,     0, kaf,   0]+slack(ncol,6),
               [  0,  0,  0,  0,kel/mu,-nig,-neaa,-nent,-nrnap,-naf,-nrp]+slack(ncol,7),
               [-mwc, 0,  0,  0,     0,   0,    0,    0,    0,   0,    0]+slack(ncol-1)+[mu]]
            )

            # constraint on the max. RNAP concentration
            if self.matrix_type == "RNAPmax":
                max_rnap = np.array(
                  [[0, 0, 0, 0, 0, 0, 0, 0, -1/mu, 0, 0]+slack(ncol-1,8)+[self.rnapmax]])
                constraints = np.concatenate((constraints, max_rnap), axis = 0)

        if "R_deg" in self.matrix_type:
            met_matrix = self.matrix_uc[0:6,:]
            ncol = 12
            slacks = np.array([slack(ncol), slack(ncol), slack(ncol),
                               slack(ncol, 0), slack(ncol, 1), slack(ncol, 2)])
            met_matrix = np.concatenate((met_matrix, slacks), axis = 1)

            krnase = self.kinetics["kexo"]/nrrna
            kdeg = self.caluclate_degradation_rate()

            constraints = np.array(
             [[-mu,0, 0, 0,     0,      0, kig,    0,    0,     0,   0,    0,   0]+slack(ncol,3),
              [ 0,-mu,0, 0,     0,      0,   0, keaa,    0,     0,   0,    0,   0]+slack(ncol,4),
              [ 0, 0,-mu,0,     0,      0,   0,    0, kent,     0,   0,    0,   0]+slack(ncol,5),
              [ 0, 0, 0,-mu,    0,      0,   0,    0,    0, krnap,   0,    0,   0]+slack(ncol,6),
              [ 0, 0, 0, 0,   -mu,      0,   0,    0,    0,     0, kaf,    0,   0]+slack(ncol,7),
              [ 0, 0, 0, 0,     0,    -mu,   0,    0,    0,     0,   0, krnase, 0]+slack(ncol,8),
              [ 0, 0, 0, 0,kel/mu,-kel/mu,-nig,-neaa,-nent,-nrnap,-naf,-nrnase,-nrp]+slack(ncol,9),
              [ 0, 0, 0, 0, -kdeg,mu+kdeg,   0,    0,    0,     0,   0,    0,   0]+slack(ncol,10),
              [-mwc,0,0, 0,     0,      0,   0,    0,    0,     0,   0,    0,   0]+slack(ncol-1)+[mu]]
            )

        if self.matrix_type == "Kostinski":
            # take the metabolite rows from the unconstrained matrix
            met_matrix = self.matrix_uc[0:5,:]

            # add slack columns
            ncol = 11
            slacks = np.array([slack(ncol), slack(ncol), slack(ncol),
                               slack(ncol, 0), slack(ncol, 1)])
            met_matrix = np.concatenate((met_matrix, slacks), axis = 1)

            # parameters from Kostinski & Reuveni
            # ribosome allocations to RNAP, rP
            phi_rnap = [0.93, 1.14, 1.35, 1.5, 1.61, 1.66][self.medium]/100  # phi^RNAP_R, in %
            phi_rp = [7.8, 9.4, 11.8, 15.3, 19.2, 23.1][self.medium]/100  # phi^rP_R, in %

            # effective translation rates for different allocation fractions
            kel_rnap = phi_rnap*kel
            kel_rp = phi_rp*kel
            kel_rest = (1-phi_rnap-phi_rp)*kel # the rest is allocated to all other proteins

            constraints = np.array(
            [[-mu,   0,   0,   0,        0, kig,    0,    0,    0,    0,   0]+slack(ncol,2),
             [  0, -mu,   0,   0,        0,   0, keaa,    0,    0,    0,   0]+slack(ncol,3),
             [  0,   0, -mu,   0,        0,   0,    0, kent,    0,    0,   0]+slack(ncol,4),
             [  0,   0,   0, -mu,        0,   0,    0,    0, krnap,   0,   0]+slack(ncol,5),
             [  0,   0,   0,   0,      -mu,   0,    0,    0,    0,  kaf,   0]+slack(ncol,6),
             [  0,   0,   0,   0,kel_rest/mu,-nig,-neaa,-nent,  0, -naf,   0]+slack(ncol,7),
             [  0,   0,   0,   0,kel_rnap/mu, 0,    0,    0,-nrnap,   0,   0]+slack(ncol,8),
             [  0,   0,   0,   0, kel_rp/mu,  0,    0,    0,    0,    0,-nrp]+slack(ncol,9),
             [-mwc,  0,   0,   0,        0,   0,    0,    0,    0,    0,   0]+slack(ncol-1)+[mu]])

        # concatenate stoichiometric matrix with constraint matrix
        self.matrix_c = np.concatenate((met_matrix, constraints), axis=0)
        self.generate_row_and_column_names()


    def generate_row_and_column_names(self):
        """add row (rows_c) and column names (columns_c) based on the shape of matrix_c"""

        if "deg" in self.matrix_type:
            rows =  self.rows_uc[:6]
        else:
            rows =  self.rows_uc[:5]

        n_slack_row = self.matrix_c.shape[0]-len(rows)
        n_slack_col = self.matrix_c.shape[1]-len(self.columns_uc)-1

        self.rows_c = rows +  [f"C{i}" for i in range(n_slack_row)]
        self.columns_c = self.columns_uc + [f"S{i}" for i in range(n_slack_col)] + ["C"]


    def run_efmtool(self):
        """run efmtool once for the matrix with constraints (matrix_c)
        Store the output in the attribute egvs (as np.array)

        Parameters
        ----------
        drop_slack: bool, if True, slack columns are dropped from the data frame

        Return
        ------
        numpy array with fluxes
        """

        options = efmtool.get_default_options()
        options["arithmetic"] = "fractional"
        options['level'] = 'WARNING'
        egvs = efmtool.calculate_efms(stoichiometry = self.matrix_c,
                                      reversibilities = [0]*self.matrix_c.shape[1],
                                      reaction_names = self.columns_c,
                                      metabolite_names = self.rows_c,
                                      options = options,
                                      jvm_options = ["--illegal-access=deny"])
        return egvs


    def update_mu(self, growth_rate):
        """
        Update growth rate and matrices in the model.

        Parameters
        ----------
        growth_rate: new growth rate value
        """
        self.growth_rate = growth_rate
        self.make_matrix_c()
        self.make_matrix_uc()


    def update_frac(self, frac):
        """
        Update protein fraction in ribosome (frac),
        rP and rRNA stoichiometries (nrp and nrrna), then matrices.

        Parameters
        ----------
        frac: protein mass fraction in ribosome
        """
        self.frac = frac
        self.stoichiometries["nrp"] = self.mol_masses["R"]*self.frac/self.mol_masses["AA"]
        self.stoichiometries["nrrna"] = self.mol_masses["R"]*(1-self.frac)/self.mol_masses["NT"]
        self.make_matrix_c()
        self.make_matrix_uc()


    def make_nice_results(self, egvs, drop_slack = True):
        """
        Convert ugly array from efmtool into a nice pandas DataFrame
        Store as Model.nice_fluxes

        Parameters
        ----------
        egvs: numpy array with the output of "run_efmtool()"
        drop_slack: bool, if True, slack columns are removed

        Returns
        -------
        pandas Data Frame with fluxes
        """

        if egvs.shape[1] == 0:
            return pd.DataFrame()

        # solve "ValueError: Big-endian buffer not supported on little-endian compiler"
        egms = egvs.byteswap().newbyteorder()

        # remove zero flux EGVs
        # subtract one because zero EGVs have 1 for a slack variable
        nonzeros = (np.count_nonzero(egms, axis=0)-1).astype(bool)
        egms = egms.T[nonzeros]

        res = pd.DataFrame(egms,
                           index = [f"EGM{(i+1)}" for i in range(egms.shape[0])],
                           columns = self.columns_c)

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


class Simulation:
    """
    Class to run simulations with efmtool and find optimal growth rate
    and protein fraciton in the ribosome

    Attributes
    ----------
    mu_list: np.array of growth rates to test
    prot_fractions: np.array of ribosome protein mass fractions (xrps) to test
    model: an instance of Model()
    fluxes: pandas DataFrame that stores growth fluxes from test_xrps()
    growth_rates: pandas DataFrame that stores growth rates from test_xrps()
    allocations: pandas DataFrame that stores ribosome allocations from test_xrps()
    mass_fractions: pandas DataFrame that stores mass fractions from test_xrps()
    """
    def __init__(self,
                 *model_pars,
                 growth_rates = np.arange(0.1, 4, 0.1),
                 prot_fractions = np.arange(0.01, 1, 0.1),
                 **kwargs):

        self.growth_rates = growth_rates
        self.prot_fractions = prot_fractions
        self.model = Model(*model_pars, **kwargs)
        self.fluxes = pd.DataFrame()
        self.max_growth_rates = pd.DataFrame()
        self.allocations = pd.DataFrame()
        self.mass_fractions = pd.DataFrame()


    def bisection_search_mu(self):
        """quickly find the biggest possible growth rate with bisection search

        Returns
        -------
        a tuple with:
            - maximum growth rate
            - pandas DataFrame of fluxes at maximum growth rate
        """

        left = 0  # The starting index of the list we have to search in
        right = len(self.growth_rates)-1  # the last index of the list we have to search in
        mid = (right + left)//2  # // means floored division
        new_mu = self.growth_rates[mid]
        last_mu = -1
        iterations = 1

        while new_mu != last_mu: # check if we already found the last mu
            self.model.update_mu(new_mu)
            egvs = self.model.run_efmtool()

            # no solution (no egvs/zero egvs in non-slack columns) -- growth rate too big
            ncol = 11
            if "R_deg" in self.model.matrix_type:
                ncol = 13
            if (egvs.shape[1] == 0) or (np.count_nonzero(egvs[:, :ncol], axis=None) == 0):
                right = mid - 1

            # otherwise growth rate is too small (or just right)
            else:
                left = mid + 1
                last_mu = new_mu

            mid = (right + left)//2
            new_mu = self.growth_rates[mid]

            iterations += 1

            if iterations > len(self.growth_rates):
                print("Optimal growth rate not found")
                last_mu = 0

        # recalculate egvs for last_mu (the last feasible solution)
        self.model.update_mu(growth_rate=last_mu)
        egvs = self.model.run_efmtool()

        return (last_mu, self.model.make_nice_results(egvs))


    def test_xrps(self, plot=True):
        """
        Vary protein fraction in ribosome (xrp) and for each xrp calculat:
        * maximum growth rate (saved in attribute "growth_rate")
        * fluxes (EGVs - elementary growth vectors) (saved in attribute "fluxes")
        * metabolite mass fractions (saved in attribute "mass_fractions")
        * ribosome allocations (saved in attribute "allocations")

        Parameters
        ----------
        plot: bool, if True, progress is printed in the terminal
        """
        last_mus = []
        all_fluxes = pd.DataFrame()

        if plot:
            initiate_progress_plot(f"{self.model.matrix_type}, {self.model.parameter_set}")

        for prog, frac in enumerate(self.prot_fractions):

            self.model.update_frac(frac)

            last_mu, fluxes = self.bisection_search_mu()
            last_mus.append(last_mu)

            fluxes["prot_fraction"] = str(round(self.model.frac, 5))
            fluxes["growth_rate"] = last_mu
            fluxes["EGVs"] = fluxes.index
            fluxes = fluxes.set_index([fluxes["prot_fraction"]+fluxes.index])
            all_fluxes = pd.concat([all_fluxes, fluxes])

            if plot:
                fill_progress_plot(self.prot_fractions[:prog], last_mus, 0.1, 4)

        # clear plot
        if plot:
            plt.clf()

        self.fluxes = all_fluxes
        self.max_growth_rates = pd.DataFrame({"prot_fraction": self.prot_fractions,
                                              "growth_rate": last_mus})
        self.calculate_mass_fractions()
        self.calculate_allocations()


    def calculate_mass_fractions(self):
        """
        Calculates metabolite mass fractions in g/g.
            mass_fractions = (Nv/mu)*mw
            N: stoichiometric matrix without constraints
            v: vector a fluxes
            mw: vector of molecular masses
        Saves pandas DataFrame as "mass_fractions" attribute
        """

        if self.fluxes.empty:
            print("No EGVs - nothing to calculate!")

        else:
            mass_fractions = pd.DataFrame(columns = self.model.rows_uc,
                                          index = self.fluxes.index)
            for egv in self.fluxes.index:
                frac = self.fluxes.loc[egv, "prot_fraction"]
                self.model.update_frac(float(frac))

                mat = pd.DataFrame(self.model.matrix_uc,
                                   index = self.model.rows_uc,
                                   columns = self.model.columns_uc)

                # molecular masses
                mwaa = self.model.mol_masses["AA"]
                mwc = self.model.mol_masses["C"]
                mwnt = self.model.mol_masses["NT"]
                mwr = self.model.mol_masses["R"]

                stoich = self.model.stoichiometries
                if mat.shape[1] == 11:  # matrix without RNase
                    mws = [mwc, mwaa, mwnt,
                           stoich["nrrna"]*mwnt, stoich["nrp"]*mwaa,
                           stoich["nig"]*mwaa, stoich["neaa"]*mwaa,
                           stoich["nent"]*mwaa, stoich["nrnap"]*mwaa,
                           stoich["naf"]*mwaa, mwr]
                if mat.shape[1] == 13:  # matrix with RNase
                    mws = [mwc, mwaa, mwnt,
                           stoich["nrrna"]*mwnt, stoich["nrp"]*mwaa, mwr,
                           stoich["nig"]*mwaa, stoich["neaa"]*mwaa,
                           stoich["nent"]*mwaa, stoich["nrnap"]*mwaa,
                           stoich["naf"]*mwaa, stoich["nrnase"]*mwaa]

                mus = self.fluxes.loc[egv, "growth_rate"]
                fluxes = self.fluxes.loc[egv][:mat.shape[1]]
                mass_fractions.loc[egv] = mat.multiply(fluxes).sum(axis=1)/float(mus)
                mass_fractions.loc[egv] = mass_fractions.loc[egv]*mws
            mass_fractions["growth_rate"] = self.fluxes.growth_rate
            mass_fractions["EGVs"] = self.fluxes.EGVs
            mass_fractions["prot_fraction"] = self.fluxes.prot_fraction

            self.mass_fractions = mass_fractions


    def calculate_allocations(self):
        """
        Calculates ribosome allocations to the different proteins in the model
        Saved in the attribute "allocations"

        For matrix without ribosome degradation:
            phi_i = (mu*n_i*w_i)/(kel*vaf)
        For matrix with ribosome degradation
            phi_i = (mu*n_i*w_i)/(kel*(vaf-vdeg))

        where
            phi_i: ribosome allocation to protein i
            mu: growth rate
            n_i: stoichiometric coefficient of protein i
            w_i: synthesis flux of protein i
            kel: translation rate
            vaf: ribosome assembly flux
            vdeg: ribosome degradation flux
        """

        if self.fluxes.empty:
            print("No EGVs - nothing to calculate!")
        else:
            protein_columns = ["wIG", "wEAA", "wENT", "wRNAP", "wAF"]
            stoich = self.model.stoichiometries
            protein_stoichiometries = [stoich["nig"], stoich["neaa"], stoich["nent"],
                                       stoich["nrnap"], stoich["naf"]]

            if "R_deg" in self.model.matrix_type:
                protein_columns = protein_columns + ["wRNase"]
                protein_stoichiometries = protein_stoichiometries + [stoich["nrnase"]]

            prot_fluxes = self.fluxes[protein_columns].multiply(protein_stoichiometries)

            # rP stoichiometry calculated from 'prot_fraction' (rP mass fraction in ribosome)
            mwr = self.model.mol_masses["R"]
            mwaa = self.model.mol_masses["AA"]
            xrps = self.fluxes["prot_fraction"]
            prot_fluxes["wrP"] = self.fluxes["wrP"]*mwr*xrps.astype(float)/mwaa

            # multiply with growth rates and divide with assembly flux and elongation rate
            kel = self.model.kinetics["kel"]
            growth_rates = self.fluxes["growth_rate"]
            v_af = self.fluxes["vAF"]

            if "R_deg" in self.model.matrix_type:
                v_deg = self.fluxes["vdeg"]
                allocations = prot_fluxes.apply(lambda x: x*growth_rates/((v_af-v_deg)*kel))
            else:
                allocations = prot_fluxes.apply(lambda x: x*growth_rates/(v_af*kel))

            # add column with protein fractions
            allocations["prot_fraction"] = xrps

            self.allocations = allocations
