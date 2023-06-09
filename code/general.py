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

def slack(n_columns, position = -100):
    """make a list with zeroes, optionally with -1 at a specified position"""
    zeros = np.zeros(n_columns, dtype = int)
    if position in range(len(zeros)):
        zeros[position] = -1
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
    plt.plot_size(70, 25)


def fill_progress_plot(xdata, ydata, ymin=0, ymax=4):
    """
    Plot data into the plot initiated with "initiate_progress_plot()"

    Parameters
    ----------
    xdata, ydata: x and y coordinates
    ymin, ymax: min and max limits of y-axis
    """
    plt.clt()
    plt.ylim(ymin, ymax)
    plt.xlim(0,1)
    plt.scatter(xdata, ydata, marker="heart", color=125)
    plt.show()


class Model:
    """
    Class to store models and results of a single run of efmtool

    Attributes
    ----------
    parameter_set: one of ["default", "activities", "activities2", 
                           "archaea", "rna_expensive", "Kostinski"]
        - changes some parameters from defaults in "set_specific_parameters()"
    matrix_type: see funcions "make_matrix_uc", "make_matrix_c"
        - "RBA": standard RBA (default)
        - "RBA_init": RBA with RNAP capacity constraint on initiation
                      instead of elongation
        - "deg", "deg_mito", "deg_hill-[n]", "deg_hill_mito-[n]": 
                matrices with rRNA degradation reaction
                deg: constant rate, deg_hill-[n] - calculated with
                Hill function where [n] is the hill factor
        - "Kostinski" - RBA with fixed ribosome and RNAP allocations
    medium: influences the values of keaa, mol_masses["C"], ktr (if parameter_set = "activities")
        - 0: succinate minimal medium
        - 1: glycerol minimal medium
        - 2: glucose minimal medium [default]
        - 3: glycerol + amino acids
        - 4: glucose + amino acids
        - 5: LB
    frac: protein mass fraction in the ribosome
    mu: growth rate
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
        - keaa - kcat of AA synthesis, depends on the value of "medium"
        - ktr: transcription rate [nt/h]
        - krnap: transcription rate [1/h]
        - kel: translation rate [aa/h]
        - kexo: exonuclease rate (RNase R) [nt/h]
    kdeg_max: max. degradation rate [1/h]
    kin_max: max. transcription initiation rate [1/h]
    matrix_uc: np.array - matrix without constraints
    matrix_c: np.array - matrix with constraints
    rows_uc, rows_c: row names for unconstrained / constrained matrices
    columns_uc, columns_c: column names for unconstrained / constrained matrices
    """
    
    # allowed attributes
    __slots__ = ["parameter_set", "matrix_type", "medium", "frac", 
                 "growth_rate", "kin_max", "kdeg_max", 
                 "stoichiometries", "mol_masses", "kinetics", 
                 "rows_c", "columns_c", "matrix_c", 
                 "rows_uc", "columns_uc", "matrix_uc"]

    def __init__(self, *model_pars, **kwargs):
        self.parameter_set = "default"
        self.matrix_type = "RBA"
        self.medium = 2
        self.frac = 0.36
        self.growth_rate = 1
        self.kin_max = 1000
        self.kdeg_max = 0

        avg_protein_length = 325
        self.stoichiometries = {"nig": 646,
                                "neaa": avg_protein_length*15,
                                "nent": avg_protein_length*15,
                                "nrnap": 3498,
                                "naf": avg_protein_length*12,
                                "nrnase": 813}

        self.mol_masses = {"R": 2300,
                           "AA": 0.109,
                           "NT": 0.3243}

        self.kinetics = {"kig": 180*TO_H,
                         "kent": 10*TO_H,
                         "kaf": 1/120*TO_H,
                         "ktr": 85*TO_H,
                         "kel": 21*TO_H,
                         "kexo": 88*TO_H}

        # if model_pars / kwargs are provided, update default parameters
        for dictionary in model_pars:
            for key, value in dictionary.items():
                setattr(self, key, value)
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.set_specific_parameters(*model_pars, **kwargs)

        # calculated stoichiometries - must be placed after all other parameters are set
        self.stoichiometries["naa"] = self.mol_masses["AA"] / self.mol_masses["C"]
        self.stoichiometries["nnt"] = ((self.mol_masses["NT"] - self.mol_masses["AA"])
                                       / self.mol_masses["C"])
        self.stoichiometries["nrp"] = self.mol_masses["R"] * self.frac / self.mol_masses["AA"]
        self.stoichiometries["nrrna"] = (self.mol_masses["R"] * (1-self.frac)
                                         / self.mol_masses["NT"])

        # initialize row and column names of matrices
        self.rows_c = []
        self.columns_c = []

        # create matrices
        self.make_matrix_uc()
        self.make_matrix_c()

        _ = self.calculate_molecular_masses()


    def set_specific_parameters(self, *model_pars, **kwargs):
        """
        change a specific set of parameters according to values
        of "Model.medium" and "Model.parameter_set"
        """
        # check if parameter set/medium are valid
        if self.parameter_set not in ["default", "activities", "activities2", 
                                      "archaea", "rna_expensive", "Kostinski"]:
            print(f"Unknown parameter set: {self.parameter_set}, using default parameters")

        if self.medium not in [0,1,2,3,4,5]:
            print(f"Unknown medium: {self.medium}, using default (2: glucose minimal)")
            self.medium = 2

        # medium specific parameters - always set
        # kcat of AA synthesis, 1/h
        self.kinetics["keaa"] = [2, 3.5, 5, 7, 8.5, 10.5][self.medium]*TO_H

        # molecular mass of the carbon source
        # 1,3: glycerol; 2,4,5: glucose; 0: succinate
        carbon_sources = {5: 0.18, 4: 0.18, 3: 0.092,
                          2: 0.18, 1: 0.092, 0: 0.118}
        self.mol_masses["C"] = carbon_sources[self.medium]

        # active fractions of R and RNAP
        active_ribosome_fraction = 0.85
        active_rnap_fraction = [13.2, 14.4, 15.0, 18.8, 24.2, 31.0][self.medium]/100

        # RNAP allocation fractions
        rnap_allocation_rrna = [0.18, 0.28, 0.42, 0.52, 0.60, 0.65][self.medium]

        if self.parameter_set == "archaea":
            self.kinetics["ktr"] = 25*TO_H
            self.kinetics["kel"] = 25/3*TO_H
            self.stoichiometries["nrnap"] = 3338
            self.mol_masses["R"] = 3040
            
            # activities from E. coli
            self.kinetics["kel"] *= active_ribosome_fraction
            self.kinetics["ktr"] *= active_rnap_fraction

        if self.parameter_set == "rna_expensive":
            self.stoichiometries["nrnap"] = 3498*15
            self.kinetics["kel"] = 21*TO_H*3
            self.kinetics["ktr"] = 85/10*TO_H

        if "activities" in self.parameter_set:
            self.kinetics["kel"] *= active_ribosome_fraction
            self.kinetics["ktr"] *= active_rnap_fraction

            if self.parameter_set == "activities2":
                self.kinetics["ktr"] *= rnap_allocation_rrna

        if self.parameter_set == "Kostinski":
            # translation rates from Kostinski & Reuveni
            kel = [12, 16.83, 21, 20.17, 21, 22.25][self.medium]*TO_H
            self.kinetics["kel"] = kel*active_ribosome_fraction
            self.kinetics["ktr"] *= active_rnap_fraction*rnap_allocation_rrna


    def calculate_molecular_masses(self):
        """generate a list of molecular masses based on matrix row names"""
        mw_aa = self.mol_masses["AA"]
        mw_nt = self.mol_masses["NT"]

        self.mol_masses["IG"] = self.stoichiometries["nig"]*mw_aa
        self.mol_masses["EAA"] = self.stoichiometries["neaa"]*mw_aa
        self.mol_masses["ENT"] = self.stoichiometries["nent"]*mw_aa
        self.mol_masses["RNAP"] = self.stoichiometries["nrnap"]*mw_aa
        self.mol_masses["RNase"] = self.stoichiometries["nrnase"]*mw_aa
        self.mol_masses["AF"] = self.stoichiometries["naf"]*mw_aa
        self.mol_masses["rP"] = self.stoichiometries["nrp"]*mw_aa
        self.mol_masses["rRNA"] = self.stoichiometries["nrrna"]*mw_nt
        
        mw_list = [self.mol_masses[metabolite] for metabolite in self.rows_uc]
        
        return mw_list

            
    def caluclate_degradation_rate(self):
        """
        Calculate rRNA degradation rate.
        The rate decreases with protein fraction in the ribosome (frac)
        * either linearly (matrix_type="deg"),
        * or according to a Hill function (matrix_type="deg_hill-[n]")
          where n is the Hill factor

        Returns
        -------
        a float with rRNA degradation rate
        """
        kdeg = self.kdeg_max * (1 - self.frac)

        if "deg_hill" in self.matrix_type:
            hill = int(self.matrix_type.split("-")[-1])
            kdeg = (self.kdeg_max
                    * (1 - (self.frac**hill / (self.frac**hill + 0.2**hill)))
                    * (1-self.frac))

        return kdeg


    def __str__(self):
        return f"""
        Parameter set: {self.parameter_set}
        Matrix type: {self.matrix_type}
        Medium: {self.medium}
        """


    def make_matrix_uc(self):
        """
        Create stoichiometric matrix without constraints (numpy array stored as "matrix_uc"),
        row names (list stored as "rows_uc") and column names (list stored as "columns_uc")
        based on the value of "matrix_type" attribute
        """
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
        self.rows_uc = ["C", "AA", "NT", "rRNA", "rP",
                        "IG", "EAA", "ENT", "RNAP", "AF", "R"]

        # matrix with ribosome degradation rate
        if "deg" in self.matrix_type:
            self.matrix_uc = np.array(
             [[1,-naa,-nnt,     0,  0,   0,   0,    0,    0,     0,    0,     0,   0],
              [0,   1,  -1,     0,  0,   0,-nig,-neaa,-nent,-nrnap,-nrnase,-naf,-nrp],
              [0,   0,   1,-nrrna,nrrna, 0,   0,    0,    0,     0,    0,     0,   0],
              [0,   0,   0,     1, -1,  -1,   0,    0,    0,     0,    0,     0,   0],
              [0,   0,   0,     0,  0,  -1,   0,    0,    0,     0,    0,     0,   1],
              [0,   0,   0,     0,  0,   1,   0,    0,    0,     0,    0,     0,   0],
              [0,   0,   0,     0,  0,   0,   1,    0,    0,     0,    0,     0,   0],
              [0,   0,   0,     0,  0,   0,   0,    1,    0,     0,    0,     0,   0],
              [0,   0,   0,     0,  0,   0,   0,    0,    1,     0,    0,     0,   0],
              [0,   0,   0,     0,  0,   0,   0,    0,    0,     1,    0,     0,   0],
              [0,   0,   0,     0,  0,   0,   0,    0,    0,     0,    1,     0,   0],
              [0,   0,   0,     0,  0,   0,   0,    0,    0,     0,    0,     1,   0]]
            )
            self.columns_uc = ["vIG", "vEAA", "vENT", "vRNAP", "vRNase", "vAF",
                               "wIG", "wEAA", "wENT", "wRNAP", "wRNase", "wAF", "wrP"]
            self.rows_uc = ["C", "AA", "NT", "rRNA", "rP", "R",
                            "IG", "EAA", "ENT", "RNAP", "RNase", "AF"]
            
        # for mitochondria, add a rP imoprt reaction
        if "mito" in self.matrix_type:
            rp_import = np.array([[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]]).T
            self.matrix_uc = np.concatenate((self.matrix_uc, rp_import), axis=1)
            self.columns_uc.append("vIrP")


    def make_matrix_c(self):
        """
        Create stoichiometric matrix without constraints (numpy array stored as "matrix_c"),
        row names (list stored as "rows_c") and column names (list stored as "columns_c")
        based on the value of "matrix_type" attribute
        """
        # get parameters
        mu = self.growth_rate
        nig, neaa, = self.stoichiometries["nig"], self.stoichiometries["neaa"]
        nent, nrnap = self.stoichiometries["nent"], self.stoichiometries["nrnap"]
        naf, nrp = self.stoichiometries["naf"], self.stoichiometries["nrp"]
        nrrna, nrnase = self.stoichiometries["nrrna"], self.stoichiometries["nrnase"]
        kig, keaa = self.kinetics["kig"], self.kinetics["keaa"]
        kent, kaf, kel = self.kinetics["kent"], self.kinetics["kaf"], self.kinetics["kel"]
        
        mwc = self.mol_masses["C"]
        krnap = self.kinetics["ktr"]/nrrna  # convert nt/h -> 1/h
        if "_init_rate" in self.matrix_type:
            krnap = self.kin_max  # if initiation rate is limiting

        # get submatrix with metabolite rows (met_matrix) and add slack columnes
        ncol = 9
        if ("deg" in self.matrix_type) or (self.matrix_type == "Kostinski"):
            ncol = 11
    
        slacks = np.array([slack(ncol), slack(ncol), slack(ncol),
                          slack(ncol,0), slack(ncol,1)])
        met_matrix = np.concatenate((self.matrix_uc[0:5, :], slacks), axis=1)

        # "RBA" matrix - enzyme capacity + dry mass constraints
        if self.matrix_type in ["RBA", "RBA_init_rate"]:
            constraints = np.array(
              [[-mu,  0,  0,  0,     0, kig,    0,    0,     0,   0,   0]+slack(ncol,2),
               [  0,-mu,  0,  0,     0,   0, keaa,    0,     0,   0,   0]+slack(ncol,3),
               [  0,  0,-mu,  0,     0,   0,    0, kent,     0,   0,   0]+slack(ncol,4),
               [  0,  0,  0,-mu,     0,   0,    0,    0, krnap,   0,   0]+slack(ncol,5),
               [  0,  0,  0,  0,   -mu,   0,    0,    0,     0, kaf,   0]+slack(ncol,6),
               [  0,  0,  0,  0,kel/mu,-nig,-neaa,-nent,-nrnap,-naf,-nrp]+slack(ncol,7),
               [-mwc, 0,  0,  0,     0,   0,    0,    0,    0,   0,    0]+slack(ncol-1)+[mu]]
            )

        # "Kostinski" matrix - like RBA but with fixed allocations of ribosome
        if self.matrix_type == "Kostinski":
            # ribosome allocation to RNAP and rP (%) from Kostinski & Reuveni
            r_allocation_rnap = [0.93, 1.14, 1.35, 1.5, 1.61, 1.66][self.medium]/100
            r_allocation_rp = [7.8, 9.4, 11.8, 15.3, 19.2, 23.1][self.medium]/100

            # translation rate "kel" is set in "set_specific_parameters"
            # the ribsome constraint is split to 3 rows, where "kel" is multiplied 
            # by the respective allocation fraction for RNAP, rP and the rest
            kel_rnap = r_allocation_rnap*kel
            kel_rp = r_allocation_rp*kel
            kel_rest = (1-r_allocation_rnap-r_allocation_rp)*kel # all other proteins

            constraints = np.array(
            [[-mu,   0,   0,   0,        0,  kig,   0,    0,    0,   0,   0]+slack(ncol,2),
             [  0, -mu,   0,   0,        0,   0,  keaa,   0,    0,   0,   0]+slack(ncol,3),
             [  0,   0, -mu,   0,        0,   0,    0,  kent,   0,   0,   0]+slack(ncol,4),
             [  0,   0,   0, -mu,        0,   0,    0,    0, krnap,  0,   0]+slack(ncol,5),
             [  0,   0,   0,   0,      -mu,   0,    0,    0,    0, kaf,   0]+slack(ncol,6),
             [  0,   0,   0,   0,kel_rest/mu,-nig,-neaa,-nent,  0,-naf,   0]+slack(ncol,7),
             [  0,   0,   0,   0,kel_rnap/mu, 0,    0,    0,-nrnap,  0,   0]+slack(ncol,8),
             [  0,   0,   0,   0, kel_rp/mu,  0,    0,    0,    0,   0,-nrp]+slack(ncol,9),
             [-mwc,  0,   0,   0,        0,   0,    0,    0,    0,   0,   0]+slack(ncol-1)+[mu]])


        # constraints like in RBA + minimum rRNA degradation rate enforced
        if "deg" in self.matrix_type:
            krnase = self.kinetics["kexo"]/nrrna
            kdeg = self.caluclate_degradation_rate()

            constraints = np.array(
             [[-mu,0, 0, 0,  0,     0, kig,    0,    0,     0,    0,    0,   0]+slack(ncol,2),
              [ 0,-mu,0, 0,  0,     0,   0, keaa,    0,     0,    0,    0,   0]+slack(ncol,3),
              [ 0, 0,-mu,0,  0,     0,   0,    0, kent,     0,    0,    0,   0]+slack(ncol,4),
              [ 0, 0, 0,-mu, 0,     0,   0,    0,    0, krnap,    0,    0,   0]+slack(ncol,5),
              [ 0, 0, 0, 0,-mu,     0,   0,    0,    0,     0, krnase,  0,   0]+slack(ncol,6),
              [ 0, 0, 0, 0,  0,   -mu,   0,    0,    0,     0,    0,   kaf,  0]+slack(ncol,7),
              [ 0, 0, 0, 0,  0,kel/mu,-nig,-neaa,-nent,-nrnap,-nrnase,-naf,-nrp]+slack(ncol,8),
              [ 0, 0, 0, 0, mu, -kdeg,   0,    0,    0,     0,    0,    0,   0]+slack(ncol,9),
              [-mwc,0,0, 0,  0,     0,   0,    0,    0,     0,    0,    0,   0]+slack(ncol-1)+[mu]]
            )
            
        # 1/3 of rP imported for free
        if "mito" in self.matrix_type:
            mwrp = self.stoichiometries["nrp"]*self.mol_masses["AA"]
            
            # add column for rP import reaction
            column =[0, 0, 0, 0, 0, 0, 0, 0, -mwrp]
            constraints = np.insert(constraints, self.matrix_uc.shape[1]-1, column, axis=1)

            # wrP = vIrP constraint
            wrP_vIrP_constraint = np.array(
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2]+slack(ncol)])
            constraints = np.concatenate((constraints, wrP_vIrP_constraint), axis = 0)

        # concatenate stoichiometric matrix with constraint matrix
        self.matrix_c = np.concatenate((met_matrix, constraints), axis=0)
        self.generate_row_and_column_names()


    def generate_row_and_column_names(self):
        """add row (rows_c) and column names (columns_c) based on the shape of matrix_c"""

        rows =  self.rows_uc[:5]
        n_slack_row = self.matrix_c.shape[0]-len(rows)
        n_slack_col = self.matrix_c.shape[1]-len(self.columns_uc)-1

        self.rows_c = rows +  [f"C{i}" for i in range(n_slack_row)]
        self.columns_c = self.columns_uc + [f"S{i}" for i in range(n_slack_col)] + ["C"]


    def run_efmtool(self):
        """
        run efmtool for matrix with constraints (matrix_c), 
        return numpy array with fluxes
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
        """Update matrices with new growth rate"""
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
    __slots__ = ["growth_rates", "prot_fractions", "model", "fluxes", 
                 "max_growth_rates", "allocations", "mass_fractions"]

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
        mid = (right + left)//2
        new_mu = self.growth_rates[mid]
        last_mu = -1
        iterations = 1

        while new_mu != last_mu: # check if we already found the last mu
            self.model.update_mu(new_mu)
            egvs = self.model.run_efmtool()

            # no solution (no egvs/zero egvs in non-slack columns) - growth rate too big
            ncol = self.model.matrix_uc.shape[1]
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
            title = (f"{self.model.matrix_type}, "
                     f"parameter set: {self.model.parameter_set}, "
                     f"medium: {self.model.medium}")
            initiate_progress_plot(title)

        for progress, frac in enumerate(self.prot_fractions):
            self.model.update_frac(frac)

            last_mu, fluxes = self.bisection_search_mu()
            last_mus.append(last_mu)

            fluxes["prot_fraction"] = str(round(self.model.frac, 5))
            fluxes["growth_rate"] = last_mu
            fluxes["EGVs"] = fluxes.index
            fluxes = fluxes.set_index([fluxes["prot_fraction"]+fluxes.index])
            all_fluxes = pd.concat([all_fluxes, fluxes])

            if plot:
                fill_progress_plot(self.prot_fractions[:progress], last_mus, 0.1, 4)

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
        Calculates metabolite mass fractions in g/g according to:
        mass_fractions = (Nv/mu)*mw
            N: stoichiometric matrix without constraints ("matrix_uc")
            v: vector a fluxes ("fluxes")
            mu: growth rate ("growth_rates")
            mw: vector of molecular masses ("molecular_masses")
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

                matrix_uc = pd.DataFrame(self.model.matrix_uc,
                                         index = self.model.rows_uc,
                                         columns = self.model.columns_uc)
                molecular_masses = self.model.calculate_molecular_masses()
                growth_rates = self.fluxes.loc[egv, "growth_rate"]
                fluxes = self.fluxes.loc[egv][:matrix_uc.shape[1]]

                # calculate mass fractions with (Nv/mu)*mw
                mass_fractions.loc[egv] = matrix_uc.multiply(fluxes).sum(axis=1)/float(growth_rates)
                mass_fractions.loc[egv] = mass_fractions.loc[egv]*molecular_masses

            mass_fractions["growth_rate"] = self.fluxes.growth_rate
            mass_fractions["EGVs"] = self.fluxes.EGVs
            mass_fractions["prot_fraction"] = self.fluxes.prot_fraction

            self.mass_fractions = mass_fractions


    def calculate_allocations(self):
        """
        Calculates ribosome allocations to the different proteins with:

        allocation_i = (mu*n_i*w_i)/(kel*vaf)
            allocation_i: ribosome allocation to protein i
            mu: growth rate
            n_i: stoichiometric coefficient of protein i
            w_i: synthesis flux of protein i
            kel: translation rate
            vaf: ribosome assembly flux
        Saved in the attribute "allocations"
        """

        if self.fluxes.empty:
            print("No EGVs - nothing to calculate!")

        else:
            protein_columns = ["wIG", "wEAA", "wENT", "wRNAP", "wAF"]
            stoich = self.model.stoichiometries
            protein_stoichiometries = [stoich["nig"], stoich["neaa"], stoich["nent"],
                                       stoich["nrnap"], stoich["naf"]]

            if "deg" in self.model.matrix_type:
                protein_columns = protein_columns + ["wRNase"]
                protein_stoichiometries = protein_stoichiometries + [stoich["nrnase"]]

            prot_fluxes = self.fluxes[protein_columns].multiply(protein_stoichiometries)

            # rP stoichiometry calculated from 'prot_fraction' (rP mass fraction in ribosome)
            mwr = self.model.mol_masses["R"]
            mwaa = self.model.mol_masses["AA"]
            prot_fractions = self.fluxes["prot_fraction"]
            prot_fluxes["wrP"] = self.fluxes["wrP"]*mwr*prot_fractions.astype(float)/mwaa

            # multiply with growth rates and divide with assembly flux and elongation rate
            kel = self.model.kinetics["kel"]
            growth_rates = self.fluxes["growth_rate"]
            v_af = self.fluxes["vAF"]
            allocations = prot_fluxes.apply(lambda x: x*growth_rates/(v_af*kel))

            allocations["prot_fraction"] = prot_fractions

            self.allocations = allocations
