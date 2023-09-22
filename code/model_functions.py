"""
General classes and functions used by other scripts
* Model definitions
* Simulations with efmtool
* Calculations of mass fractions and ribosome allocations

Author: Diana Szeliova, 20.9.2023
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
    parameter_set: identifier of the parameter set
                   if there is a substring in the form "hill-[n]" it changes
                   how degradation rate of RNA is calculated (n is an integer)
    matrix_type: see funcions "make_matrix_uc", "make_matrix_c"
        - "base": standard RBA (default)
        - "extended": matrices with rRNA degradation reaction
          "extended_mito": the same as "extended" but 1/3 rPs can be imported for free
        - "Kostinski": "base" matrix with fixed ribosome and RNAP allocations
    frac: protein mass fraction in the ribosome
    mu: growth rate
    parameters: dictionary of parameters:
        - R, AA, NT, C - molecular masses (in g/mmol) of ribosome (R),
          average amino acid (AA), average nucleotide (NT), carbon source (C)
        - nig, neaa, nent, nrnap, naf, nrnase - stoichiometric coefficients of
          transporter (nig), AA/NT synthesis enzymes (neaa/nent), RNAP (nrnap),
          ribosome assembly factors (naf) and RNase R (nrnase)
        - naa, nnt - stoich. coefficients for synthesis of AA, NT,
                     depend on the molecular mass of carbon source (parameter C)
        - nrp, nrrna - stoich. coefficients for synthesis of rP and rRNA, depend on parameter "frac"
        - kig: kcat of substrate importer [1/h]
        - kent: kcat of NT synthesis [1/h]
        - kaf: kcat of R assembly [1/h]
        - keaa - kcat of AA synthesis [1/h]
        - ktr: transcription rate [nt/h]
        - kel: translation rate [aa/h]
        - kexo: exonuclease rate (RNase R) [nt/h]
        - kdeg_max: max. rRNA degradation rate [1/h]
        - f_R - fraction of active ribosome
        - f_RNAP - fraction of active RNA polymerase
        - phi_RNAP - fraction of RNAP allocated to rRNA
        - phi_R_RNAP - fraction of ribosome allocated to RNAP production
        - phi_R_rP - fraction of ribosome allocated to rP production
    matrix_uc: np.array - matrix without constraints
    matrix_c: np.array - matrix with constraints
    rows_uc, rows_c: row names for unconstrained / constrained matrices
    columns_uc, columns_c: column names for unconstrained / constrained matrices
    """

    def __init__(self, model_pars):
        self.frac = 0.36
        self.growth_rate = 1
        self.parameters = model_pars.copy()

        # calculate effective transcription and translation rates
        self.parameters["kel"] *= self.parameters["f_R"]*TO_H
        self.parameters["ktr"] *= (self.parameters["f_RNAP"]*self.parameters["phi_RNAP"])*TO_H

        # convert rates to 1/h
        self.parameters["kig"] *= TO_H
        self.parameters["kent"] *= TO_H
        self.parameters["keaa"] *= TO_H
        self.parameters["kaf"] *= TO_H
        self.parameters["kexo"] *= TO_H

        self.calculate_degradation_rate()

        # calculated stoichiometries - must be placed after all other parameters are set
        self.parameters["naa"] = self.parameters["AA"] / self.parameters["C"]
        self.parameters["nnt"] = ((self.parameters["NT"] - self.parameters["AA"])
                                  / self.parameters["C"])
        self.parameters["nrp"] = self.parameters["R"] * self.frac / self.parameters["AA"]
        self.parameters["nrrna"] = (self.parameters["R"] * (1-self.frac)
                                 / self.parameters["NT"])

        # initialize row and column names of matrices
        self.rows_c = []
        self.columns_c = []

        # create matrices
        self.make_matrix_uc()
        self.make_matrix_c()

        _ = self.calculate_molecular_masses()


    def calculate_degradation_rate(self):
        """
        calculate degradation rate based on protein content in the ribosome (frac)
        and type of degradation function (parameter_set)
        """
        kdeg = self.parameters["kdeg_max"] * (1 - self.frac)
        if "hill" in self.parameters["parameter_set"]:
            hill = int(self.parameters["parameter_set"].split("_")[0].split("-")[-1])
            kdeg = (self.parameters["kdeg_max"]
                    * (1 - (self.frac**hill / (self.frac**hill + 0.2**hill)))
                    * (1 - self.frac))
        self.parameters["kdeg"] = kdeg


    def calculate_molecular_masses(self):
        """generate a list of molecular masses based on matrix row names"""
        mw_aa = self.parameters["AA"]

        self.parameters["IG"] = self.parameters["nig"]*mw_aa
        self.parameters["EAA"] = self.parameters["neaa"]*mw_aa
        self.parameters["ENT"] = self.parameters["nent"]*mw_aa
        self.parameters["RNAP"] = self.parameters["nrnap"]*mw_aa
        self.parameters["RNase"] = self.parameters["nrnase"]*mw_aa
        self.parameters["AF"] = self.parameters["naf"]*mw_aa
        self.parameters["rP"] = self.parameters["nrp"]*mw_aa
        self.parameters["rRNA"] = self.parameters["nrrna"]*self.parameters["NT"]

        mw_list = [self.parameters[metabolite] for metabolite in self.rows_uc]

        return mw_list


    def make_matrix_uc(self):
        """
        Create stoichiometric matrix without constraints (numpy array stored as "matrix_uc"),
        row names (list stored as "rows_uc") and column names (list stored as "columns_uc")
        based on the value of "matrix_type" attribute
        """
        naa, nnt, = self.parameters["naa"], self.parameters["nnt"]
        nig, neaa, = self.parameters["nig"], self.parameters["neaa"]
        nent, nrnap = self.parameters["nent"], self.parameters["nrnap"]
        naf, nrp = self.parameters["naf"], self.parameters["nrp"]
        nrrna, nrnase = self.parameters["nrrna"], self.parameters["nrnase"]

        # "base" matrix - without rRNA degradation rate
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

        # matrix with rRNA degradation rate
        if "extended" in self.parameters["matrix_type"]:
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
        if "mito" in self.parameters["matrix_type"]:
            rp_import = np.array([[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]]).T
            self.matrix_uc = np.concatenate((self.matrix_uc, rp_import), axis=1)
            self.columns_uc.append("vIrP")

        valid_matrix_types = ["base", "extended", "extended_mito", "Kostinski"]
        if self.parameters["matrix_type"] not in valid_matrix_types:
            raise Exception(f"Choose valid matrix type {valid_matrix_types}")


    def make_matrix_c(self):
        """
        Create stoichiometric matrix without constraints (numpy array stored as "matrix_c"),
        row names (list stored as "rows_c") and column names (list stored as "columns_c")
        based on the value of "matrix_type" attribute
        """
        # get parameters
        matrix_type = self.parameters["matrix_type"]
        mu = self.growth_rate
        nig, neaa, = self.parameters["nig"], self.parameters["neaa"]
        nent, nrnap = self.parameters["nent"], self.parameters["nrnap"]
        naf, nrp = self.parameters["naf"], self.parameters["nrp"]
        nrrna, nrnase = self.parameters["nrrna"], self.parameters["nrnase"]
        kig, keaa = self.parameters["kig"], self.parameters["keaa"]
        kent, kaf, kel = self.parameters["kent"], self.parameters["kaf"], self.parameters["kel"]
        mwc = self.parameters["C"]
        krnap = self.parameters["ktr"]/nrrna  # convert nt/h -> 1/h

        # get submatrix with metabolite rows (met_matrix) and add slack columnes
        ncol = 9
        if ("extended" in matrix_type) or (matrix_type == "Kostinski"):
            ncol = 11
        # add slack columns (if only zeros, metabolite can't accumulate, if there is a "1", it can)
        slacks = np.array([slack(ncol), slack(ncol), slack(ncol),
                          slack(ncol,0), slack(ncol,1)])
        met_matrix = np.concatenate((self.matrix_uc[0:5, :], slacks), axis=1)

        # "base" matrix - enzyme capacity + dry mass constraints
        if matrix_type in ["base"]:
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
        if matrix_type == "Kostinski":
            # the ribsome constraint is split to 3 rows, where "kel" is multiplied
            # by the respective allocation fraction for RNAP, rP and the rest
            kel_rnap = self.parameters["phi_R_RNAP"]*kel
            kel_rp = self.parameters["phi_R_rP"]*kel
            kel_rest = (1-self.parameters["phi_R_RNAP"]-self.parameters["phi_R_rP"])*kel

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


        # constraints like in "base" matrix + minimum rRNA degradation rate enforced
        if "extended" in matrix_type:
            krnase = self.parameters["kexo"]/nrrna
            kdeg = self.parameters["kdeg"]

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
        if "mito" in matrix_type:
            mwrp = self.parameters["nrp"]*self.parameters["AA"]

            # add column for rP import reaction
            column =[0, 0, 0, 0, 0, 0, 0, 0, -mwrp]
            constraints = np.insert(constraints, self.matrix_uc.shape[1]-1, column, axis=1)

            # wrP = 2*vIrP constraint - 1/3 rPs are imported
            import_constraint = np.array(
                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -2]+slack(ncol)])
            constraints = np.concatenate((constraints, import_constraint), axis = 0)

        valid_matrix_types = ["base", "extended", "extended_mito", "Kostinski"]
        if matrix_type not in valid_matrix_types:
            raise Exception(f"Choose valid matrix type {valid_matrix_types}")


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
                                      options = options)
        return egvs


    def update_mu(self, growth_rate):
        """Update matrices with new growth rate"""
        self.growth_rate = growth_rate
        self.make_matrix_c()
        self.make_matrix_uc()


    def update_frac(self, frac):
        """
        Update protein fraction in ribosome (frac),
        rP and rRNA stoichiometries (nrp and nrrna),
        degradation rate, and then matrices.

        Parameters
        ----------
        frac: protein mass fraction in ribosome
        """
        self.frac = frac
        self.parameters["nrp"] = self.parameters["R"]*self.frac/self.parameters["AA"]
        self.parameters["nrrna"] = self.parameters["R"]*(1-self.frac)/self.parameters["NT"]
        self.calculate_degradation_rate()
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
    growth_rates: np.array of growth rates to test
    prot_fractions: np.array of ribosome protein mass fractions (xrPs) to test
    model: an instance of Model()
    fluxes: pandas DataFrame that stores growth fluxes from test_xrps()
    max_growth_rates: pandas DataFrame that stores growth rates from test_xrps()
    allocations: pandas DataFrame that stores ribosome allocations from test_xrps()
    mass_fractions: pandas DataFrame that stores mass fractions from test_xrps()
    """
    __slots__ = ["growth_rates", "prot_fractions", "model", "fluxes",
                 "max_growth_rates", "allocations", "mass_fractions"]

    def __init__(self,
                 model_pars,
                 growth_rates = np.arange(0.1, 4, 0.1),
                 prot_fractions = np.arange(0.01, 1, 0.1)):

        self.growth_rates = growth_rates
        self.prot_fractions = prot_fractions
        self.model = Model(model_pars)
        self.fluxes = pd.DataFrame()
        self.max_growth_rates = pd.DataFrame()
        self.allocations = pd.DataFrame()
        self.mass_fractions = pd.DataFrame()


    def bisection_search_mu(self):
        """quickly find the highest possible growth rate with bisection search

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
        Vary protein fraction in ribosome (xrP) and for each xrP calculate:
        * maximum growth rate (saved in attribute "max_growth_rates")
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
            title = (f"{self.model.parameters['matrix_type']}, "
                     f"parameter set: {self.model.parameters['parameter_set']}, "
                     f"name: {self.model.parameters['name']}")
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
            protein_stoichiometries = [self.model.parameters["nig"],
                                       self.model.parameters["neaa"],
                                       self.model.parameters["nent"],
                                       self.model.parameters["nrnap"],
                                       self.model.parameters["naf"]]

            if "extended" in self.model.parameters["matrix_type"]:
                protein_columns = protein_columns + ["wRNase"]
                protein_stoichiometries += [self.model.parameters["nrnase"]]

            prot_fluxes = self.fluxes[protein_columns].multiply(protein_stoichiometries)

            # rP stoichiometry calculated from 'prot_fraction' (rP mass fraction in ribosome)
            mwr = self.model.parameters["R"]
            mwaa = self.model.parameters["AA"]
            prot_fractions = self.fluxes["prot_fraction"]
            prot_fluxes["wrP"] = self.fluxes["wrP"]*mwr*prot_fractions.astype(float)/mwaa

            # multiply with growth rates and divide with assembly flux and elongation rate
            kel = self.model.parameters["kel"]
            growth_rates = self.fluxes["growth_rate"]
            v_af = self.fluxes["vAF"]
            allocations = prot_fluxes.apply(lambda x: x*growth_rates/(v_af*kel))

            allocations["prot_fraction"] = prot_fractions
            self.allocations = allocations
