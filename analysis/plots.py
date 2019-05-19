# Handling ROOT files
import os
import numpy as np
import uproot
import pandas as pd
# Plotting
import matplotlib
from matplottery import Hist1D, plot_stack
import matplotlib.pyplot as plt
# Custom
from tools import GetData, CombineYears
from config import config

class Plots():

    """A set of plotting tools wrapping the matplottery library written by Nick Amin
        for use in the H --> rho/phi+gamma analysis and other related HEP analyses

    Parent to Analysis class

    Attributes
    ----------
    __colors : dict{ str: list[int, int, int] }
        dictionary of RGB values organized by sample name
    dfs : dict{ str: pandas.DataFrame }
        dictionary of pandas dataframes organized by sample
        name
    signal : str
        name of signal dataset
    analysis : str
        keyword to delineate analysis if multiple
    outDir : str, optional
        integrated luminosity (default 137.0)
    lumi : float, optional
        integrated luminosity (default 137.0)
    verbose : bool, optional
        verbosity of module (default False)

    Methods
    -------
    _initColors()
        Initialize available colors
    _initOutput()
        Initialize output directory
    MakeCut(cut)
        Make cut on a query string or list of query strings, then update dataframes
    Stacked(col, nBins, xMin, xMax, xLabel="", extra="", logY=False,
            sigWeight=None, overflow=True)
        Make plot of stacked 1D histograms of MC and data
    """

    # Internal
    self.__colors = {}

    def __init__(self, config, dfs, outDir="plots/", lumi=137.0, verbose=False):
        """
        Parameters
        ----------
        config : dict
            Configuration for analysis
        dfs : dict{ str: pandas.DataFrame }
            Dataframe used for testing BDT
        outDir : str, optional
            integrated luminosity (default 137.0)
        lumi : float, optional
            integrated luminosity (default 137.0)
        verbose : bool, optional
            verbosity of module (default False)
        """
        # Public
        self.dfs = dfs
        self.signal = config["signal"] 
        self.analysis = config["analysis"]
        self.outDir = outDir
        self.lumi = lumi
        self.verbose = verbose
        # Init functions
        self._initColors()
        self._initOutput()

    def _initColors(self):
        """Initialize available colors"""
        # List of RGB values
        rgbs=[[0.4, 0.4, 0.4], [1.0, 0.4, 1.0], [0.4, 0.0, 0.8],
              [0.4, 0.6, 1.0], [0.4, 0.8, 0.4], [0.9, 0.9, 0.9],
              [0.4, 0.4, 0.6], [0.0, 0.4, 0.0]
             ]
        # Assign relevant colors
        for i, name in enumerate(self.dfs.keys()):
            if name == self.signal:
                color = [1.0,0.0,0.0]
            elif name == "data":
                color = [1.0,1.0,1.0]
            else:
                color =  rgbs[i]
            self.__colors[name] = color

        return

    def _initOutput(self):
        """Initialize output directory"""
        out = self.outDir
        if not os.path.isdir(out): os.mkdir(out)
        if out[-1] != "/":
            self.outDir += "/"

        return

    def MakeCut(cut):
        """Make cut on a query string or list of query strings, then update dataframes

        Parameters
        ----------
        cut : str or list of str
            Query string usable in pandas.DataFrame.query()

        Raises
        ------
        TypeError
            If cut is not a string or list of strings
        """
        if type(cut) not in [str, list]:
            raise TypeError("Input BDT cut ({0}) must be {1}.".format(type(cut), [str, list]))
            return
        elif type(cut) == list:
            comp = list(np.array([type(l) for l in cut]) == np.array([str for l in cut]))
            if not np.all(comp):
                raise TypeError("All cuts in list must be query strings ({})".format(comp))
                return
            else:
                for c in cut:
                    self.MakeCut(c)
        else:
            for name, df in dfs.iteritems():
                self.dfs[name] = df.query(cut)
        return

    def Stacked(self, col, nBins, xMin, xMax, xLabel="", extra="", logY=False, 
                sigWeight=None, overflow=True):
        """Make stacked 1D histogram of MC and data for a given column name

        Parameters
        ----------
        col : str
            Column name to be plotted
        nBins, xMin, xMax : int, int/float, int/float
            Number of bins, minimum x-value, maximum x-value
        xLabel : str, optional
            Label for x-axis
        extra : str
            'Extra' string to append to plot name to prevent overwriting
        logY : bool, optional
            Toggle logarithmic y-axis (default False)
        sigWeight : float, optional
            Force weight for signal (default None)
        overflow : bool, optional
            Toggle overflow bins (default True)
        """
        bgsPlots = []
        sigPlots = []
        dataPlot = None
        for name, df in self.dfs.iteritems():
            # Get data
            x = df.loc[col].to_numpy()
            # Get scale1fb
            weights = df.loc["scale1fb"].to_numpy()
            # Manually set signal scalef1b
            if name == self.signal and np.array_equal(weights, np.ones_like(weights)):
                weights = 0.0*weights+sigWeight
            # Plot data
            thisPlot = Hist1D(x,label=name, bins=np.linspace(xMin,xMax,nBins+1), 
                              weights=weights, color=self.__colors[name], no_overflow=(not overflow))
            # Add plot to stack
            if name == self.signal: sigPlots.append(thisPlot)
            elif name == "data": dataPlot = thisPlot
            else: bgsPlots.append(thisPlot)

        if dataPlot:
            # Scale MC
            sf = dataPlot.get_integral() / sum(bgsPlots).get_integral()
            bgsPlots = [bg*sf for bg in bgsPlots]
        if not sigWeight:
            # Scale signal
            sf = float(np.amax(sum(bgsPlots).get_counts()) / np.amax(sum(sigPlots).get_counts()))
            sigPlots = [sg*sf for sg in sigPlots]

        plot_stack(bgs=bgsPlots, data=dataPlot, sigs=sigPlots, do_log=logY,
                   xlabel=xLabel, ylabel="Events",cms_type="Preliminary",
                   lumi=str(self.lumi), filename=(self.outDir+col+extra+".pdf"));
        
        if self.verbose: print("Created {} stacked plot.".format(col))

        return

if __name__ == "__main__":

    # --> Data Retrieval <-- #
    print("Loading dataframes...")
    dfs_2016 = GetData("/nfs-7/userdata/jguiang/rare-higgs/2016/v3-0-0/", verbose=False)
    dfs_2017 = GetData("/nfs-7/userdata/jguiang/rare-higgs/2017/v3-0-0/", verbose=False)
    dfs_2018 = GetData("/nfs-7/userdata/jguiang/rare-higgs/2018/v3-0-0/", verbose=False)
    print("Collecting years...")
    dataframes = CombineYears([dfs_2016, dfs_2017, dfs_2018], signal=config["signal"], verbose=True)

    # --> Plots <-- #
    print("Plotting...")
    # Plotter object
    plots = Plots(config, dataframes, verbose=True)
        # Drop rho or phi events
        if "recoRho" in col or self.analysis == "rho":
            df = df.drop(df.query("genHiggsMeson_id == 333").index)
        elif self.analysis == "phi":
            df = df.drop(df.query("genHiggsMeson_id == 113").index)
    # Cuts
    sanity = "scale1fb > -999"
    unblinded = "(recoHiggs_mass < 120 or recoHiggs_mass > 130) or isData == 0"
    ptCuts = "(recoWLepton_pt > 35 and abs(recoWLepton_id) == 11) or (recoWLepton_pt > 30 and abs(recoWLepton_id) == 13)"
    hem = "isHEM == 0"
    gold = "isGold == 1"
    filters = "passFilters == 1"
    plots.MakeCut([sanity, unblinded, ptCuts, hem, gold, filters])
    # Mass plots
    plots.Stacked("recoPhi_mass", 50,0.98,1.06, xLabel=r"$m_{K^{+}K^{-}}$ (GeV)", logY=False, overflow=False)
    plots.Stacked("recoRho_mass", 50,0.5,1.1, xLabel=r"$m_{\pi^{+}\pi^{-}}$ (GeV)", logY=False, overflow=False)
    plots.Stacked("recoHiggs_mass", 50,0,200, xLabel=r"$m_{H}$ (GeV)", logY=False)
    plots.Stacked("recoHiggs_mass", 40,100,180, xLabel=r"$m_{H}$ (GeV)", extra="_zoom", logY=False, overflow=False)
    # dR plots
    plots.Stacked("recoKpKm_dR", 50,0,0.1, xLabel=r"$dR(K^{+}, K^{-})$", logY=False)
    plots.Stacked("recoPipPim_dR", 50,0,0.1, xLabel=r"$dR(\pi^{+}, \pi^{-})$", logY=False)
    plots.Stacked("recoRhoGamma_dR", 50,0,4, xLabel=r"$dR(\rho, \gamma)$", logY=False)
    plots.Stacked("recoPhiGamma_dR", 50,0,4, xLabel=r"$dR(\phi, \gamma)$", logY=False)
    plots.Stacked("recoGammaWLepton_dR", 50,0,4, xLabel=r"$dR(\gamma, \ell)$", logY=False)
    # relIso plots
    plots.Stacked("recoPhi_relIso", 50,0,0.06, xLabel=r"$relIso_{\phi}$", logY=True)
    plots.Stacked("recoRho_relIso", 50,0,0.06, xLabel=r"$relIso_{\rho}$", logY=True)
    # pt plots
    plots.Stacked("recoPhi_pt", 50,0,200, xLabel=r"$p_{T,\phi}$", logY=True)
    plots.Stacked("recoWLepton_pt", 50,0,200, xLabel=r"$p_{T,\gamma}$", logY=True)
    plots.Stacked("recoGamma_pt", 50,0,200, xLabel=r"$p_{T,\gamma}$", logY=True)
    plots.Stacked("recoKp_pt", 50,0,200, xLabel=r"$p_{T,K^{+}}$", logY=True)
    plots.Stacked("recoKm_pt", 50,0,200, xLabel=r"$p_{T,K^{-}}$", logY=True)
    plots.Stacked("recoPip_pt", 50,0,200, xLabel=r"$p_{T,\pi^{+}}$", logY=True)
    plots.Stacked("recoPim_pt", 50,0,200, xLabel=r"$p_{T,\pi^{-}}$", logY=True)
    # eta plots
    plots.Stacked("recoPhi_eta", 50,-2.5,2.5, xLabel=r"$\eta_{T,\phi}$", logY=True)
    plots.Stacked("recoWLepton_eta", 50,-2.5,2.5, xLabel=r"$\eta_{T,\gamma}$", logY=True)
    plots.Stacked("recoGamma_eta", 50,-2.5,2.5, xLabel=r"$\eta_{T,\gamma}$", logY=True)
    plots.Stacked("recoKp_eta", 50,-2.5,2.5, xLabel=r"$\eta_{T,K^{+}}$", logY=True)
    plots.Stacked("recoKm_eta", 50,-2.5,2.5, xLabel=r"$\eta_{T,K^{-}}$", logY=True)
    plots.Stacked("recoPip_eta", 50,-2.5,2.5, xLabel=r"$\eta_{T,\pi^{+}}$", logY=True)
    plots.Stacked("recoPim_eta", 50,-2.5,2.5, xLabel=r"$\eta_{T,\pi^{-}}$", logY=True)
    # phi plots
    plots.Stacked("recoPhi_phi", 50,-4,4, xLabel=r"$\phi_{T,\phi}$", logY=True)
    plots.Stacked("recoWLepton_phi", 50,-4,4, xLabel=r"$\phi_{T,\gamma}$", logY=True)
    plots.Stacked("recoGamma_phi", 50,-4,4, xLabel=r"$phi_{\gamma}$", logY=True)
    plots.Stacked("recoKp_phi", 50,-4,4, xLabel=r"$\phi_{K^{+}}$", logY=True)
    plots.Stacked("recoKm_phi", 50,-4,4, xLabel=r"$\phi_{K^{-}}$", logY=True)
    plots.Stacked("recoPip_phi", 50,-4,4, xLabel=r"$\phi_{\pi^{+}}$", logY=True)
    plots.Stacked("recoPim_phi", 50,-4,4, xLabel=r"$\phi_{\pi^{-}}$", logY=True)
    # Magic angles
    plots.Stacked("recoMagAng_cosThetaStar", 50,-1,1, xLabel=r"$\cos(\theta^{*})$", logY=False)
    plots.Stacked("recoMagAng_cosTheta1", 50,-1,1, xLabel=r"$\cos(\theta_{1})$", logY=False)
    plots.Stacked("recoMagAng_cosTheta2", 50,-1,1, xLabel=r"$\cos(\theta_{2})$", logY=False)
    plots.Stacked("recoMagAng_Phi", 50,-4,4, xLabel=r"$\Phi$", logY=False, overflow=False)
    plots.Stacked("recoMagAng_Phi1", 50,-4,4, xLabel=r"$\Phi_{1}$", logY=False, overflow=False)
    plots.Stacked("recoMagAng_m1", 50,50,300, xLabel=r"$m_{1}$", logY=False)
    plots.Stacked("recoMagAng_m2", 50,0,200, xLabel=r"$m_{2}$", logY=False)
    print("Done.")
