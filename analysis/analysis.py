# Handling ROOT files
import json
import numpy as np
import uproot
import pandas as pd
# Machine Learning
from sklearn.metrics import roc_curve
import xgboost as xgb # BDT
# Custom
from plots import Plots
from tools import GetData, CombineYears
from config import config

class Analysis(Plots):

    """Make plots using the Plots methods and pre-trained BDT model for use in the
        H --> rho/phi+gamma analysis and other related HEP analyses

    Child to Plots class (see plots.py for shared attribute docs)

    Attributes
    ----------
    __preds : dict{ str: numpy.ndarray }
        dictionary of BDT predictions organized by sample
    bst : xgboost.core.Booster 
        pre-trained BDT
    features : list
        list of features used to train BDT
    verbose : bool
        verbosity of module (default False)

    Methods
    -------
    _initDfs()
        Prepare dataframes for BST cut
    BestBDTCut()
        Calculate optimal BDT working point from testing data       
    MakeBDTCut(bdtCut=None)
        Make cut on BDT variable, update dataframes
    """

    # Private
    self.__preds = {}

    def __init__(self, config, dfs, bst, features, **kwargs):
        """
        Parameters
        ----------
        config : dict
            Configuration for analysis
        dfs : dict{ str: pandas.DataFrame }
            Dataframe used for testing BDT
        bst : xgboost.core.Booster 
            pre-trained BDT
        features : list
            list of features used to train BDT
        """
        # Pass relevant args and remaining kwargs to parent
        super(Analysis, self).__init__(config, dfs, **kwargs)
        # Public
        self.bst = bst
        self.features = features
        self.verbose = kwargs["verbose"]
        # Init functions
        self._initDfs()

    def _initDfs(self):
        """Prepare dataframes for BST cut """
        # Get relevant *_pt columns to scale
        cols = list(set(self.dfs.columns[self.dfs.columns.str.contains("_pt")])-set(["recoWLepton_pt", "met_pt"])) # *_pt cols
        cols_s = [c+"_scaled" for c in cols] # Names of new *_pt_scaled cols
        for name, df in self.dfs.iteritems():
            # Weigh *_pt columns by reco Higgs mass
            df[cols_s] = df[cols].div(df.recoHiggs_mass, axis=0)
            # Save BST predictions
            self.__preds[name] = bst.predict(xgb.DMatrix(df[self.features]))

        return

    def BestBDTCut(x_test):
        """Calculate optimal BDT working point from testing data 
        
        Using sigma = sqrt(2*(nSig+nBg)*ln(1+nSig/nBg)-2*nSig)
        where nSig = # signal events and nBg = # background events

        Best BDT working point at maximal sigma value

        Parameters
        ----------
        x_test : pandas.DataFrame
            Dataframe used for testing BDT

        Raises
        ------
        TypeError
            If x_test is not a pandas dataframe
        """
        if type(x_test) != pd.core.frame.DataFrame:
            raise TypeError("Test data must be a pandas DataFrame")
            return
        # Signal, background counts
        counts = x_test.signal.value_counts()
        nSig = tpr*counts[True]
        nBg = fpr*counts[False]
        # Sigma
        sigma = (2*(nSig+nBg)*np.log(1+nSig/nBg)-2*nSig)**0.5
        sigma[~np.isfinite(sigma)] = 0.0
        # Optimal points
        bestEff = max(sigma)
        bestBDT = thresh[np.where(sigma == bestEff)[0][0]]

        return bestBDT

    def MakeBDTCut(bdtCut=None):
        """Make cut on BDT variable, update dataframes
        
        Parameters
        ----------
        bdtCut : int/float
            BDT working point

        Raises
        ------
        TypeError
            If BDT working poitn is not a float or integer
        ValueError
            If BDT working point not between 0.0 and 1.0
        """
        if type(bdtCut) not in [float, int]:
            raise TypeError("BDT cut must be a float or integer.")
            return
        elif bdtCut > 1.0 or bdtCut < 0.0:
            raise ValueError("BDT cut must be between 0.0 and 1.0.")
            return
        else:
            for name, df in dfs.iteritems():
                self.dfs[name] = df[self.__preds[name] > bdtCut]
        return

if __name__ == "__main__":

    # --> Data Retrieval <-- #
    print("Loading dataframes...")
    dfs_2016 = GetData("/nfs-7/userdata/jguiang/rare-higgs/2016/v3-0-0/", verbose=False)
    dfs_2017 = GetData("/nfs-7/userdata/jguiang/rare-higgs/2017/v3-0-0/", verbose=False)
    dfs_2018 = GetData("/nfs-7/userdata/jguiang/rare-higgs/2018/v3-0-0/", verbose=False)
    print("Collecting years...")
    dataframes = CombineYears([dfs_2016, dfs_2017, dfs_2018], signal=config["signal"], verbose=True)

    # --> BDT Retrieval <-- #
    # Get testing dataset
    x_test = pd.read_pickle("x_test.pkl", compression="gzip")
    # Get BDT model
    bst = pickle.load(open("bdt.pkl", "r"))
    # Get BDT features
    with open("features.json", "r") as fin:
        features = json.load(fin)

    # --> Analysis <-- #
    print("Plotting...")
    # Analysis object
    analysis = Analysis(config, dataframes, bst, features, verbose=True)
    # BDT working point
    bestBDT = analysis.BestBDTCut(x_test)
    # Make cut on data
    analysis.MakeBDTCut(bestBDT)
    # Plots
    analysis.Stacked("recoHiggs_mass", 50,0,200, xLabel=r"$m_{H}$ (GeV)", extra="_bdt", logY=False, blind=True)
    analysis.Stacked("recoHiggs_mass", 40,100,180, xLabel=r"$m_{H}$ (GeV)", extra="_bdtZoom", logY=False, blind=True, no_overflow=True)
