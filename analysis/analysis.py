# Handling ROOT files
import json
import numpy as np
import uproot
import pandas as pd
import pickle
# Machine Learning
from sklearn.metrics import roc_curve
import xgboost as xgb # BDT
# Custom
from plots import Plots
from tools import GetData, CombineYears, ScaleByCol
from config import config

class Analysis(Plots):

    """Make plots using the Plots methods and pre-trained BDT model for use in the
        H --> rho/phi+gamma analysis and other related HEP analyses

    Child to Plots class (see plots.py for shared attribute docs)

    Attributes
    ----------
    __preds : dict{ str: numpy.ndarray }
        Dictionary of dataframes indexed by sample name
    bst : xgboost.core.Booster 
        Pre-trained BDT
    features : list
        List of features used to train BDT
    verbose : bool
        Verbosity of module (default False)

    Methods
    -------
    Predict()
        Generate BDT predictions for all dataframes
    BestBDTCut(x_test, preds_test, window="")
        Calculate optimal BDT working point from testing data       
    MakeBDTCut(bdtCut)
        Make cut on BDT variable, update dataframes
    """

    # Private
    __preds = {}

    def __init__(self, config, dfs, bst, features, **kwargs):
        """
        Parameters
        ----------
        config : dict
            Configuration for analysis
        dfs : dict{ str: pandas.DataFrame }
            Dictionary of dataframes indexed by sample name
        bst : xgboost.core.Booster
            Pre-trained BDT
        features : list[str]
            List of features used to train BDT
        """
        # Pass relevant args and remaining kwargs to parent
        super(Analysis, self).__init__(config, dfs, **kwargs)
        # Public
        self.bst = bst
        self.features = features
        self.verbose = kwargs["verbose"]

    def Predict(self):
        """Generate BDT predictions for all dataframes"""
        for name, df in self.dfs.iteritems():
            self.__preds[name] = self.bst.predict(xgb.DMatrix(df[self.features]))

        return

    def BestBDTCut(self, x_test, preds_test, window=None):
        """Calculate optimal BDT working point from testing data 
        
        Using sigma = sqrt(2*(nSig+nBg)*ln(1+nSig/nBg)-2*nSig)
        where nSig = # signal events and nBg = # background events

        Best BDT working point at maximal sigma value

        Parameters
        ----------
        x_test : pandas.DataFrame
            BDT test data
        preds_test : pandas.DataFrame
            BDT test predictions
        window : numpy.ndarray, (optional)
            Numpy array of booleans removing data not intended for
            consideration in ROC curve

        Raises
        ------
        ValueError
            If window is empty (i.e. [False, False, ..., False])
        """
        if not window is None:
            if not np.any(window):
                raise ValueError("Given window is empty (all values are false).")
            else:
                x_test = x_test[window]
                preds_test = preds_test[window]
        # BDT ROC Curve
        fpr, tpr, thresh = roc_curve(x_test.signal, preds_test)
        # Signal, background counts
        counts = x_test.signal.value_counts()
        nSig = tpr*counts[True]
        nBg = fpr*counts[False]
        # Prevent divide-by-zero errors
        nonzero = (nBg != 0)
        nSig = nSig[nonzero]
        nBg = nBg[nonzero]
        # Sigma
        sigma = (2*(nSig+nBg)*np.log(1+nSig/nBg)-2*nSig)**0.5
        # Optimal points
        bestEff = max(sigma) 
        bestBDT = float(thresh[np.where(sigma == bestEff)[0][0]]) # Cast numpy.float to float

        if self.verbose:
            print("BDT working point: {}".format(bestBDT))
            print("BDT efficiency: {}".format(bestEff))

        return bestBDT

    def MakeBDTCut(self, bdtCut):
        """Make cut on BDT variable, update dataframes
        
        Parameters
        ----------
        bdtCut : int/float
            BDT working point

        Raises
        ------
        ValueError
            If BDT working point not between 0.0 and 1.0
        """
        if bdtCut > 1.0 or bdtCut < 0.0:
            raise ValueError("BDT cut must be between 0.0 and 1.0.")
        else:
            for name, df in self.dfs.iteritems():
                self.dfs[name] = df[self.__preds[name] > bdtCut]
        return
