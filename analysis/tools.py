"""A collection of uproot/numpy/pandas-based tools for general use in python
   HEP analyses

Functions
---------
CombineYears(dfs_years, signal="", verbose=False)
    Combine 2016, 2017, 2018 dataframes by sample
GetData(inDir, pklJar="pickles/", mcOnly=False, verbose=False)
    Make dictionary of Pandas dataframes
ModifyData(dfs)
    Return modified dictionary of Pandas dataframes
BDTtoJSON(fname, bst, features, labels=[])
    Return JSON of BDT model (Written by Nick Amin)
JSONtoC(fname_in, fname_out=None)
    Return C function that mimics BDT output (Written by Nick Amin)
"""

# General imports
import os, ast, json, re
# Handling ROOT files
import uproot
import pandas as pd
import numpy as np
# Custom Tools
from config import config

def CombineYears(dfs_years, signal="", verbose=False):
    """Combine 2016, 2017, 2018 dataframes by sample

    Parameters
    ----------
    dfs_years : dict{ str: list[pandas.DataFrame] }
        Dictionary of lists of pandas dataframes indexed by year
    signal : str, optional
        Name of signal sample (default "")
    verbose : bool, optional
        Verbosity of function (default False)
    """
    samples = ["ttg", "ttsl", "wjets", "wgamma", "data"]
    # Regular expressions
    exps = [re.compile("^{}.*20.*$".format(mc)) for mc in samples[:-1]]
    exps.append(re.compile("^data_Run20..._(SingleMuon|SingleElectron|EGamma)_.*$"))
    # Get dictionary of list dataframes to concatenate
    dfs = {}
    for dfs_year in dfs_years:
        for name, df in dfs_year.iteritems():
            for i, exp in enumerate(exps):
                # Signal doesn't match regexes
                if name == signal:
                    dfs[name] = df
                # Dataset-specific trigger cuts
                if "SingleMuon" in name:
                    df = df.drop(df.query("HLT_singleMu == 0").index)
                if "SingleElectron" in name or "EGamma" in name:
                    df = df.drop(df.query("HLT_singleEl == 0").index)
                # Add to dictionary
                if re.match(exp, name):
                    key = samples[i]
                    if key in dfs:
                        dfs[key].append(df)
                    else:
                        dfs[key] = [df]
    # Concatenate dataframes
    for key, df_list in dfs.iteritems():
        if key != signal:
            if verbose: print("Combining {0} dfs for {1}".format(len(df_list), key))
            dfs[key] = (pd.concat(df_list)).reset_index()
        # Add isData column
        if "isData" not in dfs[key]:
            if key == "data":
                dfs[key]["isData"] = np.ones(dfs[key].name.shape)
            else:
                dfs[key]["isData"] = np.zeros(dfs[key].name.shape)

    return dfs

def GetData(inDir, pklJar="pickles/", mcOnly=False, verbose=False):
    """Make dictionary of Pandas dataframes

    Parameters
    ----------
    inDir : str
        Path to existing directory of .root files
    pklJar : str, optional
        Path where pickled dataframes should be saved (default 'pickles/')
    mcOnly : bool, optional
        Specify if only grapping Monte Carlo samples (default False)
    verbose : bool, optional
        Verbosity of function (default False)
    """
    if not os.path.isdir(pklJar):
        os.mkdir(pklJar)
    # Check naming convention
    if inDir[-1] != "/":
        inDir += "/"
    if pklJar[-1] != "/":
        pklJar += "/"
    ver = inDir.split("/")[-2]
    ver_re = re.compile("v.-.-.")
    if not re.match(ver_re, ver):
        raise ValueError("{} is not versioned correctly.".format(inDir))
        return
    # Make dictionary of dataframes
    dfs = {}
    for output in os.listdir(inDir):
        # Extract name
        name_re = re.compile("output_*")
        if not re.match(name_re, output): continue
        name = (output.split("output_")[-1]).split(".")[0]
        # Pickle file
        pkl = "{0}/{1}/{2}.pkl".format(pklJar, ver, name)
        pickled = os.path.isfile(pkl)
        # Open file, save dataframe
        if "data" not in output or not mcOnly:
            if not pickled:
                f=uproot.open(inDir+output) # TFile
                t=f["tree"] # TTree
                dfs[name] = t.pandas.df() # dataframe
            else:
                dfs[name] = pd.read_pickle(pkl, compression="gzip")
                # Pickle dataframe
                if not pickled: dfs[name].to_pickle(pkl, compression="gzip")

    dfs = ModifyData(dfs)

    if verbose:
        print("Loaded Dataframes (Version {}):\n    ".format(ver)+"\n    ".join(dfs.keys()))

    return dfs

def ModifyData(dfs):
    """Return modified dictionary of Pandas dataframes

    Parameters
    ----------
    dfs : dict{ str: pandas.DataFrame }
        Dictionary of pandas dataframes indexed by sample name
    """
    # Add bookkeeping columns, delete bad rows
    for name, df in dfs.iteritems():
        # Add signal bool and dataset name columns
        df["name"] = name
        # Require found lepton, photon, meson
        df = df.drop(df.query("recoWLepton_pt < 0 or recoPhi_pt < 0 or recoRho_pt < 0 or recoGamma_pt < 0").index)
        # Avoid double-counting prompt photons
        genGammaMatch = ( df.genRecoGamma_isMatch == 1 )
        minGammaParton = ( df.minGammaParton_dR > 0.4 )
        if name in ["WGToLNuG", "TTGamma_SingleLeptFromT", "TTGamma_SingleLeptFromTbar"]:
            df = df.drop(df.query("not @genGammaMatch or not @minGammaParton").index)
        elif name in ["WJetsToLNu", "TTJets_SingleLeptFromT", "TTJets_SingleLeptFromTbar"]:
            df = df.drop(df.query("@genGammaMatch and @minGammaParton").index)
        # Update dataframe in dictionary
        dfs[name] = df

    return dfs

def ScaleByCol(dfs, bCols, wCol):
    """Scale all dataframe 'base' columns by one 'weight' column

    Parameters
    ----------
    dfs : dict{ str: pandas.DataFrame }
        Dictionary of pandas dataframes indexed by sample name
    bCols : list[str]
        List of base column names
    wCol : str
        Weight column name
    """
    for name, df in dfs.iteritems():
        sCols = [c+"_scaled" for c in bCols] # Names of new scaled cols
        dfs[name][sCols] = df[bCols].div(df[wCol], axis=0)

    return dfs

def BDTtoJSON(fname, bst, features, labels=[]):
    """Return JSON of BDT model (Written by Nick Amin)

    Parameters
    ----------
    fname : str
        Desired name of output file
    bst : xgboost.core.Booster
        Trained BDT
    features : list[str]
        List of features used to train BDT
    labels : list[str], optional
        Class labels (default [])
    """
    buff = "[\n"
    trees = bst.get_dump()
    ntrees = len(trees)
    for itree,tree in enumerate(trees):
        prev_depth = 0
        depth = 0
        for line in tree.splitlines():
            depth = line.count("\t")
            nascending = prev_depth - depth
            (depth == prev_depth-1)
            # print ascending, depth, prev_depth
            prev_depth = depth
            parts = line.strip().split()
            padding = "    "*depth
            for iasc in range(nascending):
                buff += "{padding}]}},\n".format(padding="    "*(depth-iasc+1))
            if len(parts) == 1:  # leaf
                nodeid = int(parts[0].split(":")[0])
                leaf = float(parts[0].split("=")[-1])
                # print "leaf: ",depth,nodeid,val
                buff += """{padding}{{ "nodeid": {nodeid}, "leaf": {leaf} }},\n""".format(
                        padding=padding,
                        nodeid=nodeid,
                        leaf=leaf,
                        )
            else:
                nodeid = int(parts[0].split(":")[0])
                split, split_condition = parts[0].split(":")[1].replace("[","").replace("]","").split("<")
                split_condition = float(split_condition)
                yes, no, missing = map(lambda x:int(x.split("=")[-1]), parts[1].split(","))
                buff += """{padding}{{ "nodeid": {nodeid}, "depth": {depth}, "split": "{split}", "split_condition": {split_condition}, "yes": {yes}, "no": {no}, "missing": {missing}, "children": [\n""".format(
                        padding=padding,
                        nodeid=nodeid,
                        depth=depth,
                        split=split,
                        split_condition=split_condition,
                        yes=yes,
                        no=no,
                        missing=missing,
                        )
        for i in range(depth):
            padding = "    "*(max(depth-1,0))
            if i == 0:
                buff += "{padding}]}}".format(padding=padding)
            else:
                buff += "\n{padding}]}}".format(padding=padding)
            depth -= 1
        if itree != len(trees)-1:
            buff += ",\n"
    buff += "\n]"
    # print buff
    to_dump = {
            "trees": list(ast.literal_eval(buff)),
            "features": features,
            "labels": map(int,np.array(labels).tolist()), # numpy array not json serializable
            }
    with open(fname, "w") as fhout:
        json.dump(to_dump,fhout,indent=2)

def JSONtoC(fname_in, fname_out=None):
    """Return C function that mimics BDT output (Written by Nick Amin)

    Parameters
    ----------
    fname_in : str
        Name of input JSON file
    fname_out : str, optional
        Desired name of output C file
    """
    with open(fname_in, "r") as fhin:
        data = json.loads(fhin.read())
        trees = data["trees"]
        features = data["features"]
        labels = data["labels"]
    def get_leaf(entry, depth=0):
        if "leaf" in entry: 
            return entry["leaf"]
        splitvar = entry["split"]
        splitval = entry["split_condition"]
        yesnode = [c for c in entry["children"] if c["nodeid"] == entry["yes"]][0]
        nonode = [c for c in entry["children"] if c["nodeid"] == entry["no"]][0]
        return "({} < {} ? {} : {})".format(splitvar, splitval, get_leaf(yesnode, depth=depth+1), get_leaf(nonode, depth=depth+1))
    buff = ""
    multi = False
    if len(labels) > 0:
        multi = True
        ntrees_per_class = len(trees) // len(labels)
        nclasses = len(labels)
    if multi:
        buff += "std::vector<float> get_prediction({}) {{\n".format(",".join(map(lambda x: "float {}".format(x), features)))
        for ic in labels:
            buff += "  float w_{} = 0.;\n".format(ic)
        for itree,j in enumerate(trees):
            iclass = int(labels[itree % nclasses])
            buff += "  w_{} += {};\n".format(iclass,get_leaf(j))
        buff += "  float w_sum = {};\n".format("+".join("exp(w_{})".format(ic) for ic in labels))
        for ic in labels:
            buff += "  w_{0} = exp(w_{0})/w_sum;\n".format(ic)
        buff += "  return {{ {} }};\n".format(",".join("w_{}".format(ic) for ic in labels))
    else:
        buff += "float get_prediction({}) {{\n".format(",".join(map(lambda x: "float {}".format(x), features)))
        buff += "  float w = 0.;\n"
        for itree,j in enumerate(trees):
            buff += "  w += {};\n".format(get_leaf(j))
        buff += "  return 1.0/(1.0+exp(-1.0*w));\n"
    buff += "}"
    if fname_out:
        with open(fname_out, "w") as fhout:
            fhout.write(buff)
    else:
        return buff

if __name__ == "__main__":
    # Tests
    GetData("/nfs-7/userdata/jguiang/rare-higgs/2017/v3-0-0", pklJar="pickles/2018", verbose=True)
