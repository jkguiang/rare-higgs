# Handling ROOT files
import os
import uproot, pandas as pd

def GetData(inDir, pklJar="pickles/", verbose=False):
    """ Return dictionary of Pandas dataframes """
    # Check if directories exist
    if not os.path.isdir(inDir):
        print("ERROR: {} does not exist.".format(inDir))
        return
    if not os.path.isdir(pklJar):
        os.mkdir(pklJar)
    # Check naming convention
    if inDir[-1] != "/":
        inDir += "/"
    if pklJar[-1] != "/":
        pklJar += "/"
    # Make dictionary of dataframes
    data = {}
    for output in os.listdir(inDir):
        # Get file name
        name = (output.split("output_")[-1]).split(".")[0]
        # Pickle file
        pkl = pklJar+name+".pkl"
        pickled = os.path.isfile(pkl)
        # Open file, save dataframe
        if not pickled:
            f=uproot.open(inDir+output) # TFile
            t=f["tree"] # TTree
            data[name] = t.pandas.df() # dataframe
        else:
            data[name] = pd.read_pickle(pkl)
        # Pickle dataframe
        if not pickled: data[name].to_pickle(pkl)

    if verbose:
        print("Loaded Dataframes:\n    "+"\n    ".join(data.keys()))

    return data

if __name__ == "__main__":
    GetData("outputs", verbose=True)
