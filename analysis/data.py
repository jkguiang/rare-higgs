# Handling ROOT files
import uproot, os

def GetData(inDir, verbose=False):
    """ Return dictionary of Pandas dataframes """

    if not os.path.isdir(inDir):
        print("ERROR: {} does not exist.")
        return

    data = {}
    for output in os.listdir(inDir):
        name = (output.split("output_")[-1]).split(".")[0]
        f=uproot.open("outputs/"+output) # TFile
        t=f["tree"] # TTree
        data[name] = t.pandas.df() # dataframe

    if verbose:
        print("Loaded Dataframes:\n    "+"\n    ".join(data.keys()))

    return data
