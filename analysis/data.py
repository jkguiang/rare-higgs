# General imports
import os, ast, json
# Handling ROOT files
import uproot, pandas as pd, numpy as np

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

def WriteJSON(fname, bst, feature_names, class_labels=[]):
    """ Return JSON of BDT (Written by Nick Amin) """
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
            "feature_names": feature_names,
            "class_labels": map(int,np.array(class_labels).tolist()), # numpy array not json serializable
            }
    with open(fname, "w") as fhout:
        json.dump(to_dump,fhout,indent=2)

if __name__ == "__main__":
    GetData("outputs", verbose=True)
