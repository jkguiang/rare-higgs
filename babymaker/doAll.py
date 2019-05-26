#! /usr/bin/env python

import ROOT, sys, os
ROOT.SetMemoryPolicy( ROOT.kMemoryStrict )

def doAll(args, debug=False):

    # Must supply sample name
    if len(args) < 1: return
    # Parse arguments
    sample = args[0]
    if not os.path.isfile(sample):
        print("{} is not a file. Running in 'debug' mode.".format(sample))
        debug = True
    # Load .so files
    ROOT.gSystem.Load("./CMS3_CORE.so")
    ROOT.gROOT.ProcessLine(".L babyTree.C++");
    ROOT.gROOT.ProcessLine(".L ScanChain.C++");
    # Make TChain
    mc = ROOT.TChain("Events")
    mc.Add(sample)
    # Run ScanChain
    if not debug:
        print("Running ScanChain...")
        ROOT.ScanChain(mc, "output.root", sample) # ScanChain(TChain, outName, sampleName, ...)
        print("Done.")

    return

if __name__ == "__main__":
    doAll(sys.argv[1:])
