# Running HiggsCombine
## Get the right environment
Install it
`>>> . setup.sh`
or source it
`>>> cd CMSSW_8_1_0`
`>>> cmsenv`
unless you want seg faults, make sure you sourced CMSSW 8.0.1
## Make test card (makes a ROOT file)
`>>> text2workspace.py testcard.txt -b -o testcard.root`
## Make a likelihood scan from r=0 to r=3 (makes a ROOT file)[1]
`>>> combine -M MultiDimFit testcard.root --algo=grid --points=100 --setParameterRanges r=0,3`
### Look at output
`>>> root higgsCombineTest.MultiDimFit.mH120.root`
`root [1] limit->Draw("2*deltaNLL:r")`
## Find best fit r value[2]
`combine -M FitDiagnostics testcard.root --robustFit=1 --saveShapes --saveWithUncertainties --saveOverallShapes`
## Set upper-bound on signal strength (r) at 95% CL (makes a ROOT file)
`combine -M AsymptoticLimits testcard.root`

# Notes
1. This can be added to any Higgs Combine command to run on an "Asimov" dataset (i.e. sets observed = expected everywhere)
`-t -1 --expectSignal=0`
Running with `--expectSignal=1` instead injects signal into dataset. When doing something you expect to have **zero** signal, you can add the flag `--rMin=-1`. This bounds the value of r(mu) beyond zero(default). If you make it *too* negative, it will yell at you. If it yells at you, check the `higgsCombine.root` file and see where the r values start to make sense.

2. In the HiggsCombine [documentation](https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideCMSDataAnalysisSchool2014HiggsCombPropertiesExercise#The_basics), the "r" variable referenced here is called mu


