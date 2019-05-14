from makeConfigs import MakeConfigs
from tqdm import tqdm

samples_2018 = {
    "data" : {
        "SingleMuon_Run2018A-17Sep2018-v2_MINIAOD_CMS4_V10-02-04"  : "data_Run2018A_SingleMuon_17Sep2018",
        "SingleMuon_Run2018B-17Sep2018-v1_MINIAOD_CMS4_V10-02-04"  : "data_Run2018B_SingleMuon_17Sep2018",
        "SingleMuon_Run2018C-17Sep2018-v1_MINIAOD_CMS4_V10-02-04"  : "data_Run2018C_SingleMuon_17Sep2018",
        "SingleMuon_Run2018D-PromptReco-v2_MINIAOD_CMS4_V10-02-04" : "data_Run2018D_SingleMuon_PromptReco-v2",
        "EGamma_Run2018A-17Sep2018-v2_MINIAOD_CMS4_V10-02-04"      : "data_Run2018A_EGamma_17Sep2018",
        "EGamma_Run2018B-17Sep2018-v1_MINIAOD_CMS4_V10-02-04"      : "data_Run2018B_EGamma_17Sep2018",
        "EGamma_Run2018C-17Sep2018-v1_MINIAOD_CMS4_V10-02-04"      : "data_Run2018C_EGamma_17Sep2018",
        "EGamma_Run2018D-PromptReco-v2_MINIAOD_CMS4_V10-02-04"     : "data_Run2018D_EGamma_PromptReco-v2",
        },
    "ttbar" : {
        "TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-04"       : "ttsl_top2018",
        "TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-04"    : "ttsl_tbar2018",
        "TTGamma_SingleLeptFromT_TuneCP5_13TeV_madgraph_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2_MINIAODSIM_CMS4_V10-02-04"    : "ttg_top2018",
        "TTGamma_SingleLeptFromTbar_TuneCP5_13TeV_madgraph_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2_MINIAODSIM_CMS4_V10-02-04" : "ttg_tbar2018",
        },
    "wjets" : {
        "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_MINIAODSIM_CMS4_V10-02-04" : "wjets2018",
        },
    "wgamma" : {
        "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-04"   : "wgamma2018",
        },
    "whiggs" : {
        "WH_HtoRhoGammaPhiGamma_privateMC_102x_MINIAOD_v1"   : "whiggs",
        },
    }

samples_2017 = {
    "data" : {
        "SingleMuon_Run2017B-31Mar2018-v1_MINIAOD_CMS4_V10-02-05"     : "data_Run2017B_SingleMuon_31Mar2018",
        "SingleMuon_Run2017C-31Mar2018-v1_MINIAOD_CMS4_V10-02-05"     : "data_Run2017C_SingleMuon_31Mar2018",
        "SingleMuon_Run2017D-31Mar2018-v1_MINIAOD_CMS4_V10-02-05"     : "data_Run2017D_SingleMuon_31Mar2018",
        "SingleMuon_Run2017E-31Mar2018-v1_MINIAOD_CMS4_V10-02-05"     : "data_Run2017E_SingleMuon_31Mar2018",
        "SingleMuon_Run2017F-31Mar2018-v1_MINIAOD_CMS4_V10-02-05"     : "data_Run2017F_SingleMuon_31Mar2018",
        "SingleElectron_Run2017B-31Mar2018-v1_MINIAOD_CMS4_V10-02-05" : "data_Run2017B_EGamma_31Mar2018",
        "SingleElectron_Run2017C-31Mar2018-v1_MINIAOD_CMS4_V10-02-05" : "data_Run2017C_EGamma_31Mar2018",
        "SingleElectron_Run2017D-31Mar2018-v1_MINIAOD_CMS4_V10-02-05" : "data_Run2017D_EGamma_31Mar2018",
        "SingleElectron_Run2017E-31Mar2018-v1_MINIAOD_CMS4_V10-02-05" : "data_Run2017E_EGamma_31Mar2018",
        "SingleElectron_Run2017F-31Mar2018-v1_MINIAOD_CMS4_V10-02-05" : "data_Run2017F_EGamma_31Mar2018",
        },
    "ttbar" : {
        "TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05"         : "ttsl_top2017",
        "TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05"      : "ttsl_tbar2017",
        "TTGamma_SingleLeptFromT_TuneCP5_PSweights_13TeV_madgraph_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05" : "ttg_top2017",
        "TTGamma_SingleLeptFromTbar_TuneCP5_PSweights_13TeV_madgraph_pythia8_RunIIFall17MiniAOD-PU2017_94X_mc2017_realistic_v11-v1_MINIAODSIM_CMS4_V10-02-05"          : "ttg_tbar2017",
        },
    "wjets" : {
        "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2_MINIAODSIM_CMS4_V10-02-05" : "wjets2017",
        },
    "wgamma" : {
        "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05"   : "wgamma2017",
        },
    }

samples_2016 = {
    "data" : {
        "SingleMuon_Run2016B-17Jul2018_ver2-v1_MINIAOD_CMS4_V10-02-05"     : "data_Run2016B_SingleMuon_17Jul2018_ver2-v1",
        "SingleMuon_Run2016C-17Jul2018-v1_MINIAOD_CMS4_V10-02-05"          : "data_Run2016C_SingleMuon_17Jul2018",
        "SingleMuon_Run2016D-17Jul2018-v1_MINIAOD_CMS4_V10-02-05"          : "data_Run2016D_SingleMuon_17Jul2018",
        "SingleMuon_Run2016E-17Jul2018-v1_MINIAOD_CMS4_V10-02-05"          : "data_Run2016E_SingleMuon_17Jul2018",
        "SingleMuon_Run2016F-17Jul2018-v1_MINIAOD_CMS4_V10-02-05"          : "data_Run2016F_SingleMuon_17Jul2018",
        "SingleMuon_Run2016G-17Jul2018-v1_MINIAOD_CMS4_V10-02-05"          : "data_Run2016G_SingleMuon_17Jul2018",
        "SingleMuon_Run2016H-17Jul2018-v1_MINIAOD_CMS4_V10-02-05"          : "data_Run2016H_SingleMuon_17Jul2018",
        "SingleElectron_Run2016B-17Jul2018_ver2-v1_MINIAOD_CMS4_V10-02-05" : "data_Run2016B_SingleElectron_17Jul2018_ver2-v1",
        "SingleElectron_Run2016C-17Jul2018-v1_MINIAOD_CMS4_V10-02-05"      : "data_Run2016C_SingleElectron_17Jul2018-v1",
        "SingleElectron_Run2016D-17Jul2018-v1_MINIAOD_CMS4_V10-02-05"      : "data_Run2016D_SingleElectron_17Jul2018-v1",
        "SingleElectron_Run2016E-17Jul2018-v1_MINIAOD_CMS4_V10-02-05"      : "data_Run2016E_SingleElectron_17Jul2018-v1",
        "SingleElectron_Run2016F-17Jul2018-v1_MINIAOD_CMS4_V10-02-05"      : "data_Run2016F_SingleElectron_17Jul2018-v1",
        "SingleElectron_Run2016G-17Jul2018-v1_MINIAOD_CMS4_V10-02-05"      : "data_Run2016G_SingleElectron_17Jul2018-v1",
        "SingleElectron_Run2016H-17Jul2018-v1_MINIAOD_CMS4_V10-02-05"      : "data_Run2016H_SingleElectron_17Jul2018-v1",
        },
    "ttbar" : {
        "TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_MINIAODSIM_CMS4_V09-04-17"         : "ttsl_top2016",
        "TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2_MINIAODSIM_CMS4_V10-02-05"    : "ttsl_top_ext2016",
        "TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_MINIAODSIM_CMS4_V09-04-17"      : "ttsl_tbar2016",
        "TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2_MINIAODSIM_CMS4_V10-02-05" : "ttsl_tbar_ext2016",
        "TTGamma_SingleLeptFromT_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_MINIAODSIM_CMS4_V10-02-05"         : "ttg_top2016",
        "TTGamma_SingleLeptFromTbar_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1_MINIAODSIM_CMS4_V10-02-05"      : "ttg_tbar2016",
        },
    "wjets" : {
        "WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2_MINIAODSIM_CMS4_V10-02-05" : "wjets2016",
        },
    "wgamma" : {
        "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1_MINIAODSIM_CMS4_V10-02-05"   : "wgamma2016",
        "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext3-v1_MINIAODSIM_CMS4_V09-04-17"   : "wgamma_ext2016",
        },
    }

if __name__ == "__main__":
    project = "rare-higgs"
    tag = "v3-0-0"
    for i in range(0,2):
        test = bool(i)
        print("Making 2018{}config files...".format(" test " if test else " "))
        base = "/hadoop/cms/store/group/snt/run2_"
        for typ, samples in tqdm(list(samples_2018.iteritems())):
            thisBase = "{0}{1}2018{2}".format(base, "data" if typ == "data" else "mc", "_private/" if typ == "whiggs" else "/")
            for sample, name in samples.iteritems():
                MakeConfigs("2018", thisBase+sample, project, name, test=test, tag=tag)
        print("Making 2017{}config files...".format(" test " if test else " "))
        for typ, samples in tqdm(list(samples_2017.iteritems())):
            thisBase = "{0}{1}2017{2}".format(base, "data" if typ == "data" else "mc", "/")
            for sample, name in samples.iteritems():
                MakeConfigs("2017", thisBase+sample, project, name, test=test, tag=tag)
        print("Making 2016{}config files...".format(" test " if test else " "))
        for typ, samples in tqdm(list(samples_2016.iteritems())):
            thisBase = "{0}{1}2016{2}".format(base, "data" if typ == "data" else "mc", "_94x/")
            for sample, name in samples.iteritems():
                MakeConfigs("2016", thisBase+sample, project, name, test=test, tag=tag)
    print("Done")
