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
        "EGamma_Run2018D-22Jan2019-v2_MINIAOD_CMS4_V10-02-04"      : "data_Run2018D_EGamma_22Jan2019",
        },
    "ttbar" : {
        "TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-04"       : "ttsl_top2018",
        "TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-04"    : "ttsl_tbar2018",
        "TTGamma_SingleLeptFromT_TuneCP5_13TeV_madgraph_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2_MINIAODSIM_CMS4_V10-02-04"    : "ttg_top2018",
        "TTGamma_SingleLeptFromTbar_TuneCP5_13TeV_madgraph_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2_MINIAODSIM_CMS4_V10-02-04" : "ttg_tbar2018",

        "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-04" : "ttg_jets2018",
        "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2_MINIAODSIM_CMS4_V10-02-04"     : "tt_jets2018",
        "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-04"      : "tt_semilep2018",
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
    "drellyan" : {
        "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-04/"   : "dyjets2018",
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

        "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05" : "ttg_jets2017",
        "TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05"                   : "tt_jets_madgraph2017",
        "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05"                  : "tt_jets_amcatnlo2017",
        "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05"      : "tt_semilep2017",
        "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05"    : "tt_semilep_psweights2017",
        },
    "wjets" : {
        "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2_MINIAODSIM_CMS4_V10-02-05" : "wjets2017",
        },
    "wgamma" : {
        "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05"   : "wgamma2017",
        },
    "drellyan" : {
        "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05/"      : "dyjets2017",
        "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1_MINIAODSIM_CMS4_V10-02-05/" : "dyjets_ext2017",
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

        "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2_MINIAODSIM_CMS4_V10-02-05" : "ttg_jets2016",
        "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1_MINIAODSIM_CMS4_V10-02-05"                : "tt_jets_matgraph2016",
        "TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2_MINIAODSIM_CMS4_V10-02-05"             : "tt_jets_amcatnlo2016",
        "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1_MINIAODSIM_CMS4_V10-02-05"                       : "tt_powheg2016",
        },
    "wjets" : {
        "WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2_MINIAODSIM_CMS4_V10-02-05" : "wjets2016",
        },
    "wgamma" : {
        "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1_MINIAODSIM_CMS4_V10-02-05"   : "wgamma2016",
        "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext3-v1_MINIAODSIM_CMS4_V09-04-17"   : "wgamma_ext2016",
        },
    "drellyan" : {
        "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2_MINIAODSIM_CMS4_V09-04-17/" : "dyjets2016",
        "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2_MINIAODSIM_CMS4_V09-04-17/" : "dyjets_ext2016",
        },
    }

if __name__ == "__main__":
    project = "rare-higgs"
    tag = "v4-0-0"
    years = ["2016", "2017", "2018"]
    samples_years = {"2016": samples_2016, "2017": samples_2017, "2018": samples_2018}
    for i in range(0,2):
        test = bool(i)
        for year, samples_year in samples_years.iteritems():
            print("Making {0}{1}config files...".format(year, " test " if test else " "))
            base = "/hadoop/cms/store/group/snt/run2_"
            for typ, samples in tqdm(list(samples_year.iteritems())):
                suffix = "/"
                if typ == "whiggs":
                    suffix = "_private/"
                elif year == "2016":
                    suffix = "_94x/"
                thisBase = "{0}{1}{2}{3}".format(base, "data" if typ == "data" else "mc", year, suffix)
                for sample, name in samples.iteritems():
                    msg = MakeConfigs(year, thisBase+sample, project, name, test=test, tag=tag)
                    if msg: tqdm.write(msg)
    print("Done")
