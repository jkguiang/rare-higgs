from metis.CMSSWTask import CMSSWTask
from metis.Sample import DirectorySample, DummySample
from metis.StatsParser import StatsParser
import time

import numpy as np

tag = "v1"
total_summary = {}
for _ in range(25):
    lhe = CMSSWTask(
            sample = DummySample(N=1, dataset="/WH_HtoRhoGammaPhiGamma/privateMC_102x/step1"),
            events_per_output = 1000,
            total_nevents = 1000000,
            pset = "gensim_cfg.py",
            cmssw_version = "CMSSW_10_2_5",
            scram_arch = "slc6_amd64_gcc700",
            tag = tag,
            split_within_files = True,
            )

    raw = CMSSWTask(
            sample = DirectorySample(
                location = lhe.get_outputdir(),
                dataset = lhe.get_sample().get_datasetname().replace("step1","GENSIM"),
                ),
            open_dataset = True,
            files_per_output = 1,
            pset = "rawsim_cfg.py",
            cmssw_version = "CMSSW_10_2_5",
            scram_arch = "slc6_amd64_gcc700",
            tag = tag,
            )

    aod = CMSSWTask(
            sample = DirectorySample(
                location = raw.get_outputdir(),
                dataset = raw.get_sample().get_datasetname().replace("GENSIM","RAWSIM"),
                ),
            open_dataset = True,
            files_per_output = 5,
            pset = "aodsim_cfg.py",
            cmssw_version = "CMSSW_10_2_5",
            scram_arch = "slc6_amd64_gcc700",
            tag = tag,
            )

    miniaod = CMSSWTask(
            sample = DirectorySample(
                location = aod.get_outputdir(),
                dataset = aod.get_sample().get_datasetname().replace("RAWSIM","AOD"),
                ),
            open_dataset = True,
            flush = True,
            files_per_output = 5,
            pset = "miniaodsim_cfg.py",
            cmssw_version = "CMSSW_10_2_5",
            scram_arch = "slc6_amd64_gcc700",
            tag = tag,
            )

    cms4 = CMSSWTask(
            sample = DirectorySample(
                location = miniaod.get_outputdir(),
                dataset = miniaod.get_sample().get_datasetname().replace("AOD","MINIAOD"),
                ),
            open_dataset = True,
            flush = True,
            files_per_output = 1,
            output_name = "merged_ntuple.root",
            pset = "/home/users/namin/2017/ProjectMetis/pset_CMS4_V10-02-04.py",
            pset_args = "data=False year=2018",
            global_tag = "102X_upgrade2018_realistic_v12",
            cmssw_version = "CMSSW_10_2_5",
            scram_arch = "slc6_amd64_gcc700",
            tag = tag,
            tarfile = "/nfs-7/userdata/libCMS3/lib_CMS4_V10-02-04_1025.tar.xz",
            )

    tasks = [lhe,raw,aod,miniaod,cms4]

    for task in tasks:
        task.process()
        summary = task.get_task_summary()
        total_summary[task.get_sample().get_datasetname()] = summary

    StatsParser(data=total_summary, webdir="~/public_html/dump/metis/").do()
    time.sleep(2.0*3600)
