import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *
generator = cms.EDFilter("Pythia8GeneratorFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.0),
    crossSection = cms.untracked.double(1.0),
    PythiaParameters = cms.PSet(
            pythia8CommonSettingsBlock,
            pythia8CP5SettingsBlock,
        processParameters = cms.vstring(
            'HiggsSM:ffbar2HW = on',
            '24:onMode = off',
            '24:onIfAny = 13 14',
            '24:onIfAny = 11 12',
            '25:oneChannel =  1 0.5 100 333 22',
            '25:addChannel =  1 0.5 100 113 22'
        ),
            parameterSets = cms.vstring('pythia8CommonSettings',
                                        'pythia8CP5Settings',
                                        'processParameters',
                                        )
    )
)
configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('\$Revision$'),
    name = cms.untracked.string('\$Source$'),
    annotation = cms.untracked.string('W,Higgs-rho,gamme, 13 TeV, TuneCP5')
)
ProductionFilterSequence = cms.Sequence(generator)
