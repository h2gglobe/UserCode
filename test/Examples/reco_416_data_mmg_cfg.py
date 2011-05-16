import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.HiggsTo2photons.h2ganalyzerOptions_cfi import options
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import pickRelValInputFiles

### Set default option values
# options.inputFiles = pickRelValInputFiles( cmsswVersion  = 'CMSSW_4_1_3',
#                                            relVal        = 'Mu',
#                                            globalTag     = 'GR_R_311_V2',
#                                            numberOfFiles = 1              )

options.outputFile = "reco_mmg_test.root"
#options.hltProcessName = "REDIGI311X"
options.isRealData = True
options.globalTag = 'GR_R_311_V2'

### Get and parse the command line arguments
options.parseArguments()


process = cms.Process("Globe")
#process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("HiggsAnalysis.HiggsTo2photons.h2ganalyzer_41X_RECO_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("PhysicsTools/PatAlgos/patSequences_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Mixing')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('HiggsAnalysis.HiggsTo2photons.ZMuSkim_cff')
process.load('HiggsAnalysis.HiggsTo2photons.photonReRecoForMMG_cfi')

#from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *

process.load("HiggsAnalysis.HiggsTo2photons.CMSSW_RelValDUMMY_cfi")
#process.skipEvents = cms.untracked.PSet(input=cms.untracked.uint32(3500))
#skipEvents = cms.untracked.uint32(3500)
process.source.fileNames = options.inputFiles

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.h2ganalyzer.RootFileName = options.outputFile
process.h2ganalyzer.Debug_Level = 0

process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)

process.h2ganalyzerPath = cms.Sequence(process.h2ganalyzer)

process.diMuonSelSeq.remove( process.ZMuHLTFilter )

process.p11 = cms.Path( process.diMuonSelSeq *
                        process.kt6PFJets *
                        process.photonReReco *
                        process.h2ganalyzerPath )

doGenSim = not options.isRealData

process.h2ganalyzer.doGenJet_algo1 = doGenSim
process.h2ganalyzer.doGenJet_algo2 = doGenSim
process.h2ganalyzer.doGenJet_algo3 = doGenSim
process.h2ganalyzer.doGenParticles = doGenSim
process.h2ganalyzer.doGenVertices  = doGenSim
process.h2ganalyzer.doGenerator    = doGenSim

process.h2ganalyzer.doReducedGen = doGenSim

process.h2ganalyzer.doSimHits   = doGenSim
process.h2ganalyzer.doSimTracks = doGenSim
process.h2ganalyzer.doSimTrackPlusSimVertex = doGenSim

process.h2ganalyzer.doAodSim = doGenSim
process.h2ganalyzer.doPileup = doGenSim

process.h2ganalyzer.doL1 = True
process.h2ganalyzer.doHLT = True
process.GlobalTag.globaltag = options.globalTag + "::All"
process.h2ganalyzer.HLTParameters.PrimaryTriggerResultsTag = cms.InputTag(
    "TriggerResults", "", options.hltProcessName
)
process.h2ganalyzer.HLTParameters.useSecondaryTrigger = cms.bool(False)
process.h2ganalyzer.HLTParameters.TriggerResultsTag = cms.InputTag(
    "hltTriggerSummaryAOD", "", options.hltProcessName
)

## For tab-completion and history during interactive inspection
if __name__ == "__main__" : import user
