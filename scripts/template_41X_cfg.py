import FWCore.ParameterSet.Config as cms

#DATA TYPE
flagData = 'OFF'
flagMC = 'OFF'

#SKIM TYPE
flagSkimDiphoton = 'OFF'
flagNoSkim = 'OFF'
flagMMgSkim = 'OFF'

#ADDITIONAL OPTIONS
flagAOD = 'ON'
jobMaker = 'jobmaker unknown'

if (not((flagNoSkim is 'ON') ^ (flagSkimDiphoton is 'ON') ^ (flagMMgSkim is 'ON'))):
  print "You must skim or not skim... these are your options"
  exit(-1)


process = cms.Process("Globe") 
process.load("Configuration.StandardSequences.GeometryDB_cff") 
if flagAOD is 'ON':
  process.load("HiggsAnalysis.HiggsTo2photons.h2ganalyzer_41X_AOD_cfi")
else:
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

process.load("HiggsAnalysis.HiggsTo2photons.CMSSW_RelValDUMMY_cfi")
#process.skipEvents = cms.untracked.PSet(input=cms.untracked.uint32(3500))
#skipEvents = cms.untracked.uint32(3500)

hltLabel = "HLT"

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500


process.superClusterMerger =  cms.EDProducer(
  "EgammaSuperClusterMerger",
  src = cms.VInputTag(cms.InputTag('correctedHybridSuperClusters'),
                      cms.InputTag('correctedMulti5x5SuperClustersWithPreshower'))
  )


process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(15),
                                           maxd0 = cms.double(2)
                                           )

process.goodPhotonsLowPtCut = cms.EDFilter("PhotonSelector",
                                      src = cms.InputTag("photons"),
                                      cut = cms.string(
                                      "abs(superCluster.eta) < 2.5"
                                      " && superCluster.energy*sin(superCluster.position.theta) > 20."
                                      " && hadronicOverEm < 0.5 "
                                      " && trkSumPtHollowConeDR03 < 2.0*(3.5 + 0.001*superCluster.energy*sin(superCluster.position.theta))"
                                      " && ecalRecHitSumEtConeDR03 < 2.0*(4.2 + 0.006*superCluster.energy*sin(superCluster.position.theta))"
                                      " && hcalTowerSumEtConeDR03 < 2.0*(2.2 + 0.0025*superCluster.energy*sin(superCluster.position.theta))"
                                                      )
                                   )

process.TwoPhotonsLowPtCut = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("goodPhotonsLowPtCut"),
                                       minNumber = cms.uint32(2)
                                     )

process.goodPhotonsHighPtCut = cms.EDFilter("PhotonSelector",
                                      src = cms.InputTag("photons"),
                                      cut = cms.string(
                                      "abs(superCluster.eta) < 2.5"
                                      " && superCluster.energy*sin(superCluster.position.theta) > 30."
                                      " && hadronicOverEm < 0.5 "
                                      " && trkSumPtHollowConeDR03 < 2.0*(3.5 + 0.001*superCluster.energy*sin(superCluster.position.theta))"
                                      " && ecalRecHitSumEtConeDR03 < 2.0*(4.2 + 0.006*superCluster.energy*sin(superCluster.position.theta))"
                                      " && hcalTowerSumEtConeDR03 < 2.0*(2.2 + 0.0025*superCluster.energy*sin(superCluster.position.theta))"
                                                      )
                                   )

process.OnePhotonsHighPtCut = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("goodPhotonsHighPtCut"),
                                       minNumber = cms.uint32(1)
                                     )


process.dummySelector = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("gsfElectrons"),
                                     minNumber = cms.uint32(0)
                                     )


#process.goodEvents = cms.Sequence(process.noScraping * process.primaryVertexFilter)
#process.pathToCheck = cms.Sequence(process.L10and34 *process.noScraping*process.primaryVertexFilter*process.HFCoincidence*process.L140or41)
#process.pathToCheck2 = cms.Sequence(process.L10and34 *process.noScraping*process.primaryVertexFilter*process.L140or41*process.HFCoincidence)

process.diMuonSelSeq.remove(process.ZMuHLTFilter)

if flagSkimDiphoton == 'ON':
  process.eventFilter1 = cms.Sequence(process.superClusterMerger*process.goodPhotonsLowPtCut*process.TwoPhotonsLowPtCut) # for bkg
elif flagNoSkim == 'ON':    
  process.eventFilter1 = cms.Sequence(process.dummySelector)   #for signal MC
elif flagMMgSkim == 'ON':
  process.eventFilter1 = cms.Sequence(process.diMuonSelSeq*process.photonReReco)


process.h2ganalyzer.RootFileName = 'aod_mc_test.root'
process.h2ganalyzer.Debug_Level = 0

process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)

process.h2ganalyzerPath = cms.Sequence(process.h2ganalyzer)
process.p11 = cms.Path(process.eventFilter1*process.kt6PFJets*process.h2ganalyzerPath)

process.h2ganalyzer.JobMaker = jobMaker

if flagMC is 'ON':
  process.h2ganalyzer.doGenJet_algo1 = True
  process.h2ganalyzer.doGenJet_algo2 = True
  process.h2ganalyzer.doGenJet_algo3 = True
  process.h2ganalyzer.doGenParticles = True
  process.h2ganalyzer.doReducedGen = True
elif flagData is 'ON':
  process.h2ganalyzer.doPileup = False
  process.h2ganalyzer.doGenJet_algo1 = False
  process.h2ganalyzer.doGenJet_algo2 = False
  process.h2ganalyzer.doGenJet_algo3 = False
  process.h2ganalyzer.doGenParticles = False
  process.h2ganalyzer.doGenVertices = False
  process.h2ganalyzer.doReducedGen = False

if flagMC is 'ON' and flagAOD is 'OFF':
  process.h2ganalyzer.doSimTracks = True
  process.h2ganalyzer.doSimTrackPlusSimVertex = False

if flagAOD is 'ON':
  process.h2ganalyzer.doAodSim = True
  process.h2ganalyzer.doHcal = False
  process.h2ganalyzer.doHFHcal = False
  process.h2ganalyzer.doPreshowerHits = False
else:
  process.h2ganalyzer.doAodSim = False
  process.h2ganalyzer.doHcal = True
  process.h2ganalyzer.doHFHcal = True
  process.h2ganalyzer.doPreshowerHits = True
  

process.h2ganalyzer.doL1 = True
process.h2ganalyzer.doHLT = True

process.GlobalTag.globaltag = "START39_v8::All"
process.h2ganalyzer.HLTParameters.PrimaryTriggerResultsTag = cms.InputTag("TriggerResults","", hltLabel)
process.h2ganalyzer.HLTParameters.useSecondaryTrigger = cms.bool(False)
process.h2ganalyzer.HLTParameters.TriggerResultsTag = cms.InputTag("hltTriggerSummaryAOD","", hltLabel)
