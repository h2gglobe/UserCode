#include "GlobeAnalyzer.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Tools.h"

#include <stdio.h>
#include <time.h>

#define TIMINGDEBUG 0

static string memory_usage() {
  ostringstream mem;
  //PP("hi");
  ifstream proc("/proc/self/status");
  string s;
  while(getline(proc, s), !proc.fail()) {
    if(s.substr(0, 6) == "VmSize") {
      mem << s;
      return mem.str();
    }
  }
  return mem.str();
}

double diffclock1(clock_t clock1,clock_t clock2) {
  double diffticks=clock1-clock2;
  double diffms=(diffticks*10)/CLOCKS_PER_SEC;
  return diffms;
}

GlobeAnalyzer::GlobeAnalyzer(const edm::ParameterSet& iConfig) {

  fileName = iConfig.getParameter<std::string>("RootFileName");
  jobmaker = iConfig.getParameter<std::string>("JobMaker");

  globalCountersNames = iConfig.getParameter<std::vector<std::string> >("globalCounters");

  doElectronStd = iConfig.getParameter<bool>("doElectron_std");
  doMuon = iConfig.getParameter<bool>("doMuon");
  doJetAlgo1 = iConfig.getParameter<bool>("doJet_algo1");
  doJetAlgo2 = iConfig.getParameter<bool>("doJet_algo2");
  doJetAlgo3 = iConfig.getParameter<bool>("doJet_algo3");
  doJetAlgoPF1 = iConfig.getParameter<bool>("doJet_algoPF1");
  doJetAlgoPF2 = iConfig.getParameter<bool>("doJet_algoPF2");
  doJetAlgoPF3 = iConfig.getParameter<bool>("doJet_algoPF3");
  doGenJetAlgo1 = iConfig.getParameter<bool>("doGenJet_algo1");
  doGenJetAlgo2 = iConfig.getParameter<bool>("doGenJet_algo2");
  doGenJetAlgo3 = iConfig.getParameter<bool>("doGenJet_algo3");
  doGenerator = iConfig.getParameter<bool>("doGenerator");
  doGenParticles = iConfig.getParameter<bool>("doGenParticles");
  doGenVertices = iConfig.getParameter<bool>("doGenVertices");
  doPhoton = iConfig.getParameter<bool>("doPhoton"); 
  doAllConversions = iConfig.getParameter<bool>("doAllConversions"); 
  doCaloTower = iConfig.getParameter<bool>("doCaloTower"); 
  doHcal = iConfig.getParameter<bool>("doHcal"); 
  doEcal = iConfig.getParameter<bool>("doEcal"); 
  doL1 = iConfig.getParameter<bool>("doL1");
  doHLT = iConfig.getParameter<bool>("doHLT");
  doVertices_std = iConfig.getParameter<bool>("doVertices_std"); 
  doVertices_nobs = iConfig.getParameter<bool>("doVertices_nobs"); 
  doMet = iConfig.getParameter<bool>("doMet"); 
  dotcMet = iConfig.getParameter<bool>("dotcMet");
  doPFMet = iConfig.getParameter<bool>("doPFMet");
  doSimHits = iConfig.getParameter<bool>("doSimHits"); 
  doSimTracks = iConfig.getParameter<bool>("doSimTracks"); 
  doTracks = iConfig.getParameter<bool>("doTracks"); 
  doGsfTracks = iConfig.getParameter<bool>("doGsfTracks"); 
  doTrackingParticles = iConfig.getParameter<bool>("doTrackingParticles"); 
  doEcalRecHits = iConfig.getParameter<bool>("doEcalRecHits");
  doReducedGen = iConfig.getParameter<bool>("doReducedGen");
  doLeptons = iConfig.getParameter<bool>("doLeptons");
  //doHt = iConfig.getParameter<bool>("doHt");
  doPFCandidates = iConfig.getParameter<bool>("doPFCandidates");
  //doPAT = iConfig.getParameter<bool>("doPAT");
  doRho = iConfig.getParameter<bool>("doRho");
  doPileup = iConfig.getParameter<bool>("doPileup");
  doPdfWeight = iConfig.getParameter<bool>("doPdfWeight");
  debug_level = iConfig.getParameter<int>("Debug_Level");

  common = new GlobeCommon(iConfig);

  if(doGenerator && doGenParticles) {
    std::cout << "doGenerator and doGenParticles cannot be true at the same time." << std::endl;
    throw;
  }

  if (doGenerator)
    gen = new GlobeGenerator(iConfig);

  if (doGenParticles)
    genP = new GlobeGenParticles(iConfig);
  
  if (doGenVertices)
    genV = new GlobeGenVertices(iConfig);

  if (doSimHits)
    simhits = new GlobeSimHits(iConfig);

  if (doSimTracks)
    simtracks = new GlobeSimTracks(iConfig);

  if (doL1)
    level1   = new GlobeL1(iConfig); 

  if (doHLT)
    hlt = new GlobeHLT(iConfig); 
   
  if (doEcal)
    ecalclusters = new GlobeEcalClusters(iConfig);

  if (doCaloTower)
    calotowers = new GlobeCaloTowers(iConfig);

  if (doTracks)
    tracks = new GlobeTracks(iConfig);

  if (doGsfTracks)
    gsfTracks = new GlobeGsfTracks(iConfig);

  if (doTrackingParticles)
    trackingParticles = new GlobeTrackingParticles(iConfig);

  if (doVertices_std) 
    vertex_std   = new GlobeVertex(iConfig, "std");

  if (doVertices_nobs) 
    vertex_nobs   = new GlobeVertex(iConfig, "nobs");

  if (doPhoton)
    photons = new GlobePhotons(iConfig);
  else
    photons = 0;

  if (doAllConversions)
    allConversions = new GlobeConversions(iConfig, "std");

  if (doElectronStd)
    std_electrons = new GlobeElectrons(iConfig, "std");
  else
    std_electrons = 0;

  if (doMuon)
    muons = new GlobeMuons(iConfig);
  else
    muons = 0;

  if (doMet)
    met = new GlobeMET(iConfig, "met"); 

  if (dotcMet)
    tcmet = new GlobeMET(iConfig, "tcmet");

  if (doPFMet)
    pfmet = new GlobeMET(iConfig, "pfmet");

  if (doGenJetAlgo1)
    algo1_genJets = new GlobeGenJets(iConfig, "algo1");
  if (doGenJetAlgo2)
    algo2_genJets = new GlobeGenJets(iConfig, "algo2");
  if (doGenJetAlgo3)
    algo3_genJets = new GlobeGenJets(iConfig, "algo3");

  if (doJetAlgo1)
    algo1_jets = new GlobeJets(iConfig, "algo1");
  if (doJetAlgo2)
    algo2_jets = new GlobeJets(iConfig, "algo2");
  if (doJetAlgo3)
    algo3_jets = new GlobeJets(iConfig, "algo3");
  if (doJetAlgoPF1)
    algoPF1_jets = new GlobeJets(iConfig, "algoPF1");
  if (doJetAlgoPF2)
    algoPF2_jets = new GlobeJets(iConfig, "algoPF2");
  if (doJetAlgoPF3)
    algoPF3_jets = new GlobeJets(iConfig, "algoPF3");

  selector = new GlobeSelector(iConfig);

  if (doReducedGen)
    reducedgen = new GlobeReducedGen(iConfig);

  if (doPFCandidates)
    pfCandidates = new GlobePFCandidates(iConfig);

  if (doRho) {
    rho1 = new GlobeRho(iConfig, "algo1");
    rho2 = new GlobeRho(iConfig, "algo2");
    rho3 = new GlobeRho(iConfig, "algo3");
  }

  if (doPileup)
    pileup = new GlobePileup(iConfig);

  if (doPdfWeight)
    pdfweights = new GlobePdfWeights(iConfig);

  if (doLeptons)
    leptons = new GlobeLeptons();
  if (!doMuon or !doElectronStd or !doPhoton) {
    std::cout << "WARNING: doLeptons needs doMuons, doElectronStd and doPhoton true" << std::endl;
    doLeptons = false;
  }
  
  if (doEcalRecHits)
    ecalrechits = new GlobeEcalHits(iConfig); 
  if (!doElectronStd and !doMuon and !doPhoton) {
    std::cout << "WARNING: EcalRecHits needs Electrons, Muons or Photons." << std::endl;
    doEcalRecHits = false;
  }
  
  if (doHcal)
    hcalhits = new GlobeHcal(iConfig);
  if (!doElectronStd and !doMuon and !doPhoton) {
    std::cout << "WARNING: HcalRecHits needs Electrons, Muons and Photons." << std::endl;
    doHcal = false;
  }

  if (!doLeptons) {
    std::cout << "WARNING: doReducedGen needs doLeptons true" << std::endl;
    doReducedGen = false;
  }
  
  readConfiguration(iConfig);
  nProcessedEvents = 0;
}

GlobeAnalyzer::~GlobeAnalyzer() {}

void GlobeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  clock_t begin, t;
  if (TIMINGDEBUG)
    begin=clock();
  
  if(debug_level > 9) std::cout << "GlobeAnalyzer: Begin" << std::endl;
  
  tot_events++;
  
  common->analyze(iEvent, iSetup);

  
  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T1: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
    std::cout << "Mem1: " << memory_usage() << std::endl;
  }
  

  //PHOTONS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: photons" << std::endl;  
  if (doPhoton)
    photons->analyze(iEvent, iSetup);

  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T1.2: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
    std::cout << "Mem2: " << memory_usage() << std::endl;
  }
  
  //ALL CONVERSIONS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: conversions" << std::endl;
  if (doAllConversions)
    allConversions->analyze(iEvent, iSetup);

  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T1: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
  }

  //ELECTRONS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: std_electrons" << std::endl;
  if (doElectronStd)
    std_electrons->analyze(iEvent, iSetup);
  
  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T2: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
    std::cout << "Mem3: " << memory_usage() << std::endl;
  }
  
  //MUONS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: muons" << std::endl;
  if (doMuon)
    muons->analyze(iEvent, iSetup);
  
  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T3: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
    std::cout << "Mem4: " << memory_usage() << std::endl;
  }

  //Leptons
  if (doLeptons)
    leptons->Zero();
  
  if(debug_level > 2) 
    std::cout << "GlobeAnalyzer: leptons->addmuons" << std::endl;
  if(doMuon && doLeptons)
    leptons->addMuons(muons);
  
  if(debug_level > 2) 
    std::cout << "GlobeAnalyzer: leptons->addelectrons" << std::endl;
  if(doElectronStd && doLeptons)
    leptons->addElectrons(std_electrons);
  
  if(debug_level > 2) std::cout << "GlobeAnalyzer: leptons->addphotons" << std::endl;
  if(doPhoton && doLeptons)
    leptons->addPhotons(photons);  
  
  //GENERATOR
  if(debug_level > 2) std::cout << "GlobeAnalyzer: gen" << std::endl;
  if (doGenerator)
    gen->analyze(iEvent, iSetup);
  
  //GENPARTICLES
  if(debug_level > 2) std::cout << "GlobeAnalyzer: genP" << std::endl;
  if (doGenParticles)
    genP->analyze(iEvent, iSetup);

  //GENVERTICES 
  if(debug_level > 2) std::cout << "GlobeAnalyzer: genV" << std::endl;
  if (doGenVertices)
    genV->analyze(iEvent, iSetup);

  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T4: " << double(diffclock1(t,begin)) << " ms"<< std::endl;    
    std::cout << "Mem5: " << memory_usage() << std::endl;
  }

  //SIMHITS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: simhits" << std::endl;
  if (doSimHits)
    simhits->analyze(iEvent, iSetup);
  
  //SIMTRACKS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: simtracks" << std::endl;
  if (doSimTracks)
    simtracks->analyze(iEvent, iSetup);  
  
  //LEVEL 1 TRIGGER
  if(debug_level > 2) std::cout << "GlobeAnalyzer: level1" << std::endl;
  if (doL1)
    level1->analyze(iEvent, iSetup);
  
  //HLTRIGGER
  if(debug_level > 2) std::cout << "GlobeAnalyzer: hlt" << std::endl;
  if (doHLT)
    hlt->analyze(iEvent, iSetup);
  
  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T5: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
    std::cout << "Mem6: " << memory_usage() << std::endl;
  }

  //TRACKS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: tracks" << std::endl;
  if (doTracks)
    tracks->analyze(iEvent, iSetup);
  
  if(debug_level > 2) std::cout << "GlobeAnalyzer: gsftracks" << std::endl;
  if (doGsfTracks)
    gsfTracks->analyze(iEvent, iSetup);
  
  //TRACKING PARTICLES
  if(debug_level > 2) std::cout << "GlobeAnalyzer: trackingparticles" << std::endl;
  if (doTrackingParticles) {
    trackingParticles->analyze(iEvent, iSetup,tracks);
    tracks->GetAssociatedTrackingParticleIndex(iEvent,iSetup,trackingParticles);
  }

  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T6: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
    std::cout << "Mem7: " << memory_usage() << std::endl;
  }

  //ECAL REC HITS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: ecalrechits" << std::endl;
  if (doEcalRecHits)
    ecalrechits->analyze(iEvent, iSetup, std_electrons, muons, photons);
  
  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T7: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
    std::cout << "Mem7: " << memory_usage() << std::endl;
  }
  
  //ECAL CLUSTERS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: ecalclusters" << std::endl;
  if (doEcal)
    ecalclusters->analyze(iEvent, iSetup);
  
  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T8: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
    std::cout << "Mem8: " << memory_usage() << std::endl;
  }
  
  //HCAL HITS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: hcalhits" << std::endl;
  if (doHcal)
    hcalhits->analyze(iEvent, iSetup, std_electrons, muons, photons);
  
  //CALO TOWERS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: calotowers" << std::endl;
  if (doCaloTower)
    calotowers->analyze(iEvent, iSetup);
  
  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T9: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
    std::cout << "Mem9: " << memory_usage() << std::endl;
  }

  //VERTEX
  if(debug_level > 2) std::cout << "GlobeAnalyzer: vertex_std" << std::endl;
  if (doVertices_std) 
    vertex_std->analyze(iEvent, iSetup);
  
  if(debug_level > 2) std::cout << "GlobeAnalyzer: vertex_nobs" << std::endl;
  if (doVertices_nobs) 
    vertex_nobs->analyze(iEvent, iSetup);

  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T10: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
    std::cout << "Mem10: " << memory_usage() << std::endl;
  }

  //MET
  if(debug_level > 2) std::cout << "GlobeAnalyzer: met" << std::endl;
  if (doMet)
    met->analyze(iEvent, iSetup);
  if (dotcMet)
    tcmet->analyze(iEvent, iSetup);
  if (doPFMet)
    pfmet->analyze(iEvent, iSetup);

  //GENJETS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: algo1genjet" << std::endl;
  if (doGenJetAlgo1)
    algo1_genJets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: algo2genjet" << std::endl;
  if (doGenJetAlgo2)
    algo2_genJets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: algo3genjet" << std::endl;
  if (doGenJetAlgo3)
    algo3_genJets->analyze(iEvent, iSetup);

  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T11: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
    std::cout << "Mem11: " << memory_usage() << std::endl;
  }

  //JETS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: algo1_jets" << std::endl;
  if (doJetAlgo1)
    algo1_jets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: algo2_jets" << std::endl;
  if (doJetAlgo2)
    algo2_jets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: algo3_jets" << std::endl;
  if (doJetAlgo3)
    algo3_jets->analyze(iEvent, iSetup);
    
  if(debug_level > 2) std::cout << "GlobeAnalyzer: algopf1_jets" << std::endl;
  if (doJetAlgoPF1)
    algoPF1_jets->analyze(iEvent, iSetup);
  if(debug_level > 2) std::cout << "GlobeAnalyzer: algopf2_jets" << std::endl;
  if (doJetAlgoPF2)
    algoPF2_jets->analyze(iEvent, iSetup);
  if(debug_level > 2) std::cout << "GlobeAnalyzer: algopf3_jets" << std::endl;
  if (doJetAlgoPF3)
    algoPF3_jets->analyze(iEvent, iSetup);

  if (TIMINGDEBUG) {
    t=clock();
    std::cout << "Time elapsed T12: " << double(diffclock1(t,begin)) << " ms"<< std::endl;
    std::cout << "Mem12: " << memory_usage() << std::endl;
  }

  //HT
  //if (doHt) {
  //  bool doLeptonHT=doJetAlgo1 && doMet && doPhoton && doLeptons;
  //  if(debug_level > 2) std::cout << "GlobeAnalyzer: leptonHT" << std::endl;
  //  if(doLeptonHT)
  //    ht->fillLeptonHT(algo1_jets, met, leptons);
  //  
  //  bool doCaloTowerHT=doJetAlgo1 && doMet && doPhoton && doLeptons;
  //  if(debug_level > 2) std::cout << "GlobeAnalyzer: caloHT" << std::endl;
  //  if(doCaloTowerHT)
  //    ht->fillCaloTowerHT(met, calotowers);

  //  if(debug_level > 2) std::cout << "GlobeAnalyzer: trackHT" << std::endl;
  //  if(doTracks)
  //    ht->fillTrackHT(iEvent);
  //}

  if(doReducedGen && doGenerator && doLeptons) {
    if(debug_level > 2) 
      std::cout << "GlobeAnalyzer: reducedgen" << std::endl;
    reducedgen->fillRedGenList(gen, leptons);
  }

  if(doReducedGen && doGenParticles && doLeptons) {
    if(debug_level > 2) 
      std::cout << "GlobeAnalyzer: reducedgenparticles" << std::endl;
    reducedgen->fillRedGenList(genP, leptons);
  }

  // PF CANDIDATES
  if (doPFCandidates)
    pfCandidates->analyze(iEvent, iSetup, tracks, muons, photons);

  //PAT
  //if (doPAT) {
  //  if(debug_level > 2) std::cout << "GlobeAnalyzer: PAT" << std::endl;
  //  pat->analyze(iEvent, iSetup, std_electrons, photons, algo1_jets);
  //}

  if (doRho) {
    rho1->analyze(iEvent, iSetup);
    rho2->analyze(iEvent, iSetup);
    rho3->analyze(iEvent, iSetup);
  }

  if (doPileup)
    pileup->analyze(iEvent, iSetup);

  if (doPdfWeight)
    pdfweights->analyze(iEvent, iSetup);

  if(debug_level > 2) 
    std::cout << "GlobeAnalyzer: selectorbits" << std::endl;

  if (doGenParticles || doGenerator)
    selector_bits = selector->select(std_electrons, muons, photons, gen, leptons, reducedgen).to_ulong();
  else
    selector_bits = selector->select(std_electrons, muons, photons).to_ulong();
  
  if(debug_level > 2) 
    std::cout << "GlobeAnalyzer: selectorbits = " << selector_bits << std::endl;
  
  if (selector_bits > 0) {
    sel_events++;
    
    if(debug_level > 2) 
      std::cout << "GlobeAnalyzer: fill my tree!" << std::endl;
    
    // fill the tree
    tree->Fill();
    nProcessedEvents++;
  }

  if (nProcessedEvents % 500 == 0) {
    file->cd();
    tree->AutoSave();
    nProcessedEvents = 0;
  }

  if(debug_level > 9) std::cout << "GlobeAnalyzer: End" << std::endl;
}

void GlobeAnalyzer::beginJob() { 

  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("event", "Event data");
  tree2 = new TTree("global_variables", "Global Parameters"); // need a different tree to fill once per job
  lumitree = new TTree("lumi", "Processed lumi sections");

  common->defineBranch(tree);
  common->defineLumiBranch(lumitree);
  
  if (doPhoton)
    photons->defineBranch(tree);

  if (doAllConversions)
    allConversions->defineBranch(tree);  

  if (doEcalRecHits)
    ecalrechits->defineBranch(tree);

  if (doEcal)
    ecalclusters->defineBranch(tree);

  if (doCaloTower)
    calotowers->defineBranch(tree);

  if (doHcal)
    hcalhits->defineBranch(tree);

  if (doL1)
    level1->defineBranch(tree);

  if (doHLT)
    hlt->defineBranch(tree);

  if (doVertices_std)
    vertex_std->defineBranch(tree);

  if (doVertices_nobs)
    vertex_nobs->defineBranch(tree);

  if (doMet)
    met->defineBranch(tree);

  if(dotcMet)
    tcmet->defineBranch(tree);

  if(doPFMet)
    pfmet->defineBranch(tree);

  if (doSimHits)
    simhits->defineBranch(tree);

  if (doSimTracks)
    simtracks->defineBranch(tree);

  if (doTracks)
    tracks->defineBranch(tree);

  if (doGsfTracks)
    gsfTracks->defineBranch(tree);

  if (doTrackingParticles)
    trackingParticles->defineBranch(tree);
 
  if (doElectronStd)
    std_electrons->defineBranch(tree);
 
  if (doMuon)
    muons->defineBranch(tree);

  if (doJetAlgo1)
    algo1_jets->defineBranch(tree);
  if (doJetAlgo2)
    algo2_jets->defineBranch(tree);
  if (doJetAlgo3)
    algo3_jets->defineBranch(tree);
  if (doJetAlgoPF1)
    algoPF1_jets->defineBranch(tree);
  if (doJetAlgoPF2)
    algoPF2_jets->defineBranch(tree);
  if (doJetAlgoPF3)
    algoPF3_jets->defineBranch(tree);

  if (doGenerator)
    gen->defineBranch(tree);

  if (doGenParticles)
    genP->defineBranch(tree);

  if (doGenVertices)
    genV->defineBranch(tree);

  if (doGenJetAlgo1)
    algo1_genJets->defineBranch(tree);
  if (doGenJetAlgo2)
    algo2_genJets->defineBranch(tree);
  if (doGenJetAlgo3)
    algo3_genJets->defineBranch(tree);

  if (doLeptons)
    leptons->defineBranch(tree);
  
  if (doReducedGen) 
    reducedgen->defineBranch(tree);

  if (doPFCandidates)
     pfCandidates->defineBranch(tree);

  if (doRho) {
    rho1->defineBranch(tree);
    rho2->defineBranch(tree);
    rho3->defineBranch(tree);
  }

  if (doPileup)
     pileup->defineBranch(tree);
  
  if (doPdfWeight)
    pdfweights->defineBranch(tree);

  defineBranch();
  
  tot_events = 0;
  sel_events = 0;
}

void GlobeAnalyzer::endJob() { 

  fillTree();
  file->cd();
  if (doPileup) {
    TH1D* h = pileup->getHisto(); 
    Int_t last_bin = h->GetNbinsX();
    h->SetBinContent(last_bin-1, h->GetBinContent(last_bin)+h->GetBinContent(last_bin-1));
    h->Write();  
    
    h = pileup->getHistoTrue(); 
    last_bin = h->GetNbinsX();
    h->SetBinContent(last_bin-1, h->GetBinContent(last_bin)+h->GetBinContent(last_bin-1));
    h->Write();  
  }

  tree->Write(0, TObject::kWriteDelete);
  tree2->Write(0, TObject::kWriteDelete);
  lumitree->Write();

  file->Close();
}

void GlobeAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & l, const edm::EventSetup & es) {
  common->endLumiBlock(l,es);
  lumitree->Fill();
  for(size_t ii=0; ii< globalCounters.size(); ++ii ) {
	  edm::Handle<edm::MergeableCounter> ctrHandle;
	  l.getByLabel(globalCountersNames[ii], ctrHandle);
	  globalCounters[ii] += ctrHandle->value;
  }
}


void GlobeAnalyzer::readConfiguration(const edm::ParameterSet& iConfig) {
  parameters = new std::vector<std::string>;
  parameters->push_back(iConfig.dump());
}

void GlobeAnalyzer::defineBranch() {

  tree->Branch("selector_bits", &selector_bits, "selector_bits/I");

  tree2->Branch("version", &version, "version/I");
  tree2->Branch("type", &type, "type/I");
  tree2->Branch("tot_events", &tot_events, "tot_events/I");
  tree2->Branch("sel_events", &sel_events, "sel_events/I");
  tree2->Branch("parameters", "std::vector<std::string>", &parameters); 
  tree2->Branch("jobmaker", "std::string", &jobmaker); 
  
  globalCounters.clear();
  globalCounters.resize(globalCountersNames.size(),0);
  for(size_t ii=0; ii< globalCounters.size(); ++ii ) {
	  tree2->Branch( globalCountersNames[ii].c_str(), &globalCounters[ii], (globalCountersNames[ii]+"/I").c_str() );
  }
  
}


void GlobeAnalyzer::fillTree() {
  char cmd[500];
  char* descr = getenv("CMSSW_BASE");
  sprintf(cmd, "cat %s/src/HiggsAnalysis/HiggsTo2photons/data/CVS/Tag", descr);
  ExecCommand exec(cmd);
  version = exec.getTag();
  type = 0;

  tree2->Fill();
}

DEFINE_FWK_MODULE(GlobeAnalyzer);
