import FWCore.ParameterSet.Config as cms

## Import the photon re-reco sequence
from HiggsAnalysis.HiggsTo2photons.photonReReco_cff import *

## Loosen the re-reco cuts for the Z -> mu mu gamma tag-and-probe
##+ WARNING: These are extremely loose settings!
photonCore.minSCEt = 2.0
photons.minSCEtBarrel = 2.0
photons.minSCEtEndcap = 2.0
photons.maxHoverEBarrel = 10.0
photons.maxHoverEEndcap = 10.0

