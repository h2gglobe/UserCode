import FWCore.ParameterSet.Config as cms

## Sequence to re-reco photons in >= 39x, see
##+ https://hypernews.cern.ch/HyperNews/CMS/get/egamma/960.html

## Some of these are loaded in the *_cfg.py; therefore they are 
##+ commented out here
from Configuration.StandardSequences.Services_cff import *
from Configuration.StandardSequences.MagneticField_38T_cff import *
from Configuration.StandardSequences.GeometryDB_cff import *
from Configuration.StandardSequences.Reconstruction_cff import *
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *
from RecoEgamma.EgammaPhotonProducers.conversionTracks_cff import *

photonReReco = cms.Sequence( ckfTracksFromConversions *
                             conversionSequence *
                             photonSequence *
                             photonIDSequence            )
