import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

### Setup 'analysis'  options
options = VarParsing ('analysis')

### Register custom options
options.register( "globalTag",
    "GR10_P_V7::All", # default value is latest prompt reco (August 2010)
    VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.varType.string,         # bool, string, int, or float
    "Global tag to be used."
)

options.register( "isRealData",
    False, # default value
    VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.varType.bool,           # bool, string, int, or float
    "Is this real data?"
)

options.register( "hltProcessName",
    "HLT",                             # default value
    VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.varType.string,         # bool, string, int, or float
    "Name of the Process that produced the HLT information."
)

