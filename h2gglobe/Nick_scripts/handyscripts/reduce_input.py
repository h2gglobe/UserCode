import ROOT
from python.configProducer import *

import sys

ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("../../libLoopAll.so");

ROOT.gBenchmark.Start("Reduction")

ut = ROOT.Util();
cfg = configProducer(ut,sys.argv[1],1)
ut.LoopAndFillHistos()

ROOT.gBenchmark.Show("Reduction");


