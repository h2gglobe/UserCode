import ROOT
from python.configProducer import *

ROOT.gROOT.SetStyle('Plain')
ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("../libLoopAll.so");


ROOT.gBenchmark.Start("Analysis");
from ROOT import Util

ut = Util();
#cfg = configProducer(ut,"inputfiles_mc.dat",2)
"""  
ut.LoopAndFillHistos();
ut.WriteHist();  
ut.WriteCounters();  

ROOT.gBenchmark.Show("Analysis");
"""
