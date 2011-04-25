import ROOT
from python.configProducer import *

ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("../libLoopAll.so");

ROOT.gBenchmark.Start("Analysis");

ut = ROOT.Util();
cfg = configProducer(ut,"inputfiles.dat",2)
  
ut.LoopAndFillHistos();
ut.WriteHist();  
ut.WriteCounters();  

ROOT.gBenchmark.Show("Analysis");

