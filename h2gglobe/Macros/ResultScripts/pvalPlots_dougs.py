print "Loading Root..."

#import pdb; pdb.set_trace()
from ROOT import *
import os
import sys
import array

gStyle.SetOptStat(000000)
gStyle.SetCanvasBorderMode(0);
gStyle.SetCanvasColor(kWhite);
gROOT.SetBatch(1)
#ROOT.gROOT.ProcessLine(".x ../tdrstyle.cc")

print "Setting Initial Parameters."
can = TCanvas("Plots","Plots")
can.SetLogy(True)
can.SetGrid(True)
leg = TLegend(0.65, 0.74, 0.89, 0.89)
leg.SetFillColor(0)
leg.SetBorderSize(1)
mytext = TLatex()
mytext.SetTextSize(0.04)
mytext.SetNDC()
Minimum=0.0001
Maximum=10
Sigmas=3
#intlumi=[5.1,3.8]
#Energy=[7,8]
intlumi=[12.2]
Energy=[8]
legands=[]
legands=["8TeV Observed","8TeV Expected"]
#legands=["Combined Observed", "Combined Category 1", "Combined Category 2", "Combined Category 3", "Combined Category 4", "Combined dijet Categories"]
FrequantistPValues=[]
#FrequantistPValues=["higgsCombinePValue.HybridNew.mH125.0.root"]
#FrequantistPValues=["higgsCombinePValue.HybridNew.mH123.5.root","higgsCombinePValue.HybridNew.mH124.0.root","higgsCombinePValue.HybridNew.mH124.5.root","higgsCombinePValue.HybridNew.mH125.0.root","higgsCombinePValue.HybridNew.mH125.5.root","higgsCombinePValue.HybridNew.mH126.0.root"]
def getPValue(mass,tree):
  for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    if tree.mh==mass:
      return tree.limit,tree.limitErr
  print "Mass "+str(mass)+" not found exiting"
  exit(1)

Masses=array.array("f",[x * 0.1 for x in range(1100,1501,10)])
files = sys.argv[1:-1]

if len(files)==2 and len(legands)==0: legands=["Observed","Expected"]
elif len(legands)==0:
  for i in range(len(files)+1):
    if i==0:legands.append("Observed")
    else: legands.append("Category "+str(i))

#colors = [ kBlack, kRed, kOrange, kGreen, kCyan, kBlue, kMagenta]
colors = [ kBlack, kRed, kBlue, kBlack]

multigraph = TMultiGraph()
multigraph.SetTitle(";M_{H} (GeV/c^{2});Local P-Value")
multigraph.SetMinimum(Minimum)
multigraph.SetMaximum(Maximum)

for file,legand,color in zip(files,legands,colors):
  print file
  PValues=[]
  PValuesError=[]
  tree = TFile.Open(file).Get("limit")
  
  for Mass in Masses:
    pval,perr=getPValue(Mass,tree)
    PValues.append(pval)
    PValuesError.append(perr)

  graph=TGraphErrors(len(PValues),Masses,array.array("f",PValues),array.array("f",[0]*len(Masses)),array.array("f",PValuesError))
  graph.SetLineColor(color)
  graph.SetLineWidth(2)
  #if legand.find("Expected")!=-1: graph.SetLineStyle(7)
  graph.SetFillColor(kWhite)
  leg.AddEntry(graph,legand,'L')
  multigraph.Add(graph)

if len(FrequantistPValues)>0:
  FrequentistGraph = TGraphErrors(len(FrequantistPValues))
  for i in range(len(FrequantistPValues)):
    FrequentistTree = TFile.Open(FrequantistPValues[i]).Get("limit")
    FrequentistTree.GetEntry(0)
    FrequentistGraph.SetPoint(i,FrequentistTree.mh,FrequentistTree.limit)
    FrequentistGraph.SetPointError(i,0,FrequentistTree.limitErr)
  FrequentistGraph.SetLineColor(kRed)
  FrequentistGraph.SetLineWidth(1)
  FrequentistGraph.SetMarkerStyle(8)
  FrequentistGraph.SetMarkerColor(kBlue)
  FrequentistGraph.SetFillColor(kWhite)
  #multigraph.Add(FrequentistGraph)
  leg.AddEntry(FrequentistGraph,"Assembly of Toys",'lep')

multigraph.Draw("AC")
multigraph.GetXaxis().SetRangeUser(110,150)
if len(FrequantistPValues)>0: FrequentistGraph.Draw("EP SAME")
leg.Draw("")

line=[]
label=[]
for i in range(Sigmas):
  y = RooStats.SignificanceToPValue(i+1)
  line.append(TLine(110, y, 150, y))
  line[i].SetLineWidth(2)
  line[i].SetLineStyle(2)
  line[i].SetLineColor(kRed)
  line[i].Draw()
  label.append(TLatex(110 + 2, y * 1.1, "%d #sigma" % (i+1)))
  label[i].SetTextAlign(11);
  label[i].Draw()

if len(intlumi)==2: mytext.DrawLatex(0.18,0.82,"#splitline{CMS preliminary}{#splitline{#sqrt{s} = %i TeV L = %.1f fb^{-1}}{#sqrt{s} = %i TeV L = %.1f fb^{-1}}}" %(int(Energy[0]),float(intlumi[0]),int(Energy[1]),float(intlumi[1])))
else: mytext.DrawLatex(0.18,0.8,"#splitline{CMS preliminary}{#sqrt{s} = %i TeV L = %.1f fb^{-1}}" %(int(Energy[0]),float(intlumi[0])))

can.SaveAs(sys.argv[-1])
print "Done!"
