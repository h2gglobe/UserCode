#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cstring>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TList.h"
#include "TString.h"
#include "Normalization.C"

using namespace std;

int davidCalls=0;
int fracCalls=0;
int linCalls=0;
int systCalls=0;
int bkgCalls=0;

TH1F* Interpolate(double massLow, TH1F* low, double massHigh, TH1F* high, double massInt){

  if (low->GetNbinsX()!=high->GetNbinsX()) std::cout << "Cannot interpolate differently binned histograms" << std::endl;
  assert(low->GetNbinsX()==high->GetNbinsX());

  TH1F *interp = (TH1F*)low->Clone();
  for (int i=0; i<interp->GetNbinsX(); i++){
    double OutLow = low->GetBinContent(i+1);
    double OutHigh = high->GetBinContent(i+1);
    double OutInt = (OutHigh*(massInt-massLow)-OutLow*(massInt-massHigh))/(massHigh-massLow);
    interp->SetBinContent(i+1,OutInt);
  }

  std::string name = low->GetName();
  std::string strLow = Form("%3.1f",massLow);
  std::string strInt = Form("%3.1f",massInt);
  int ind = name.rfind(strLow);
  name.replace(ind-6,11,strInt);
  interp->SetName(name.c_str());

  return interp;
   
}

std::pair<int,int> findNearest(double mass){
  
  double bestdiff=1000.0;
  double nbestdiff=1000.0;
  int nearest;
  int nnearest;

  for (int i=115; i<155; i+=5){
    if (i==145) continue;
    double diff = std::fabs((mass+0.001)-double(i));
    if (diff<bestdiff) {
      bestdiff=diff;
      nearest = i;
    }
  }
  for (int i=115; i<155; i+=5){
    if (i==145) continue;
    double diff = std::fabs((mass+0.001)-double(i));
    if (diff>bestdiff && diff<nbestdiff){
      nbestdiff=diff;
      nnearest=i;
    }
  }
  std::pair<int,int> nearestPair(nearest,nnearest);
  return nearestPair;
}

int getIndex(int mass){

  int index;
  if (mass==115) index=0;
  if (mass==120) index=1;
  if (mass==125) index=2;
  if (mass==130) index=3;
  if (mass==135) index=4;
  if (mass==140) index=5;
  if (mass==150) index=6;

  return index;
}

TH1F* linearBin(TH1F *hist){
  assert(hist->GetEntries()!=0);
  assert(hist->GetNbinsX()>2);
  int nBins = hist->GetNbinsX();
  TH1F *temp = new TH1F(Form("lin_%s_%d",hist->GetName(),linCalls),Form("lin_%s_%d",hist->GetName(),linCalls),nBins,0,nBins);
  for (int i=0; i<nBins; i++){
    temp->SetBinContent(i+1,hist->GetBinContent(i+1));
    temp->SetBinError(i+1,hist->GetBinError(i+1));
  }
  linCalls++;
  return temp;
}

void plotBkgModel(TList* HistList, std::string name){
  
  gROOT->SetBatch();
  system("mkdir -p plots/ada/bkgMod");
  system("mkdir -p plots/grad/bkgMod");

  std::string bdt;
  TString str = HistList->At(0)->GetName();
  if (str.Contains("ada")) bdt="ada";
  else if (str.Contains("grad")) bdt="grad";
  else std::cout << "Error find BDT type" << std::endl;
  assert (str.Contains("ada") || str.Contains("grad"));
  
  gStyle->SetOptStat(0);
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  int color[4] = {kBlue,kMagenta,kRed,kGreen+2};

  TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.15,0.6,0.55,0.85);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  TPaveText *txt = new TPaveText(0.6,0.6,0.8,0.85,"NDC");
  txt->SetFillColor(0);
  txt->SetLineColor(0);
  txt->AddText("#int L = 1.66 fb^{-1}");

  for (int i=1; i<HistList->GetEntries(); i++){
    if (((TH1F*)HistList->At(i))->GetNbinsX()!=((TH1F*)HistList->At(0))->GetNbinsX()) std::cout << "Plot problem: calling plot for histograms with different number of bins" << std::endl;
    assert (((TH1F*)HistList->At(i))->GetNbinsX()==((TH1F*)HistList->At(0))->GetNbinsX());
    TH1F *temp = linearBin((TH1F*)HistList->At(i));
    temp->Scale(((TH1F*)HistList->At(0))->Integral()/temp->Integral());
    temp->SetLineColor(color[i-1]);
    temp->SetMarkerStyle(20);
    temp->SetMarkerColor(color[i-1]);
    temp->SetTitle(Form("Background model in signal region from sidebands %s %s",bdt.c_str(),name.c_str()));
    temp->GetXaxis()->SetTitle("BDT Output Bin");
    temp->GetYaxis()->SetRangeUser(0.0,2.*(((TH1F*)HistList->At(0))->GetMaximum()));
    if (i==1) temp->Draw("p");
    else temp->Draw("same p");
    if (i==1) leg->AddEntry(temp,"Low 2 sideband","lep");
    if (i==2) leg->AddEntry(temp,"Low 1 sideband","lep");
    if (i==3) leg->AddEntry(temp,"High 1 sideband","lep");
    if (i==4) leg->AddEntry(temp,"High 2 sideband","lep");
  }
  leg->Draw("same");
  txt->Draw("same");

  canv->Print(("plots/"+bdt+"/bkgMod/"+name+".png").c_str(),"png");
  
  delete canv;
  delete txt;
  delete leg;

  bkgCalls++;
}


void plotDavid(TH1F* bkgT, TH1F* sigT, TH1F* dataT, std::string name){
  
  gROOT->SetBatch();
  system("mkdir -p plots/ada/david");
  system("mkdir -p plots/grad/david");
  system("mkdir -p plots/ada/diff");
  system("mkdir -p plots/grad/diff");

  std::string bdt;
  TString str = dataT->GetName();
  if (str.Contains("ada")) bdt="ada";
  else if (str.Contains("grad")) bdt="grad";
  else std::cout << "Error find BDT type" << std::endl;
  assert (str.Contains("ada") || str.Contains("grad"));

  if (bkgT->GetNbinsX() != sigT->GetNbinsX() || sigT->GetNbinsX() != dataT->GetNbinsX()) std::cout << "Plot problem: calling plot for histograms with different number of bins" << std::endl;
  assert(bkgT->GetNbinsX() == sigT->GetNbinsX() || sigT->GetNbinsX() == dataT->GetNbinsX());
  
  TH1F *bkg = linearBin(bkgT);
  TH1F *sig2 = linearBin(sigT);
  TH1F *sig5 = (TH1F*)sig2->Clone();
  TH1F *sig10 = (TH1F*)sig2->Clone();
  TH1F *data = linearBin(dataT);
  TH1F *diff = (TH1F*)data->Clone();
  diff->Add(bkg,-1);

  gStyle->SetOptStat(0);
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();

  TCanvas *canv = new TCanvas();
  bkg->SetLineColor(kBlue);
  bkg->SetFillColor(kBlue-9);
  sig2->SetLineColor(kRed);
  sig2->SetLineStyle(3);
  sig2->Scale(2.);
  sig5->SetLineColor(kRed);
  sig5->SetLineStyle(7);
  sig5->Scale(5.);
  sig10->SetLineColor(kRed);
  sig10->Scale(10.);
  data->SetMarkerStyle(20);
  bkg->SetTitle(("Background, data and signal distributions for "+bdt+" "+name).c_str());
  sig2->SetTitle(("Background, data and signal distributions for "+bdt+" "+name).c_str());
  data->SetTitle(("Background, data and signal distributions for "+bdt+" "+name).c_str());
  bkg->GetXaxis()->SetTitle("BDT Output Bin");
  sig2->GetXaxis()->SetTitle("BDT Output Bin");
  data->GetXaxis()->SetTitle("BDT Output Bin");
  bkg->GetYaxis()->SetTitle("Events");
  sig2->GetYaxis()->SetTitle("Events");
  data->GetYaxis()->SetTitle("Events");

  TLegend *leg = new TLegend(0.15,0.6,0.55,0.85);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->AddEntry(bkg,"Background","f");
  leg->AddEntry(data,"Data","lep");
  leg->AddEntry(sig10,"Signal (2, 5, 10#~#times SM)","l");
  TPaveText *txt = new TPaveText(0.6,0.6,0.8,0.85,"NDC");
  txt->SetFillColor(0);
  txt->SetLineColor(0);
  txt->AddText("#int L = 1.66 fb^{-1}");

  bkg->GetYaxis()->SetRangeUser(0.0,2.*(data->GetMaximum()));

  bkg->Draw("e2");
  sig2->Draw("same hist");
  sig5->Draw("same hist");
  sig10->Draw("same hist");
  data->Draw("same e");
  leg->Draw("same");
  txt->Draw("same");

  canv->Print(("plots/"+bdt+"/david/"+name+".png").c_str(),"png");
  canv->Clear();
  TLegend *leg2 = new TLegend(0.15,0.6,0.55,0.85);
  leg2->SetFillColor(0);
  leg2->SetLineColor(0);
  sig10->GetYaxis()->SetRangeUser((-1*sig10->GetMaximum())+10,sig10->GetMaximum()+20);
  sig10->GetXaxis()->SetTitle("BDT Output Bin");
  sig10->GetYaxis()->SetTitle("Events");
  sig10->SetTitle(Form("Data, background difference compared to signal %s %s",bdt.c_str(),name.c_str()));
  diff->SetMarkerStyle(20);
  diff->GetXaxis()->SetTitle("BDT Output Bin");
  diff->GetYaxis()->SetTitle("Events");
  leg2->AddEntry(diff,"Data - background model","lep");
  leg2->AddEntry(sig10,"Signal (2, 5, 10 #times SM)","l");
  TF1 *line = new TF1("line","0.0",0.0,sig10->GetNbinsX()+1);
  line->SetLineColor(kBlue);
  sig10->Draw("hist");
  sig5->Draw("same hist");
  sig2->Draw("same hist");
  line->Draw("same");
  diff->Draw("p same");
  leg2->Draw("same");
  txt->Draw("same");
  canv->Print(("plots/"+bdt+"/diff/"+name+".png").c_str(),"png");

  delete canv;
  delete txt;
  delete leg;
  
  davidCalls++;
}

void plotFrac(TList* HistList, TH1F* compT, std::string name, bool fracAll){

  gROOT->SetBatch();
  system("mkdir -p plots/ada/fracs");
  system("mkdir -p plots/grad/fracs");

  std::string bdt;
  TString str = compT->GetName();
  if (str.Contains("ada")) bdt="ada";
  if (str.Contains("grad")) bdt="grad";

  int nHists=HistList->GetEntries();
  TH1F *comp = linearBin(compT);
  TH1F *uandd[nHists];
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  TCanvas *canv = new TCanvas("","",700,700);
  
  canv->Divide(1,2,0.1,0.1);
  canv->cd(1);
  gPad->SetPad(0.01,0.3,0.99,0.99);
  gPad->SetBottomMargin(0.0);
  gPad->SetLeftMargin(0.1);
  comp->SetLineColor(1);
  comp->SetFillColor(kGray);
  comp->Draw("hist");
  comp->SetTitle(("Up, down and interpolated templates for "+bdt+" "+name).c_str());
  comp->GetYaxis()->SetRangeUser(0.0,((TH1F*)HistList->At(0))->GetMaximum()+0.5);
  comp->GetYaxis()->SetTitle("Events / bin");
  comp->GetXaxis()->SetTitle("BDT output bin");
  int mass;
  if (name=="syst120") mass=120;
  if (name=="syst135") mass=135;
  TLegend *leg = new TLegend(0.2,0.2,0.6,0.8);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  if (name=="syst120" || name=="syst135"){
    leg->AddEntry(comp,Form("True %d signal binned at %d",mass,mass),"f");
  }
  
  canv->cd(2);
  gPad->SetPad(0.01,0.01,0.99,0.3);
  gPad->SetTopMargin(0.0);
  gPad->SetBottomMargin(0.2);
  gPad->SetLeftMargin(0.1);
  TF1 *line = new TF1("line","0.0",0.,comp->GetBinLowEdge(comp->GetNbinsX()+1));
  line->SetLineColor(kBlack);
  
  for (int i=0; i<nHists; i++){
    TH1F *temp = linearBin((TH1F*)HistList->At(i));
    temp->SetLineColor((i+1)*2);
    if (!fracAll && i==nHists-1) temp->SetLineColor(kBlack); 
    canv->cd(1);
    temp->DrawCopy("same e");
    uandd[i]= (TH1F*)temp->Clone();
    uandd[i]->Add(comp,-1);
    uandd[i]->Divide(temp);
    uandd[i]->Scale(-100.);
    //uandd[i]->Divide(comp);
    if (name=="syst120" || name=="syst135") {
      if (i==0) leg->AddEntry(uandd[i],Form("True %d signal binned at %d",mass+5,mass),"lep");
      if (i==1) leg->AddEntry(uandd[i],Form("True %d signal binned at %d",mass-5,mass),"lep");
      if (i==nHists-1) leg->AddEntry(uandd[i],Form("Interpolated %d signal",mass),"lep");
    }
    //uandd[i]->Scale(100.0);
    canv->cd(2);
    if (fracAll){
      if (i==0) uandd[i]->Draw("e");
      else {
        uandd[i]->Draw("same e");
        line->Draw("same");
      }
    }
    else {
      if (i==nHists-1) {
        uandd[i]->Draw("e");
        line->Draw("same");
      }
    }
    uandd[i]->SetTitle("");
    uandd[i]->GetYaxis()->SetLabelSize(0.08);
    uandd[i]->GetYaxis()->SetRangeUser(-50.,50.);
    uandd[i]->GetYaxis()->SetTitle("#frac{#Delta_{int}}{int} %");
    uandd[i]->GetYaxis()->SetTitleOffset(0.4);
    uandd[i]->GetYaxis()->SetTitleSize(0.08);
    uandd[i]->GetXaxis()->SetLabelSize(0.08);
    uandd[i]->GetXaxis()->SetTitle("BDT output");
    uandd[i]->GetXaxis()->SetTitleSize(0.08);
    canv->cd(1);
    if (name=="syst120" || name=="syst135") leg->Draw("same");
  }
  

  canv->Print(("plots/"+bdt+"/fracs/"+name+".png").c_str(),"png");
  delete canv;
  
  fracCalls++;
  
}

void plotSystFracs(TList* HistList, TH1F* compT, std::string name){

  gROOT->SetBatch();
  system("mkdir -p plots/ada/systs");
  system("mkdir -p plots/grad/systs");

  std::string bdt;
  TString str = compT->GetName();
  if (str.Contains("ada")) bdt="ada";
  if (str.Contains("grad")) bdt="grad";

  int nHists=HistList->GetEntries();
  TH1F *comp = linearBin(compT);
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  TCanvas *canv = new TCanvas("","",700,700);
  
  TLegend *leg = new TLegend(0.6,0.6,0.88,0.88);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  
  TF1 *line = new TF1("line","0.0",0.,comp->GetBinLowEdge(comp->GetNbinsX()+1));
  line->SetLineColor(kBlack);
  TF1 *line1 = new TF1("line","10.0",0.,comp->GetBinLowEdge(comp->GetNbinsX()+1));
  line1->SetLineColor(kGray+2);
  line1->SetLineStyle(2);
  TF1 *line2 = new TF1("line","-10.0",0.,comp->GetBinLowEdge(comp->GetNbinsX()+1));
  line2->SetLineColor(kGray+2);
  line2->SetLineStyle(2);
  
  int colors[10] = {kBlue,kMagenta,kGreen,kCyan,kRed,kBlue+3,kOrange+1,kSpring-1,kMagenta+3,kGreen+3};
  int color=0;
  for (int i=0; i<nHists; i++){
    TH1F *systHist = linearBin((TH1F*)HistList->At(i));
    systHist->Add(comp,-1);
    systHist->Divide(comp);
    systHist->Scale(100.);
    std::string systStr = systHist->GetName();
    int ind = systStr.rfind("cat0");
    std::string systName = systStr.substr(ind+5,systStr.size());
    systHist->SetLineColor(colors[color]);
    systHist->SetTitle(Form("%s",name.c_str()));
    systHist->GetYaxis()->SetTitle("Difference over nominal %");
    systHist->GetYaxis()->SetTitleOffset(1.4);
    systHist->GetXaxis()->SetTitle("BDT output");
    systHist->GetYaxis()->SetRangeUser(-100.,250);
    if (int(systName.find("Up"))>0){
      systHist->SetLineStyle(1);
      systName = systName.substr(0,systName.rfind("Up"));
      leg->AddEntry(systHist,systName.c_str(),"l");
      color++;
    }
    else if (int(systName.find("Down"))>0){
      systHist->SetLineStyle(2);
      systName = systName.substr(0,systName.rfind("Down"));
    }
    if (i==0) systHist->DrawCopy("hist");
    else systHist->DrawCopy("same hist");
  }
  leg->Draw("same");
  line->Draw("same");
  line1->Draw("same");
  line2->Draw("same");
  canv->Print(("plots/"+bdt+"/systs/"+name+".png").c_str(),"png");
  
  delete canv;
  
  systCalls++;
  
}

std::string getTime(){
  
  std::string timeStr;
  time_t now = time(0);
  struct tm* tm = localtime(&now);
  timeStr = Form("%02d:%02d_%02d-%02d-%02d",tm->tm_hour,tm->tm_min,tm->tm_mday,tm->tm_mon+1,tm->tm_year-100);

  return timeStr;

}

std::pair<TH1F*,TH1F*> GetPMSyst(TH1F* trueHist, TH1F* intHist, std::string name){
  TH1F *diff = (TH1F*)trueHist->Clone();
  diff->Add(intHist,-1);
  diff->Divide(intHist);
  diff->Scale(3.);
  TH1F *up = (TH1F*)intHist->Clone();
  up->Add(diff,1);
  up->SetName((name+"Up01_sigma").c_str());
  TH1F *down = (TH1F*)intHist->Clone();
  down->Add(diff,-1);
  down->SetName((name+"Down01_sigma").c_str());

  std::pair<TH1F*,TH1F*> result(up,down);
  return result; 
}

int BDTInterpolation(std::string inFileName,bool Diagnose=false, bool doNorm=true, bool doSidebands=false){

  std::cout << getTime() << std::endl;

  system("rm -r plots");
  system("mkdir plots");

  //gROOT->SetBatch();
  gStyle->SetOptStat(0);

  // input flags
  bool all=0;
 
  if (Diagnose) {
    std::cout << "Diagnostics turned on \n Output plots available in \"plots\" directory \n Diagnostic log available in \"plots/BDTInterpolationDiagnostics.txt\"" << std::endl;
    all = true;
  }
  else std::cout << "Diagnostics turned off" << std::endl;
  if (doNorm) std::cout << "Normalization turned on" << std::endl;
  else std::cout << "Normalization turned off" << std::endl;
  if (doSidebands) std::cout << "Background model from sidebands turned on" << std::endl;
  else std::cout << "Background model from sidebands turned off" << std::endl;

  TFile *inFile = new TFile(inFileName.c_str());
  //TFile *inFile = new TFile("/vols/cms02/nw709/hgg/src_cvs/oct13/CMSSW_4_2_8/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/Macros/CMS-HGG_1658pb_mva.root");
  //TFile *inFile = new TFile("RefWorkspaces/CMS-HGG_1658pb_mva.root");
  TFile *outFile = new TFile(Form("%s_interpolated.root",inFileName.c_str()),"RECREATE");
  TFile *systTestF = new TFile("systTest.root","RECREATE");
  ofstream diagFile("plots/BDTInterpolationDiagnostics.txt");

  const int nBDTs=2;
  const int nMasses=7;
  std::string BDTtype[nBDTs] = {"ada","grad"};
  std::string BDTmasses[nMasses] = {"115.0","120.0","125.0","130.0","135.0","140.0","150.0"};

// ----- else can just get rebinned histograms straight out of workspace

  std::cout << "Extracting histograms from workspace........." << std::endl;

  diagFile << "Diagnostics for Signal Interpolation run at " << getTime() << std::endl;
  diagFile << "---------------------------------------------------------------------" << std::endl;
  diagFile << "Following orginal histograms rewritten into new workspace: " << std::endl;
  
  // ------ stuff for interpolation systematic -------
  TH1F *systHists120[2][3];
  TH1F *systHists135[2][3];

  std::string syst120mass[3] = {"120.0","115.0","125.0"};
  std::string syst135mass[3] = {"135.0","130.0","140.0"};
  // ---------------------------------------------------
  
  // make plots of background model from sidebands
  if (doSidebands){
    for (int bdt=0; bdt<nBDTs; bdt++){
      for (double mass=115.; mass<150.5; mass+=0.5){
        TList *bkgModelList = new TList();
        TH1F *bkgModel[5];
        bkgModel[0] = (TH1F*)inFile->Get(Form("th1f_bkg_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
        bkgModel[1] = (TH1F*)inFile->Get(Form("th1f_bkg_2low_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
        bkgModel[2] = (TH1F*)inFile->Get(Form("th1f_bkg_1low_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
        bkgModel[3] = (TH1F*)inFile->Get(Form("th1f_bkg_1high_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
        bkgModel[4] = (TH1F*)inFile->Get(Form("th1f_bkg_2high_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
        for (int i=0; i<5; i++) {
          bkgModelList->Add(bkgModel[i]);
        }
        std::string name(Form("%3.1f",mass));
        plotBkgModel(bkgModelList,name);
      }
    }
  }
  
  // write original histograms in out file
  TList *HistList = inFile->GetListOfKeys();
  for (int j=0; j<HistList->GetSize(); j++){
    TH1F *temp = (TH1F*)inFile->Get(HistList->At(j)->GetName());
    TString name = temp->GetName();

    // store stuff for interpolation systematic
    for (int bdt=0; bdt<2; bdt++){
      for (int syst=0; syst<3; syst++){
        if (name.Contains(("th1f_sig_"+BDTtype[bdt]+"_"+syst120mass[0]+"_"+syst120mass[syst]).c_str()) && !name.Contains("sigma")){
          systHists120[bdt][syst] = (TH1F*)inFile->Get(name.Data());
          systHists120[bdt][syst]->SetLineColor(syst+2);
        }
        if (name.Contains(("th1f_sig_"+BDTtype[bdt]+"_"+syst135mass[0]+"_"+syst135mass[syst]).c_str()) && !name.Contains("sigma")){
          systHists135[bdt][syst] = (TH1F*)inFile->Get(name.Data());
          systHists135[bdt][syst]->SetLineColor(syst+2);
        }
      }
    }
    // ---------------------------------------------
    std::string tName = temp->GetName();
    for (int i=0; i<nMasses; i++){
      int ind = tName.find(BDTmasses[i]+"_"+BDTmasses[i]);
      if (ind>0) tName.replace(ind,11,BDTmasses[i]);
    }
    temp->SetName(tName.c_str());
      
    outFile->cd();
    temp->Write();
    diagFile << "Histo written: " << temp->GetName() << std::endl;
  }

  // ----------- Do stuff for interpolation systematic ------------
  std::cout << "Creating interpolation systematic templates......." << std::endl;

  diagFile << "---------------------------------------------------------------------" << std::endl;
  diagFile << "Writing following interpolation systematic templates: " << std::endl;
  
  for (int bdt=0; bdt<2; bdt++){
    double norm120 = GetNorm(115.0,systHists120[bdt][1],125.0,systHists120[bdt][2],120.0);
    double norm135 = GetNorm(130.0,systHists135[bdt][1],140.0,systHists135[bdt][2],135.0);
    if (doNorm) {
      systHists120[bdt][1]->Scale(1./(GetXsection(115.0)*GetBR(115.0)));
      systHists120[bdt][2]->Scale(1./(GetXsection(125.0)*GetBR(125.0)));
      systHists135[bdt][1]->Scale(1./(GetXsection(130.0)*GetBR(130.0)));
      systHists135[bdt][2]->Scale(1./(GetXsection(140.0)*GetBR(140.0)));
    }
    TH1F *int120 = Interpolate(115.0,systHists120[bdt][1],125.0,systHists120[bdt][2],120.0);
    TH1F *int135 = Interpolate(130.0,systHists135[bdt][1],140.0,systHists135[bdt][2],135.0);
    if (doNorm) {
      systHists120[bdt][1]->Scale(GetXsection(115.0)*GetBR(115.0));
      systHists120[bdt][2]->Scale(GetXsection(125.0)*GetBR(125.0));
      systHists135[bdt][1]->Scale(GetXsection(130.0)*GetBR(130.0));
      systHists135[bdt][2]->Scale(GetXsection(140.0)*GetBR(140.0));
    }
    if (doNorm) {
      int120->Scale(norm120/int120->Integral());
      int135->Scale(norm135/int135->Integral());
    }
    
    TList *int120list = new TList();
    int120list->Add(systHists120[bdt][1]);
    int120list->Add(systHists120[bdt][2]);
    int120list->Add(int120);
    plotFrac(int120list,systHists120[bdt][0],"syst120",false);
    
    TList *int135list = new TList();
    int135list->Add(systHists135[bdt][1]);
    int135list->Add(systHists135[bdt][2]);
    int135list->Add(int135);
    plotFrac(int135list,systHists135[bdt][0],"syst135",false);
    TH1F *frac=(TH1F*)systHists135[bdt][0]->Clone();
    frac->Add(int135,-1);
    frac->Divide(int135);
    TH1F *linFrac = (TH1F*)linearBin(frac);
    linFrac->GetYaxis()->SetRangeUser(-0.4,0.4);
    TCanvas *c = new TCanvas();
    systTestF->cd();
    linFrac->Write();
    linFrac->Draw();
    TF1 *fit = new TF1("fit","[0]*(x-36.)^2",0.,36.);
    linFrac->Fit(fit,"q");
    fit->DrawCopy("same");
    fit->SetLineColor(kGray);
    fit->SetLineStyle(2);
    double var = fit->GetParameter(0);
    double err = fit->GetParError(0);
    fit->SetParameter(0,var+err);
    fit->DrawCopy("same");
    fit->SetParameter(0,var-err);
    fit->DrawCopy("same");

    c->Print(Form("plots/%s/fracs/systTest135.png",BDTtype[bdt].c_str()),"png");

    // get plus and minus templates
    std::pair<TH1F*,TH1F*> result135 = GetPMSyst(systHists135[bdt][0],int135,"th1f_sig_"+BDTtype[bdt]+"_135.0_cat0_sigInt");
    TH1F* sigIntDown = (TH1F*)linearBin(result135.first);
    TH1F* sigIntUp = (TH1F*)linearBin(result135.second);
    TH1F* sigIntCent = (TH1F*)linearBin(int135);
    diagFile << (result135.first)->GetName() << std::endl;
    diagFile << (result135.second)->GetName() << std::endl;
    
    TCanvas *test = new TCanvas();
    sigIntDown->SetLineColor(2);
    sigIntUp->SetLineColor(4);
    sigIntDown->Draw();
    sigIntUp->Draw("same");
    sigIntCent->Draw("same'");
    test->Print(("plots/"+BDTtype[bdt]+"/fracs/PMsysts135.png").c_str(),"png");
   
    outFile->cd();
    (result135.first)->Write();
    (result135.second)->Write();
   
  }

    // -------------------- systematic stuff done ----------------------

  // get lists of middle, upper and lower templates for each mass
  TList *orgHistList[2][7];
  TList *orgHistListBelow[2][7];
  TList *orgHistListAbove[2][7];
  for (int i=0; i<2; i++) {
    for (int j=0; j<7; j++) {
      orgHistList[i][j]=new TList();
      orgHistListBelow[i][j]=new TList();
      orgHistListAbove[i][j]=new TList();
    }
  }

  for (int j=0; j<HistList->GetSize(); j++){
    TString HistName(HistList->At(j)->GetName());
    for (int bdt=0; bdt<nBDTs; bdt++){
      for (int bdtmass=0; bdtmass<nMasses; bdtmass++){
        for (int k=-1; k<2; k++){
          if ((bdtmass==0 && k==-1) || (bdtmass==nMasses-1 && k==1)) continue;
          if (HistName.Contains(("sig_"+BDTtype[bdt]+"_"+BDTmasses[bdtmass]+"_"+BDTmasses[bdtmass+k]).c_str())){
            TH1F *temp = (TH1F*)inFile->Get(HistName.Data());
            if (k==-1) orgHistListBelow[bdt][bdtmass]->Add(temp);
            if (k==0) orgHistList[bdt][bdtmass]->Add(temp);
            if (k==1) orgHistListAbove[bdt][bdtmass]->Add(temp);
          }
        }
      }
    }
  }
  diagFile << "---------------------------------------------------------------------" << std::endl;
  diagFile << "Following histo's being used for interpolation: " << std::endl;
  for (int bdt=0; bdt<2; bdt++){
    for (int mass=0; mass<nMasses; mass++){
      diagFile << "BDT: " << BDTtype[bdt] << std::endl;
      diagFile << "Mass: " << BDTmasses[mass] << std::endl;
      for (int syst=0; syst<orgHistList[bdt][mass]->GetSize(); syst++){
        diagFile << "   Central: " << orgHistList[bdt][mass]->At(syst)->GetName() << std::endl;
        if (mass!=0) diagFile << "   Lower:   " << orgHistListBelow[bdt][mass]->At(syst)->GetName() << std::endl;
        if (mass!=nMasses-1) diagFile << "   Upper:   " << orgHistListAbove[bdt][mass]->At(syst)->GetName() << std::endl;
      }
    }
  }

  // now have orgHistList for each BDT and above and below signals.
  // the contain 17 histos - 1 for actual sig and 8 up and down for each systematic

  // ---------- at this points have signal hists for each mass points as well as +- models all rebinned at mass points --------

  std::cout << "Calculating interpolated signal templates........" << std::endl;
  diagFile << "---------------------------------------------------------------------" << std::endl;
  diagFile << "Interpolating intermediate signals from following lower and upper templates" << std::endl;

  TList *orgHistListInt[2][71];
  for (int i=0; i<2; i++) for (int j=0; j<71; j++) orgHistListInt[i][j] = new TList();
  TH1F *systUp, *systDown, *systTrue;
  TH1F *background, *data, *signal;

  int i=0;
  // loop over mass points etc.
  for (double mass=115.; mass<150.5; mass+=0.5){
    if (int(mass)%2==0) std::cout << Form("%3.0f",((mass-115.)/35.)*100.) << "% done" << std::endl;
    //points we have signal for
    if (int(mass*2)%10==0 && mass!=145.0) {
      for (int bdt=0; bdt<2; bdt++){  
        background = (TH1F*)inFile->Get(Form("th1f_bkg_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
        data = (TH1F*)inFile->Get(Form("th1f_data_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
        signal = (TH1F*)inFile->Get(Form("th1f_sig_%s_%3.1f_%3.1f_cat0",BDTtype[bdt].c_str(),mass,mass));
        std::string name = Form("%3.1f",mass);
        if (Diagnose) {
          plotDavid(background,signal,data,name);
          TList *systList = new TList();
          TH1F *central;
          int bdtmass = getIndex((findNearest(mass).first));
          for (int syst=0; syst<orgHistList[bdt][bdtmass]->GetSize(); syst++){
            if (syst==0) central = (TH1F*)orgHistList[bdt][bdtmass]->At(syst);
            else {
              TH1F *tempSig = (TH1F*)orgHistList[bdt][bdtmass]->At(syst);
              systList->Add(tempSig);
            }
          }
          plotSystFracs(systList,central,Form("%3.1f_systFracs",mass));
        }
        // store some systematic histograms
        if (int(mass)>115 && int(mass)<150){
          int bdtmass = getIndex(int(mass));
          systTrue = (TH1F*)orgHistList[bdt][bdtmass]->At(0)->Clone();
          systUp = (TH1F*)orgHistListAbove[bdt][bdtmass]->At(0)->Clone();
          systDown = (TH1F*)orgHistListBelow[bdt][bdtmass]->At(0)->Clone();

          double trueNorm = systTrue->Integral();
          systUp->Scale(trueNorm/systUp->Integral());
          systDown->Scale(trueNorm/systDown->Integral());
          systUp->SetName(Form("th1f_sig_%s_%3.1f_cat0_sigIntNewUp01_sigma",BDTtype[bdt].c_str(),mass));
          systDown->SetName(Form("th1f_sig_%s_%3.1f_cat0_sigIntNewDown01_sigma",BDTtype[bdt].c_str(),mass));
          
          outFile->cd();
          systUp->Write();
          systDown->Write();
          TCanvas *lil = new TCanvas();
          systUp->SetLineColor(kRed);
          systDown->SetLineColor(kBlue);
          systDown->Draw("e");
          systUp->Draw("same e");
          systTrue->Draw("same e");
          lil->Print(Form("plots/%s/systNew_%3.1f.png",BDTtype[bdt].c_str(),mass),"png");
          delete lil;
        }
      }
      continue;
    }
    std::pair<int,int> nearestPair = findNearest(mass);
    bool above;
    int nearest = nearestPair.first;
    int nextNear =nearestPair.second;
    if (nearest-nextNear < 0) above = true;
    else above = false;
    
    //std::cout << mass << " bracketed by: " << nearestPair.first << " " << nearestPair.second << " " << above << std::endl;

    int bdtmass = getIndex(nearest); // gives index of bdt to use

    // loop bdt type
    for (int bdt=0; bdt<2; bdt++){
      diagFile << "Mass: " << mass << std::endl;
      diagFile << "BDT: " << BDTtype[bdt] << std::endl;
      // loop different histos in list (signal and systematics)
      background = (TH1F*)inFile->Get(Form("th1f_bkg_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
      data = (TH1F*)inFile->Get(Form("th1f_data_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
      TList *systList = new TList();
      TH1F *central;
      for (int syst=0; syst<orgHistList[bdt][bdtmass]->GetSize(); syst++){
        TH1F *tempSig = (TH1F*)orgHistList[bdt][bdtmass]->At(syst);
        central = (TH1F*)tempSig->Clone();
        TH1F *tempAbove = (TH1F*)orgHistListAbove[bdt][bdtmass]->At(syst);
        TH1F *tempBelow = (TH1F*)orgHistListBelow[bdt][bdtmass]->At(syst);
        TH1F* tempInt;
        TList* plotList = new TList();
        if (above){
          if (doNorm) tempSig->Scale(1./(GetXsection(double(nearest))*GetBR(double(nearest))));
          if (doNorm) tempAbove->Scale(1./(GetXsection(double(nextNear))*GetBR(double(nextNear))));
          tempInt = Interpolate(double(nearest),tempSig,double(nextNear),tempAbove,mass);
          if (doNorm) tempSig->Scale(GetXsection(double(nearest))*GetBR(double(nearest)));
          if (doNorm) tempAbove->Scale(GetXsection(double(nextNear))*GetBR(double(nextNear)));
          double norm = GetNorm(double(nearest),tempSig,double(nextNear),tempAbove,mass);
          if (doNorm) tempInt->Scale(norm/tempInt->Integral());
          // diagnostic stuff
          diagFile << "   Interpolated:     " << tempInt->GetName() << std::endl;
          diagFile << "   from:             " << tempSig->GetName() << std::endl;
          diagFile << "   and:              " << tempAbove->GetName() << std::endl;
          std::string histName = tempSig->GetName();
          std::string systName = histName.substr(histName.rfind("cat0"),histName.size());
          plotList->Add(tempSig);
          plotList->Add(tempAbove);
          if (Diagnose){
            if (syst==0) plotFrac(plotList,tempInt,Form("%3.1f_%s",mass,systName.c_str()),true);
            else systList->Add(tempInt);
          }
        }
        else{
          if (doNorm) tempSig->Scale(1./(GetXsection(double(nearest))*GetBR(double(nearest))));
          if (doNorm) tempBelow->Scale(1./(GetXsection(double(nextNear))*GetBR(double(nextNear))));
          tempInt = Interpolate(double(nextNear),tempBelow,double(nearest),tempSig,mass);
          if (doNorm) tempSig->Scale(GetXsection(double(nearest))*GetBR(double(nearest)));
          if (doNorm) tempBelow->Scale(GetXsection(double(nextNear))*GetBR(double(nextNear)));
          double norm = GetNorm(double(nextNear),tempBelow,double(nearest),tempSig,mass);
          if (doNorm) tempInt->Scale(norm/tempInt->Integral());
          // diagnostic stuff
          diagFile << "   Interpolated:     " << tempInt->GetName() << std::endl;
          diagFile << "   from:             " << tempBelow->GetName() << std::endl;
          diagFile << "   and:              " << tempSig->GetName() << std::endl;
          std::string histName = tempSig->GetName();
          std::string systName = histName.substr(histName.rfind("cat0"),histName.size());
          plotList->Add(tempBelow);
          plotList->Add(tempSig);
          if (Diagnose){
            if (syst==0) plotFrac(plotList,tempInt,Form("%3.1f_%s",mass,systName.c_str()),true);
            else systList->Add(tempInt);
          }
        }
        orgHistListInt[bdt][i]->Add(tempInt);
        std::string name = Form("%3.1f",mass);
        if (syst==0 && Diagnose) plotDavid(background,tempInt,data,name);
        delete plotList;
      }
      if (Diagnose) plotSystFracs(systList,central,Form("%3.1f_systFracs",mass));
    }
    i++;
  }
 
  outFile->cd();

  diagFile << "---------------------------------------------------------------------" << std::endl;
  diagFile << "Writing following interpolated histograms to file: " << std::endl;
  for (int l=0; l<2; l++) for (int j=0; j<71; j++) for (int k=0; k<orgHistListInt[l][j]->GetSize(); k++) {
    orgHistListInt[l][j]->At(k)->Write();
    diagFile << orgHistListInt[l][j]->At(k)->GetName() << std::endl;
  }
  
  TList* endList = outFile->GetListOfKeys();
  for (double mass=115.0; mass<150.5; mass+=0.5){
    diagFile << mass << std::endl;
    std::pair<int,int> nearestPair = findNearest(mass);
    double nearest = nearestPair.first;
    double nextNear = nearestPair.second;
    for (int k=0; k<endList->GetSize(); k++){
      TString hName = endList->At(k)->GetName();
      if (hName.Contains("sigma") || !hName.Contains("sig")) continue;
      if (hName.Contains(Form("a_%3.1f_cat0",nearest)) || hName.Contains(Form("d_%3.1f_cat0",nearest)) || hName.Contains(Form("%3.1f_%3.1f_cat0",nearest,nextNear)) || hName.Contains(Form("%3.1f_cat0",mass))) {
        TH1F *temp = (TH1F*)outFile->Get(hName);
        diagFile << "    " << setw(30) << temp->GetName();
        diagFile << "    " << setw(8) << temp->GetEntries();
        diagFile << "    " << setw(8) << temp->Integral() << std::endl;
      }
    }
  }
    
  outFile->Close();

  std::cout << "Diagnostics log written to \"BDTInterpolationDiagnostics.txt\"" << std::endl;
  if (davidCalls>0) std::cout << davidCalls << " plots written to: plots/ada/david/ \n                    plots/grad/david/ " << std::endl;
  if (davidCalls>0) std::cout << davidCalls << " plots written to: plots/ada/diff/ \n                    plots/grad/diff/ " << std::endl;
  if (fracCalls>0) std::cout << fracCalls << " plots written to: plots/ada/frac/ \n                      plots/grad/frac/" << std::endl;
  if (systCalls>0) std::cout << systCalls << " plots written to: plots/ada/syst/ \n                    plots/grad/syst/ " << std::endl;
  if (bkgCalls>0) std::cout << bkgCalls << " plots written to: plots/ada/bkgMod/ \n                    plots/grad/bkgMod/ " << std::endl;

  std::string TIME_DATE = getTime();
  if (Diagnose){
    system("whoami > temp.txt");
    std::ifstream temp("temp.txt");
    std::string result;
    temp >> result;
    temp.close();
    system("rm temp.txt");
    system(("python make_html.py "+result+" "+TIME_DATE+" "+inFileName).c_str());
    system(("mkdir -p ~/public_html/h2g/MVA/SigInt/Diagnostics/"+TIME_DATE+"/").c_str());
    system(("cp -r plots ~/public_html/h2g/MVA/SigInt/Diagnostics/"+TIME_DATE+"/").c_str());
    system(("rm ~/public_html/h2g/MVA/SigInt/Diagnostics/Current"));
    system(("ln -s ~/public_html/h2g/MVA/SigInt/Diagnostics/"+TIME_DATE+" ~/public_html/h2g/MVA/SigInt/Diagnostics/Current").c_str());
    std::cout << ("Plots avaiable to view in ~/public_html/h2g/MVA/SigInt/Diagnostics/"+TIME_DATE+"/").c_str() << std::endl;
    std::cout << "If working on /vols/ at IC plots avaliable to view at www.hep.ph.ic.ac.uk/~"+result+"/h2g/MVA/SigInt/Diagnostics/Current/plots/plots.html" << std::endl;
  }
  std::cout << "Checking all relevant histograms have been written to workspace......" << std::endl;
  //system(("python checkOK.py "+string(outFile->GetName())).c_str());
  std::cout << "New workspace written to " << outFile->GetName() << std::endl;
  return 0;
}

