#include <iostream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGenericPdf.h"
#include "RooPlot.h"

#include "Fit.C"
#include "Rebin.C"

using namespace std;
using namespace RooFit;

const double mLow=100.;
const double mHigh=180.;
const double sidebandWidth=0.02;
const int nMasses=9;

int getIndex(int mass){
  int index=-1;
  if (mass==105) index=0;
  else if (mass==110) index=1;
  else if (mass==115) index=2;
  else if (mass==120) index=3;
  else if (mass==125) index=4;
  else if (mass==130) index=5;
  else if (mass==135) index=6;
  else if (mass==140) index=7;
  else if (mass==150) index=8;
  else {
    cout << mass << " is an invalid mass. Bailing out." << endl;
    exit(1);
  }
  return index;
}

int getMass(int index){
  
  int mass=-1;
  if (index==0) mass=105;
  else if (index==1) mass=110;
  else if (index==2) mass=115;
  else if (index==3) mass=120;
  else if (index==4) mass=125;
  else if (index==5) mass=130;
  else if (index==6) mass=135;
  else if (index==7) mass=140;
  else if (index==8) mass=150;
  else {
    cout << "index " << index << " not found. Bailing out." << endl;
    exit(1);
  }
  return mass;
}


void mergeHistograms(std::string nameHist, TH1F* hist1, TH1F* hist2){
   
   //Get Bin Low edges of histogram 1
   int nbins1 = hist1->GetNbinsX();
   int nbins2 = hist2->GetNbinsX();
   int nbinsTot = nbins1+nbins2;

   double *arrBins1 = new double[nbins1];
   double *arrBins2 = new double[nbins2+1];
   double *arrBinsTot = new double[nbinsTot+1];

   for (int i=1;i<=nbins1;i++){
     arrBinsTot[i-1]=hist1->GetBinLowEdge(i);
   }
   for (int i=1;i<=nbins2+1;i++){	// Include upper edge in 2nd hist
     arrBinsTot[i+nbins1-1]=hist2->GetBinLowEdge(i);
   }

   const char *histoname = hist1->GetName();
   const char *histotitle = hist1->GetTitle();

   TH1F *newHist = new TH1F(Form("NUMPTYNAME%s",hist1->GetName()),histotitle,nbinsTot,arrBinsTot);
   newHist->SetName(hist1->GetName());
   for (int i=1;i<=nbins1;i++){
	newHist->SetBinContent(i,hist1->GetBinContent(i));
	newHist->SetBinError(i,hist1->GetBinError(i));
   } 
   for (int i=1;i<=nbins2;i++){
	newHist->SetBinContent(i+nbins1,hist2->GetBinContent(i));
	newHist->SetBinError(i+nbins1,hist2->GetBinError(i));
   } 

   std::cout << "RooContainer::MergeHistograms -- Replacing histogram - " 
	     << newHist->GetName()
	     << std::endl;
   // Now the dangerous part!
   *hist1 = *newHist;
}

void makeSignalOutputHistogram(TFile *inFile, string new_name, string old_name, vector<double> binEdges, string vbf_name, vector<double> binEdgesVBF, bool includeVBF){

  TH1F *old = (TH1F*)inFile->Get(old_name.c_str());
  TH1F *rebinned = rebinBinnedDataset(new_name,old,binEdges);
  
  string new_vbf_name = vbf_name;
  new_vbf_name.replace(new_vbf_name.find("VBF"),3,"vbf");

  if (includeVBF){
    TH1F *vbf = (TH1F*)inFile->Get(vbf_name.c_str());
    TH1F *vbf_rebinned = rebinBinnedDataset(new_vbf_name,vbf,binEdgesVBF);
    mergeHistograms(new_name,rebinned,vbf_rebinned);
  }
  rebinned->Write();

  for (int i=0; i<nSyst; i++){
    
    TH1F *down = (TH1F*)inFile->Get(Form("%s_%sDown01_sigma",old_name.c_str(),systematics[i].c_str()));
    TH1F *up = (TH1F*)inFile->Get(Form("%s_%sUp01_sigma",old_name.c_str(),systematics[i].c_str()));
    TH1F *down_rebinned = rebinBinnedDataset(Form("%s_%sDown01_sigma",new_name.c_str(),systematics[i].c_str()),down,binEdges);
    TH1F *up_rebinned = rebinBinnedDataset(Form("%s_%sUp01_sigma",new_name.c_str(),systematics[i].c_str()),up,binEdges);
    
    if (includeVBF){
      TH1F *down_vbf = (TH1F*)inFile->Get(Form("%s_%sDown01_sigma",vbf_name.c_str(),systematics[i].c_str()));
      TH1F *up_vbf = (TH1F*)inFile->Get(Form("%s_%sUp01_sigma",vbf_name.c_str(),systematics[i].c_str()));
      TH1F *down_vbf_rebinned = rebinBinnedDataset(Form("%s_%sDown01_sigma",new_vbf_name.c_str()),down_vbf,binEdgesVBF);
      TH1F *up_vbf_rebinned = rebinBinnedDataset(Form("%s_%sUp01_sigma",new_vbf_name.c_str()),up_vbf,binEdgesVBF);
      
      mergeHistograms(Form("%s_%sDown01_sigma",new_name.c_str(),systematics[i].c_str()),down_rebinned,down_vbf_rebinned);
      mergeHistograms(Form("%s_%sUp01_sigma",new_name.c_str(),systematics[i].c_str()),up_rebinned,up_vbf_rebinned);
    }
    down_rebinned->Write();
    up_rebinned->Write();
    
  }

}

void makeOutputHistogram(TFile *inFile, string new_name, string old_name, vector<double> binEdges, string vbf_name, vector<double> binEdgesVBF, bool includeVBF){
  TH1F *old = (TH1F*)inFile->Get(old_name.c_str());
  TH1F *rebinned = rebinBinnedDataset(new_name,old,binEdges);
  
  if (includeVBF){
    TH1F *vbf = (TH1F*)inFile->Get(vbf_name.c_str());
    TH1F *vbf_rebinned = rebinBinnedDataset(Form("vbf_%s",new_name.c_str()),vbf,binEdgesVBF);
    mergeHistograms(new_name,rebinned,vbf_rebinned);
  }
  rebinned->Write();

}

void checkMass(int mass_hyp){
  if (mass_hyp!=110  && mass_hyp!=115 && mass_hyp!=120 && mass_hyp!=125 && mass_hyp !=130 && mass_hyp!=135 && mass_hyp!=140 && mass_hyp!=150){
    cout << Form("Invalid mass: %d",mass_hyp) << endl;
    exit(1);
  }
}

pair<int,int> getUpDownHyp(int mass){
  checkMass(mass);
  int down = mass-5;
  int up = mass+5;
  if (mass==140) up = 150;
  if (mass==150) {
    down = 140;
    up   = 150;
  }
  if (mass==110) down = 100;
}

pair<double,double> getUpDown(int mass){
  double mass_h_low=2.5;
  double mass_h_high=2.1;
  if (mass==140){ 
    mass_h_low  =-2.5;
    mass_h_high =4.6;
  } 
  else if (mass==150){  
    mass_h_low  =-5.0;
    mass_h_high =0.1;
  } 
  else {
    mass_h_low  =-2.5;
    mass_h_high =2.1;
  }

  return pair<double,double>(mass_h_low,mass_h_high);
}

TH1F *getSignalComb(int mass_hyp,string type="grad",int cat){
  
  if (type=="grad") type="BDT_grad";
  else if (type=="ada") type="BDT_ada";
  else if (type=="VBF") type="VBF";
  else {
    cout << type << " is not a valid signal combination type. Bailing out." << endl;
    exit(1);
  }
  TH1F *sig = ((TH1F*)inFile->Get(Form("sig_%s_ggh_%3d.0_cat%d",type.c_str(),mass_hyp,cat)))->Clone();;
  sig->Add((TH1F*)inFile->Get(Form("sig_%s_vbf_%3d.0_cat%d",type.c_str(),mass_hyp,cat)));
  sig->Add((TH1F*)inFile->Get(Form("sig_%s_wzh_%3d.0_cat%d",type.c_str(),mass_hyp,cat)));
  sig->Add((TH1F*)inFile->Get(Form("sig_%s_tth_%3d.0_cat%d",type.c_str(),mass_hyp,cat)));
  sig->SetName(Form("sig_%s_all_%3d.0_cat%d",type.c_str(),mass_hyp,cat));
  return sig;
}

void FitAndRebin(string filename,int mass_hyp=120, int cat=0){

  checkMass(mass_hyp);

  TFile *inFile = TFile::Open(filename.c_str(),"UPDATE");

  // various stuff that needs to be grabbed hold of
  bool includeVBF=true;
  int numberOfSidebands=6;
  int nSyst=8;

  // define stuff for fit
  RooWorkspace *inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
  RooRealVar *massVar    = (RooRealVar*)inWS->var("CMS_hgg_mass");
  RooDataSet *data    = (RooDataSet*)inWS->data(Form("data_mass_cat%d",cat));
  double fitNormalisation = -1.;
  double fitNormError = -1.;

  // get histograms for binning algorithm
  TH1F *bkgForBinningGRAD = (TH1F*)inFile->Get(Form("bkg_BDT_grad_%3d.0_cat%d",binning_mass,cat));
  TH1F *sigForBinningGRAD = getSignalComb(binning_mass,"grad",cat);
  
  TH1F *bkgForBinningADA = (TH1F*)inFile->Get(Form("bkg_BDT_ada_%3d.0_cat%d",binning_mass,cat));
  TH1F *sigForBinningADA = getSignalComb(binning_mass,"ada",cat);

  TH1F *bkgForBinningVBF = (TH1F*)inFile->Get(Form("bkg_VBF_%3d.0_cat%d",binning_mass,cat));
  TH1F *sigForBinningVBF = getSignalComb(binning_mass,"vbf",cat);
  
  // calculate bin edges for central mass and masses above and below
  vector<double> binEdgesGrad = significanceOptimizedBinning(sigForBinningGRAD,bkgForBinningGRAD,10);
  vector<double> binEdgesAda = significanceOptimizedBinning(sigForBinningADA,bkgForBinningADA,10);
  vector<double> binEdgesVBF;
  binEdgesVBF.push_back(1.);
  binEdgesVBF.push_back(1.04);

  // get mass iterator for background and data
  pair<double,double> massHL = getUpDown(mass_hyp);
  double mass_h_low = massHL.first;
  double mass_h_high = massHL.second;

  // background and data
  for (double mass_diff=mass_h_low; mass_diff<mass_h_high; mass_diff+=0.5){
    double mass = mass_hyp+mass_diff;
    
    // fit with power law
    pair<double,double> fitRes = FitPow(data,massVar,mass,true);
    fitNormalisation = fitRes.first;
    fitNormError = fitRes.second;

    // rebin bkg and data and add VBF on the end 
    makeOutputHistogram(inFile,Form("bkg_mc_grad_%3.1f_cat%d",mass,cat),Form("bkg_BDT_grad_%3.1f_cat%d",mass,cat),binEdgesGrad,Form("bkg_VBF_%3.1f_cat%d",mass,cat),binEdgesVBF,includeVBF);
    makeOutputHistogram(inFile,Form("data_grad_%3.1f_cat%d",mass,cat),Form("data_BDT_grad_%3.1f_cat%d",mass,cat),binEdgesGrad,Form("data_VBF_%3.1f_cat%d",mass,cat),binEdgesVBF,includeVBF);
    makeOutputHistogram(inFile,Form("bkg_mc_ada_%3.1f_cat%d",mass,cat),Form("bkg_BDT_ada_%3.1f_cat%d",mass,cat),binEdgesGrad,Form("bkg_VBF_%3.1f_cat%d",mass,cat),binEdgesVBF,includeVBF);
    makeOutputHistogram(inFile,Form("data_ada_%3.1f_cat%d",mass,cat),Form("data_BDT_ada_%3.1f_cat%d",mass,cat),binEdgesGrad,Form("data_VBF_%3.1f_cat%d",mass,cat),binEdgesVBF,includeVBF);
   
    // rebin sidebands for bkg mc and bkg model (from data)
    for (int sideband_i=1; sideband_i<=numberOfSidebands; sideband_i++){
      makeOutputHistogram(inFile,Form("bkg_mc_%dhigh_grad_%3.1f_cat%d",sideband_i,mass,cat),Form("bkg_%dhigh_BDT_grad_%3.1f_cat%d",sideband_i,mass,cat),binEdgesGrad,Form("bkg_%dhigh_VBF_%3.1f_cat%d",sideband_i,mass,cat),binEdgesVBF,includeVBF);
      makeOutputHistogram(inFile,Form("bkg_mc_%dlow_grad_%3.1f_cat%d",sideband_i,mass,cat),Form("bkg_%dlow_BDT_grad_%3.1f_cat%d",sideband_i,mass,cat),binEdgesGrad,Form("bkg_%dlow_VBF_%3.1f_cat%d",sideband_i,mass,cat),binEdgesVBF,includeVBF);
      makeOutputHistogram(inFile,Form("bkg_%dhigh_grad_%3.1f_cat%d",sideband_i,mass,cat),Form("data_%dhigh_BDT_grad_%3.1f_cat%d",sideband_i,mass,cat),binEdgesGrad,Form("data_%dhigh_VBF_%3.1f_cat%d",sideband_i,mass,cat),binEdgesVBF,includeVBF);
      makeOutputHistogram(inFile,Form("bkg_%dlow_grad_%3.1f_cat%d",sideband_i,mass,cat),Form("data_%dlow_BDT_grad_%3.1f_cat%d",sideband_i,mass,cat),binEdgesGrad,Form("data_%dlow_VBF_%3.1f_cat%d",sideband_i,mass,cat),binEdgesVBF,includeVBF);
      makeOutputHistogram(inFile,Form("bkg_mc_%dhigh_ada_%3.1f_cat%d",sideband_i,mass,cat),Form("bkg_%dhigh_BDT_ada_%3.1f_cat%d",sideband_i,mass,cat),binEdgesGrad,Form("bkg_%dhigh_VBF_%3.1f_cat%d",sideband_i,mass,cat),binEdgesVBF,includeVBF);
      makeOutputHistogram(inFile,Form("bkg_mc_%dlow_ada_%3.1f_cat%d",sideband_i,mass,cat),Form("bkg_%dlow_BDT_ada_%3.1f_cat%d",sideband_i,mass,cat),binEdgesGrad,Form("bkg_%dlow_VBF_%3.1f_cat%d",sideband_i,mass,cat),binEdgesVBF,includeVBF);
      makeOutputHistogram(inFile,Form("bkg_%dhigh_ada_%3.1f_cat%d",sideband_i,mass,cat),Form("data_%dhigh_BDT_ada_%3.1f_cat%d",sideband_i,mass,cat),binEdgesGrad,Form("data_%dhigh_VBF_%3.1f_cat%d",sideband_i,mass,cat),binEdgesVBF,includeVBF);
      makeOutputHistogram(inFile,Form("bkg_%dlow_ada_%3.1f_cat%d",sideband_i,mass,cat),Form("data_%dlow_BDT_ada_%3.1f_cat%d",sideband_i,mass,cat),binEdgesGrad,Form("data_%dlow_VBF_%3.1f_cat%d",sideband_i,mass,cat),binEdgesVBF,includeVBF);
    }
  }

  // now do signal (which has to include systematics)
  int i==getIndex(mass_hyp);
  for (int j=-1; j<2; j++){
    if ((i==1 && j==-1) || (i==nMasses-1 && j==1)) continue;
    int mass=getMass(i+j);

    // loop on signal processes
    // rebin signal (and systematics)
    makeSignalOutputHistogram(inFile,Form("sig_grad_ggh_%d.0_%d.0_cat%d",mass_hyp,mass,cat),Form("sig_BDT_grad_ggh_%d.0_cat%d",mass,cat),binEdgesGrad,Form("sig_VBF_ggh_%d.0_cat%d",mass,cat),binEdgesVBF,includeVBF);
    makeSignalOutputHistogram(inFile,Form("sig_grad_vbf_%d.0_%d.0_cat%d",mass_hyp,mass,cat),Form("sig_BDT_grad_vbf_%d.0_cat%d",mass,cat),binEdgesGrad,Form("sig_VBF_vbf_%d.0_cat%d",mass,cat),binEdgesVBF,includeVBF);
    makeSignalOutputHistogram(inFile,Form("sig_grad_wzh_%d.0_%d.0_cat%d",mass_hyp,mass,cat),Form("sig_BDT_grad_wzh_%d.0_cat%d",mass,cat),binEdgesGrad,Form("sig_VBF_wzh_%d.0_cat%d",mass,cat),binEdgesVBF,includeVBF);
    makeSignalOutputHistogram(inFile,Form("sig_grad_tth_%d.0_%d.0_cat%d",mass_hyp,mass,cat),Form("sig_BDT_grad_tth_%d.0_cat%d",mass,cat),binEdgesGrad,Form("sig_VBF_tth_%d.0_cat%d",mass,cat),binEdgesVBF,includeVBF);
    makeSignalOutputHistogram(inFile,Form("sig_ada_ggh_%d.0_%d.0_cat%d",mass_hyp,mass,cat),Form("sig_BDT_ada_ggh_%d.0_cat%d",mass,cat),binEdgesGrad,Form("sig_VBF_ggh_%d.0_cat%d",mass,cat),binEdgesVBF,includeVBF);
    makeSignalOutputHistogram(inFile,Form("sig_ada_vbf_%d.0_%d.0_cat%d",mass_hyp,mass,cat),Form("sig_BDT_ada_vbf_%d.0_cat%d",mass,cat),binEdgesGrad,Form("sig_VBF_vbf_%d.0_cat%d",mass,cat),binEdgesVBF,includeVBF);
    makeSignalOutputHistogram(inFile,Form("sig_ada_wzh_%d.0_%d.0_cat%d",mass_hyp,mass,cat),Form("sig_BDT_ada_wzh_%d.0_cat%d",mass,cat),binEdgesGrad,Form("sig_VBF_wzh_%d.0_cat%d",mass,cat),binEdgesVBF,includeVBF);
    makeSignalOutputHistogram(inFile,Form("sig_ada_tth_%d.0_%d.0_cat%d",mass_hyp,mass,cat),Form("sig_BDT_ada_tth_%d.0_cat%d",mass,cat),binEdgesGrad,Form("sig_VBF_tth_%d.0_cat%d",mass,cat),binEdgesVBF,includeVBF);
  }



  inFile->Close();
}
