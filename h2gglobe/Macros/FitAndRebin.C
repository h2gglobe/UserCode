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

using namespace std;
using namespace RooFit;

const double mLow=100.;
const double mHigh=180.;
const double sidebandWidth=0.02;

pair<double,double> FitPow(RooDataSet* data, RooRealVar* var, double mass,bool diagnose=true){

  if (mass<110. || mass>150.) {
    cout << Form("%3.1f invalid mass. Not in range [110,150]",mass) << endl;
    exit(1);
  }

  RooRealVar *r1 = new RooRealVar(Form("r1_%3.1f",mass),Form("r1_%3.1f",mass),-8.,-20.,0.);
  RooRealVar *r2 = new RooRealVar(Form("r2_%3.1f",mass),Form("r2_%3.1f",mass),-4.,-20.,0.);
  RooRealVar *f1 = new RooRealVar(Form("f1_%3.1f",mass),Form("f1_%3.1f",mass),0.5,0.,1.);

  RooGenericPdf *fit  = new RooGenericPdf(Form("data_pow_model_m%3.1f",mass),Form("data_pow_model_m%3.1f",mass),"(1-@3)*TMath::Power(@0,@1) + @3*TMath::Power(@0,@2)",RooArgList(*var,*r1,*r2,*f1));

  double sidebL = mass*(1-sidebandWidth);
  double sidebH = mass*(1+sidebandWidth);
  var->setRange(Form("rangeLow_m%3.1f",mass),mLow,sidebL);
  var->setRange(Form("rangeHig_m%3.1f",mass),sidebH,mHigh);
  var->setRange(Form("sigReg_m%3.1f",mass),sidebL,sidebH);

  RooFitResult *fitRes = fit->fitTo(*data,Range(Form("rangeLow_m%3.1f",mass),Form("rangeHig_m%3.1f",mass)),Save(true),Strategy(1));

  if (diagnose){
    RooPlot *frame = var->frame();
    data->plotOn(frame,Binning(160));
    fit->plotOn(frame,Range(Form("rangeLow_m%3.1f",mass),Form("rangeHig_m%3.1f",mass)));
    TCanvas *c1 = new TCanvas();
    frame->Draw();
    c1->SaveAs(Form("fit_m%3.1f.pdf",mass));
  }
  // integral in sig region
  RooAbsReal *integral = fit->createIntegral(*var,NormSet(*var),Range(Form("sigReg_m%3.1f",mass)));
  // comb integral in two sideband regions 
  RooAbsReal *sidebandInt = fit->createIntegral(*var,NormSet(*var),Range(Form("rangeLow_m%3.1f",mass),Form("rangeHig_m%3.1f",mass)));
  double tempEv = data->sumEntries(Form("CMS_hgg_mass>=%3.1f && CMS_hgg_mass<%3.1f",mLow,sidebL))+data->sumEntries(Form("CMS_hgg_mass>%3.1f && CMS_hgg_mass<=%3.1f",sidebH,mHigh));
  RooConstVar *sidebandNevents = new RooConstVar(Form("sbEvents_m%3.1f",mass),Form("sbEvents_m%3.1f",mass),tempEv);
  RooFormulaVar normIntVar(Form("normIntVar_m%3.1f",mass),Form("normIntVar_m%3.1f",mass),"@0*@1/@2",RooArgSet(*sidebandNevents,*integral,*sidebandInt));

  double result     = normIntVar.getVal();
  double fullError  = normIntVar.getPropagatedError(*fitRes);

  return pair<double,double>(result,fullError);
  
}

void FitAndRebin(string filename,int mass_hyp=120){

  if (mass_hyp!=110  && mass_hyp!=115 && mass_hyp!=120 && mass_hyp!=125 && mass_hyp !=130 && mass_hyp!=135 && mass_hyp!=140 && mass_hyp!=150){
    cout << Form("Invalid mass: %d",mass_hyp) << endl;
    exit(1);
  }

  TFile *inFile = TFile::Open(filename.c_str());
  RooWorkspace *inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
  RooRealVar *mass    = (RooRealVar*)inWS->var("CMS_hgg_mass");
  RooRealVar *intLumi = (RooRealVar*)inWS->var("IntLumi");
  RooRealVar *bdt     = (RooRealVar*)inWS->var("BDT");
  RooRealVar *vbf     = (RooRealVar*)inWS->var("VBF");
  RooDataSet *data    = (RooDataSet*)inWS->data("data_mass_cat0");

  pair<double,double> fitRes = FitPow(data,mass,120.0,true);

  /*
  RooContainer *rooContainer = new RooContainer();
  //rooContainer->SetNCategories(1);
  rooContainer->nsigmas = 1;
  rooContainer->sigmaRange = 3;
  rooContainer->SaveRooDataHists(false);
  rooContainer->Verbose(false);
  rooContainer->AddObservable("CMS_hgg_mass",*mass);
  rooContainer->AddConstant("IntLumi",*intLumi);
  rooContainer->AddObservable("BDT",*bdt);
  rooContainer->AddObservable("VBF",*vbf);

  

  RooContainer::soverBOptimizedBinning
  */
  inFile->Close();
}
