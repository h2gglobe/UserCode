#ifndef ROOCONTAINER
#define ROOCONTAINER

// ROOT includes
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"

// RooFit includes
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooGlobalFunc.h"

// RooStats includes
#include "RooWorkspace.h"

// standard includes
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

class RooContainer {

  public:

   RooContainer(int ncat=1);
   ~RooContainer(){};

   void AddRealVar(std::string,float,float);
   void AddRealVar(std::string,float,float,float);
   void AddGenericPdf(std::string,std::string,
		      std::vector<std::string> &, 
		      int form,
		      double norm_guess=10);
   void ComposePdf(std::string, std::string
			     ,std::vector<std::string> &);

   void CreateDataSet(std::string,int nbins=-1); 
   void MakeSystematics(std::string,std::string);

   void FitToData(std::string,std::string 
	     ,double x1=-999,double x2=-999,double x3=-999,double x4=-999);
   void FitToSystematicSet(std::string,std::string,std::string
	     ,double x1=-999,double x2=-999,double x3=-999,double x4=-999);
   void InputDataPoint(std::string,int,float,float w=1.);
   void InputSystematicSet(std::string,std::string,int
		,std::vector<float>,float w=1.);
   void Save();
   
  private:

   void addRealVar(std::string,float,float);
   void addRealVar(std::string,float,float,float);
   void addGenericPdf(std::string,std::string,
		      std::vector<std::string> &,
		      int form, 
		      double norm_guess=10);
   void composePdf(std::string , std::string 
			     ,std::vector<std::string> &);

   void createDataSet(std::string,std::string,int);
   void makeSystematics(std::string,std::string);
   
   void fitToData(std::string,std::string,std::string
	         ,double,double,double,double);
   void fitToSystematicSet(std::string,std::string,std::string
	     ,double,double,double,double);
   void writeRooDataHist(std::string, TH1F *);
   void writeRooPlot(RooPlot *);
   int ncat;
   int nsigmas;

   std::string getcatName(std::string,int);   
   std::string getsysName(std::string,std::string);   
   std::string getsysindexName(std::string,std::string
			 ,int,int);
   
   std::vector<RooAbsPdf*> v_gen_;

   std::map<std::string, RooRealVar> m_real_var_;
   std::map<std::string, RooExtendPdf> m_exp_;
   std::map<std::string, RooAddPdf> m_pdf_;
   std::map<std::string,TH1F> m_th1f_;

   std::map<std::string, float> m_var_min_;
   std::map<std::string, float> m_var_max_;

   std::map<std::string,RooDataSet> data_; 

   //std::map<std::string,std::vector<RooRealVar*> > m_vars_up_; 
   //std::map<std::string,std::vector<RooRealVar*> > m_vars_dn_;

   std::map<std::string,std::vector<RooDataSet*> > data_up_; 
   std::map<std::string,std::vector<RooDataSet*> > data_dn_;

   std::map<std::string,std::vector<TH1F*> > m_th1f_up_; 
   std::map<std::string,std::vector<TH1F*> > m_th1f_dn_;
 
   std::map<std::string,int> bins_;
   std::map<std::string,float> inits_;
   std::map<RooPlot*,RooFitResult*> fit_res_;

   RooWorkspace ws;

};


#endif
