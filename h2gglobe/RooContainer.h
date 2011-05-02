#ifndef ROOCONTAINER
#define ROOCONTAINER

// RooFit includes
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"

// RooStats includes
#include "RooWorkspace.h"

#include <map>
#include <vector>
// For now only conatins expnentials and real vars and datasets.
// Abstract to addPdf soon

class RooContainer {

  public:

   RooContainer(int ncat=1);
   ~RooContainer(){};

   void AddRealVar(std::string,float,float);
   void AddRealVar(std::string,float,float,float);
   void AddGenericPdf(std::string,std::string,
		      std::vector<std::string> &, 
		      double norm_guess=10);
   void ComposePdf(std::string, std::string
			     ,std::vector<std::string> &);

   void CreateDataSet(std::string);
   
   void FitToData(std::string,std::string,int bins=-1);
   void FitToData(std::string,std::string 
	         ,double, double, double, double,int bins=-1);
   void SetRealVar(std::string, int, float, float w=1.);
   void Save();
   
  private:

   void addRealVar(std::string,float,float);
   void addRealVar(std::string,float,float,float);
   void addGenericPdf(std::string,std::string,
		      std::vector<std::string> &, 
		      double norm_guess=10);
   void composePdf(std::string , std::string 
			     ,std::vector<std::string> &);

   void createDataSet(std::string);
   
   void fitToData(std::string,std::string,int bins=-1);
   void fitToData(std::string,std::string
	         ,double, double, double, double,int bins=-1);
   void writeRooDataHist(std::string, RooDataSet *);
   void writeRooPlot(RooPlot *);
   int ncat;
   std::string getcatName(std::string, int);   

   std::map<std::string, RooRealVar> m_real_var_;
   std::map<std::string, RooExtendPdf> m_exp_;
   std::map<std::string, RooGenericPdf> m_gen_;
   std::map<std::string, RooAddPdf> m_pdf_;

   std::map<std::string, float> m_var_min_;
   std::map<std::string, float> m_var_max_;

   std::map<std::string,RooDataSet*> data_; 
   std::map<RooPlot*,RooFitResult*> fit_res_;

   RooWorkspace ws;

};


#endif
