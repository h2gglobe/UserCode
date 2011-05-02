//ROOT includes
#include "TROOT.h"
#include "TCanvas.h"

//RooFit includes
#include "RooDataHist.h"
#include "RooContainer.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooGlobalFunc.h"

//standard includes
#include <cmath>
#include <iostream>

using namespace RooFit;

RooContainer::RooContainer(int n):ncat(n){}

void RooContainer::AddRealVar(std::string name,float xmin,float xmax){
  for (int cat=0;cat<ncat;cat++){
    addRealVar(getcatName(name,cat),xmin,xmax);
  }
}

void RooContainer::AddRealVar(std::string name,float init, float vmin, float vmax){
  for (int cat=0;cat<ncat;cat++){
    addRealVar(getcatName(name,cat),init,vmin,vmax);
  }
}

void RooContainer::AddGenericPdf(std::string name,std::string formula
				,std::vector<std::string> & var
				,double norm_guess ){
  for (int cat=0;cat<ncat;cat++){
    std::vector<std::string> cat_var;
    for (std::vector<std::string>::iterator it=var.begin()
	;it!=var.end()
	;++it){
      cat_var.push_back(getcatName(*it,cat));
    }  
    addGenericPdf(getcatName(name,cat),formula,cat_var,norm_guess);
  }
}

void RooContainer::ComposePdf(std::string name, std::string  composition
			     ,std::vector<std::string> & formula){
  for (int cat=0;cat<ncat;cat++){
    std::vector<std::string> cat_formula;
    for (std::vector<std::string>::iterator it=formula.begin()
	;it!=formula.end()
	;++it){
      cat_formula.push_back(getcatName(*it,cat));
    }  
    composePdf(getcatName(name,cat),composition,cat_formula);
  }
}

void RooContainer::CreateDataSet(std::string name){
  for (int cat=0;cat<ncat;cat++){
    createDataSet(getcatName(name,cat));
  }
}

void RooContainer::FitToData(std::string name_func, std::string  name_var, int bins){
  for (int cat=0;cat<ncat;cat++){
    fitToData(getcatName(name_func,cat),getcatName(name_var,cat),bins);
  }
}
void RooContainer::FitToData(std::string name_func, std::string  name_var
			    ,double x1, double x2, double x3, double x4, int bins){
  for (int cat=0;cat<ncat;cat++){
    fitToData(getcatName(name_func,cat),getcatName(name_var,cat)
	     ,x1,x2,x3,x4,bins);
  }
}

void RooContainer::Save(){

  std::cout << "RooContainer::Save -- Saving To File "
            << std::endl;

  std::map<RooPlot*,RooFitResult*>::iterator it;
  for(it  = fit_res_.begin()
     ;it != fit_res_.end()
     ;it++ ){
      
       writeRooPlot((it->first));
       //it->first->Write();
       //it->second->Write();
  }
  
  std::map<std::string,RooDataSet*>::iterator it_d = data_.begin();
  std::map<std::string,RooDataSet*>::iterator it_e = data_.end();

  for(;it_d != it_e; it_d++){
     writeRooDataHist((*it_d).first,it_d->second);
  }
  ws.SetName("cms-hgg-workspace");
  ws.Write();
}

void RooContainer::SetRealVar(std::string var_name, int cat, float x, float w){
 
  if (cat>-1 && cat<ncat){
    std::string name = getcatName(var_name,cat);
    std::map<std::string, RooRealVar>::const_iterator it_var  = m_real_var_.find(name);
    if (it_var == m_real_var_.end()) 
      std::cout << "Warning, No DataSet named "<< name << std::endl;
    else{
      float min_x = m_var_min_[name];
      float max_x = m_var_max_[name];

      if (x > min_x && x < max_x){
        m_real_var_[name] = x;
        data_[name]->add(RooArgSet(m_real_var_[name]),w);
      }
    }
  }
}
// ---------------------------------------------------------------------------------
std::string RooContainer::getcatName(std::string name, int i){
  char* char_string = Form("%s_cat%d",name.c_str(),i);
  std::string output(char_string);
  return output;
}
void RooContainer::addRealVar(std::string name ,float xmin,float xmax){

  RooRealVar temp(name.c_str(),name.c_str(),xmin,xmax);

  m_real_var_.insert(pair<std::string,RooRealVar >(name,temp));
  m_var_min_.insert(pair<std::string, float >(name,xmin));
  m_var_max_.insert(pair<std::string, float >(name,xmax));

  std::cout << "RooContainer::AddRealVar -- Appended the variable " 
	    << name <<std::endl;
 
}

void RooContainer::addRealVar(std::string name ,float init, float vmin, float vmax){
  RooRealVar temp(name.c_str(),name.c_str(),init,vmin,vmax);
  m_real_var_.insert(pair<std::string,RooRealVar>(name,temp));
  std::cout << "RooContainer::AddRealVar -- Appended the variable " 
	    << name <<std::endl;
  
}

void RooContainer::addGenericPdf(std::string name,std::string formula
				,std::vector<std::string> & var
				,double norm_guess ){

    RooArgList roo_args;
    for (std::vector<std::string>::iterator it_var = var.begin()
	;it_var != var.end()
	;it_var++
	){
	  std::cout << "RooContainer::AddGenericPdf -- Adding Parameter " 
		    << *it_var << std::endl;
	  roo_args.add(m_real_var_[*it_var]);
	}

    std::cout << "RooContainer::AddGenericPdf -- Added all variables" 
	      << std::endl;
    RooGenericPdf temp_1(Form("comp_%s",name.c_str()),name.c_str(),formula.c_str(),roo_args);	
    m_gen_.insert(std::pair<std::string , RooGenericPdf>(name,temp_1));

    RooRealVar temp_var(Form("norm_%s",name.c_str()),name.c_str(),norm_guess,0.0,10000);
    m_real_var_.insert(pair<std::string,RooRealVar>(name,temp_var));

    RooExtendPdf  temp(name.c_str(),name.c_str(),m_gen_[name],m_real_var_[name]);

    std::cout << "RooContainer::AddGenericPdf -- Made extended PDF " 
	      << name << std::endl;			       
    m_exp_.insert(pair<std::string,RooExtendPdf>(name,temp));
}

void RooContainer::composePdf(std::string name, std::string  composition
			     ,std::vector<std::string> & formula){

    RooArgList roo_args;
    RooArgList roo_funs;

    for (std::vector<std::string>::iterator it_fun = formula.begin()
	;it_fun != formula.end()
	;it_fun++
	){
	  std::cout << "RooContainer::ComposePdf -- Including Function " 
		    << *it_fun << std::endl;
	  roo_funs.add(m_exp_[*it_fun]);
	  roo_args.add(m_real_var_[(*it_fun)]);
	}

    RooAddPdf temp(name.c_str(),composition.c_str(),roo_funs,roo_args);
    //RooAddPdf temp(name,composition,roo_funs);
    std::cout << "RooContainer::ComposePdf -- Created Composed PDF " 
	      << name << std::endl;			       
    m_pdf_.insert(pair<std::string,RooAddPdf>(name,temp));
}


void RooContainer::createDataSet(std::string name){

    std::map<std::string,RooRealVar>::const_iterator test=m_real_var_.find(name);
    if (test != m_real_var_.end()){ 
      data_[name] = new RooDataSet(name.c_str(),name.c_str(),RooArgSet(m_real_var_[name]) );
      cout << "RooContainer::CreateDataSet -- Created RooDataSet from " << name << endl;
    } 
    else {
      std::cout << "WARNING -- RooContainer::CreateDataSet -- No RealVar found Named "
                << name
		<< " Expect a segmentation fault!!! -- WARNING"
		<< std::endl;
    }	

}

void RooContainer::fitToData(std::string name_func, std::string  name_var, int bins){

    bool use_composed_pdf = false;
    double chi_square;
    std::cout << "RooContainer::FitToData -- Fitting function " 
	      << name_func 
	      << " To data " 
	      << name_var
	      << std::endl; 

    RooFitResult *fit_result; // -> Save this to the Workspace?
    // Look in the composed pdf before checking the standards
    std::map<std::string ,RooAddPdf>::const_iterator it_pdf_ = m_pdf_.find(name_func);
    if (it_pdf_ != m_pdf_.end()){
      fit_result = m_pdf_[name_func].fitTo(*(data_[name_var]));
      use_composed_pdf = true;
    }
    else {
      fit_result = m_exp_[name_func].fitTo(*(data_[name_var]));
    }

    float x_min = m_var_min_[name_var];
    float x_max = m_var_max_[name_var];
 
    RooPlot *xframe = m_real_var_[name_var].frame(x_min,x_max);

    if (bins > 0)  data_[name_var]->plotOn(xframe,Binning(bins));
    else  data_[name_var]->plotOn(xframe);

    if (use_composed_pdf){
      m_pdf_[name_func].plotOn(xframe,LineColor(4));
      int npdfs = m_pdf_[name_func].pdfList().getSize();
    
      for (int i=0;i<npdfs;i++){
	m_pdf_[name_func].plotOn(xframe
				,Components(*(m_pdf_[name_func].pdfList().at(i)))
				,LineColor(i+1)
				,LineStyle(kDashed));
	m_pdf_[name_func].paramOn(xframe);
      } 
    }
    else {
	m_exp_[name_func].plotOn(xframe,LineColor(4));
	m_exp_[name_func].paramOn(xframe);
    }

    xframe->SetName(Form("%s_%s",name_func.c_str(),name_var.c_str()));
    fit_res_.insert(std::pair<RooPlot*,RooFitResult*>(xframe,fit_result));
}

void RooContainer::fitToData(std::string name_func, std::string  name_var
			    ,double x1, double x2, double x3, double x4, int bins){

    
    float x_min = m_var_min_[name_var];
    float x_max = m_var_max_[name_var];

    bool use_composed_pdf = false;
    double chi_square;
    std::cout << "RooContainer::FitToData -- Fitting function " 
	      << name_func 
	      << " To data " 
	      << name_var
	      << std::endl; 

    RooFitResult *fit_result;
    // Look in the composed pdf before checking the standards
    std::map<std::string ,RooAddPdf>::const_iterator it_pdf_ = m_pdf_.find(name_func);
    if (it_pdf_ != m_pdf_.end()){

      if (x1 < x_min || x4 > x_max){
        std::cout << "RooContainer::FitToData -- WARNING!! Ranges outside of DataSet Range!" 
		  << std::endl;
        fit_result = m_pdf_[name_func].fitTo(*(data_[name_var]));
        std::cout << " Fitted To Full Range -- WARNING!!" 
		  << std::endl;
      } else {
        m_real_var_[name_var].setRange("rnge1",x1,x2);
        m_real_var_[name_var].setRange("rnge2",x3,x4);
        fit_result = m_pdf_[name_func].fitTo(*(data_[name_var]),Range("rnge1,rnge2"));
      }
      use_composed_pdf = true;
    }
    else {
     if (x1 < x_min || x4 > x_max){
        std::cout << "RooContianer::FitToData -- WARNING!! Ranges outside of DataSet Range!" 
		  << std::endl;
        fit_result = m_exp_[name_func].fitTo(*(data_[name_var]));
        std::cout << " Fitted To Full Range -- WARNING!!" 
		  << std::endl;
      } else {
        m_real_var_[name_var].setRange("rnge1",x1,x2);
        m_real_var_[name_var].setRange("rnge2",x3,x4);
        fit_result = m_exp_[name_func].fitTo(*(data_[name_var]),Range("rnge1,rnge2"));
      }
    }

    RooPlot *xframe = m_real_var_[name_var].frame(x_min,x_max);

    if (bins > 0)  data_[name_var]->plotOn(xframe,Binning(bins));
    else  data_[name_var]->plotOn(xframe);

    if (use_composed_pdf){
      m_pdf_[name_func].plotOn(xframe,LineColor(4));
      int npdfs = m_pdf_[name_func].pdfList().getSize();
    
      for (int i=0;i<npdfs;i++){
	m_pdf_[name_func].plotOn(xframe
			,Components(*(m_pdf_[name_func].pdfList().at(i)))
			,LineColor(i+1)
			,LineStyle(kDashed));
	m_pdf_[name_func].paramOn(xframe);
      }
    }

    else {
	m_exp_[name_func].plotOn(xframe,LineColor(4));
	m_exp_[name_func].paramOn(xframe);
    }

    xframe->SetName(Form("%s_%s",name_func.c_str(),name_var.c_str()));

    fit_res_.insert(std::pair<RooPlot*,RooFitResult*>(xframe,fit_result));
}


void RooContainer::writeRooDataHist(std::string name, RooDataSet* data){
  RooDataHist tmp(Form("hist_%s",name.c_str()),name.c_str(),RooArgList(m_real_var_[name]),*data);
  ws.import(tmp);
  //tmp.Write();
}


void RooContainer::writeRooPlot(RooPlot *plot){
  // Dont want to let the plots pop up!
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  // ---------------------------------
  TCanvas *can = new TCanvas(Form("can_%s",plot->GetName())
			 ,plot->GetName(),900,600) ;    
  can->cd(); plot->Draw();
  can->Write();
  delete can;
}

