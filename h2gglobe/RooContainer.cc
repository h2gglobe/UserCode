/****************************************************** 
   RooContainer.cc  
   Original Author - Nicholas Wardle - Imperial College
*******************************************************/
   
#include "RooContainer.h"

using namespace RooFit;

RooContainer::RooContainer(int n):ncat(n),nsigmas(30){}

// ----------------------------------------------------------------------------------------------------
void RooContainer::AddRealVar(std::string name,float xmin,float xmax){
  for (int cat=0;cat<ncat;cat++){
    addRealVar(getcatName(name,cat),xmin,xmax);
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::AddRealVar(std::string name,float init, float vmin, float vmax){
  for (int cat=0;cat<ncat;cat++){
    addRealVar(getcatName(name,cat),init,vmin,vmax);
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::AddGenericPdf(std::string name,std::string formula
				,std::vector<std::string> & var, int form
				,double norm_guess ){
  for (int cat=0;cat<ncat;cat++){
    std::vector<std::string> cat_var;
    for (std::vector<std::string>::iterator it=var.begin()
	;it!=var.end()
	;++it){
      cat_var.push_back(getcatName(*it,cat));
    }  
    addGenericPdf(getcatName(name,cat),formula,cat_var,form,norm_guess);
  }
}

// ----------------------------------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------------------------------
void RooContainer::CreateDataSet(std::string name,int nbins){
  for (int cat=0;cat<ncat;cat++){
    std::string cat_name = getcatName(name,cat);
    createDataSet(cat_name,cat_name,nbins);
    
    float xmin = m_var_min_[cat_name];
    float xmax = m_var_max_[cat_name];

    if (nbins >0) {
        m_real_var_[cat_name].setBins(nbins);
      }
      else {
        m_real_var_[cat_name].setBins((int)(xmax-xmin));
      }
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::MakeSystematics(std::string s_name, std::string sys_name){
  for (int cat=0;cat<ncat;cat++){
    makeSystematics(getcatName(s_name,cat),sys_name);
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::FitToData(std::string name_func, std::string  name_var
			    ,double x1, double x2, double x3, double x4){
  for (int cat=0;cat<ncat;cat++){
    fitToData(getcatName(name_func,cat),getcatName(name_var,cat),getcatName(name_var,cat)
	     ,x1,x2,x3,x4);
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::FitToSystematicSet(std::string name_func,std::string name_var
	     ,std::string sys_name
	     ,double x1,double x2,double x3,double x4){
  // Assuming here that dont care about previous fits to the other data-sets
  for (int cat=0;cat<ncat;cat++) {
    fitToSystematicSet(getcatName(name_func,cat),getcatName(name_var,cat)
		      ,sys_name,x1,x2,x3,x4);
  }    
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::Save(){

  std::cout << "RooContainer::Save -- Saving To File "
            << std::endl;

  std::map<RooPlot*,RooFitResult*>::iterator it;
  for(it  = fit_res_.begin()
     ;it != fit_res_.end()
     ;it++ ){
      
       writeRooPlot((it->first));
  }
  
  std::map<std::string,TH1F>::iterator it_d = m_th1f_.begin();
  std::map<std::string,TH1F>::iterator it_e = m_th1f_.end();

  for(;it_d != it_e; it_d++){
     writeRooDataHist((*it_d).first,&(it_d->second));
  }
  ws.SetName("cms-hgg-workspace");
  ws.Write();
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::InputDataPoint(std::string var_name, int cat, float x, float w){
 
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
        data_[name].add(RooArgSet(m_real_var_[name]),w);
        m_th1f_[name].Fill(x,w);
      }
    }
  }

  else {
    std::cout << "No Category Number " << cat 
              << ", category must be from 0 to " << ncat-1
	      << std::endl;
  }
}


// ----------------------------------------------------------------------------------------------------
void RooContainer::InputSystematicSet(std::string s_name, std::string sys_name,int cat
		  ,std::vector<float> x, float w){

  if (x.size() != 2*nsigmas){
    std::cout << "WARNING -- RooContainer::InputSystematicSet -- Size of vector must be equal to "
	      << 2*nsigmas << " -- WARNING" << std::endl;	
  }
 
 else{

   if (cat>-1 && cat<ncat){

    std::string cat_name = getcatName(s_name,cat);
    std::string name = getsysName(cat_name,sys_name);

    std::map<std::string,RooRealVar>::iterator it_var  = m_real_var_.find(cat_name);

    if (it_var == m_real_var_.end()) 
      std::cout << "WARNING -- RooContainer::InpusSystematicSet -- No DataSet named "<< cat_name << std::endl;
  
    else {

      std::vector<RooDataSet*>::iterator data_set_up = data_up_[name].begin();
      std::vector<RooDataSet*>::iterator data_set_dn = data_dn_[name].begin();
      std::vector<TH1F*>::iterator th1f_set_up = m_th1f_up_[name].begin();
      std::vector<TH1F*>::iterator th1f_set_dn = m_th1f_dn_[name].begin();

      std::vector<float>::iterator val = x.begin();

      float min_x = m_var_min_[cat_name];
      float max_x = m_var_max_[cat_name];

      // Loop over the first nsigmas elements as the -1 -> dn sys
      for (;data_set_dn != data_dn_[name].end()
	   ;data_set_dn++,th1f_set_dn++,val++){

        if (*val > min_x && *val < max_x){
           it_var->second = *val;
           (*data_set_dn)->add(RooArgSet(it_var->second),w);
           (*th1f_set_dn)->Fill(*val,w);
        }
      }

      // Loop over the second nsigmas elements as the +1 -> up sys
      for (;data_set_up != data_up_[name].end()
	   ;data_set_up++,th1f_set_up++,val++){

        if (*val > min_x && *val < max_x){
           it_var->second = *val;
           (*data_set_up)->add(RooArgSet(it_var->second),w);
           (*th1f_set_up)->Fill(*val,w);
        }
      }
    }
   }
   else {
    std::cout << "No Category Number " << cat 
              << ", category must be from 0 to " << ncat-1
	      << std::endl;
   }
  }
}

// ----------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
std::string RooContainer::getcatName(std::string name, int i){
  char* char_string = Form("%s_cat%d",name.c_str(),i);
  std::string output(char_string);
  return output;
}

// ----------------------------------------------------------------------------------------------------
std::string RooContainer::getsysindexName(std::string name,std::string sys_name
				    ,int sys,int direction){
  char* char_string;
  if (direction >= 0)
    char_string = Form("%s_%s-up%.2d-sigma",name.c_str(),sys_name.c_str(),sys);
  else
    char_string = Form("%s_%s-down%.2d-sigma",name.c_str(),sys_name.c_str(),sys);
  std::string output(char_string);
  return output;
}

// ----------------------------------------------------------------------------------------------------
std::string RooContainer::getsysName(std::string name,std::string sys_name){
  char* char_string = Form("%s_%s",name.c_str(),sys_name.c_str());
  std::string output(char_string);
  return output;
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::addRealVar(std::string name ,float xmin,float xmax){

  RooRealVar temp(name.c_str(),name.c_str(),xmin,xmax);

  m_real_var_.insert(pair<std::string,RooRealVar >(name,temp));
  m_var_min_.insert(pair<std::string, float >(name,xmin));
  m_var_max_.insert(pair<std::string, float >(name,xmax));

  std::cout << "RooContainer::AddRealVar -- Appended the variable " 
	    << name <<std::endl; 
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::addRealVar(std::string name ,float init, float vmin, float vmax){
  RooRealVar temp(name.c_str(),name.c_str(),init,vmin,vmax);
  m_real_var_.insert(pair<std::string,RooRealVar>(name,temp));
  std::cout << "RooContainer::AddRealVar -- Appended the variable " 
	    << name <<std::endl;
  
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::addGenericPdf(std::string name,std::string formula
				,std::vector<std::string> & var, int form
				,double norm_guess ){

    RooAbsPdf *temp_1;

    if (form==0){
       RooArgList roo_args;
      for (std::vector<std::string>::iterator it_var = var.begin()
	   ;it_var != var.end()
	   ;it_var++
	   ){
	     std::cout << "RooContainer::AddGenericPdf -- Adding Parameter " 
		       << *it_var << std::endl;

	     std::map<std::string,RooRealVar>::iterator real_var =  m_real_var_.find(*it_var);
	     if (real_var != m_real_var_.end())
	       roo_args.add((*real_var).second);
	     else 
      		std::cout << "WARNING -- RooContainer::AddGenericPdf -- No Variable Found Named " 
	        	  << *it_var << std::endl;
  	   }

      std::cout << "RooContainer::AddGenericPdf -- Added all variables" 
	        << std::endl;

      temp_1 = new RooGenericPdf(Form("comp_%s",name.c_str()),name.c_str(),formula.c_str(),roo_args);	
      v_gen_.push_back(temp_1);

    } else if (form == 1) { //RooExponential  - x,slope
	if (var.size() == 2){
	  temp_1 = new RooExponential(Form("comp_%s",name.c_str()),name.c_str(),m_real_var_[var[0]],m_real_var_[var[1]]);
          v_gen_.push_back(temp_1);
	} else {
		
          std::cout << "WARNING -- RooContainer::AddGenericPdf -- Need 2 arguments for RooExponential, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
    } else if (form == 2) { //RooGaussian - x,mean,sigma
	if (var.size() == 3){
	  temp_1 = new RooGaussian(Form("comp_%s",name.c_str()),name.c_str(),m_real_var_[var[0]],m_real_var_[var[1]],m_real_var_[var[2]]);
          v_gen_.push_back(temp_1);
	} else {
		
          std::cout << "WARNING -- RooContainer::AddGenericPdf -- Need 3 arguments for RooGaussian, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
    } else if (form == 3) { //RooBreitWigner - x,centre,width
	if (var.size() == 3){
	  temp_1 = new RooBreitWigner(Form("comp_%s",name.c_str()),name.c_str(),m_real_var_[var[0]],m_real_var_[var[1]],m_real_var_[var[2]]);
          v_gen_.push_back(temp_1);
	} else {
		
          std::cout << "WARNING -- RooContainer::AddGenericPdf -- Need 3 arguments for RooBreitWigner, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
    } else if (form == 4) { //RooCBShape - x,centre,sigma,slope,n
	if (var.size() == 5){
	  temp_1 = new RooCBShape(Form("comp_%s",name.c_str()),name.c_str(),m_real_var_[var[0]],m_real_var_[var[1]],m_real_var_[var[2]]
				 ,m_real_var_[var[3]],m_real_var_[var[4]]);
          v_gen_.push_back(temp_1);
	} else {
		
          std::cout << "WARNING -- RooContainer::AddGenericPdf -- Need 5 arguments for RooCBShape, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
    } else if (form == 5) { //RooVoigtian - x,centre,width,sigma
	if (var.size() == 4){
	  temp_1 = new RooVoigtian(Form("comp_%s",name.c_str()),name.c_str(),m_real_var_[var[0]],m_real_var_[var[1]],m_real_var_[var[2]]
				  ,m_real_var_[var[3]]);
          v_gen_.push_back(temp_1);
	} else {
		
          std::cout << "WARNING -- RooContainer::AddGenericPdf -- Need 4 arguments for RooVoigtian, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
    }

    else {
	
       std::cout << "WARNING -- RooContainer::AddGenericPdf -- No Mode << form "
		 << "Understood -- WARNING"
	         << std::endl;
    }
		       
    RooRealVar temp_var(Form("norm_%s",name.c_str()),name.c_str(),norm_guess,0.0,10000);
    m_real_var_.insert(pair<std::string,RooRealVar>(name,temp_var));

    RooExtendPdf temp(name.c_str(),name.c_str(),*temp_1,m_real_var_[name]);
    m_exp_.insert(pair<std::string,RooExtendPdf>(name,temp));

    std::cout << "RooContainer::AddGenericPdf -- Made extended PDF " 
	        << name << std::endl;	
}

// ----------------------------------------------------------------------------------------------------
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

	  std::map<std::string,RooExtendPdf>::iterator exp = m_exp_.find(*it_fun);

	  if (exp != m_exp_.end()){
	    roo_funs.add((*exp).second);
	    roo_args.add(m_real_var_[(*it_fun)]);

	  } else {
      	    std::cout << "WARNING -- RooContainer::ComposePdf -- No Pdf Found Named " 
	              << *it_fun << std::endl;
	  }
	}

    RooAddPdf temp(name.c_str(),composition.c_str(),roo_funs,roo_args);

    std::cout << "RooContainer::ComposePdf -- Created Composed PDF " 
	      << name << std::endl;			       

    m_pdf_.insert(pair<std::string,RooAddPdf>(name,temp));
}


// ----------------------------------------------------------------------------------------------------
void RooContainer::createDataSet(std::string name,std::string data_name,int nbins){

    std::map<std::string,RooRealVar>::const_iterator test=m_real_var_.find(name);
    if (test != m_real_var_.end()){ 

      RooDataSet data_tmp(data_name.c_str(),data_name.c_str(),RooArgSet(m_real_var_[name]) );
      data_.insert(std::pair<std::string,RooDataSet>(data_name,data_tmp));
      bins_[name] = nbins;

      float xmin = m_var_min_[name];
      float xmax = m_var_max_[name];

      if (nbins >0) {
	m_th1f_[data_name] = TH1F(Form("th1f-%s",data_name.c_str()),Form("th1f-%s",data_name.c_str()),nbins,xmin,xmax);
      }
      else {
	m_th1f_[data_name] = TH1F(Form("th1f-%s",data_name.c_str()),Form("th1f-%s",data_name.c_str()),(int)(xmax-xmin),xmin,xmax);
      }

      cout << "RooContainer::CreateDataSet -- Created RooDataSet from " << name 
	   << " with name " << data_name <<endl;
      cout << "RooContainer::CreateDataSet -- Created TH1F from " << name << endl;
    } 

    else {
      std::cout << "WARNING -- RooContainer::CreateDataSet -- No RealVar found Named "
                << name
		<< " CRASH Expected soon!!! -- WARNING"
		<< std::endl;
    }	
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::fitToData(std::string name_func, std::string name_data, std::string  name_var
			    ,double x1, double x2, double x3, double x4){

    bool use_composed_pdf = false;
    double chi_square;
    std::cout << "RooContainer::FitToData -- Fitting function " 
	      << name_func 
	      << " To data " 
	      << name_var
	      << std::endl; 

    std::map<std::string,RooDataSet>::iterator it_data = data_.find(name_data);
    if (it_data == data_.end()){
      std::cout << "WARNING -- RooContainer::FitToData -- No DataSet Found Named "
 	 	<< name_var << std::endl;
      return;
    }

    RooFitResult *fit_result;
    // Look in the composed pdf before checking the standards
    std::map<std::string ,RooAddPdf>::iterator it_pdf = m_pdf_.find(name_func);
    std::map<std::string ,RooExtendPdf>::iterator it_exp = m_exp_.find(name_func);

    int bins = bins_[name_var];
    float x_min = m_var_min_[name_var];
    float x_max = m_var_max_[name_var];

    RooRealVar *real_var = &m_real_var_[name_var];

    if (it_pdf != m_pdf_.end()){
      use_composed_pdf = true;

     if (x1 < -990. && x2 < -990.  && x3 < -990. && x4 < -990.){
    std::cout<<"DEBUG1"<<std::endl;
        fit_result = (it_pdf->second).fitTo(it_data->second);
    std::cout<<"DEBUG2"<<std::endl;
     } else {
      if (x1 < x_min || x4 > x_max){
        std::cout << "RooContainer::FitToData -- WARNING!! Ranges outside of DataSet Range!" 
		  << std::endl;
        fit_result = (it_pdf->second).fitTo(it_data->second);
        std::cout << " Fitted To Full Range -- WARNING!!" 
		  << std::endl;
      } else {
        real_var->setRange("rnge1",x1,x2);
        real_var->setRange("rnge2",x3,x4);
    std::cout<<"DEBUG1"<<std::endl;
        fit_result = (it_pdf->second).fitTo((it_data->second),Range("rnge1,rnge2"));
    std::cout<<"DEBUG2"<<std::endl;
      }
     }
    }

    else {
      if (it_exp != m_exp_.end()){
       if (x1 < -990. && x2 < -990.  && x3 < -990. && x4 < -990.){
    std::cout<<"DEBUG1 -- YES HERE"<<std::endl;
        fit_result = (it_exp->second).fitTo(it_data->second);
    std::cout<<"DEBUG2"<<std::endl;
       } else {
        if (x1 < x_min || x4 > x_max){
          std::cout << "RooContianer::FitToData -- WARNING!! Ranges outside of DataSet Range!" 
		    << std::endl;
          fit_result = (it_exp->second).fitTo((it_data->second));
          std::cout << " Fitted To Full Range -- WARNING!!" 
		    << std::endl;
        } else {
          real_var->setRange("rnge1",x1,x2);
          real_var->setRange("rnge2",x3,x4);
          fit_result = (it_exp->second).fitTo((it_data->second),Range("rnge1,rnge2"));
        }
       }
      }
      else
	std::cout << "WARNING -- RooContainer::FitToData -- No Pdf Found Named "
	 	  << name_func << std::endl;
    }

    RooPlot *xframe = (*real_var).frame(x_min,x_max);

    if (bins > 0) (it_data->second).plotOn(xframe,Binning(bins));
    else  (it_data->second).plotOn(xframe);

    if (use_composed_pdf){
      (it_pdf->second).plotOn(xframe,LineColor(4));
      int npdfs = (it_pdf->second).pdfList().getSize();
    
      for (int i=0;i<npdfs;i++){
	(it_pdf->second).plotOn(xframe
			,Components(*((it_pdf->second).pdfList().at(i)))
			,LineColor(i+1)
			,LineStyle(kDashed));
	(it_pdf->second).paramOn(xframe);
      }
    }

    else {
	(it_exp->second).plotOn(xframe,LineColor(4));
	(it_exp->second).paramOn(xframe);
    }

    xframe->SetName(Form("%s_%s",name_func.c_str(),name_data.c_str()));

    fit_res_.insert(std::pair<RooPlot*,RooFitResult*>(xframe,fit_result));
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::fitToSystematicSet(std::string name_func,std::string name_var
	     ,std::string sys_name
	     ,double x1,double x2,double x3,double x4){
  for (int sys=1;sys<=nsigmas;sys++){
    fitToData(name_func,getsysindexName(name_var,sys_name,sys,-1),name_var,x1,x2,x3,x4);
    fitToData(name_func,getsysindexName(name_var,sys_name,sys,1),name_var,x1,x2,x3,x4);
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::writeRooDataHist(std::string name, TH1F* hist){

  RooRealVar tmp_var(name.c_str(),name.c_str(),hist->GetBinLowEdge(1)
		,hist->GetBinLowEdge(hist->GetNbinsX()+1));
  tmp_var.setBins(hist->GetNbinsX());

  RooDataHist tmp(Form("roohist-%s",name.c_str()),name.c_str(),RooArgList(tmp_var),hist);

  ws.import(tmp);
  hist->Write();
}


// ----------------------------------------------------------------------------------------------------
void RooContainer::writeRooPlot(RooPlot *plot){
  // Dont want to let the plots pop up!
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  // ---------------------------------
  TCanvas *can = new TCanvas(Form("plot-%s",plot->GetName()),plot->GetName(),900,600) ;    
  can->cd(); plot->Draw();
  can->Write();
  delete can;
}


// ----------------------------------------------------------------------------------------------------
void RooContainer::makeSystematics(std::string s_name, std::string sys_name){

  std::map<std::string,RooDataSet>::iterator test_it = data_.find(s_name);

  if (test_it != data_.end()){

    std::vector<RooDataSet*> v_sys_up;
    std::vector<RooDataSet*> v_sys_dn;

    std::vector<TH1F*> v_th1f_up;
    std::vector<TH1F*> v_th1f_dn;

    float min  = m_var_min_[s_name];
    float max  = m_var_max_[s_name];
    int bins   = bins_[s_name];
  
    for (int sys=1;sys<=nsigmas;sys++){

      std::string name_up = getsysindexName(s_name,sys_name,sys,1);
      std::string name_dn = getsysindexName(s_name,sys_name,sys,-1);
    
      createDataSet(s_name,name_up,bins);
      createDataSet(s_name,name_dn,bins);

      v_sys_up.push_back(&(data_[name_up]));
      v_sys_dn.push_back(&(data_[name_dn]));

      v_th1f_up.push_back(&(m_th1f_[name_up]));
      v_th1f_dn.push_back(&(m_th1f_[name_dn]));
    }

  std::string map_name = getsysName(s_name,sys_name);

  data_up_.insert(std::pair<std::string,std::vector<RooDataSet*> >(map_name,v_sys_up));
  data_dn_.insert(std::pair<std::string,std::vector<RooDataSet*> >(map_name,v_sys_dn));

  m_th1f_up_.insert(std::pair<std::string,std::vector<TH1F*> >(map_name,v_th1f_up));
  m_th1f_dn_.insert(std::pair<std::string,std::vector<TH1F*> >(map_name,v_th1f_dn));
  
  }
  else
    std::cout << "RooContainer::MakeSystematics -- Cannot Create Systematics Study, "
	      << "No DATASET Found Named " << s_name << std::endl;
}
// ----------------------------------------------------------------------------------------------------
//EOF

