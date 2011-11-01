#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include <string>
#include <iostream>

// Macro to produce templates reflecting the statistical error bin to bin from the background model

void bkgUncorrelatedErrors(std::string FileName){
      int nSys = 10;
      TFile *file = (TFile*) TFile::Open(FileName.c_str(),"UPDATE");
      std::string types[2]={"ada","grad"};
      for (int i = 0;i<2;i++){
	for (float mass=115.0;mass<=150.0;mass+=0.5){
		
		TH1F* hbkg=(TH1F*) file->Get(Form("th1f_bkg_%s_%3.1f_cat0",types[i].c_str(),mass));
		int nBins=hbkg->GetNbinsX();	
		int sys_n=0;
		for (int j=nBins;(j>0&&sys_n<nSys);j--){
			
			double Integral = hbkg->Integral();
			double BinVal = hbkg->GetBinContent(j);
			double BinErr = hbkg->GetBinError(j);
			TH1F *hup = (TH1F*)hbkg->Clone();
			TH1F *hdn = (TH1F*)hbkg->Clone();

			for (int bin=1;bin<=nBins;bin++){
				if (bin==j) {
					std::cout << "FFS!!!!!!!!!! " << bin <<std::endl;
					hup->SetBinContent(bin,BinVal+BinErr);
					hdn->SetBinContent(bin,BinVal-BinErr);					
				}
			}
			
			hup->Scale(Integral/hup->Integral());
			hdn->Scale(Integral/hdn->Integral());
			hup->SetName(Form("th1f_bkg_%s_%3.1f_cat0_bkgUCorr%dUp01_sigma",types[i].c_str(),mass,sys_n));
			hdn->SetName(Form("th1f_bkg_%s_%3.1f_cat0_bkgUCorr%dDown01_sigma",types[i].c_str(),mass,sys_n));
			file->cd();
			hup->Write();
			hdn->Write();
			sys_n++;
		}
	}

      }

}
