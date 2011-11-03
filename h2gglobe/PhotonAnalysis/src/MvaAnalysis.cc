#include "../interface/MvaAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
MvaAnalysis::MvaAnalysis()  : 
    name_("MvaAnalysis"),
    vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams)
{

    systRange  = 3.; // in units of sigma
    nSystSteps = 1;    
    nMasses  = 9;
}

// ----------------------------------------------------------------------------------------------------
MvaAnalysis::~MvaAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::Term(LoopAll& l) 
{
    if (doTraining){
        for (int i = 0; i<nMasses;i++){
            mvaFile_->cd();
            if (0<signalTrainTree_[i]->GetEntries()    ){
                signalTrainTree_[i]->Write(("sig_train"+names[i]).c_str());
            }
            if (0<signalTestTree_[i]->GetEntries()    ){
                signalTestTree_[i]->Write(("sig_test"+names[i]).c_str());
            }
            if (0<backgroundTrainTree_[i]->GetEntries()    ){
                backgroundTrainTree_[i]->Write(("bkg_train"+names[i]).c_str());
            }
            if (0<backgroundTestTree_[i]->GetEntries()    ){
                backgroundTestTree_[i]->Write(("bkg_test"+names[i]).c_str());
            }
        }
        if (0<backgroundTrainTree_all->GetEntries() ){
            backgroundTrainTree_all->Write("bkg_train_all");
        }
        if (0<backgroundTestTree_all->GetEntries() ){
            backgroundTestTree_all->Write("bkg_test_all");
        }
        mvaFile_->Close();
    }
    else{
        // loop hypothesis
        for (double mass=115.0; mass<150.2; mass+=0.5){
	    if (mass < masses[2] || mass > masses[nMasses-1] ) continue;	

            // define hypothesis masses for the sidebands
            float mass_hypothesis = mass;

            // define the signal Region
            float sideband_boundaries[2];

            sideband_boundaries[0] = mass_hypothesis*(1.-sidebandWidth);
            sideband_boundaries[1] = mass_hypothesis*(1.+sidebandWidth);

            l.rooContainer->FitToData("data_pow_model", "data_mass",massMin,sideband_boundaries[0],sideband_boundaries[1],massMax);
            std::vector<std::pair<double,double> > N_sigErr = l.rooContainer->GetFitNormalisationsAndErrors("data_pow_model","data_mass",sideband_boundaries[0],sideband_boundaries[1]);
            
            l.rooContainer->AddNormalisationSystematics(Form("bkg_norm_%3.1f",mass),N_sigErr, 1); // 1 means it effect the background only

            // Calculate weights to apply to the sidebands
            std::vector<double> N_sig = l.rooContainer->GetFitNormalisations("data_pow_model","data_mass",sideband_boundaries[0],sideband_boundaries[1]);

            bool scale = true;//true;
	    std::vector<string> ada_bkgsets;
	    std::vector<string> grad_bkgsets;
	    std::vector<string> ada_datasets;
	    std::vector<string> grad_datasets;

	    for (int sideband_i=1;sideband_i<=numberOfSidebands;sideband_i++) {
		ada_bkgsets.push_back(Form("bkg_%dlow_BDT_ada_%3.1f",sideband_i,mass));
		ada_bkgsets.push_back(Form("bkg_%dhigh_BDT_ada_%3.1f",sideband_i,mass));
		ada_datasets.push_back(Form("data_%dlow_BDT_ada_%3.1f",sideband_i,mass));
		ada_datasets.push_back(Form("data_%dhigh_BDT_ada_%3.1f",sideband_i,mass));

		grad_bkgsets.push_back(Form("bkg_%dlow_BDT_grad_%3.1f",sideband_i,mass));
		grad_bkgsets.push_back(Form("bkg_%dhigh_BDT_grad_%3.1f",sideband_i,mass));
		grad_datasets.push_back(Form("data_%dlow_BDT_grad_%3.1f",sideband_i,mass));
		grad_datasets.push_back(Form("data_%dhigh_BDT_grad_%3.1f",sideband_i,mass));
		
	    }
	    // Also for the bkg MC we want to throw in the signal region for the binning algos
            ada_bkgsets.push_back(Form("bkg_BDT_ada_%3.1f",mass)) ;
            grad_bkgsets.push_back(Form("bkg_BDT_grad_%3.1f",mass)) ;

	    
            l.rooContainer->SumMultiBinnedDatasets(Form("bkg_BDT_ada_all_%3.1f",mass),ada_bkgsets, N_sig, scale);
            l.rooContainer->SumMultiBinnedDatasets(Form("bkg_BDT_grad_all_%3.1f",mass),grad_bkgsets, N_sig, scale);
		
            l.rooContainer->SumMultiBinnedDatasets(Form("data_BDT_sideband_ada_%3.1f",mass),ada_datasets, N_sig, scale);
            l.rooContainer->SumMultiBinnedDatasets(Form("data_BDT_sideband_grad_%3.1f",mass),grad_datasets, N_sig, scale);

      }
      for (int i=2; i<nMasses; i++){
        // Need to binning on MC!
        std::vector <std::vector<double> > optimizedGradBins =  l.rooContainer->OptimizedBinning("bkg_BDT_grad_all"+names[i],50,false,true,-1);
        std::vector<std::vector <double> > optimizedAdaBins =  l.rooContainer->OptimizedBinning("bkg_BDT_ada_all"+names[i],50,false,true,-1);
        
	      double mass_h_low;      
              double mass_h_high;
 
	      if (i==7){ // the 140 case
	       mass_h_low  =-2.5;
	       mass_h_high =4.6;
	      } else if (i==8){  // the 150 case
	       mass_h_low  =-5.0;
	       mass_h_high =0.1;
	      } else {
	       mass_h_low  =-2.5;
	       mass_h_high =2.1;
	      }

        for (double mass=mass_h_low; mass<mass_h_high; mass+=0.5){

          double mass_hypothesis = masses[i]+mass;
	  if (mass_hypothesis < 115 || mass_hypothesis > 150) continue;
      
          l.rooContainer->RebinBinnedDataset(Form("bkg_grad_%3.1f",mass_hypothesis),Form("data_BDT_sideband_grad_%3.1f",mass_hypothesis),optimizedGradBins,false);
          l.rooContainer->RebinBinnedDataset(Form("bkg_mc_grad_%3.1f",mass_hypothesis),Form("bkg_BDT_grad_%3.1f",mass_hypothesis),optimizedGradBins,false);
          l.rooContainer->RebinBinnedDataset(Form("bkg_mc_balg_grad_%3.1f",mass_hypothesis),Form("bkg_BDT_grad_all_%3.1f",mass_hypothesis),optimizedGradBins,false); 
          l.rooContainer->RebinBinnedDataset(Form("data_grad_%3.1f",mass_hypothesis),Form("data_BDT_grad_%3.1f",mass_hypothesis),optimizedGradBins,false);

          l.rooContainer->RebinBinnedDataset(Form("bkg_ada_%3.1f",mass_hypothesis),Form("data_BDT_sideband_ada_%3.1f",mass_hypothesis),optimizedAdaBins,false);
          l.rooContainer->RebinBinnedDataset(Form("bkg_mc_ada_%3.1f",mass_hypothesis),Form("bkg_BDT_ada_%3.1f",mass_hypothesis),optimizedAdaBins,false);
          l.rooContainer->RebinBinnedDataset(Form("bkg_mc_balg_ada_%3.1f",mass_hypothesis),Form("bkg_BDT_ada_all_%3.1f",mass_hypothesis),optimizedAdaBins,false);
          l.rooContainer->RebinBinnedDataset(Form("data_ada_%3.1f",mass_hypothesis),Form("data_BDT_ada_%3.1f",mass_hypothesis),optimizedAdaBins,false);

	  for (int sideband_i=1;sideband_i<=numberOfSidebands;sideband_i++){
           l.rooContainer->RebinBinnedDataset(Form("bkg_%dhigh_grad_%3.1f",sideband_i,mass_hypothesis),Form("data_%dhigh_BDT_grad_%3.1f",sideband_i,mass_hypothesis),optimizedGradBins,false);
           l.rooContainer->RebinBinnedDataset(Form("bkg_%dlow_grad_%3.1f",sideband_i,mass_hypothesis),Form("data_%dlow_BDT_grad_%3.1f",sideband_i,mass_hypothesis),optimizedGradBins,false);
           l.rooContainer->RebinBinnedDataset(Form("bkg_mc_%dhigh_grad_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dhigh_BDT_grad_%3.1f",sideband_i,mass_hypothesis),optimizedGradBins,false);
           l.rooContainer->RebinBinnedDataset(Form("bkg_mc_%dlow_grad_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dlow_BDT_grad_%3.1f",sideband_i,mass_hypothesis),optimizedGradBins,false);

           l.rooContainer->RebinBinnedDataset(Form("bkg_%dhigh_ada_%3.1f",sideband_i,mass_hypothesis),Form("data_%dhigh_BDT_ada_%3.1f",sideband_i,mass_hypothesis),optimizedAdaBins,false);
           l.rooContainer->RebinBinnedDataset(Form("bkg_%dlow_ada_%3.1f",sideband_i,mass_hypothesis),Form("data_%dlow_BDT_ada_%3.1f",sideband_i,mass_hypothesis),optimizedAdaBins,false);
           l.rooContainer->RebinBinnedDataset(Form("bkg_mc_%dhigh_ada_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dhigh_BDT_ada_%3.1f",sideband_i,mass_hypothesis),optimizedAdaBins,false);
           l.rooContainer->RebinBinnedDataset(Form("bkg_mc_%dlow_ada_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dlow_BDT_ada_%3.1f",sideband_i,mass_hypothesis),optimizedAdaBins,false);
	  }
          
        }

        for (int j=-1; j<2; j++){
        if ((i==2 && j==-1) || (i==nMasses-1 && j==1)) continue;
          l.rooContainer->RebinBinnedDataset("sig_grad"+names[i]+names[i+j],"sig_BDT_grad"+names[i]+names[i+j],optimizedGradBins,true);
          l.rooContainer->RebinBinnedDataset("sig_ada"+names[i]+names[i+j],"sig_BDT_ada"+names[i]+names[i+j],optimizedAdaBins,true);
        } 
      }
    }
    PhotonAnalysis::Term(l);
}

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::Init(LoopAll& l) 
{
    if(PADEBUG) 
	cout << "InitRealMvaAnalysis START"<<endl;

    nevents=0., sumwei=0.; 
    sumaccept=0., sumsmear=0., sumev=0.;
    
    names[0]="_105.0";
    BDTnames[0]="_105";
    masses[0] = 105.;

    names[1]="_110.0";
    BDTnames[1]="_110";
    masses[1] = 110.;

    names[2]="_115.0";
    BDTnames[2]="_115";
    masses[2] = 115.;

    names[3]="_120.0";
    BDTnames[3]="_120";
    masses[3] = 120.;

    names[4]="_125.0";
    BDTnames[4]="_125";
    masses[4] = 125.;

    names[5]="_130.0";
    BDTnames[5]="_130";
    masses[5] = 130.;

    names[6]="_135.0";
    BDTnames[6]="_135";
    masses[6] = 135.;

    names[7]="_140.0";
    BDTnames[7]="_140";
    masses[7] = 140.;

    names[8]="_150.0";
    BDTnames[8]="_150";
    masses[8] = 150.;
/*
    names[9]="_121.0";
    BDTnames[9]="_121";
    masses[9] = 121.;

    names[10]="_123.0";
    BDTnames[10]="_123";
    masses[10] = 123.;
*/
    for (float mass = 90.;
               mass<150.;
               mass*= (1+signalRegionWidth)/(1-signalRegionWidth)){
        bkg_masses.push_back(mass); 
        cout<<"background mass hypothesis: "<< mass << endl;
    }
    std::string outputfilename = (std::string) l.histFileName;
    //
    // These parameters are set in the configuration file
    std::cout
	<< "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< "MvaAnalysis " << "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< "leadEtCut "<< leadEtCut << "\n"
	<< "subleadEtCut "<< subleadEtCut << "\n"
	<< "doTriggerSelection "<< doTriggerSelection << "\n"
	<< "nEtaCategories "<< nEtaCategories << "\n"
	<< "nR9Categories "<< nR9Categories << "\n"		
	<< "nPtCategories "<< nPtCategories << "\n"		
	<< "doEscaleSyst "<< doEscaleSyst << "\n"
	<< "doEresolSyst "<< doEresolSyst << "\n"
	<< "doEcorrectionSyst "<< doEcorrectionSyst << "\n"
	<< "efficiencyFile " << efficiencyFile << "\n"
	<< "doPhotonIdEffSyst "<< doPhotonIdEffSyst << "\n"
	<< "doR9Syst "<< doR9Syst << "\n"
	<< "doVtxEffSyst "<< doVtxEffSyst << "\n"
	<< "doTriggerEffSyst "<< doTriggerEffSyst << "\n"
	<< "doKFactorSyst "<< doKFactorSyst << "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< std::endl;

    // call the base class initializer
    PhotonAnalysis::Init(l);

    // Avoid reweighing from histo conainer
    for(size_t ind=0; ind<l.histoContainer.size(); ind++) {
	l.histoContainer[ind].setScale(1.);
    }
    
    diPhoCounter_ = l.countersred.size();
    l.countersred.resize(diPhoCounter_+1);

    // initialize the analysis variables
    nCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nCategories_ *= nR9Categories;
    if( nPtCategories != 0 ) nCategories_ *= nPtCategories;

    nPhotonCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nPhotonCategories_ *= nR9Categories;
    

    effSmearPars.categoryType = "2CatR9_EBEE";
    effSmearPars.n_categories = 4;
    effSmearPars.efficiency_file = efficiencyFile;

    diPhoEffSmearPars.n_categories = 8;
    diPhoEffSmearPars.efficiency_file = efficiencyFile;

    if( doEcorrectionSmear ) {
        // instance of this smearer done in PhotonAnalysis
        photonSmearers_.push_back(eCorrSmearer);
    }
    if( doEscaleSmear ) {
        photonSmearers_.push_back(eScaleSmearer);
    }
    if( doEresolSmear ) {
	// energy resolution smearing
	std::cerr << __LINE__ << std::endl; 
	eResolSmearer = new EnergySmearer( eSmearPars );
	eResolSmearer->name("E_res");
	eResolSmearer->doEnergy(false);
	eResolSmearer->scaleOrSmear(false);
	photonSmearers_.push_back(eResolSmearer);
    }
    if( doPhotonIdEffSmear ) {
	// photon ID efficiency 
	std::cerr << __LINE__ << std::endl; 
	idEffSmearer = new EfficiencySmearer( effSmearPars );
	idEffSmearer->name("idEff");
	idEffSmearer->setEffName("ratioTP");
	idEffSmearer->init();
	idEffSmearer->doPhoId(true);
	photonSmearers_.push_back(idEffSmearer);
    }
    if( doR9Smear ) {
	// R9 re-weighting
	r9Smearer = new EfficiencySmearer( effSmearPars );
	r9Smearer->name("r9Eff");
	r9Smearer->setEffName("ratioR9");
	r9Smearer->init();
	r9Smearer->doR9(true);
	photonSmearers_.push_back(r9Smearer);
    }
    if( doVtxEffSmear ) {
	// Vertex ID
	std::cerr << __LINE__ << std::endl; 
	vtxEffSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );   // triplicate TF1's here
	vtxEffSmearer->name("vtxEff");
	vtxEffSmearer->setEffName("ratioVertex");
	vtxEffSmearer->doVtxEff(true);
	vtxEffSmearer->init();
	diPhotonSmearers_.push_back(vtxEffSmearer);
    }
    if( doTriggerEffSmear ) {
	// trigger efficiency
	std::cerr << __LINE__ << std::endl; 
	triggerEffSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );
	triggerEffSmearer->name("triggerEff");
	triggerEffSmearer->setEffName("effL1HLT");
	triggerEffSmearer->doVtxEff(false);
	triggerEffSmearer->init();
	diPhotonSmearers_.push_back(triggerEffSmearer);
    }
    if(doKFactorSmear) {
	// kFactor efficiency
	std::cerr << __LINE__ << std::endl; 
	kFactorSmearer = new KFactorSmearer( kfacHist );
	kFactorSmearer->name("kFactor");
	kFactorSmearer->init();
	genLevelSmearers_.push_back(kFactorSmearer);
    }


    // RooContainer stuff
    // FIXME move these params to config file
    //l.rooContainer->SetNCategories(nCategories_);
    l.rooContainer->SetNCategories(1);
    l.rooContainer->nsigmas = nSystSteps;
    l.rooContainer->sigmaRange = systRange;
    // RooContainer does not support steps different from 1 sigma
    //assert( ((float)nSystSteps) == systRange );
    if( doEcorrectionSmear && doEcorrectionSyst ) {
        // instance of this smearer done in PhotonAnalysis
        systPhotonSmearers_.push_back(eCorrSmearer);
	std::vector<std::string> sys(1,eCorrSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEscaleSmear && doEscaleSyst ) {
	systPhotonSmearers_.push_back( eScaleSmearer );
	std::vector<std::string> sys(1,eScaleSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEresolSmear && doEresolSyst ) {
	systPhotonSmearers_.push_back( eResolSmearer );
	std::vector<std::string> sys(1,eResolSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doPhotonIdEffSmear && doPhotonIdEffSyst ) {
	systPhotonSmearers_.push_back( idEffSmearer );
	std::vector<std::string> sys(1,idEffSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doR9Smear && doR9Syst ) {
	systPhotonSmearers_.push_back( r9Smearer );
	std::vector<std::string> sys(1,r9Smearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doVtxEffSmear && doVtxEffSyst ) {
	systDiPhotonSmearers_.push_back( vtxEffSmearer );
	std::vector<std::string> sys(1,vtxEffSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doTriggerEffSmear && doTriggerEffSyst ) {
	systDiPhotonSmearers_.push_back( triggerEffSmearer );
	std::vector<std::string> sys(1,triggerEffSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if(doKFactorSmear && doKFactorSyst) {
	systGenLevelSmearers_.push_back(kFactorSmearer);
	std::vector<std::string> sys(1,kFactorSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
	
    // ----------------------------------------------------
    // ----------------------------------------------------
    // Global systematics - Lumi
    l.rooContainer->AddGlobalSystematic("lumi",1.06,1.00);
    // ----------------------------------------------------


    l.rooContainer->AddObservable("CMS_hgg_mass",massMin,massMax);

    l.rooContainer->AddConstant("IntLumi",l.intlumi_);
 
    l.rooContainer->AddRealVar("r1",-4,-10.,0);
    l.rooContainer->AddRealVar("r2",-1,-10.,0);
    l.rooContainer->AddRealVar("f2",0.2,0.,1.);

    // Power law
    std::vector<std::string> data_pow_pars(3,"p");	 
    data_pow_pars[0] = "r1";
    data_pow_pars[1] = "r2";
    data_pow_pars[2] = "f2";
    l.rooContainer->AddGenericPdf("data_pow_model", "(1-@3)*TMath::Power(@0,@1) + @3*TMath::Power(@0,@2)","CMS_hgg_mass",data_pow_pars,0);
        
    if (doTraining){
	TString outfileName( "TMVA_input_" + (std::string) l.histFileName);
        mvaFile_ = TFile::Open( outfileName, "RECREATE" );
        mvaFile_->cd();
        for (int i = 0; i<nMasses;i++){
            signalTrainTree_[i] = new TTree(("sig_train"+names[i]).c_str(),("SignalTrainTree"+names[i]).c_str());
            SetBDTInputTree(signalTrainTree_[i]);
            signalTestTree_[i]  = new TTree(("sig_test"+names[i]).c_str(), ("SignalTestTree"+names[i]).c_str());
            SetBDTInputTree(signalTestTree_[i]);

            backgroundTrainTree_[i] = new TTree(("bkg_train"+names[i]).c_str(),("BackgroundTrainTree"+names[i]).c_str());
            SetBDTInputTree(backgroundTrainTree_[i]);
            backgroundTestTree_[i]  = new TTree(("bkg_test"+names[i]).c_str(), ("BackgroundTestTree"+names[i]).c_str());
            SetBDTInputTree(backgroundTestTree_[i]);
        }
        backgroundTrainTree_all = new TTree("bkg_train_all" ,"BackgroundTrainTree_all");
        SetBDTInputTree(backgroundTrainTree_all);
        backgroundTestTree_all  = new TTree("bkg_test_all"   ,"BackgroundTestTree_all");
        SetBDTInputTree(backgroundTestTree_all);
    }
    else{

        l.rooContainer->AddObservable("BDT" ,-1.,1.);

        //Set up TMVA reader
        tmvaReader_= new TMVA::Reader();
 		tmvaReader_->AddVariable("H_ptOverM", &_H_ptOverM);
 		tmvaReader_->AddVariable("H_eta", &_H_eta);
 		tmvaReader_->AddVariable("d_phi", &_d_phi);
 		tmvaReader_->AddVariable("max_eta", &_max_eta);
 		tmvaReader_->AddVariable("min_r9", &_min_r9);
 		tmvaReader_->AddVariable("pho1_eta", &_pho1_eta);
 		tmvaReader_->AddVariable("pho2_eta", &_pho2_eta);
 		tmvaReader_->AddVariable("pho1_ptOverM", &_pho1_ptOverM);
 		tmvaReader_->AddVariable("pho2_ptOverM", &_pho2_ptOverM);
 	    tmvaReader_->AddVariable("deltaMOverM", &_deltaMOverM);
 		tmvaReader_->AddVariable("sigmaMOverM", &_sigmaMOverM);

        //Invariant Mass Spectra
        l.rooContainer->CreateDataSet("CMS_hgg_mass","data_mass",nDataBins);
        l.rooContainer->CreateDataSet("CMS_hgg_mass","bkg_mass" ,nDataBins);    	 
		
	int nBDTbins = 5000;

        // Usual datasets

        for (double mass=115.0; mass<150.5; mass+=0.5){

            //Adaptive Boost
            l.rooContainer->CreateDataSet("BDT",Form("data_BDT_ada_%3.1f",mass)	     ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT",Form("bkg_BDT_ada_%3.1f",mass)       ,nBDTbins);

            l.rooContainer->CreateDataSet("BDT",Form("data_BDT_grad_%3.1f",mass)     ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT",Form("bkg_BDT_grad_%3.1f",mass)      ,nBDTbins);

	    for (int sideband_i=1;sideband_i<=numberOfSidebands;sideband_i++){

              l.rooContainer->CreateDataSet("BDT",Form("bkg_%dlow_BDT_ada_%3.1f",sideband_i,mass)   ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("bkg_%dhigh_BDT_ada_%3.1f",sideband_i,mass)  ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("data_%dlow_BDT_ada_%3.1f",sideband_i,mass)  ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("data_%dhigh_BDT_ada_%3.1f",sideband_i,mass) ,nBDTbins);
		
              l.rooContainer->CreateDataSet("BDT",Form("bkg_%dlow_BDT_grad_%3.1f",sideband_i,mass)  ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("bkg_%dhigh_BDT_grad_%3.1f",sideband_i,mass) ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("data_%dlow_BDT_grad_%3.1f",sideband_i,mass) ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("data_%dhigh_BDT_grad_%3.1f",sideband_i,mass),nBDTbins);

	    }

        }

        // loop signal mass points signal datasets
        for (int i = 2; i<nMasses;i++){  // We are ignoring masses 105 and 110 for now
          for (int j=-1; j<2; j++){
            if ((i==2 && j==-1) || (i==nMasses-1 && j==1)) continue;
            l.rooContainer->CreateDataSet("BDT","sig_BDT_ada"+names[i]+names[i+j]       ,nBDTbins);    
            l.rooContainer->CreateDataSet("BDT","sig_BDT_grad"+names[i]+names[i+j]      ,nBDTbins);  
        
            // Make the signal Systematic Sets
            l.rooContainer->MakeSystematics("BDT","sig_BDT_grad"+names[i]+names[i+j] ,-1)	;
            l.rooContainer->MakeSystematics("BDT","sig_BDT_ada"+names[i]+names[i+j]  ,-1)	;
          }

        }
        //TMVA Reader
        tmvaReader_->BookMVA("BDT_ada_123", mvaWeightsFolder+"/TMVAClassification_BDT_ada_123_all.weights.xml");
        tmvaReader_->BookMVA("BDT_grad_123",mvaWeightsFolder+"/TMVAClassification_BDT_grad_123_all.weights.xml");
    }

    if(PADEBUG) 
	cout << "InitRealMvaAnalysis END"<<endl;
}

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{

    if(PADEBUG) 
	cout << "Analysis START; cur_type is: " << l.itype[l.current] <<endl;
    int cur_type = l.itype[l.current];
    float weight = l.sampleContainer[l.current_sample_index].weight;

    if (!doTraining){
	if (splitSignalSample  && cur_type < 0) {
	  if (jentry%2!=0) return;
	  else weight*=2;
	}

	if (splitBackgroundSample && cur_type > 0) {
	  if (jentry%2!=0) return;
	  else weight*=2;
	}
    }
	
 
    l.FillCounter( "Processed", 1. );
    assert( weight > 0. );  
    l.FillCounter( "XSWeighted", weight );
    nevents+=1.;

    //PU reweighting
    unsigned int n_pu = l.pu_n;
    if ( cur_type !=0 && puHist != "") {
	bool hasSpecificWeight = weights.find( cur_type ) != weights.end() ; 
	if( cur_type < 0 && !hasSpecificWeight && jentry == 1 ) {
	    std::cerr  << "WARNING no pu weights specific for sample " << cur_type << std::endl;
	}
	std::vector<double> & puweights = hasSpecificWeight ? weights[ cur_type ] : weights[0]; 
	if(n_pu<puweights.size()){
	    weight *= puweights[n_pu];
	    sumwei+=puweights[n_pu]; 
	}    
	else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	    cout <<"n_pu ("<< n_pu<<") too big ("<<puweights.size()<<") ["<< l.itype[l.current]<<"], event will not be reweighted for pileup"<<endl;
	}
    }
    
    assert( weight >= 0. );  
    l.FillCounter( "PUWeighted", weight );
    
    if( jentry % 10000 ==  0 ) {
	    std::cout << " nevents " <<  nevents << " sumpuweights " << sumwei << " ratio " << sumwei / nevents 
		      << " equiv events " << sumev << " accepted " << sumaccept << " smeared " << sumsmear << " "  
		      <<  sumaccept / sumev << " " << sumsmear / sumaccept
		      << std::endl;
    }
    // ------------------------------------------------------------
    //PT-H K-factors
    double gPT = 0;
    TLorentzVector gP4(0,0,0,0);
    if (cur_type<0){            // if background sample, gP4 remains 4vect(0)
	for (int gi=0;gi<l.gp_n;gi++){
	    if (l.gp_pdgid[gi]==25){
		gP4 = *((TLorentzVector*)l.gp_p4->At(gi));
		gPT = gP4.Pt();
		break;
	    }
	}
    }

    // ------------------------------------------------------------
    // smear all of the photons!
    std::pair<int,int> diphoton_index;
   
    // do gen-level dependent first (e.g. k-factor); only for signal
    double genLevWeight=1; 
    if(cur_type<0){
	for(std::vector<BaseGenLevelSmearer*>::iterator si=genLevelSmearers_.begin(); si!=genLevelSmearers_.end(); si++){
	    float genWeight=1;
	    (*si)->smearEvent( genWeight,gP4, l.pu_n, cur_type, 0. );
	    if( genWeight < 0. ) {
		std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
		assert(0);
	    }
	    genLevWeight*=genWeight;
	}
    }

    // Nominal smearing
    std::vector<float> smeared_pho_energy(l.pho_n,0.); 
    std::vector<float> smeared_pho_r9(l.pho_n,0.); 
    std::vector<float> smeared_pho_weight(l.pho_n,1.);
   
    for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
	std::vector<std::vector<bool> > p;
	PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), 
				    // *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
				    ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), 
				    energyCorrected[ipho],
				    l.pho_isEB[ipho], l.pho_r9[ipho],
				    l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_),
				    (energyCorrectedError!=0?energyCorrectedError[ipho]:0)
				    );
	float pweight = 1.;
	// smear MC. But apply energy shift to data 
	if( cur_type != 0 && doMCSmearing ) { // if it's MC
	    for(std::vector<BaseSmearer *>::iterator si=photonSmearers_.begin(); si!= photonSmearers_.end(); ++si ) {
		float sweight = 1.;
		(*si)->smearPhoton(phoInfo,sweight,l.run,0.);	   
		if( sweight < 0. ) {
			std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
			assert(0);
		}
		pweight *= sweight;
	    }
	} else if(cur_type == 0 ) {          // if it's data
	    float sweight = 1.; 
	    if (doEcorrectionSmear){
	      eCorrSmearer->smearPhoton(phoInfo,sweight,l.run,0.);     // This Smearer is the same as for MC so can just re-use it
	    }
	    eScaleDataSmearer->smearPhoton(phoInfo,sweight,l.run,0.);
	    pweight *= sweight;	
	}

	smeared_pho_energy[ipho] = phoInfo.energy();
	smeared_pho_r9[ipho] = phoInfo.r9();
	smeared_pho_weight[ipho] = pweight;
    }
   
    sumev += weight;
    // FIXME pass smeared R9
    int diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
    /// std::cerr << "Selected pair " << l.dipho_n << " " << diphoton_id << std::endl;
    if (diphoton_id > -1 ) {

	diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
    	// bring all the weights together: lumi & Xsection, single gammas, pt kfactor
	float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;

	l.countersred[diPhoCounter_]++;

	TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
	float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
	TLorentzVector Higgs = lead_p4 + sublead_p4; 	
	TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);

	bool CorrectVertex;
	// FIXME pass smeared R9
	int category = 0;//l.gge Society says thDiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
	int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
	if( cur_type != 0 && doMCSmearing ) {
	    float pth = Higgs.Pt();
	    for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
		float rewei=1.;
		(*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
		if( rewei < 0. ) {
		    std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
		    assert(0);
		}
		evweight *= rewei;
	    }
	    CorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
	}
	float mass    = Higgs.M();
	float ptHiggs = Higgs.Pt();
      
	assert( evweight >= 0. ); 


	// Mass Resolution of the Event
	massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);

	double massResolution = massResolutionCalculator->massResolution();


        if (doTraining ){//train on EVEN!!!!! EVEN!!!!
            if (cur_type > 0 ){// Background 
                for (int i = 0; i<bkg_masses.size();i++) {//Nicks alternative background method
                    float mass_hypothesis = bkg_masses[i];
                    float boundaries[2];
                    boundaries[0] = mass_hypothesis*(1-signalRegionWidth);
                    boundaries[1] = mass_hypothesis*(1+signalRegionWidth);
                    if( mass>boundaries[0] && mass<boundaries[1] ){//Signal mass window cut
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis,evweight,category);
                        _sideband = i;
                        if (jentry%2==0) backgroundTrainTree_all->Fill();
                        else if (jentry%2==1) backgroundTestTree_all->Fill();
                    }
                }
                for (int i = 0; i<nMasses;i++) {
                    // define hypothesis masses for the sidebands
                    float mass_hypothesis = masses[i];
                    float mass_hypothesis_low = mass_hypothesis*(1-signalRegionWidth)/(1+signalRegionWidth);
                    float mass_hypothesis_high = mass_hypothesis*(1+signalRegionWidth)/(1-signalRegionWidth);
                    // define the sidebands
                    float sideband_boundaries[4];
                    sideband_boundaries[0] = mass_hypothesis_low*(1-signalRegionWidth);
                    sideband_boundaries[1] = mass_hypothesis*(1-signalRegionWidth);
                    sideband_boundaries[2] = mass_hypothesis*(1+signalRegionWidth);
                    sideband_boundaries[3] = mass_hypothesis_high*(1+signalRegionWidth);

                    if( mass>sideband_boundaries[1] && mass<sideband_boundaries[2] ){//Signal mass window cut
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis,evweight,category);
                        _sideband = 0;
                        if (jentry%2==0) backgroundTrainTree_[i]->Fill();
                        else if (jentry%2==1) backgroundTestTree_[i]->Fill();
                    }
                    else if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1] ){//lower sideband mass window cut
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis_low,evweight,category);
                        _sideband = -1;
                        if (jentry%2==0) backgroundTrainTree_[i]->Fill();
                        else if (jentry%2==1) backgroundTestTree_[i]->Fill();
                    }
                    else if( mass>sideband_boundaries[2] && mass<sideband_boundaries[3] ){//upper sideband mass window cut
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis_high,evweight,category);
                        _sideband = 1;
                        if (jentry%2==0) backgroundTrainTree_[i]->Fill();
                        else if (jentry%2==1) backgroundTestTree_[i]->Fill();
                    }
                } 
            }
            else { //Signal 
                int i0 = SignalType(cur_type); 
                if (i0<0){
                    cout<<"CAN'T TRAIN ON DATA!\n";
                    return;
                }
                float mass_hypothesis = masses[i0];
                float sideband_boundaries[2];
                sideband_boundaries[0] = mass_hypothesis*(1-signalRegionWidth);
                sideband_boundaries[1] = mass_hypothesis*(1+signalRegionWidth);
                 
                if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1] ){//Signal mass window cut
                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,masses[i0],evweight,category);
                    _sideband = 0;
                    if (jentry%2==0) signalTrainTree_[i0]->Fill();
                    else if (jentry%2==1) signalTestTree_[i0]->Fill();
                }
            }
        }
        else{

        // --- Fill invariant mass spectrum -------
          if (cur_type==0){  // Data
            l.rooContainer->InputDataPoint("data_mass",category,mass);
          } else if (cur_type>0){ // Background MC
            l.rooContainer->InputDataPoint("bkg_mass",category,mass,evweight);
          }
           
          // ------ Deal with Signal MC first
            if (cur_type<0){ // signal MC
            // Iterate over each BDT mass. 
            for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
              if (SignalType(cur_type)<i-1 || SignalType(cur_type) > i+1) continue;
                // loop over +- signal masses
              for (int j=-1; j<2; j++){
                if ((i==2 && j==-1) || (i==nMasses-1 && j==1)) continue;
                // define hypothesis masses for the sidebands
                float mass_hypothesis = masses[i+j];

                // define the sidebands
                float sideband_boundaries[2];
                sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
                sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);

                //Signal Window
                if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis,evweight);

                    float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                    float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

                    if (SignalType(cur_type)==i+j){
                        l.rooContainer->InputDataPoint("sig_BDT_ada"+names[i]+names[i+j],category,bdt_ada,evweight);
                        l.rooContainer->InputDataPoint("sig_BDT_grad"+names[i]+names[i+j] ,category,bdt_grad,evweight);
                    }
                 }
               } // end loop over sig masses
             } // end loop over BDTs
           }
           // ---- Now deal with background MC and data
           else {

            // Iterate over each BDT mass. 
            for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now

	      double mass_h_low;      
              double mass_h_high;
 
	      if (i==7){ // the 140 case because currently 145 is missing
	       mass_h_low  =-2.5;
	       mass_h_high =4.6;
	      } else if (i==8){  // the 150 case
	       mass_h_low  =-5.0;
	       mass_h_high =0.1;
	      } else {
	       mass_h_low  =-2.5;
	       mass_h_high =2.1;
	      }

              for (double mass_h=mass_h_low; mass_h<mass_h_high; mass_h+=0.5){
           
                float mass_hypothesis = masses[i]+mass_h;
	        if (mass_hypothesis < masses[2] || mass_hypothesis > masses[nMasses-1]) continue;
              
                float mass_hypothesis_low = mass_hypothesis*(1-signalRegionWidth)/(1+signalRegionWidth) -sidebandShift;
                float mass_hypothesis_high = mass_hypothesis*(1+signalRegionWidth)/(1-signalRegionWidth)+sidebandShift;

                float sideband_boundaries[2];
                sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
                sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);
              
                //Signal Window
                if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                  //std::cout << "sig region" << std::endl;
                  SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis,evweight);
                    
                  float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                  float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

                  //std::cout << "ada: " << bdt_ada << "  grad: " << bdt_grad << std::endl;

                  if (cur_type == 0 ){//data
                      l.rooContainer->InputDataPoint(Form("data_BDT_ada_%3.1f",mass_hypothesis),category,bdt_ada,evweight);
                      l.rooContainer->InputDataPoint(Form("data_BDT_grad_%3.1f",mass_hypothesis),category,bdt_grad,evweight);
                  }
                  if (cur_type > 0 ){// background MC
                      l.rooContainer->InputDataPoint(Form("bkg_BDT_ada_%3.1f",mass_hypothesis),category,bdt_ada,evweight);
                      l.rooContainer->InputDataPoint(Form("bkg_BDT_grad_%3.1f",mass_hypothesis) ,category,bdt_grad,evweight);
                  }
                }

		else {
		  // Loop over N lower sidebands
		  for (int sideband_i = 1 ; sideband_i <= numberOfSidebands ; sideband_i++){

		   double hypothesisModifier = (1.-sidebandWidth)/(1+sidebandWidth) - sidebandShift;
		   mass_hypothesis_low     = (mass_hypothesis*(1.-signalRegionWidth)/(1.+signalRegionWidth)-sidebandShift)*(TMath::Power(hypothesisModifier,sideband_i-1)) - sidebandShift;
		   double sideband_boundaries_low = mass_hypothesis_low*(1.-sidebandWidth);
		   double sideband_boundaries_high= mass_hypothesis_low*(1.+sidebandWidth);

		   if ( mass>sideband_boundaries_low && mass<sideband_boundaries_high){
                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis_low,evweight);
                    float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                    float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

                    if (cur_type == 0 ){//data
                         l.rooContainer->InputDataPoint(Form("data_%dlow_BDT_ada_%3.1f",sideband_i,mass_hypothesis),category,bdt_ada,evweight);
                         l.rooContainer->InputDataPoint(Form("data_%dlow_BDT_grad_%3.1f",sideband_i,mass_hypothesis) ,category,bdt_grad,evweight);
                    }
                    else if (cur_type > 0 ){// background MC
                        l.rooContainer->InputDataPoint(Form("bkg_%dlow_BDT_ada_%3.1f",sideband_i,mass_hypothesis),category,bdt_ada,evweight);
                        l.rooContainer->InputDataPoint(Form("bkg_%dlow_BDT_grad_%3.1f",sideband_i,mass_hypothesis),category,bdt_grad,evweight);
                    }
		   }
		  }

		  // Loop over N higher sidebands
		  for (int sideband_i = 1 ; sideband_i <= numberOfSidebands ; sideband_i++){

		   double hypothesisModifier = (1.+sidebandWidth)/(1-sidebandWidth) + sidebandShift;
		   mass_hypothesis_high     = (mass_hypothesis*(1.+signalRegionWidth)/(1.-signalRegionWidth)+sidebandShift)*(TMath::Power(hypothesisModifier,sideband_i-1)) + sidebandShift;
		   double sideband_boundaries_low = mass_hypothesis_high*(1.-sidebandWidth);
		   double sideband_boundaries_high= mass_hypothesis_high*(1.+sidebandWidth);

		   if ( mass>sideband_boundaries_low && mass<sideband_boundaries_high){
                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis_high,evweight);
                    float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                    float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

                    if (cur_type == 0 ){//data
                         l.rooContainer->InputDataPoint(Form("data_%dhigh_BDT_ada_%3.1f",sideband_i,mass_hypothesis),category,bdt_ada,evweight);
                         l.rooContainer->InputDataPoint(Form("data_%dhigh_BDT_grad_%3.1f",sideband_i,mass_hypothesis) ,category,bdt_grad,evweight);
                    }
                    else if (cur_type > 0 ){// background MC
                        l.rooContainer->InputDataPoint(Form("bkg_%dhigh_BDT_ada_%3.1f",sideband_i,mass_hypothesis),category,bdt_ada,evweight);
                        l.rooContainer->InputDataPoint(Form("bkg_%dhigh_BDT_grad_%3.1f",sideband_i,mass_hypothesis),category,bdt_grad,evweight);
                    }
		   }
		  }
              } 
            } // end loop over mass hypotheses
        } // end loop over BDT
      } // Background MC or Data

      l.FillCounter( "Accepted", weight );
      l.FillCounter( "Smeared", evweight );
      sumaccept += weight;
      sumsmear += evweight;

     }
    }
    if(PADEBUG) 
	cout<<"myFillHistRed END"<<endl;


    // Now do some MC Systematic studies - For MVA Analysis, assume that we only care about Signal Systematics Here, ie make the check cur_type < 0 to increase speed
    // --------------------------------------------------------------------------------------------------------------------------------------------------------------
    if( cur_type < 0 && doMCSmearing && !doTraining){
	// fill steps for syst uncertainty study
	float systStep = systRange / (float)nSystSteps;
	// di-photon smearers systematics
	if (diphoton_id > -1 ) {
	       
	    float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
	    float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
	    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	    TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
	 
	    for(std::vector<BaseGenLevelSmearer*>::iterator si=systGenLevelSmearers_.begin(); si!=systGenLevelSmearers_.end(); si++){
		std::vector<double> bdt_grad_errors[3];
		std::vector<double> bdt_ada_errors[3];
		std::vector<double> weights[3];
		std::vector<int>    categories[3];
	   
		for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		    if( syst_shift == 0. ) { continue; } // skip the central value
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
	     
		    int category = 0;  // Category always 0 fro MVA analysis
		    double genLevWeightSyst=1; 
	     
		    for(std::vector<BaseGenLevelSmearer *>::iterator sj=genLevelSmearers_.begin(); sj!= genLevelSmearers_.end(); ++sj ) {
			float swei=1.;
			if( *si == *sj ) { 
			    (*si)->smearEvent(swei, gP4, l.pu_n, cur_type, syst_shift );
			} else {
			    (*sj)->smearEvent(swei, gP4, l.pu_n, cur_type, 0. );
			}
			genLevWeightSyst *= swei;
		    }

		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeightSyst;
		    float mass = Higgs.M();
		    float ptHiggs = Higgs.Pt();

		    // Mass Resolution of the Event
		    massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);

		    double massResolution = massResolutionCalculator->massResolution();

                    // Iterate over each mass point. 
	            for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
              		if (SignalType(cur_type) < i-1 || SignalType(cur_type) > i+1) continue;
              		  for (int j=-1; j<2; j++){ // loop over +- signal masses
                	    if ((i==2 && j==-1) || (i==nMasses-1 && j==1)) continue;

               		    // define hypothesis masses for the sidebands
                	    float mass_hypothesis = masses[i+j];
                	    // define the sidebands
                	    float sideband_boundaries[2];
                	    sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
                	    sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);

	     	    	    if ((SignalType(cur_type)==i+j)){
                     	    //Signal Window
                		if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                    		SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis,evweight);
                    		float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                    		float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
			
                    		categories[j+1].push_back(category);
                    		bdt_ada_errors[j+1].push_back(bdt_ada);
                    		bdt_grad_errors[j+1].push_back(bdt_grad);
                    		weights[j+1].push_back(evweight);

                		} else {

                    		categories[j+1].push_back(-1);
                    		bdt_ada_errors[j+1].push_back(-100.);
                    		bdt_grad_errors[j+1].push_back(-100.);
                    		weights[j+1].push_back(0.);
               			}
	     		   }
            		 }
        	     }

	    	}// end loop on systematics steps
		// Fill In the Corect Systematic Set ---------------------------------------------------------------------------//
                // Iterate over each mass point. 
                for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
                  if (SignalType(cur_type)<i-1 || SignalType(cur_type) > i+1) continue;
                    // loop over +- signal masses
                  for (int j=-1; j<2; j++){
                    if ((i==2 && j==-1) || (i==nMasses-1 && j==1)) continue;
	     	      if ((SignalType(cur_type)!=i+j)) continue;
                      l.rooContainer->InputSystematicSet("sig_BDT_ada"+names[i]+names[i+j],(*si)->name(),categories[j+1],bdt_ada_errors[j+1],weights[j+1]);
                      l.rooContainer->InputSystematicSet("sig_BDT_grad"+names[i]+names[i+j],(*si)->name(),categories[j+1],bdt_grad_errors[j+1],weights[j+1]);
                 }
	        }
		// -------------------------------------------------------------------------------------------------------------//
 
	    } // end loop on smearers 
		 

	    for(std::vector<BaseDiPhotonSmearer *>::iterator si=systDiPhotonSmearers_.begin(); si!= systDiPhotonSmearers_.end(); ++si ) {
		std::vector<double> bdt_grad_errors[3];
		std::vector<double> bdt_ada_errors[3];
		std::vector<double> weights[3];
		std::vector<int> categories[3];
		       
		for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		    if( syst_shift == 0. ) { continue; } // skip the central value
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
			       
		    // restart with 'fresh' wait for this round of systematics
		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
			       
		    // FIXME pass smeared R9 and di-photon
		    int category = 0; 
		    int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
		    for(std::vector<BaseDiPhotonSmearer *>::iterator sj=diPhotonSmearers_.begin(); sj!= diPhotonSmearers_.end(); ++sj ) {
			float swei=1.;
			float pth = Higgs.Pt();
			if( *si == *sj ) { 
			    (*si)->smearDiPhoton( Higgs, *vtx, swei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), syst_shift );
			} else { 
			    (*sj)->smearDiPhoton( Higgs, *vtx, swei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
			}
			evweight *= swei;
		    }
		    float mass = Higgs.M();
		    float ptHiggs = Higgs.Pt();
		    // Mass Resolution of the Event
		    massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);

		    double massResolution = massResolutionCalculator->massResolution();
                   // Iterate over each mass point. 
        	   for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
          	     if (SignalType(cur_type) < i-1 || SignalType(cur_type) > i+1) continue;
         	     for (int j=-1; j<2; j++){ // loop over +- signal masses
           		if ((i==2 && j==-1) || (i==nMasses-1 && j==1)) continue;

           		// define hypothesis masses for the sidebands
           		float mass_hypothesis = masses[i+j];
		        // define the sidebands
           		float sideband_boundaries[2];
           		sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
           		sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);

	 		if ((SignalType(cur_type)==i+j)){
           		//Signal Window
           		 if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
             		  SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis,evweight);
             		  float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
             		  float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

             		  categories[j+1].push_back(category);
             		  bdt_ada_errors[j+1].push_back(bdt_ada);
             		  bdt_grad_errors[j+1].push_back(bdt_grad);
             		  weights[j+1].push_back(evweight);

		    	 } else {

		       	  categories[j+1].push_back(-1);
		       	  bdt_ada_errors[j+1].push_back(-100.);
		       	  bdt_grad_errors[j+1].push_back(-100.);
		       	  weights[j+1].push_back(0.);
		         }
	     		}
		     }	
    		  }
    		}// end loop on systematics steps

		// Fill In the Corect Systematic Set ---------------------------------------------------------------------------//
                // Iterate over each mass point. 
                for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
                  if (SignalType(cur_type)<i-1 || SignalType(cur_type) > i+1) continue;
                    // loop over +- signal masses
                  for (int j=-1; j<2; j++){
                    if ((i==2 && j==-1) || (i==nMasses-1 && j==1)) continue;
	     	    if ((SignalType(cur_type)!=i+j)) continue;
                      l.rooContainer->InputSystematicSet("sig_BDT_ada"+names[i]+names[i+j],(*si)->name(),categories[j+1],bdt_ada_errors[j+1],weights[j+1]);
                      l.rooContainer->InputSystematicSet("sig_BDT_grad"+names[i]+names[i+j],(*si)->name(),categories[j+1],bdt_grad_errors[j+1],weights[j+1]);
                 }
	        }
		// -------------------------------------------------------------------------------------------------------------//
 
	    } // end loop on smearers 
		 
	} // Close If on CiC Selection
	// Now the systematic Steps which can effect the CiC Selection

       
	// loop over the smearers included in the systematics study
	for(std::vector<BaseSmearer *>::iterator  si=systPhotonSmearers_.begin(); si!= systPhotonSmearers_.end(); ++si ) {
	    std::vector<double> bdt_grad_errors[3];
	    std::vector<double> bdt_ada_errors[3];
	    std::vector<double> weights[3];
	    std::vector<int> categories[3];
	   
	    // loop over syst shift
	    for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		if( syst_shift == 0. ) { continue; } // skip the central value
		// smear the photons 
		for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
		    std::vector<std::vector<bool> > p;
		    PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), 
						/// *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
						((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), 
						energyCorrected[ipho],
						l.pho_isEB[ipho], l.pho_r9[ipho],
						l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_),
			    		        (energyCorrectedError!=0?energyCorrectedError[ipho]:0)
				    	      );
		  
		    float pweight = 1.;
		    for(std::vector<BaseSmearer *>::iterator  sj=photonSmearers_.begin(); sj!= photonSmearers_.end(); ++sj ) {
			float sweight = 1.;
			if( *si == *sj ) {
			    // move the smearer under study by syst_shift
			    (*si)->smearPhoton(phoInfo,sweight,l.run,syst_shift);
			} else {
			    // for the other use the nominal points
			    (*sj)->smearPhoton(phoInfo,sweight,l.run,0.);
			}
			pweight *= sweight;
		    }
		    smeared_pho_energy[ipho] = phoInfo.energy();
		    smeared_pho_r9[ipho] = phoInfo.r9();
		    smeared_pho_weight[ipho] = pweight;
		}
	       
		// analyze the event
		// FIXME pass smeared R9
		int diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
	       
		if (diphoton_id > -1 ) {
		   
		    diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] *genLevWeight;
		   
		    float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
		    float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
		    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
		    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
		    TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
		   
		    int category = 0; 
		    int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
		    if( cur_type != 0 && doMCSmearing ) {
			for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
			    float rewei=1.;
			    float pth = Higgs.Pt();
			    (*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
			    evweight *= rewei;
			}
		    }
		    float mass = Higgs.M();
			
		    float ptHiggs = Higgs.Pt();
		    // Mass Resolution of the Event
		    massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);

		    double massResolution = massResolutionCalculator->massResolution();

                   // Iterate over each mass point. 
            	   for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
                     if (SignalType(cur_type) < i-1 || SignalType(cur_type) > i+1) continue;
              	     for (int j=-1; j<2; j++){ // loop over +- signal masses
                	if ((i==2 && j==-1) || (i==nMasses-1 && j==1)) continue;

              		// define hypothesis masses for the sidebands
               		float mass_hypothesis = masses[i+j];
               		// define the sidebands
           		float sideband_boundaries[2];
           		sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
           		sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);

	 		if ((SignalType(cur_type)==i+j)){
               		  //Signal Window
               		  if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                 		SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis,evweight);
                 		float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                 		float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                 
                 		categories[j+1].push_back(category);
                 		bdt_ada_errors[j+1].push_back(bdt_ada);
                 		bdt_grad_errors[j+1].push_back(bdt_grad);
                 		weights[j+1].push_back(evweight);

		    	  } else {

		        	categories[j+1].push_back(-1);
		       		bdt_ada_errors[j+1].push_back(-100.);
		       		bdt_grad_errors[j+1].push_back(-100.);
		       		weights[j+1].push_back(0.);
		    	  }
			}
		     }
    		   }

		   } else { // In case CiC selection fails now
      		     for (int j=-1; j<2; j++){ // loop over +- signal masses
         		categories[j+1].push_back(-1);
         		bdt_ada_errors[j+1].push_back(-100.);
         		bdt_grad_errors[j+1].push_back(-100.);
         		weights[j+1].push_back(0.);
      		     }
		   }

	  }// end loop on systematics steps

	  // Fill In the Corect Systematic Set ---------------------------------------------------------------------------//
          // Iterate over each mass point. 
          for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
              if (SignalType(cur_type)<i-1 || SignalType(cur_type) > i+1) continue;
                // loop over +- signal masses
              for (int j=-1; j<2; j++){
                if ((i==2 && j==-1) || (i==(nMasses-1) && j==1)) continue;

	     	  if ((SignalType(cur_type)!=i+j)) continue;

                  l.rooContainer->InputSystematicSet("sig_BDT_ada"+names[i]+names[i+j],(*si)->name(),categories[j+1],bdt_ada_errors[j+1],weights[j+1]);
                  l.rooContainer->InputSystematicSet("sig_BDT_grad"+names[i]+names[i+j],(*si)->name(),categories[j+1],bdt_grad_errors[j+1],weights[j+1]);
              }
	      
          }
	  // -------------------------------------------------------------------------------------------------------------//
 
	} // Close on Smearers
    } // End Signal Systematic Study 
    // --------------------------------------------------------------------------------------------------------------------------------------------------------------
}

// ----------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
    vtxAna_.setBranchAdresses(t,"vtx_std_");
    vtxAna_.getBranches(t,"vtx_std_",s);
}

// ----------------------------------------------------------------------------------------------------
bool MvaAnalysis::SelectEvents(LoopAll& l, int jentry) 
{
    return true;
}
// ----------------------------------------------------------------------------------------------------

int MvaAnalysis::SignalType(int cur_type){
    int i0 = -1;
    if (cur_type == -13 || cur_type == -14 ||  cur_type == -15 || cur_type == -16){//105 
        i0 = 0;}
    else if (cur_type == -17 || cur_type == -18 ||  cur_type == -19 || cur_type == -20){//110 
        i0 = 1;}
    else if (cur_type == -21 || cur_type == -22 ||  cur_type == -23 || cur_type == -24){//115 
        i0 = 2;}
    else if (cur_type == -25 || cur_type == -26 ||  cur_type == -27 || cur_type == -28){//120 
        i0 = 3;}
    else if (cur_type == -37 || cur_type == -38 ||  cur_type == -39 || cur_type == -40){//125 
        i0 = 4;}
    else if (cur_type == -29 || cur_type == -30 ||  cur_type == -31 || cur_type == -32){//130 
        i0 = 5;}
    else if (cur_type == -41 || cur_type == -42 ||  cur_type == -43 || cur_type == -44){//135 
        i0 = 6;}
    else if (cur_type == -33 || cur_type == -34 ||  cur_type == -35 || cur_type == -36){//140 
        i0 = 7;}
//    else if (cur_type == -45 || cur_type == -46 ||  cur_type == -47 || cur_type == -48){//145 
//        i0 = 8;}
    else if (cur_type == -49 || cur_type == -50 ||  cur_type == -51 || cur_type == -52){//150 
        i0 = 8;}
//    else if (cur_type == -53 || cur_type == -54 ||  cur_type == -55 || cur_type == -56){//121 
//        i0 = 10;}
//    else if (cur_type == -57 || cur_type == -58 ||  cur_type == -59 || cur_type == -60){//123 
//        i0 = 11;}
    return i0;
}


void MvaAnalysis::SetBDTInputVariables(TLorentzVector *lead_p4, TLorentzVector *sublead_p4, 
                                       double lead_r9, double sublead_r9, 
                                       double massResolution, double mass_hypothesis, double evweight, int cat){

	TLorentzVector Higgs = *lead_p4 + *sublead_p4; 	
	float mass    = Higgs.M();

    _log_H_pt =  log10( Higgs.Pt());
    _H_eta = fabs(Higgs.Eta());
    _d_phi = fabs(lead_p4->DeltaPhi(*sublead_p4));
    _max_eta = max(fabs(lead_p4->Eta()),fabs(sublead_p4->Eta()));
    _min_r9  = min(lead_r9,sublead_r9);
    _pho1_eta = lead_p4->Eta();
    _pho2_eta = sublead_p4->Eta();

    _mgg = mass;
    _pho1_phi = lead_p4->Phi();
    _pho1_pt = lead_p4->Pt();
    _pho1_r9 = lead_r9;

    _pho2_phi = sublead_p4->Phi();
    _pho2_pt = sublead_p4->Pt();
    _pho2_r9 = sublead_r9;

    _H_pt = Higgs.Pt();
    _Ht = lead_p4->Pt()+sublead_p4->Pt();

    _d_eta = lead_p4->Eta()-sublead_p4->Eta();
    _mod_d_eta = fabs(lead_p4->Eta()-sublead_p4->Eta());
    _cos_theta_star = fabs(lead_p4->E()-sublead_p4->E())/Higgs.P();

    _wt= evweight;

    _pho1_ptOverM = lead_p4->Pt()/mass_hypothesis;
    _pho2_ptOverM = sublead_p4->Pt()/mass_hypothesis;
    _deltaMOverM = (mass-mass_hypothesis)/mass_hypothesis;
    _deltaMOverSigmaM = (mass-mass_hypothesis)/massResolution;
    _sigmaMOverM = massResolution/mass;
    _H_ptOverM    = Higgs.Pt()/mass_hypothesis;
    _cat        = cat;

}

void MvaAnalysis::SetBDTInputTree(TTree *tree){
    mvaFile_->cd();
    TBranch *b_log_H_pt             = tree->Branch("log_H_pt", &_log_H_pt , "log_H_pt/F");;
    TBranch *b_H_ptOverM            = tree->Branch("H_ptOverM", &_H_ptOverM , "_H_ptOverM/F");;
    TBranch *b_H_eta                = tree->Branch("H_eta", &_H_eta , "H_eta/F");;
    TBranch *b_d_phi                = tree->Branch("d_phi", &_d_phi,"d_phi/F");
    TBranch *b_max_eta              = tree->Branch("max_eta", &_max_eta , "max_eta/F");;
    TBranch *b_min_r9               = tree->Branch("min_r9", &_min_r9 , "min_r9/F");;
    TBranch *b_pho1_eta             = tree->Branch("pho1_eta", &_pho1_eta , "pho1_eta/F");
    TBranch *b_pho2_eta             = tree->Branch("pho2_eta", &_pho2_eta , "pho2_eta/F");
    TBranch *b_pho1_ptOverM         = tree->Branch("pho1_ptOverM", &_pho1_ptOverM , "pho1_ptOverM/F");;
    TBranch *b_pho2_ptOverM         = tree->Branch("pho2_ptOverM", &_pho2_ptOverM , "pho2_ptOverM/F");;
    TBranch *b_deltaMOverM          = tree->Branch("deltaMOverM", &_deltaMOverM,"deltaMOverM/F");
    TBranch *b_deltaMOverSigmaM     = tree->Branch("deltaMOverSigmaM", &_deltaMOverSigmaM,"deltaMOverSigmaM/F");
    TBranch *b_sigmaMOverM          = tree->Branch("sigmaMOverM", &_sigmaMOverM,"sigmaMOverM/F");
    TBranch *b_mgg                  = tree->Branch("mgg", &_mgg, "mgg/F");
    TBranch *b_pho1_phi             = tree->Branch("pho1_phi", &_pho1_phi , "pho1_phi/F");
    TBranch *b_pho1_pt              = tree->Branch("pho1_pt", &_pho1_pt , "pho1_pt/F");;
    TBranch *b_pho1_r9              = tree->Branch("pho1_r9", &_pho1_r9 , "_pho1_r9/F");;
    TBranch *b_pho2_phi             = tree->Branch("pho2_phi", &_pho2_phi , "pho2_phi/F");
    TBranch *b_pho2_pt              = tree->Branch("pho2_pt", &_pho2_pt , "pho2_pt/F");;
    TBranch *b_pho2_r9              = tree->Branch("pho2_r9", &_pho2_r9 , "_pho2_r9/F");;
    TBranch *b_H_pt                 = tree->Branch("H_pt", &_H_pt , "H_pt/F");;
    TBranch *b_Ht                   = tree->Branch("Ht", &_Ht , "Ht/F");;
    TBranch *b_d_eta                = tree->Branch("d_eta", &_d_eta,"d_eta/F");
    TBranch *b_mod_d_eta            = tree->Branch("mod_d_eta", &_mod_d_eta,"mod_d_eta/F");
    TBranch *b_cos_theta_star       = tree->Branch("cos_theta_star", &_cos_theta_star , "cos_theta_star/F");;
    TBranch *b_wt                   = tree->Branch("wt", &_wt, "wt/F");
    TBranch *b_cat                  = tree->Branch("cat", &_cat, "cat/I");
    TBranch *b_sideband             = tree->Branch("sideband", &_sideband, "sideband/I");
}

void MvaAnalysis::ResetAnalysis(){
    // Reset Random Variable on the EnergyResolution Smearer
    eResolSmearer->resetRandom();
}
