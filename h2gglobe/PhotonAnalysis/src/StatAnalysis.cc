#include "../interface/StatAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
StatAnalysis::StatAnalysis()  : 
    name_("StatAnalysis"),
    vtxAna_(vtxAlgoParams),
//    vtxConv_(vtxAlgoParams),
	useNVert(3)
{

    systRange  = 3.; // in units of sigma
    nSystSteps = 1;    

    nVtxCategories = 0;
}

// ----------------------------------------------------------------------------------------------------
StatAnalysis::~StatAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::Term(LoopAll& l) 
{

    std::string outputfilename = (std::string) l.histFileName;
    // Make Fits to the data-sets and systematic sets
    l.rooContainer->FitToData("data_pol_model","data_mass");  // Fit to full range of dataset
  
    l.rooContainer->WriteSpecificCategoryDataCards(outputfilename,"data_mass","sig_mass","data_pol_model");
    l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass","data_pol_model");
    // mode 0 as above, 1 if want to bin in sub range from fit,

    // Write the data-card for the Combinations Code, needs the output filename, makes binned analysis DataCard
    // Assumes the signal datasets will be called signal_name+"_mXXX"
//    l.rooContainer->GenerateBinnedPdf("bkg_mass_rebinned","data_pol_model","data_mass",1,50,1); // 1 means systematics from the fit effect only the backgroundi. last digit mode = 1 means this is an internal constraint fit 
//    l.rooContainer->WriteDataCard(outputfilename,"data_mass","sig_mass","bkg_mass_rebinned");

    eventListText.close();

    std::cout << " nevents " <<  nevents << " " << sumwei << std::endl;

//	kfacFile->Close();
//	PhotonAnalysis::Term(l);
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::Init(LoopAll& l) 
{
    if(PADEBUG) 
	cout << "InitRealStatAnalysis START"<<endl;

    nevents=0., sumwei=0.; 
    sumaccept=0., sumsmear=0., sumev=0.;
    
    std::string outputfilename = (std::string) l.histFileName;
    eventListText.open(Form("%s_ascii_events.txt",outputfilename.c_str()));
    //
    // These parameters are set in the configuration file
    std::cout
	<< "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< "StatAnalysis " << "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< "leadEtCut "<< leadEtCut << "\n"
	<< "subleadEtCut "<< subleadEtCut << "\n"
	<< "doTriggerSelection "<< doTriggerSelection << "\n"
	<< "nEtaCategories "<< nEtaCategories << "\n"
	<< "nR9Categories "<< nR9Categories << "\n"		
	<< "nPtCategories "<< nPtCategories << "\n"		
	<< "nVtxCategories "<< nVtxCategories << "\n"		
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

    // avoid recalculated the CIC ID every time
    l.runCiC = reRunCiC;
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
    if( nVtxCategories != 0 ) nCategories_ *= nVtxCategories;
    
    nPhotonCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nPhotonCategories_ *= nR9Categories;
    
    //// This is done in PhotonAnalysis now GF
    //eSmearPars.categoryType = "2CatR9_EBEE";
    //eSmearPars.byRun = false;
    //eSmearPars.n_categories = 4;
    //
    //// E scale is shifted for data, NOT for MC 
    //eSmearPars.scale_offset["EBHighR9"] = 0.;
    //eSmearPars.scale_offset["EBLowR9"]  = 0.;
    //eSmearPars.scale_offset["EEHighR9"] = 0.;
    //eSmearPars.scale_offset["EELowR9"]  = 0.;
    //// E scale systematics are applied to MC, NOT to data
    //eSmearPars.scale_offset_error["EBHighR9"] = scale_offset_error_EBHighR9;
    //eSmearPars.scale_offset_error["EBLowR9"]  = scale_offset_error_EBLowR9;
    //eSmearPars.scale_offset_error["EEHighR9"] = scale_offset_error_EEHighR9;
    //eSmearPars.scale_offset_error["EELowR9"]  = scale_offset_error_EELowR9;
    //// E resolution smearing applied to MC 
    //eSmearPars.smearing_sigma["EBHighR9"] = smearing_sigma_EBHighR9;
    //eSmearPars.smearing_sigma["EBLowR9"]  = smearing_sigma_EBLowR9;
    //eSmearPars.smearing_sigma["EEHighR9"] = smearing_sigma_EEHighR9;
    //eSmearPars.smearing_sigma["EELowR9"]  = smearing_sigma_EELowR9;
    //// E resolution systematics applied to MC 
    //eSmearPars.smearing_sigma_error["EBHighR9"] = smearing_sigma_error_EBHighR9;
    //eSmearPars.smearing_sigma_error["EBLowR9"]  = smearing_sigma_error_EBLowR9;
    //eSmearPars.smearing_sigma_error["EEHighR9"] = smearing_sigma_error_EEHighR9;
    //eSmearPars.smearing_sigma_error["EELowR9"]  = smearing_sigma_error_EELowR9;
    // MC would need Paul's corrections, of its own GF


    // This is done in PhotonAnalysis now
    //// eSmearDataPars.categoryType = "2CatR9_EBEE";
    //// eSmearDataPars.n_categories = 4;
    //// 
    //// // initialize smearer specific to energy shifts in DATA; use opposite of energy scale shift
    //// eSmearDataPars.scale_offset["EBHighR9"] = -1*scale_offset_EBHighR9;
    //// eSmearDataPars.scale_offset["EBLowR9"]  = -1*scale_offset_EBLowR9;
    //// eSmearDataPars.scale_offset["EEHighR9"] = -1*scale_offset_EEHighR9;
    //// eSmearDataPars.scale_offset["EELowR9"]  = -1*scale_offset_EELowR9;
    //// // no energy scale systematics applied to data
    //// eSmearDataPars.scale_offset_error["EBHighR9"] = 0.;
    //// eSmearDataPars.scale_offset_error["EBLowR9"]  = 0.;
    //// eSmearDataPars.scale_offset_error["EEHighR9"] = 0.;
    //// eSmearDataPars.scale_offset_error["EELowR9"]  = 0.;
    //// // E resolution smearing NOT applied to data 
    //// eSmearDataPars.smearing_sigma["EBHighR9"] = 0.;
    //// eSmearDataPars.smearing_sigma["EBLowR9"]  = 0.;
    //// eSmearDataPars.smearing_sigma["EEHighR9"] = 0.;
    //// eSmearDataPars.smearing_sigma["EELowR9"]  = 0.;
    //// // E resolution systematics NOT applied to data 
    //// eSmearDataPars.smearing_sigma_error["EBHighR9"] = 0.;
    //// eSmearDataPars.smearing_sigma_error["EBLowR9"]  = 0.;
    //// eSmearDataPars.smearing_sigma_error["EEHighR9"] = 0.;
    //// eSmearDataPars.smearing_sigma_error["EELowR9"]  = 0.;
    // DATA would need Paul's corrections, of its own (different from eSmearPars which is for MC) GF


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
        // Moved to PhotonAnalysis GF 
	//// energy scale systematics to MC
        //eScaleSmearer = new EnergySmearer( eSmearPars );
	//eScaleSmearer->name("E_scale");
	//eScaleSmearer->doEnergy(true);
	//eScaleSmearer->scaleOrSmear(true);
        photonSmearers_.push_back(eScaleSmearer);

	//// Moved to PhotonAnalysis PM
	//// // energy scale corrections to Data
	//// eScaleDataSmearer = new EnergySmearer( eSmearDataPars );
	//// eScaleDataSmearer->name("E_scale_data");
	//// eScaleDataSmearer->doEnergy(true);
	//// eScaleDataSmearer->scaleOrSmear(true);
	//photonDataSmearers_.push_back(eScaleDataSmearer); // must not be included among MC smearers; will be singled out upon need // GF questions?
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

    // Define the number of categories for the statistical analysis and
    // the systematic sets to be formed

    // FIXME move these params to config file
    l.rooContainer->SetNCategories(nCategories_);
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

    // Create observables for shape-analysis with ranges
    // l.rooContainer->AddObservable("mass" ,100.,150.);
    l.rooContainer->AddObservable("CMS_hgg_mass" ,massMin,massMax);

    l.rooContainer->AddConstant("IntLumi",l.intlumi_);

    // SM Model
    l.rooContainer->AddConstant("XSBR_150",0.01428+0.001307+0.000641+0.000066  );
    l.rooContainer->AddConstant("XSBR_145",0.018820+0.001676+0.000891+0.000090 );
    l.rooContainer->AddConstant("XSBR_140",0.0234109+0.00203036+0.001163597+0.000117189);
    l.rooContainer->AddConstant("XSBR_135",0.0278604+0.002343+0.001457559+0.000145053  );
    l.rooContainer->AddConstant("XSBR_130",0.0319112+0.00260804+0.001759636+0.000173070);
    l.rooContainer->AddConstant("XSBR_125",0.0350599+0.00277319+0.002035123+0.000197718);
    l.rooContainer->AddConstant("XSBR_120",0.0374175+0.00285525+0.002285775+0.00021951 );
    l.rooContainer->AddConstant("XSBR_123",0.0369736+0.00281352+0.00213681+0.00020663);
    l.rooContainer->AddConstant("XSBR_121",0.0369736+0.00284082+0.00223491+0.00021510);
    l.rooContainer->AddConstant("XSBR_115",0.0386169+0.00283716+0.002482089+0.000235578);
    l.rooContainer->AddConstant("XSBR_110",0.0390848+0.00275406+0.002654575+0.000247629);
    l.rooContainer->AddConstant("XSBR_105",0.0387684+0.00262016+0.002781962+0.000255074);

    // FF model	
    l.rooContainer->AddConstant("ff_XSBR_150",0.00259659+0.00127278);
    l.rooContainer->AddConstant("ff_XSBR_145",0.00387544+0.00205969);
    l.rooContainer->AddConstant("ff_XSBR_140",0.00565976+0.003243602);
    l.rooContainer->AddConstant("ff_XSBR_135",0.00825+0.00513225   );
    l.rooContainer->AddConstant("ff_XSBR_130",0.0122324+0.00825316 );
    l.rooContainer->AddConstant("ff_XSBR_125",0.0186494+0.01368598 );
    l.rooContainer->AddConstant("ff_XSBR_123",0.022212+0.0168696  );
    l.rooContainer->AddConstant("ff_XSBR_121",0.0266484+0.0209646 );
    l.rooContainer->AddConstant("ff_XSBR_120",0.0293139+0.02346729 );
    l.rooContainer->AddConstant("ff_XSBR_115",0.0482184+0.04218386 );
    l.rooContainer->AddConstant("ff_XSBR_110",0.083181+0.08017625  );
    l.rooContainer->AddConstant("ff_XSBR_105",0.151616+0.1609787   );

    // Background modeling 
    l.rooContainer->AddRealVar("pol0",-0.01,-1.5,1.5);
    l.rooContainer->AddRealVar("pol1",-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("modpol0","@0*@0","pol0");
    l.rooContainer->AddFormulaVar("modpol1","@0*@0","pol1");

    std::vector<std::string> data_pol_pars(2,"p");	 
    data_pol_pars[0] = "modpol0";
    data_pol_pars[1] = "modpol1";
    l.rooContainer->AddGenericPdf("data_pol_model",
	  "0","CMS_hgg_mass",data_pol_pars,72);	// >= 71 means RooBernstein of order >= 1
        
    // -----------------------------------------------------
    // Make some data sets from the observables to fill in the event loop		  
    // Binning is for histograms (will also produce unbinned data sets)
    l.rooContainer->CreateDataSet("CMS_hgg_mass","data_mass"    ,nDataBins); // (100,110,150) -> for a window, else full obs range is taken 
    l.rooContainer->CreateDataSet("CMS_hgg_mass","bkg_mass"     ,nDataBins);    	  	
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_m105",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_m110",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_m115",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_m120",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_m121",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_m123",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_m125",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_m130",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_m135",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_m140",nDataBins);    
//    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_m145",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_m150",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_rv_m105",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_rv_m110",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_rv_m115",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_rv_m120",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_rv_m121",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_rv_m123",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_rv_m125",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_rv_m130",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_rv_m135",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_rv_m140",nDataBins);    
//    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_rv_m145",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_rv_m150",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_wv_m105",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_wv_m110",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_wv_m115",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_wv_m120",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_wv_m121",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_wv_m123",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_wv_m125",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_wv_m130",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_wv_m135",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_wv_m140",nDataBins);    
//    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_wv_m145",nDataBins);    
    l.rooContainer->CreateDataSet("CMS_hgg_mass","sig_mass_wv_m150",nDataBins);    

    // Make more data sets to represent systematic shitfs , 
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_mass_m105",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_mass_m110",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_mass_m115",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_mass_m120",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_mass_m121",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_mass_m123",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_mass_m125",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_mass_m130",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_mass_m135",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_mass_m140",-1);	
//    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_mass_m145",-1);	
    l.rooContainer->MakeSystematics("CMS_hgg_mass","sig_mass_m150",-1);	
	
    /* -----------------------------------------------------------------------------------------
       KFactors Reweighting
       ------------------------------------------------------------------------------------------- */
    if(PADEBUG) 
	cout << "InitRealStatAnalysis END"<<endl;
	

    //vertex-related TMVAs
    //// tmvaPerEvtReader_ = new TMVA::Reader( "!Color:!Silent" );
    //// MVA_.resize(useNVert);
    //// dZ_.resize(useNVert);
    //// tmvaPerEvtReader_->AddVariable( "diphoPt0", &diphoRelPt_ );
    //// tmvaPerEvtReader_->AddVariable( "nVert"	 , &nVert_ 		);
    //// tmvaPerEvtReader_->AddVariable( "MVA0" 	 , &MVA_[0]		);
    //// tmvaPerEvtReader_->AddVariable( "MVA1"    , &MVA_[1]		);
    //// tmvaPerEvtReader_->AddVariable( "dZ1"     , &dZ_[1]		);
    //// tmvaPerEvtReader_->AddVariable( "MVA2"    , &MVA_[2]		);
    //// tmvaPerEvtReader_->AddVariable( "dZ2"     , &dZ_[2]		);
    //// /// tmvaPerEvtReader_->AddVariable( "nConv"     , &nConv_		);
    //// tmvaPerEvtReader_->BookMVA( tmvaPerEvtMethod, tmvaPerEvtWeights );

    // FIXME book of additional variables
}

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
    if(PADEBUG) 
	cout << "Analysis START; cur_type is: " << l.itype[l.current] <<endl;
   
    int cur_type = l.itype[l.current];
    float weight = l.sampleContainer[l.current_sample_index].weight;
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
    if(cur_type!=0){
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
	} else if( cur_type == 0 ) {          // if it's data
	    float sweight = 1.;
	    if( doEcorrectionSmear )  { 
	      eCorrSmearer->smearPhoton(phoInfo,sweight,l.run,0.); 
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

        // MicroAnalysis-specific stuff
        // TODO dipho selection should match the Higgs truth (not just the best dipho)
        vtxAna_.setPairID(diphoton_id);
        nVert_ = l.vtx_std_n;

        //put stupid values in vectors
        MVA_.assign(useNVert,-10);
        dZ_.assign(useNVert,-100);

        //get list of vertices as ranked for the selected diphoton
        vector<int> & rankedVtxs = (*l.vtx_std_ranked_list)[diphoton_id];
        /// vtxAna_.preselection(rankedVtxs);
	/// if( mvaVertexSelection ) {
	/// 	rankedVtxs = vtxAna_.rank(*tmvaPerVtxReader_,tmvaPerVtxMethod);
	/// } else {
	/// 	vtxAna_.evaluate(*tmvaPerVtxReader_,tmvaPerVtxMethod);
	/// }
        for (size_t vi=0;vi<rankedVtxs.size();vi++) {
        	if(vi>=useNVert) break;
        	MVA_[vi] = vtxAna_.mva(rankedVtxs[vi]);
        	dZ_[vi] = ( *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[vi]) - *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[0])).Z();

        	TLorentzVector lead_pho = l.get_pho_p4( l.dipho_leadind[diphoton_id], rankedVtxs[vi], &smeared_pho_energy[0]);
        	TLorentzVector sublead_pho = l.get_pho_p4( l.dipho_subleadind[diphoton_id], rankedVtxs[vi], &smeared_pho_energy[0]);
        	TLorentzVector dipho = lead_pho+sublead_pho;
        	/// if(vi==0) diphoRelPt_ = dipho.Pt()/dipho.M();
        	if(vi==0) diphoRelPt_ = dipho.Pt();
        }
        // VtxEvtMVA_ = tmvaPerEvtReader_->EvaluateMVA(tmvaPerEvtMethod);
	/// VtxEvtMVA_ = vtxAna_.perEventMva( *tmvaPerEvtReader_, tmvaPerEvtMethod, rankedVtxs );
	VtxProb_ = vtxAna_.vertexProbability( VtxEvtMVA_ ); 
	
        //TODO microanalysis imported stuff
	l.countersred[diPhoCounter_]++;

	TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
	float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
	TLorentzVector Higgs = lead_p4 + sublead_p4; 	
	TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);

	bool CorrectVertex;
	int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories,nVtxCategories,VtxEvtMVA_);
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
	
	// Monitor vertex per-event MVA inputs
        for (size_t vi=0;vi<rankedVtxs.size();vi++) {
        	if(vi>=useNVert) break;
		
		l.FillHist(Form("vtx_mva_%d",vi),0,MVA_[vi],evweight);
		l.FillHist(Form("vtx_mva_%d",vi),category+1,MVA_[vi],evweight);
		if( vi>0 ) { 
			l.FillHist(Form("vtx_dz_%d",vi),0,dZ_[vi],evweight); 
			l.FillHist(Form("vtx_dz_%d",vi),category+1,dZ_[vi],evweight); 
		}
        }
	l.FillHist("vtx_n",0,nVert_,evweight);
	l.FillHist("vtx_n",category+1,nVert_,evweight);
	l.FillHist("vtx_evt_mva",0,VtxEvtMVA_,evweight);
	l.FillHist("vtx_evt_mva",category+1,VtxEvtMVA_,evweight);
	l.FillHist2D("vtx_prob_vs_pt",0,ptHiggs,VtxProb_,evweight);
	l.FillHist2D("vtx_prob_vtx_n",category+1,nVert_,VtxProb_,evweight);

	// Monitor correlations 
	std::vector<std::pair<std::string,float> > mva_input_corr;
	mva_input_corr.push_back( std::make_pair("vtx_n",nVert_) );
	mva_input_corr.push_back( std::make_pair("dipho_pt",ptHiggs) );
	for (size_t vi=0;vi<rankedVtxs.size();vi++) {
		mva_input_corr.push_back( std::make_pair(Form("vtx_mva_%d",vi),MVA_[vi]) );
	}
	for (size_t vi=1;vi<rankedVtxs.size();vi++) {
		mva_input_corr.push_back( std::make_pair(Form("vtx_dz_%d",vi),dZ_[vi]) );
	}
	for( size_t ivar=0; ivar<mva_input_corr.size(); ++ivar ) {
		for( size_t jvar=ivar+1; jvar<mva_input_corr.size(); ++jvar ) {
			std::string name = mva_input_corr[ivar].first+"_vs_"+mva_input_corr[jvar].first;
			float & ival = mva_input_corr[ivar].second;
			float & jval = mva_input_corr[jvar].second;
			l.FillHist2D( name, 0,          ival, jval, evweight );
			l.FillHist2D( name, category+1, ival, jval, evweight );
		}
	}
	
	l.FillCounter( "Accepted", weight );
	l.FillCounter( "Smeared", evweight );
	sumaccept += weight;
 	sumsmear += evweight;

	// control plots 
	l.FillHist("all_mass",0, Higgs.M(), evweight);
	l.FillHist("all_mass",category+1, Higgs.M(), evweight);
	if( mass>=massMin && mass<=massMax  ) {
		l.FillHist("mass",0, Higgs.M(), evweight);
		l.FillHist("pt",0, Higgs.Pt(), evweight);
		l.FillHist("eta",0, Higgs.Eta(), evweight);
		
		l.FillHist("pho_pt",0,lead_p4.Pt(), evweight);
		l.FillHist("pho1_pt",0,lead_p4.Pt(), evweight);
		l.FillHist("pho_eta",0,lead_p4.Eta(), evweight);
		l.FillHist("pho1_eta",0,lead_p4.Eta(), evweight);
		l.FillHist("pho_r9",0, lead_r9, evweight);
		l.FillHist("pho1_r9",0, lead_r9, evweight);
		
		l.FillHist("pho_pt",0,sublead_p4.Pt(), evweight);
		l.FillHist("pho2_pt",0,sublead_p4.Pt(), evweight);
		l.FillHist("pho_eta",0,sublead_p4.Eta(), evweight);
		l.FillHist("pho2_eta",0,sublead_p4.Eta(), evweight);
		l.FillHist("pho_r9",0, sublead_r9, evweight);
		l.FillHist("pho1_r9",0, sublead_r9, evweight);
		
		l.FillHist("mass",category+1, Higgs.M(), evweight);
		l.FillHist("pt",category+1, Higgs.Pt(), evweight);
		l.FillHist("eta",category+1, Higgs.Eta(), evweight);
		
		l.FillHist("pho_pt",category+1,lead_p4.Pt(), evweight);
		l.FillHist("pho1_pt",category+1,lead_p4.Pt(), evweight);
		l.FillHist("pho_eta",category+1,lead_p4.Eta(), evweight);
		l.FillHist("pho1_eta",category+1,lead_p4.Eta(), evweight);
		l.FillHist("pho_r9",category+1, lead_r9, evweight);
		l.FillHist("pho1_r9",category+1, lead_r9, evweight);
		
		l.FillHist("pho_pt",category+1,sublead_p4.Pt(), evweight);
		l.FillHist("pho2_pt",category+1,sublead_p4.Pt(), evweight);
		l.FillHist("pho_eta",category+1,sublead_p4.Eta(), evweight);
		l.FillHist("pho2_eta",category+1,sublead_p4.Eta(), evweight);
		l.FillHist("pho_r9",category+1, sublead_r9, evweight);
		l.FillHist("pho1_r9",category+1, sublead_r9, evweight);
		
		l.FillHist("pho_n",category+1,l.pho_n, evweight);
	}

	if (cur_type==0){
	  eventListText << setprecision(4) <<"Type = "<< cur_type <<  "Run = " << l.run << "  LS = " << l.lumis << "  Event = " << l.event << "  SelVtx = " << l.vtx_std_sel << "  CAT4 = " << category % 4 << "  ggM = " << mass << " gg_Pt =  " << ptHiggs;
	  eventListText << endl;
	}
       
	// --------------------------------------------------------------------------------------------- 
	if (cur_type == 0 ){
	    l.rooContainer->InputDataPoint("data_mass",category,mass);
	}
	if (cur_type > 0 && cur_type != 3 && cur_type != 4)
	    l.rooContainer->InputDataPoint("bkg_mass",category,mass,evweight);
	else if (cur_type == -13|| cur_type == -14 || cur_type == -15|| cur_type == -16){
	    l.rooContainer->InputDataPoint("sig_mass_m105",category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m105",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_mass_wv_m105",category,mass,evweight);
	}
	else if (cur_type == -17 || cur_type == -18 || cur_type == -19|| cur_type == -20){
	    l.rooContainer->InputDataPoint("sig_mass_m110",category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m110",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_mass_wv_m110",category,mass,evweight);
	}
	else if (cur_type == -21 || cur_type == -22 || cur_type == -23|| cur_type == -24){
	    l.rooContainer->InputDataPoint("sig_mass_m115",category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m115",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_mass_wv_m115",category,mass,evweight);
	}
	else if (cur_type == -25 || cur_type == -26 || cur_type == -27|| cur_type == -28){
	    l.rooContainer->InputDataPoint("sig_mass_m120",category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m120",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_mass_wv_m120",category,mass,evweight);
	}
	else if (cur_type == -29 || cur_type == -30 || cur_type == -31|| cur_type == -32){
	    l.rooContainer->InputDataPoint("sig_mass_m130",category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m130",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_mass_wv_m130",category,mass,evweight);
	}
	else if (cur_type == -33 || cur_type == -34 || cur_type == -35|| cur_type == -36){
	    l.rooContainer->InputDataPoint("sig_mass_m140",category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m140",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_mass_wv_m140",category,mass,evweight);
	}
	else if (cur_type == -37 || cur_type == -38 || cur_type == -39|| cur_type == -40){
	    l.rooContainer->InputDataPoint("sig_mass_m125",category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m125",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_mass_wv_m125",category,mass,evweight);
	}
	else if (cur_type == -41 || cur_type == -42 || cur_type == -43|| cur_type == -44){
	    l.rooContainer->InputDataPoint("sig_mass_m135",category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m135",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_mass_wv_m135",category,mass,evweight);
	}


//	else if (cur_type == -45 || cur_type == -46 || cur_type == -47|| cur_type == -48){
//	    l.rooContainer->InputDataPoint("sig_mass_m145",category,mass,evweight);
//	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m145",category,mass,evweight);
//	    else l.rooContainer->InputDataPoint("sig_mass_wv_m145",category,mass,evweight);
//	}
	else if (cur_type == -49 || cur_type == -50 || cur_type == -51|| cur_type == -52){
	    l.rooContainer->InputDataPoint("sig_mass_m150",category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m150",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_mass_wv_m150",category,mass,evweight);
	}
	else if (cur_type == -53 || cur_type == -54 || cur_type == -55|| cur_type == -56){
	    l.rooContainer->InputDataPoint("sig_mass_m121",category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m121",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_mass_wv_m121",category,mass,evweight);
	}
	else if (cur_type == -57 || cur_type == -58 || cur_type == -59|| cur_type == -60){
	    l.rooContainer->InputDataPoint("sig_mass_m123",category,mass,evweight);
	    if (CorrectVertex) l.rooContainer->InputDataPoint("sig_mass_rv_m123",category,mass,evweight);
	    else l.rooContainer->InputDataPoint("sig_mass_wv_m123",category,mass,evweight);
	}
       
    }
   
   
    // Systematics
    if( cur_type != 0 && doMCSmearing ) { 
	// fill steps for syst uncertainty study
	float systStep = systRange / (float)nSystSteps;
	// di-photon smearers systematics
	if (diphoton_id > -1 ) {
	       
	    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	    TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
	 
	    for(std::vector<BaseGenLevelSmearer*>::iterator si=systGenLevelSmearers_.begin(); si!=systGenLevelSmearers_.end(); si++){
		std::vector<double> mass_errors;
		std::vector<double> weights;
		std::vector<int>    categories;
	   
		for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		    if( syst_shift == 0. ) { continue; } // skip the central value
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
	     
		    int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories,nVtxCategories,VtxEvtMVA_);
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
		    categories.push_back(category);
		    mass_errors.push_back(mass);
		    weights.push_back(evweight);
		}// end loop on systematics steps
	   
		if (cur_type == -13|| cur_type == -14 || cur_type == -15|| cur_type == -16)
		    l.rooContainer->InputSystematicSet("sig_mass_m105",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -17|| cur_type == -18 || cur_type == -19|| cur_type == -20)
		    l.rooContainer->InputSystematicSet("sig_mass_m110",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -21|| cur_type == -22 || cur_type == -23|| cur_type == -24)
		    l.rooContainer->InputSystematicSet("sig_mass_m115",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -25|| cur_type == -26 || cur_type == -27|| cur_type == -28)
		    l.rooContainer->InputSystematicSet("sig_mass_m120",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -29|| cur_type == -30 || cur_type == -31|| cur_type == -32)
		    l.rooContainer->InputSystematicSet("sig_mass_m130",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -33|| cur_type == -34 || cur_type == -35|| cur_type == -36)
		    l.rooContainer->InputSystematicSet("sig_mass_m140",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -37|| cur_type == -38 || cur_type == -39|| cur_type == -40)
		    l.rooContainer->InputSystematicSet("sig_mass_m125",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -41|| cur_type == -42 || cur_type == -43|| cur_type == -44)
		    l.rooContainer->InputSystematicSet("sig_mass_m135",(*si)->name(),categories,mass_errors,weights);

	//	else if (cur_type == -45|| cur_type == -46 || cur_type == -47|| cur_type == -48)
	//	    l.rooContainer->InputSystematicSet("sig_mass_m145",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -49|| cur_type == -50 || cur_type == -51|| cur_type == -52)
		    l.rooContainer->InputSystematicSet("sig_mass_m150",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -53|| cur_type == -54 || cur_type == -55|| cur_type == -56)
		    l.rooContainer->InputSystematicSet("sig_mass_m121",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -57|| cur_type == -58 || cur_type == -59|| cur_type == -60)
		    l.rooContainer->InputSystematicSet("sig_mass_m123",(*si)->name(),categories,mass_errors,weights);
	    }// end loop on smearers 
		 

	    for(std::vector<BaseDiPhotonSmearer *>::iterator si=systDiPhotonSmearers_.begin(); si!= systDiPhotonSmearers_.end(); ++si ) {
		std::vector<double> mass_errors;
		std::vector<double> weights;
		std::vector<int> categories;
		       
		for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		    if( syst_shift == 0. ) { continue; } // skip the central value
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
			       
		    // restart with 'fresh' wait for this round of systematics
		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
			       
		    // FIXME pass smeared R9 and di-photon
		    int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories,nVtxCategories,VtxEvtMVA_);
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
		    categories.push_back(category);
		    mass_errors.push_back(mass);
		    weights.push_back(evweight);
		}
		if (cur_type == -13|| cur_type == -14 || cur_type == -15|| cur_type == -16)
		    l.rooContainer->InputSystematicSet("sig_mass_m105",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -17|| cur_type == -18 || cur_type == -19|| cur_type == -20)
		    l.rooContainer->InputSystematicSet("sig_mass_m110",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -21|| cur_type == -22 || cur_type == -23|| cur_type == -24)
		    l.rooContainer->InputSystematicSet("sig_mass_m115",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -25|| cur_type == -26 || cur_type == -27|| cur_type == -28)
		    l.rooContainer->InputSystematicSet("sig_mass_m120",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -29|| cur_type == -30 || cur_type == -31|| cur_type == -32)
		    l.rooContainer->InputSystematicSet("sig_mass_m130",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -33|| cur_type == -34 || cur_type == -35|| cur_type == -36)
		    l.rooContainer->InputSystematicSet("sig_mass_m140",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -37|| cur_type == -38 || cur_type == -39|| cur_type == -40)
		    l.rooContainer->InputSystematicSet("sig_mass_m125",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -41|| cur_type == -42 || cur_type == -43|| cur_type == -44)
		    l.rooContainer->InputSystematicSet("sig_mass_m135",(*si)->name(),categories,mass_errors,weights);
	//	else if (cur_type == -45|| cur_type == -46 || cur_type == -47|| cur_type == -48)
	//	    l.rooContainer->InputSystematicSet("sig_mass_m145",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -49|| cur_type == -50 || cur_type == -51|| cur_type == -52)
		    l.rooContainer->InputSystematicSet("sig_mass_m150",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -53|| cur_type == -54 || cur_type == -55|| cur_type == -56)
		    l.rooContainer->InputSystematicSet("sig_mass_m121",(*si)->name(),categories,mass_errors,weights);
		else if (cur_type == -57|| cur_type == -58 || cur_type == -59|| cur_type == -60)
		    l.rooContainer->InputSystematicSet("sig_mass_m123",(*si)->name(),categories,mass_errors,weights);
	    }

	}
       
	// loop over the smearers included in the systematics study
	for(std::vector<BaseSmearer *>::iterator  si=systPhotonSmearers_.begin(); si!= systPhotonSmearers_.end(); ++si ) {
	    std::vector<double> mass_errors;
	    std::vector<double> weights;
	    std::vector<int> categories;
	   
	    // loop over syst shift
	    for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		if( syst_shift == 0. ) { continue; } // skip the central value
		// smear the photons 
		for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
		    std::vector<std::vector<bool> > p;
		    //std::cout << "GF check: " <<  l.pho_residCorrEnergy[ipho] << "  " << l.pho_residCorrResn[ipho] << std::endl;
		    PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), 
						/// *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
						((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), 
						energyCorrected[ipho],
						l.pho_isEB[ipho], l.pho_r9[ipho],
						l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_));
		   
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
		   
		    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
		    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
		    TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
		   
		    int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories,nVtxCategories,VtxEvtMVA_);
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
		   
	       	    categories.push_back(category);
	            mass_errors.push_back(mass);
	            weights.push_back(evweight);

		} else {
		    mass_errors.push_back(0.);   
		    weights.push_back(0.);   
		    categories.push_back(-1);
		}
	       
	    }
	    if (cur_type == -13|| cur_type == -14 || cur_type == -15|| cur_type == -16)
		l.rooContainer->InputSystematicSet("sig_mass_m105",(*si)->name(),categories,mass_errors,weights);
	    else if (cur_type == -17|| cur_type == -18 || cur_type == -19|| cur_type == -20)
		l.rooContainer->InputSystematicSet("sig_mass_m110",(*si)->name(),categories,mass_errors,weights);
	    else if (cur_type == -21|| cur_type == -22 || cur_type == -23|| cur_type == -24)
		l.rooContainer->InputSystematicSet("sig_mass_m115",(*si)->name(),categories,mass_errors,weights);
	    else if (cur_type == -25|| cur_type == -26 || cur_type == -27|| cur_type == -28)
		l.rooContainer->InputSystematicSet("sig_mass_m120",(*si)->name(),categories,mass_errors,weights);
	    else if (cur_type == -29|| cur_type == -30 || cur_type == -31|| cur_type == -32)
		l.rooContainer->InputSystematicSet("sig_mass_m130",(*si)->name(),categories,mass_errors,weights);
	    else if (cur_type == -33|| cur_type == -34 || cur_type == -35|| cur_type == -36)
		l.rooContainer->InputSystematicSet("sig_mass_m140",(*si)->name(),categories,mass_errors,weights);
	    else if (cur_type == -37|| cur_type == -38 || cur_type == -39|| cur_type == -40)
	        l.rooContainer->InputSystematicSet("sig_mass_m125",(*si)->name(),categories,mass_errors,weights);
	    else if (cur_type == -41|| cur_type == -42 || cur_type == -43|| cur_type == -44)
	        l.rooContainer->InputSystematicSet("sig_mass_m135",(*si)->name(),categories,mass_errors,weights);
	    //else if (cur_type == -45|| cur_type == -46 || cur_type == -47|| cur_type == -48)
	    //    l.rooContainer->InputSystematicSet("sig_mass_m145",(*si)->name(),categories,mass_errors,weights);
	    else if (cur_type == -49|| cur_type == -50 || cur_type == -51|| cur_type == -52)
	        l.rooContainer->InputSystematicSet("sig_mass_m150",(*si)->name(),categories,mass_errors,weights);
	    else if (cur_type == -53|| cur_type == -54 || cur_type == -55|| cur_type == -56)
	        l.rooContainer->InputSystematicSet("sig_mass_m121",(*si)->name(),categories,mass_errors,weights);
	    else if (cur_type == -57|| cur_type == -58 || cur_type == -59|| cur_type == -60)
	        l.rooContainer->InputSystematicSet("sig_mass_m123",(*si)->name(),categories,mass_errors,weights);
       
	}
       
       
    }
   
    if(PADEBUG) 
	cout<<"myFillHistRed END"<<endl;
}

// ----------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------
void StatAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
    vtxAna_.setBranchAdresses(t,"vtx_std_");
    vtxAna_.getBranches(t,"vtx_std_",s);
}

// ----------------------------------------------------------------------------------------------------
bool StatAnalysis::SelectEvents(LoopAll& l, int jentry) 
{
    return true;
}
// ----------------------------------------------------------------------------------------------------
double StatAnalysis::GetDifferentialKfactor(double gPT, int Mass)
{

/*  
    if (Mass <=110 ) return thm110->GetBinContent(thm110->FindFixBin(gPT));
    else if (Mass ==120 ) return thm120->GetBinContent(thm120->FindFixBin(gPT));
    else if (Mass ==130 ) return thm130->GetBinContent(thm130->FindFixBin(gPT));
    else if (Mass ==140 ) return thm140->GetBinContent(thm140->FindFixBin(gPT));
    else if (Mass ==115 ) return (0.5*thm110->GetBinContent(thm110->FindFixBin(gPT)) +0.5*thm120->GetBinContent(thm120->FindFixBin(gPT)));
*/
    return 1.0;
/*
  int  genMasses[4] = {110,120,130,140};
  if (Mass<=genMasses[0] ) return kfactorHistograms[0]->GetBinContent(kfactorHistograms[0]->FindBin(gPT));
  else if (Mass<genMasses[nMasses-1]) {

  TH1D *hm1,*hm2;
  double m1=0,m2=0;
  for (int m=0;m<nMasses;m++){
  if (Mass<genMasses[m+1]){
  hm1=kfactorHistograms[m];
  hm2=kfactorHistograms[m+1];
  m1 = genMasses[m];
  m2 = genMasses[m+1];
  //	cout << "Gen Mass: "<< Mass << " Using "<<m1<< " " << m2<< " Hist name check " << hm1->GetName()<<" " <<hm2->GetName()<<endl;
  break;
  }
  }
  if ((int)Mass == (int)m1 ){
  //cout << "Found the appropriate historgam "<<hm1->GetName()<<endl;
  return hm1->GetBinContent(hm1->FindBin(gPT));
  } else {

  TH1D *hm = (TH1D*) hm1->Clone("hm");
  double alpha = ((float) (Mass-m1))/(m2-m1); // make sure ms are not integers
  hm->Add(hm1,hm2,alpha,(1-alpha));
  return hm->GetBinContent(hm->GetBinContent(hm->FindBin(gPT)));
  }

  }
  else return kfactorHistograms[nMasses-1]->GetBinContent(kfactorHistograms[nMasses-1]->FindBin(gPT));
*/
}


