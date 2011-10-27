#include "../interface/MicroAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
MicroAnalysis::MicroAnalysis()  : 
    uFileName("uTree.root"),
	tmvaMethod("Likelihood"),
	tmvaWeights(""),
	storeNVert(3),
	name_("MicroAnalysis"),
    pho1_(0), pho2_(0), dipho_(0), vtxVars_(100)
{
	cout << "Constructing MicroAnalysis" << endl;
	useMvaRanking = false;
}

// ----------------------------------------------------------------------------------------------------
MicroAnalysis::~MicroAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void MicroAnalysis::Term(LoopAll& l) 
{
	if( rooFile_ ) {
		rooFile_->Close();
	} 
	cout<<endl;
	cout << "Closing-up MicroAnalysis" << endl;
	//TODO close up
	uFile_->cd();
	uTree_->Write(0,TObject::kWriteDelete);
	evTree_->Write(0,TObject::kWriteDelete);
	uFile_->Close();

	dtEBEB_.erase(dtEBEB_.begin(),dtEBEB_.end());
	dtEBEE_.erase(dtEBEE_.begin(),dtEBEE_.end());
	dtEEEE_.erase(dtEEEE_.begin(),dtEEEE_.end());
	rooFileTime_->Close(); 
}

// ----------------------------------------------------------------------------------------------------
void MicroAnalysis::Init(LoopAll& l) 
{
	StatAnalysis::Init(l);

	cout << "Initializing MicroAnalysis" << endl;
	if(PADEBUG) cout << "InitRealMicroAnalysis START"<<endl;
	uFile_ = TFile::Open(uFileName,"recreate");
	uFile_->cd(); 

	// per vertex tree
	uTree_ = new TTree("utree","MicroAnalysis per Vertex Tree");
	//LINK http://root.cern.ch/root/html/TTree.html
	uTree_->Branch("dZToGen",&dZToGen_);
	/// uTree_->Branch("dZToClosest",&dZtoClosest_);
	uTree_->Branch("nVert",&nVert_);
	uTree_->Branch("nPU",&nPU_);
	uTree_->Branch("evWeight",&evWeight_);
	uTree_->Branch("pho1",&pho1_,32000,0);
	uTree_->Branch("pho2",&pho2_,32000,0);
	uTree_->Branch("dipho",&dipho_,32000,0);
	/// uTree_->Branch("ptasym",&ptasym_);
	/// uTree_->Branch("ptbal",&ptbal_);
	/// uTree_->Branch("logsumpt2",&logsumpt2_);
	uTree_->Branch("absDeltaEta",&absDeltaEta_);
	uTree_->Branch("tofCorrTdiff",&tofCorrTdiff_);
	uTree_->Branch("isClosestToGen",&isClosestToGen_);
	uTree_->Branch("nConv",&nConv_);
	/// uTree_->Branch("dZToConv",&dZToConv_);
	uTree_->Branch("pullToConv",&pullToConv_);
	/// uTree_->Branch("convCompat",&convCompat_);
	uTree_->Branch("rprod",&rprod_);
	uTree_->Branch("convRprod",&convRprod_);

	vtxVarNames_.push_back("ptvtx"), vtxVarNames_.push_back("ptasym"), vtxVarNames_.push_back("ptratio"), vtxVarNames_.push_back("ptbal"), vtxVarNames_.push_back("logsumpt2"), vtxVarNames_.push_back("ptmax3"), vtxVarNames_.push_back("ptmax"), vtxVarNames_.push_back("nchthr"), vtxVarNames_.push_back("sumtwd"); 
	assert( vtxVars_.size() >= vtxVarNames_.size());
	for( size_t iv=0; iv<vtxVarNames_.size(); ++iv ){
		uTree_->Branch(vtxVarNames_[iv].c_str(),&vtxVars_[iv]);
	}
	uTree_->Print();
	//todo add delta eta


	//TMVA
	tmvaVariables_.push_back("ptbal"), tmvaVariables_.push_back("ptasym"), tmvaVariables_.push_back("logsumpt2");
	rankVariables_ = tmvaVariables_;
	if( addConversionToMva ) {
		tmvaVariables_.push_back("limPullToConv");
		tmvaVariables_.push_back("nConv");
	}
	tmvaReader_ = new TMVA::Reader( "!Color:!Silent" );
	HggVertexAnalyzer::bookVariables( *tmvaReader_, tmvaVariables_ );
	tmvaReader_->BookMVA( tmvaMethod, tmvaWeights );

	MVA_.resize(storeNVert);
	rankprod_.resize(storeNVert);
	convRankprod_.resize(storeNVert);
	alpha_.resize(storeNVert);
	prob_.resize(storeNVert);
	corrProb_.resize(storeNVert);
	wrongProb_.resize(storeNVert);
	zRMSn_.resize(storeNVert);
	dZ_.resize(storeNVert);
	diphoM_.resize(storeNVert);
	diphoCosTheta_.resize(storeNVert);
	diphoCosDeltaPhi_.resize(storeNVert);
	diphoPt_.resize(storeNVert);
	dZToConvV_.resize(storeNVert);
	pullToConvV_.resize(storeNVert);

	evTree_ = new TTree("evtree","MicroAnalysis per Event Tree");
	evTree_->Branch("dZTrue",&dZTrue_);
	evTree_->Branch("alphaTrue",&alphaTrue_);
	evTree_->Branch("zTrue",&zTrue_);
	evTree_->Branch("rTrue",&rTrue_);
	evTree_->Branch("zRMS",&zRMS_);
	evTree_->Branch("category",&category_,"category/I");
	evTree_->Branch("mTrue",&mTrue_);
	evTree_->Branch("mTrueVtx",&mTrueVtx_);
	evTree_->Branch("nVert",&nVert_);
	evTree_->Branch("evWeight",&evWeight_);
	evTree_->Branch("nConv",&nConv_);
	evTree_->Branch("ConvCompat",&convCompat_);
	for(int i=0;i<storeNVert; i++){
		evTree_->Branch(Form("MVA%d",i),&MVA_[i]);
		evTree_->Branch(Form("rankprod%d",i),&rankprod_[i]);
		evTree_->Branch(Form("convRankprod%d",i),&convRankprod_[i]);
		evTree_->Branch(Form("alpha%d",i),&alpha_[i]);
		/// evTree_->Branch(Form("prob%d",i),&prob_[i]);
		/// evTree_->Branch(Form("corrProb%d",i),&corrProb_[i]);
		/// evTree_->Branch(Form("wrongProb%d",i),&wrongProb_[i]);
		evTree_->Branch(Form("zRMSn%d",i),&zRMSn_[i]);
		evTree_->Branch(Form("dZ%d",i),&dZ_[i]);
		evTree_->Branch(Form("diphoM%d",i),&diphoM_[i]);
		evTree_->Branch(Form("diphoCosTheta%d",i),&diphoCosTheta_[i]);
		evTree_->Branch(Form("diphoCosDeltaPhi%d",i),&diphoCosDeltaPhi_[i]);
		evTree_->Branch(Form("diphoPt%d",i),&diphoPt_[i]);
		evTree_->Branch(Form("dZToConv%d",i),&dZToConvV_[i]);
		evTree_->Branch(Form("pullToConv%d",i),&pullToConvV_[i]);
	}
	
	vtxProb_ = 0;
	rooFile_     = TFile::Open(roofitProbs);
	if( rooFile_ ) {
		RooWorkspace * w = (RooWorkspace*)rooFile_->Get("w");
		assert(w);
		RooRealVar * vtxLike = w->var("vtxLike");

		RooAbsPdf * correctVtxUnbinnedPdf = w->pdf("correctVtxUnbinnedPdf");
		RooAbsPdf * wrongVtxUnbinnedPdf = w->pdf("wrongVtxUnbinnedPdf");

		RooAbsReal * correctVtxCdf = correctVtxUnbinnedPdf->createCdf(RooArgSet(*vtxLike));
		RooAbsReal *  wrongVtxCdf   = wrongVtxUnbinnedPdf->createCdf(RooArgSet(*vtxLike));

		RooFormulaVar * wrongVtxInvCdf = new RooFormulaVar("wrongVtxInvCdf","1.-@0",RooArgList(*wrongVtxCdf) );
			
		vtxProb_ = new RooGenFunction(*(w->function("vtxProb")),RooArgList(*(w->var("vtxLike")),*(w->var("nVtx"))),RooArgList()); 
		sigProb_ = new RooGenFunction(*correctVtxCdf,RooArgList(*vtxLike),RooArgList());
		bkgProb_ = new RooGenFunction(*wrongVtxInvCdf,RooArgList(*vtxLike),RooArgList());
	}

    // getting di-cluster time resolution from file, and making 1d projection as a function of absDeta
    rooFileTime_ = TFile::Open(timePerfHist);
    dtVSdEtaEBEB_ = (TH2F*) rooFileTime_->Get("EBEB/TOF-corr cluster time difference VS #Delta#eta RightVertex");
    dtVSdEtaEBEE_ = (TH2F*) rooFileTime_->Get("EBEE/TOF-corr cluster time difference VS #Delta#eta RightVertex");
    dtVSdEtaEEEE_ = (TH2F*) rooFileTime_->Get("EBEE/TOF-corr cluster time difference VS #Delta#eta RightVertex");
    dtVSdEtaEEEE_ = (TH2F*) dtVSdEtaEEEE_->Clone("EEEE/TOF-corr cluster time difference VS #Delta#eta RightVertex"); // temporarily, make a clone
    cout << "histo: " << dtVSdEtaEBEB_->GetName() << ", " << dtVSdEtaEBEB_->GetName() << " and " << dtVSdEtaEEEE_->GetName() << " gotten from the file; timeOptimismFactor is: " << timeOptimismFactor <<endl;


    char tmpBuffer[10];
    std::vector<TH2F*>       dtVSdEtaS_;     dtVSdEtaS_.push_back(dtVSdEtaEBEB_); dtVSdEtaS_.push_back(dtVSdEtaEBEE_); dtVSdEtaS_.push_back(dtVSdEtaEEEE_);
    std::vector<std::string> dtVSdEtaNameS_; dtVSdEtaNameS_.push_back("EBEB");    dtVSdEtaNameS_.push_back("EBEE");    dtVSdEtaNameS_.push_back("EEEE");
    // project the 2d histo tof-corr-deltaT_VS_absDeltaEta into 1d slices at absDeltaEta==const
    for(int u=0; u<dtVSdEtaS_.size(); u++){
      TH2F* the2dHistoRW    = dtVSdEtaS_[u];
      int   numDeltaEtaBins = the2dHistoRW->GetNbinsX(); 
      int   numDeltaTBins   = the2dHistoRW->GetNbinsY();
      
      for(int theBinX=1; theBinX<=numDeltaEtaBins; theBinX++){ 
	sprintf (tmpBuffer, "%d", theBinX);
	std::string tmp1dHistoName = dtVSdEtaNameS_[u] + std::string(" TOF-corr cluster #Deltat bin ") + std::string(tmpBuffer); 
	TH1F theHisto(tmp1dHistoName.c_str(),tmp1dHistoName.c_str(),the2dHistoRW->GetNbinsY(),the2dHistoRW->GetYaxis()->GetBinLowEdge(1),the2dHistoRW->GetYaxis()->GetBinUpEdge(numDeltaTBins)); 
	
	for(int theBinY=1; theBinY<=numDeltaTBins; theBinY++){ 
	  float tmpValue = the2dHistoRW->GetBinContent(theBinX,theBinY); 
	  theHisto.SetBinContent(theBinY,tmpValue); 
	} // loop on ye bins 
	
	if     (dtVSdEtaEBEB_==the2dHistoRW) dtEBEB_.push_back(theHisto);
	else if(dtVSdEtaEBEE_==the2dHistoRW) dtEBEE_.push_back(theHisto);
	else if(dtVSdEtaEEEE_==the2dHistoRW) dtEEEE_.push_back(theHisto);
	else  { std::cout << "the2dHistoRW " << the2dHistoRW << " does not match any of the foreseen. Bailing out." << std::endl; assert(0); }

      }// loop on absDeltaEta==const bins
    }// loop on type of events
    if(PADEBUG)	cout << "dtEBEB_ has " << dtEBEB_.size() << " elements, dtEBEE_ has " << dtEBEE_.size() << " and dtEEEE_ has " << dtEEEE_.size() << " elements." <<endl;

    if(PADEBUG)	cout << "InitRealMicroAnalysis END"<<endl;
}

// ----------------------------------------------------------------------------------------------------
void MicroAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
	static std::vector<HggVertexAnalyzer::getter_t> varMeths_(0);
	if( varMeths_.empty() ) {
		for( size_t iv=0; iv<vtxVarNames_.size(); ++iv ){
			varMeths_.push_back(HggVertexAnalyzer::dictionary()[vtxVarNames_[iv]].first);
		}
	}

    if(PADEBUG)	cout << "Analysis START; cur_type is: " << l.itype[l.current] <<endl;
   
    int cur_type = l.itype[l.current];
    float weight = l.sampleContainer[l.current_sample_index].weight;

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
    
    //flush tree
    if( jentry % 1000 ==  0 ) {
	    //
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
						((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), l.pho_residCorrEnergy[ipho],
						l.pho_isEB[ipho], l.pho_r9[ipho],
						l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_) );
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
		} else if( doEscaleSmear && cur_type == 0 ) {          // if it's data
			float sweight = 1.;
			eScaleDataSmearer->smearPhoton(phoInfo,sweight,l.run,0.);
			pweight *= sweight;
		}
		smeared_pho_energy[ipho] = phoInfo.energy();
		smeared_pho_r9[ipho] = phoInfo.r9();
		smeared_pho_weight[ipho] = pweight;
    }
   
    sumev += weight;
    // FIXME pass smeared R9
    int diphoton_id = l.DiphotonCiCSelection(l.phoNOCUTS, l.phoNOCUTS, leadEtCut, subleadEtCut, 4,false, &smeared_pho_energy[0] );
    // Better do this without applying the ID
    /// int diphoton_id = 0;
    // std::cerr << "Selected pair " << l.dipho_n << " " << diphoton_id << std::endl;
    if (diphoton_id <= -1 || l.dipho_n < 1 ) return;
    
    
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
    int category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
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
    
    
    // MicroAnalysis-specific stuff
    // TODO dipho selection should match the Higgs truth (not just the best dipho)
    vtxAna_.setPairID(diphoton_id);
    nVert_ = l.vtx_std_n;
    nPU_ = l.pu_n;
    evWeight_ = evweight;
    pho1_ = &lead_p4;
    pho2_ = &sublead_p4;
    dipho_ = &Higgs;

    TVector3 *genVtx = (TVector3*)l.gv_pos->At(0);
    int closest_idx = -1;
    float closest_dist = 1e6;
    vector<float> closest_reco(nVert_,1e6);
    for (int vi=0;vi<l.vtx_std_n;vi++){
    	TVector3 *curVtx = (TVector3*)l.vtx_std_xyz->At(vi);
    	float dist = (*curVtx-*genVtx).Mag();
    	if (dist < closest_dist){
//    		cout << dist << " ";
    		closest_idx = vi;
    		closest_dist = dist;
    	}

        for (int vj=0;vj<l.vtx_std_n;vj++){
        	if (vj==vi) continue;
        	TVector3 *otherVtx = (TVector3*)l.vtx_std_xyz->At(vj);
        	float dist = (*curVtx-*otherVtx).Mag();
        	if (dist < closest_reco[vi]){
        		closest_reco[vi] = dist;
        	}
        }
    }
//    cout<<endl;
    assert(closest_idx!=-1);

    // Look for conversions and combine conversion information
    nConv_=0;
    convCompat_=-1;
    Float_t zconv=0;
    Float_t szconv=0;
//	cout << "filling photon infos" <<endl;
//	cout << l.dipho_leadind[diphoton_id] << " " << l.dipho_subleadind[diphoton_id] << endl;
	PhotonInfo pho1Info=l.fillPhotonInfos(l.dipho_leadind[diphoton_id]		,vtxAlgoParams.useAllConversions);
//	cout << "lead done" <<endl;
	PhotonInfo pho2Info=l.fillPhotonInfos(l.dipho_subleadind[diphoton_id]	,vtxAlgoParams.useAllConversions);
//	cout << "sublead done" <<endl;
    if ( (pho1Info.isAConversion() || pho2Info.isAConversion() ) )  {
    	nConv_=1;
    	if (pho1Info.isAConversion()  && !pho2Info.isAConversion() ){
    		zconv  = vtxConv_.vtxZ (pho1Info);
    		szconv = vtxConv_.vtxdZ(pho1Info);
    	}
    	if (pho2Info.isAConversion() && !pho1Info.isAConversion()){
    		zconv  = vtxConv_.vtxZ (pho2Info);
    		szconv = vtxConv_.vtxdZ(pho2Info);
    	}

    	if (pho1Info.isAConversion() && pho2Info.isAConversion()){
    		nConv_=2;
    		float z1  = vtxConv_.vtxZ (pho1Info);
    		float sz1 = vtxConv_.vtxdZ(pho1Info);

    		float z2  = vtxConv_.vtxZ (pho2Info);
    		float sz2 = vtxConv_.vtxdZ(pho2Info);

    		zconv  = (z1/sz1/sz1 + z2/sz2/sz2)/(1./sz1/sz1 + 1./sz2/sz2 );  // weighted average
    		szconv = sqrt( 1./(1./sz1/sz1 + 1./sz2/sz2)) ;

    		convCompat_=abs(z1-z2)/sqrt(sz1*sz1+sz2*sz2);
    	}
    }

    vector<int> & rankedVtxs = (*l.vtx_std_ranked_list)[diphoton_id];
    vtxAna_.preselection(rankedVtxs);
    vector<int> noConvRank = vtxAna_.rankprod(rankVariables_);
    vtxAna_.setNConv(nConv_);
    for (int vi=0;vi<l.vtx_std_n;vi++){
    	ptasym_ = vtxAna_.ptasym(vi);
    	ptbal_ = vtxAna_.ptbal(vi);
    	logsumpt2_ = vtxAna_.logsumpt2(vi);
	rprod_ = vtxAna_.rcomb(vi); 
	convRprod_ = rprod_;
	if( vi == rankedVtxs[0] && vi != noConvRank[0] ) {
		convRprod_  = 0.;
	}
    	TVector3 *curVtx = (TVector3*)l.vtx_std_xyz->At(vi);
    	dZToGen_ = (*curVtx-*genVtx).Mag();
    	dZtoClosest_ = closest_reco[vi];

    	isClosestToGen_ = (vi == closest_idx);
	
	// abs value of delta eta difference
	absDeltaEta_ = fabs( pho1_->PseudoRapidity() - pho2_->PseudoRapidity() );
	
	// tofCorrTdiff_ is the difference (t_cluster_lead - t_cluster_sub), corrected for the tof difference assuming the current vertex
	TVector3 caloPosLead    = ( * (TVector3*) l.pho_calopos->At(  l.dipho_leadind[diphoton_id] ) ) ;
	TVector3 caloPosSubLead = ( * (TVector3*) l.pho_calopos->At(  l.dipho_subleadind[diphoton_id] ) ) ;
	TVector3 currentVertex  = ( * (TVector3*)l.vtx_std_xyz->At(vi) );
	TVector3 closestVertex  = ( * (TVector3*)l.vtx_std_xyz->At(closest_idx) );
	tofCorrTdiff_           = getDeltaTof(caloPosLead, caloPosSubLead, closestVertex); // "true" time difference, given SC and true_vtx positions
	tofCorrTdiff_          -= getDeltaTof(caloPosLead, caloPosSubLead, currentVertex); //  tof correction: for current vtx;
                                                                                           //  if current vtx is not true => residual, which is the discriminating variable
	// at last adding smearing mimiking resolution effects; use "tof-corrected (t_cluster_lead - t_cluster_sub)" from Z data
	tofCorrTdiff_          += timeOptimismFactor * getTimeResol(absDeltaEta_, l.pho_isEB[l.dipho_leadind[diphoton_id]] , l.pho_isEB[l.dipho_subleadind[diphoton_id]] );
	//std::cout << "position: " << caloPosLead.PseudoRapidity() << " eta: " << pho1_->PseudoRapidity() << " sublead: " << caloPosSubLead.PseudoRapidity() << " eta2: " << pho2_->PseudoRapidity()  << " vertex: " << currentVertex.Z() << " tofCorrTdiff_: " << tofCorrTdiff_ << std::endl; // DEBUG

    	dZToConv_=-1;
    	pullToConv_=-1;
        if(nConv_>0){
        	dZToConv_ = abs( curVtx->Z() - zconv );
        	pullToConv_ = dZToConv_/szconv;
        }

	vtxAna_.setPullToConv( vi, pullToConv_ );
	/// vtxAna_.setTofCorrTdiff( vi, tofCorrTdiff_ );
	for( size_t ivar=0; ivar<vtxVarNames_.size(); ++ivar ){
		vtxVars_[ivar] = (vtxAna_.*(varMeths_[ivar]))(vi);
	}

    	uTree_->Fill();
	if( jentry % 1000 == 0 ) {
		uFile_->cd();
		uTree_->Write(0,TObject::kWriteDelete);
	}
    }

    //per event tree

    //put stupid values in vectors
    MVA_.assign(storeNVert,-10);
    rankprod_.assign(storeNVert,-10);
    convRankprod_.assign(storeNVert,-10);
    prob_.assign(storeNVert,0.);
    corrProb_.assign(storeNVert,0.);
    wrongProb_.assign(storeNVert,0.);
    zRMSn_.assign(storeNVert,0.);
    dZ_.assign(storeNVert,-100);
    diphoM_.assign(storeNVert,-2);
    diphoCosTheta_.assign(storeNVert,-2);
    diphoCosDeltaPhi_.assign(storeNVert,-2);
    diphoPt_.assign(storeNVert,-2);
    dZToConvV_.assign(storeNVert,-1);
    pullToConvV_.assign(storeNVert,-1);

    TLorentzVector truevtx_lead_pho = l.get_pho_p4( l.dipho_leadind[diphoton_id], genVtx, &smeared_pho_energy[0]);
    TLorentzVector truevtx_sublead_pho = l.get_pho_p4( l.dipho_subleadind[diphoton_id], genVtx, &smeared_pho_energy[0]);
    TLorentzVector truevtx_dipho = truevtx_lead_pho+truevtx_sublead_pho;
    mTrue_ = gP4.M(); 
    mTrueVtx_ = truevtx_dipho.M(); 
    zTrue_ = genVtx->Z();
    rTrue_ = genVtx->Mag();
    category_ = category;
    alphaTrue_ = truevtx_lead_pho.Angle(truevtx_sublead_pho.Vect());

    //get list of vertices as ranked for the selected diphoton
    if( useMvaRanking ) {
	    rankedVtxs = vtxAna_.rank(*tmvaReader_,tmvaMethod);
    } else {
	    vtxAna_.evaluate(*tmvaReader_,tmvaMethod);
    }
    dZTrue_ = ( *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[0]) - *genVtx).Z();
    Double_t vtxProbInputs[2];
    float wtot = 0.;
    zRMS_ = 0.;
    for (size_t vi=0;vi<rankedVtxs.size();vi++) {
    	if(vi>=storeNVert) break;
    	MVA_[vi] = vtxAna_.mva(rankedVtxs[vi]);
    	rankprod_[vi] = vtxAna_.rcomb(rankedVtxs[vi]);
	convRankprod_[vi] = vi == 0 && rankedVtxs[0] != noConvRank[0] ? 0. : rankprod_[vi];
    	dZ_[vi] = ( *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[vi]) - *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[0])).Z();
	
    	TLorentzVector lead_pho = l.get_pho_p4( l.dipho_leadind[diphoton_id], rankedVtxs[vi], &smeared_pho_energy[0]);
    	TLorentzVector sublead_pho = l.get_pho_p4( l.dipho_subleadind[diphoton_id], rankedVtxs[vi], &smeared_pho_energy[0]);
    	TLorentzVector dipho = lead_pho+sublead_pho;
    	diphoM_[vi]				= dipho.M();
		diphoCosTheta_[vi]		= TMath::TanH(0.5*(lead_pho.Rapidity()-sublead_pho.Rapidity()));
		diphoCosDeltaPhi_[vi]	= TMath::Cos(lead_pho.Phi()-sublead_pho.Phi());
		diphoPt_[vi]			= dipho.Pt();
		alpha_[vi] = lead_pho.Angle(sublead_pho.Vect());
	
		if( vtxProb_ ) {
			vtxProbInputs[0] = MVA_[vi];
			vtxProbInputs[1] = nVert_;
			prob_[vi] = (*vtxProb_)( vtxProbInputs );
			corrProb_[vi] = (*sigProb_)( MVA_[vi] );
			wrongProb_[vi] = (*bkgProb_)( MVA_[vi] );
		}
		TVector3 *curVtx = (TVector3*)l.vtx_std_xyz->At(rankedVtxs[vi]);
		dZToConvV_[vi]	= abs( curVtx->Z() - zconv );
		pullToConvV_[vi]= dZToConvV_[vi]/szconv;

		if( vi > 0 ) {
			zRMSn_[vi] = zRMS_ + dZ_[vi]*dZ_[vi]*prob_[vi];
			zRMS_ = zRMSn_[vi];
			wtot  += prob_[vi];
			zRMSn_[vi] /= wtot;
		}
    }
    zRMS_ /= wtot;
    evTree_->Fill();
    
    vtxAna_.clear();
}

//http://root.cern.ch/root/html/TH1.html#TH1:GetRandom
float MicroAnalysis::getTimeResol(float absDeltaEta, bool iseb1, bool iseb2){

  if (iseb1 && iseb2) {
    int theBin = dtVSdEtaEBEB_->GetXaxis()->FindBin(absDeltaEta_);

    if(theBin>dtVSdEtaEBEB_->GetNbinsX()) theBin=dtVSdEtaEBEB_->GetNbinsX(); // relocate overflows into last real bin

    // protect against 1d histograms which are empty and canont give random number 
    while (dtEBEB_[theBin-1].Integral()<5 && theBin>0 && theBin<dtVSdEtaEBEB_->GetNbinsX() ) {
      theBin--;}

    return dtEBEB_[theBin-1].GetRandom();
  }
  else if ( (!iseb1 && iseb2) || (iseb1 && !iseb2) ) {
    int theBin = dtVSdEtaEBEE_->GetXaxis()->FindBin(absDeltaEta_);

    if(theBin>dtVSdEtaEBEE_->GetNbinsX()) theBin=dtVSdEtaEBEE_->GetNbinsX(); // relocate overflows into last real bin

    // protect against 1d histograms which are empty and canont give random number 
    while (dtEBEE_[theBin-1].Integral()<5 && theBin>0 && theBin<dtVSdEtaEBEB_->GetNbinsX() ) {
      theBin--;}
    
    return dtEBEE_[theBin-1].GetRandom();
  }
  else {
    int theBin = dtVSdEtaEEEE_->GetXaxis()->FindBin(absDeltaEta_);
    if(theBin>dtVSdEtaEEEE_->GetNbinsX()) theBin=dtVSdEtaEEEE_->GetNbinsX(); // relocate overflows into last real bin
    while (dtEEEE_[theBin-1].Integral()<5 && theBin>0 && theBin<dtVSdEtaEEEE_->GetNbinsX() ) {
      theBin++;}
    return dtEEEE_[theBin-1].GetRandom();
  }
}


float MicroAnalysis::getDeltaTof(TVector3 &posLead, TVector3 &posSubLead, TVector3 &posVertex){

  //std::cout << "\n getExtraTravelTime(posLead,posVertex): " << ( getExtraTravelTime(posLead,posVertex) - getExtraTravelTime(posSubLead,posVertex) ) << std::endl; 
  return getExtraTravelTime(posLead,posVertex) - getExtraTravelTime(posSubLead,posVertex);
}


float MicroAnalysis::getExtraTravelTime(TVector3 &posSC, TVector3 &posVertex){
  float travelled = sqrt( pow(posSC.X()-posVertex.X(), 2) + 
			  pow(posSC.Y()-posVertex.Y(), 2) + 
			  pow(posSC.Z()-posVertex.Z(), 2)  );
  float nominal   = sqrt( pow(posSC.X(), 2) + 
			  pow(posSC.Y(), 2) + 
			  pow(posSC.Z(), 2)  );
  
  //std::cout << "posSC.X(): " << posSC.X() << " posVertex.X(): " << posVertex.X() << " posSC.Z(): " << posSC.Z() << " posVertex.Z(): " << posVertex.Z() << " travelled: " << travelled << " nominal: " << nominal << " \t return: " << ((travelled-nominal)/100./speedOfLight*1e9 ) << std::endl; // DEBUG
  
  return (travelled-nominal)/100./speedOfLight*1.e9;
}
