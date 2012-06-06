#include "../interface/VertexOptimizationAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
VertexOptimizationAnalysis::VertexOptimizationAnalysis()  
{
    name_ = "VertexOptimizationAnalysis";

    minBiasRefName = "aux/minBiasRef.root";
}

// ----------------------------------------------------------------------------------------------------
VertexOptimizationAnalysis::~VertexOptimizationAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void VertexOptimizationAnalysis::Term(LoopAll& l) 
{
    l.outputFile->cd();
    uTree_->Write();
    hMinBiasSpecturm_->Write();
    hHiggsSpecturm_->Write();
    uTree_->SetDirectory(0);
    hMinBiasSpecturm_->SetDirectory(0);
    hHiggsSpecturm_->SetDirectory(0);
}

// ----------------------------------------------------------------------------------------------------
void VertexOptimizationAnalysis::Init(LoopAll& l) 
{
    doSystematics = false;

    TFile * mbRef = TFile::Open(minBiasRefName);
    hMinBiasRef_ = (TH1*)mbRef->Get("minBiasSpecturm")->Clone("hMinBiasRef");
    hMinBiasRef_->SetDirectory(0);
    mbRef->Close();
    
    StatAnalysis::Init(l);
}

void VertexOptimizationAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
    // per vertex tree
    pho1_=0, pho2_=0, dipho_=0;
    TDirectory * pwd = gDirectory;
    outputTree->GetDirectory()->cd();
    uTree_ = new TTree("vtxOptTree","Vertex optimization tree");
    uTree_->Branch("nVert",   &nVert_   );
    uTree_->Branch("nPU",     &nPU_     );
    uTree_->Branch("evWeight",&evWeight_);
    uTree_->Branch("ksprob",&ksprob_);
    uTree_->Branch("isClosestToGen",&isClosestToGen_);
    uTree_->Branch("passCiC",&passCiC_);
    uTree_->Branch("pho1",&pho1_,32000,0);
    uTree_->Branch("pho2",&pho2_,32000,0);
    uTree_->Branch("dipho",&dipho_,32000,0);
    
    vtxVarNames_.push_back("ptvtx"), vtxVarNames_.push_back("ptasym"), vtxVarNames_.push_back("ptratio"), 
	vtxVarNames_.push_back("ptbal"), vtxVarNames_.push_back("logsumpt2"), vtxVarNames_.push_back("ptmax3"), 
	vtxVarNames_.push_back("ptmax"), vtxVarNames_.push_back("nchthr"), vtxVarNames_.push_back("sumtwd"),
	vtxVarNames_.push_back("pulltoconv"), vtxVarNames_.push_back("limpulltoconv"), 
	vtxVarNames_.push_back("nch"), vtxVarNames_.push_back("nconv"), vtxVarNames_.push_back("nlegs"), 
	vtxVarNames_.push_back("mva"); 
    vtxVars_.resize( vtxVarNames_.size() );
    for( size_t iv=0; iv<vtxVarNames_.size(); ++iv ) {
	uTree_->Branch(vtxVarNames_[iv].c_str(),&vtxVars_[iv]);
    }

    hMinBiasSpecturm_ = new TH1F("minBiasSpecturm","minBiasSpecturm;;p_{T} (GeV/c)",200,0,20);
    hHiggsSpecturm_   = new TH1F("higgsSpecturm","higgsSpecturm;;p_{T} (GeV/c)",200,0,20);
    
    //////////// // per event tree
    //////////// evTree_ = new TTree("evtree","MicroAnalysis per Event Tree");
    //////////// outputTree->Branch("dZTrue",&dZTrue_);
    //////////// outputTree->Branch("alphaTrue",&alphaTrue_);
    //////////// outputTree->Branch("zTrue",&zTrue_);
    //////////// outputTree->Branch("rTrue",&rTrue_);
    //////////// outputTree->Branch("zRMS",&zRMS_);
    //////////// outputTree->Branch("category",&category_,"category/I");
    //////////// outputTree->Branch("mTrue",&mTrue_);
    //////////// outputTree->Branch("mTrueVtx",&mTrueVtx_);
    //////////// outputTree->Branch("nVert",&nVert_);
    //////////// outputTree->Branch("evWeight",&evWeight_);
    //////////// outputTree->Branch("nConv",&nConv_);
    //////////// outputTree->Branch("ConvCompat",&convCompat_);
    //////////// for(int i=0;i<storeNVert; i++){
    //////////// 	evTree_->Branch(Form("MVA%d",i),&MVA_[i]);
    //////////// 	evTree_->Branch(Form("alpha%d",i),&alpha_[i]);
    //////////// 	evTree_->Branch(Form("zRMSn%d",i),&zRMSn_[i]);
    //////////// 	evTree_->Branch(Form("dZ%d",i),&dZ_[i]);
    //////////// }
    
    pwd->cd();
    PhotonAnalysis::ReducedOutputTree(l,0);
}

// ----------------------------------------------------------------------------------------------------
bool VertexOptimizationAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, 
					      float & mass, float & evweight, int & category, int & diphoton_id,
					      bool & isCorrectVertex, float &kinematic_bdtout,
					      bool isSyst, 
					      float syst_shift, bool skipSelection,
					      BaseGenLevelSmearer *genSys, 
					      BaseSmearer *phoSys, BaseDiPhotonSmearer * diPhoSys)
{
    static std::vector<HggVertexAnalyzer::getter_t> varMeths_(0);
    if( varMeths_.empty() ) {
	for( size_t iv=0; iv<vtxVarNames_.size(); ++iv ){
	    varMeths_.push_back(HggVertexAnalyzer::dictionary()[vtxVarNames_[iv]].first);
	}
    }
    
    assert( ! isSyst );

    int cur_type = l.itype[l.current];
    float sampleweight = l.sampleContainer[l.current_sample_index].weight;
    /// diphoton_id = -1;
    
    std::pair<int,int> diphoton_index;
   
    // do gen-level dependent first (e.g. k-factor); only for signal
    genLevWeight=1.;
    if(cur_type!=0 ) {
	applyGenLevelSmearings(genLevWeight,gP4,l.pu_n,cur_type,genSys,syst_shift);
    }

    // first apply corrections and smearing on the single photons 
    smeared_pho_energy.clear(); smeared_pho_energy.resize(l.pho_n,0.); 
    smeared_pho_r9.clear();     smeared_pho_r9.resize(l.pho_n,0.); 
    smeared_pho_weight.clear(); smeared_pho_weight.resize(l.pho_n,1.);
    applySinglePhotonSmearings(smeared_pho_energy, smeared_pho_r9, smeared_pho_weight, cur_type, l, energyCorrected, energyCorrectedError,
			       phoSys, syst_shift);
    
    diphoton_id = l.DiphotonCiCSelection(l.phoNOCUTS, l.phoNOCUTS, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 

    if (diphoton_id > -1 ) {
	int closest_id = -1;
	float minDist = 999.;
	for(int vi=0; vi<l.vtx_std_n; ++vi) {
	    float dist =  fabs( ( *((TVector3*)l.vtx_std_xyz->At(vi)) - *((TVector3*)l.gv_pos->At(0)) ).Z() );
	    if( dist < minDist ) {
		closest_id = vi;
		minDist = dist;
	    }
	}
	
        diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
	evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;

	std::vector<std::vector<bool> > ph_passcut;
	passCiC_ = ( l.PhotonCiCSelectionLevel(diphoton_index.first,  closest_id, ph_passcut, 4, 0, &smeared_pho_energy[0] ) >= l.phoSUPERTIGHT && 
		     l.PhotonCiCSelectionLevel(diphoton_index.second, closest_id, ph_passcut, 4, 1, &smeared_pho_energy[0] ) >= l.phoSUPERTIGHT );
	
        TLorentzVector lead_p4, sublead_p4, Higgs;
        float lead_r9, sublead_r9;
        TVector3 * vtx;
	fillDiphoton(lead_p4, sublead_p4, Higgs, lead_r9, sublead_r9, vtx, &smeared_pho_energy[0], l, diphoton_id);  
      
        // FIXME pass smeared R9
	category = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
	mass     = Higgs.M();

	// apply di-photon level smearings and corrections
        int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
        if( cur_type != 0 && doMCSmearing ) {
	    applyDiPhotonSmearings(Higgs, *vtx, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), evweight, zero_, zero_,
				   diPhoSys, syst_shift);
        }
	isCorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
        float ptHiggs = Higgs.Pt();
	
	// fill optimization tree
	vtxAna_.setPairID(diphoton_id);
	nVert_    = l.vtx_std_n;
	nPU_      = l.pu_n;
	evWeight_ = evweight;
	pho1_ = &lead_p4;
	pho2_ = &sublead_p4;
	dipho_ = &Higgs;
	
	for(int vi=0; vi<l.vtx_std_n; ++vi) {
	    TH1 * h = (TH1*)hMinBiasRef_->Clone("h");
	    h->Reset("ICE");
	    isClosestToGen_ = (vi == closest_id);
	    
	    for( size_t ivar=0; ivar<vtxVarNames_.size(); ++ivar ) {
		vtxVars_[ivar] = (vtxAna_.*(varMeths_[ivar]))(vi);
	    }
	    
	    int ntks = l.vtx_std_ntks[vi];
	    for(int ti=0; ti<ntks; ++ti) {
		int tkind = (*l.vtx_std_tkind)[vi][ti];
		if( tkind >= l.tk_n ) { continue; }
		TLorentzVector* tkp4 = (TLorentzVector*)l.tk_p4->At(tkind);
		if( ! isClosestToGen_ ) {
		    hMinBiasSpecturm_->Fill(tkp4->Pt());
		} else {
		    hHiggsSpecturm_->Fill(tkp4->Pt());
		}
		h->Fill(tkp4->Pt());
	    }
	    
	    ksprob_ = hMinBiasRef_->KolmogorovTest(h);
	    delete h;
	    
	    uTree_->Fill();
	}
    }
    
    return false;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
