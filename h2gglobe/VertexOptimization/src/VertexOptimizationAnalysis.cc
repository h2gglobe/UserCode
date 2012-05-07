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
}

// ----------------------------------------------------------------------------------------------------
VertexOptimizationAnalysis::~VertexOptimizationAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void VertexOptimizationAnalysis::Term(LoopAll& l) 
{
}

// ----------------------------------------------------------------------------------------------------
void VertexOptimizationAnalysis::Init(LoopAll& l) 
{
    doVtxEffSmear = false;
    doTriggerEffSmear = false;
    doSystematics = false;

    StatAnalysis::Init(l);

    // per vertex tree
    uFile_ = TFile::Open(uFileName,"recreate");
    uFile_->cd(); 

    uTree_ = new TTree("utree","MicroAnalysis per Vertex Tree");
    uTree_->Branch("nVert",   &nVert_   );
    uTree_->Branch("nPU",     &nPU_     );
    uTree_->Branch("evWeight",&evWeight_);
    
    vtxVarNames_.push_back("ptvtx"), vtxVarNames_.push_back("ptasym"), vtxVarNames_.push_back("ptratio"), 
	vtxVarNames_.push_back("ptbal"), vtxVarNames_.push_back("logsumpt2"), vtxVarNames_.push_back("ptmax3"), 
	vtxVarNames_.push_back("ptmax"), vtxVarNames_.push_back("nchthr"), vtxVarNames_.push_back("sumtwd"),
	vtxVarNames_.push_back("pulltoconv"), vtxVarNames_.push_back("nch"); 
    vtxVars_.resize( vtxVarNames_.size() );
    for( size_t iv=0; iv<vtxVarNames_.size(); ++iv ) {
	uTree_->Branch(vtxVarNames_[iv].c_str(),&vtxVars_[iv]);
    }
    
    //////////// // per event tree
    //////////// evTree_ = new TTree("evtree","MicroAnalysis per Event Tree");
    //////////// evTree_->Branch("dZTrue",&dZTrue_);
    //////////// evTree_->Branch("alphaTrue",&alphaTrue_);
    //////////// evTree_->Branch("zTrue",&zTrue_);
    //////////// evTree_->Branch("rTrue",&rTrue_);
    //////////// evTree_->Branch("zRMS",&zRMS_);
    //////////// evTree_->Branch("category",&category_,"category/I");
    //////////// evTree_->Branch("mTrue",&mTrue_);
    //////////// evTree_->Branch("mTrueVtx",&mTrueVtx_);
    //////////// evTree_->Branch("nVert",&nVert_);
    //////////// evTree_->Branch("evWeight",&evWeight_);
    //////////// evTree_->Branch("nConv",&nConv_);
    //////////// evTree_->Branch("ConvCompat",&convCompat_);
    //////////// for(int i=0;i<storeNVert; i++){
    //////////// 	evTree_->Branch(Form("MVA%d",i),&MVA_[i]);
    //////////// 	evTree_->Branch(Form("alpha%d",i),&alpha_[i]);
    //////////// 	evTree_->Branch(Form("zRMSn%d",i),&zRMSn_[i]);
    //////////// 	evTree_->Branch(Form("dZ%d",i),&dZ_[i]);
    //////////// }
    
}

// ----------------------------------------------------------------------------------------------------
bool VertexOptimizationAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, int & category, int & diphoton_id,
			      bool & isCorrectVertex,
			      bool isSyst=false, 
			      float syst_shift=0., bool skipSelection=false,
			      BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0)
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
        diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
	evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;

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
            isCorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
        }
        float ptHiggs = Higgs.Pt();
	
	vtxAna_.setPairID(diphoton_id);
	nVert_    = l.vtx_std_n;
	nPU_      = l.pu_n;
	evWeight_ = evweight;
	for(int vi=0; vi<l.vtx_std_n; ++vi) {
	    
	    for( size_t ivar=0; ivar<vtxVarNames_.size(); ++ivar ){
		vtxVars_[ivar] = (vtxAna_.*(varMeths_[ivar]))(vi);
	    }
	    
	}
    }	
}

// ----------------------------------------------------------------------------------------------------
bool VertexOptimizationAnalysis::SelectEvents(LoopAll&, int)
{
    return true;
}

// ----------------------------------------------------------------------------------------------------
void VertexOptimizationAnalysis::FillReductionVariables(LoopAll& l, int jentry)
{
}
   
// ----------------------------------------------------------------------------------------------------
bool VertexOptimizationAnalysis::SelectEventsReduction(LoopAll&, int)
{
    return true;
}

// ----------------------------------------------------------------------------------------------------
bool VertexOptimizationAnalysis::SkimEvents(LoopAll&, int)
{
    return true;
}

// ----------------------------------------------------------------------------------------------------
void VertexOptimizationAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
    vtxAna_.branches(outputTree,"vtx_std_");    
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
