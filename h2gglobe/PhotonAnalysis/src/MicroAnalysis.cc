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
	pho1_(0), pho2_(0), dipho_(0)
{
	cout << "Constructing MicroAnalysis" << endl;
}

// ----------------------------------------------------------------------------------------------------
MicroAnalysis::~MicroAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void MicroAnalysis::Term(LoopAll& l) 
{
	cout<<endl;
	cout << "Closing-up MicroAnalysis" << endl;
	//TODO close up
	uFile_->cd();
	uTree_->Write(0,TObject::kWriteDelete);
	evTree_->Write(0,TObject::kWriteDelete);
	uFile_->Close();
}

// ----------------------------------------------------------------------------------------------------
void MicroAnalysis::Init(LoopAll& l) 
{
	StatAnalysis::Init(l);

	cout << "Initializing MicroAnalysis" << endl;
	if(PADEBUG) cout << "InitRealMicroAnalysis START"<<endl;
	uFile_ = TFile::Open(uFileName,"recreate");

	// per vertex tree
	uTree_ = new TTree("utree","MicroAnalysis per Vertex Tree");
	//LINK http://root.cern.ch/root/html/TTree.html
	uTree_->Branch("dZToGen",&dZToGen_);
	uTree_->Branch("dZToClosest",&dZtoClosest_);
	uTree_->Branch("nVert",&nVert_);
	uTree_->Branch("nPU",&nPU_);
	uTree_->Branch("evWeight",&evWeight_);
	uTree_->Branch("pho1",&pho1_,32000,0);
	uTree_->Branch("pho2",&pho2_,32000,0);
	uTree_->Branch("dipho",&dipho_,32000,0);
	uTree_->Branch("ptasym",&ptasym_);
	uTree_->Branch("ptbal",&ptbal_);
	uTree_->Branch("logsumpt2",&logsumpt2_);
	uTree_->Branch("isClosestToGen",&isClosestToGen_);


	//TMVA
	tmvaVariables_.push_back("ptbal"), tmvaVariables_.push_back("ptasym"), tmvaVariables_.push_back("logsumpt2");
	tmvaReader_ = new TMVA::Reader( "!Color:!Silent" );
	HggVertexAnalyzer::bookVariables( *tmvaReader_, tmvaVariables_ );
	tmvaReader_->BookMVA( tmvaMethod, tmvaWeights );

	MVA_.resize(storeNVert);
	dZ_.resize(storeNVert);
	diphoM_.resize(storeNVert);
	diphoCosTheta_.resize(storeNVert);
	diphoDeltaPhi_.resize(storeNVert);
	diphoPt_.resize(storeNVert);

	evTree_ = new TTree("evtree","MicroAnalysis per Event Tree");
	evTree_->Branch("dZTrue",&dZTrue_);
	evTree_->Branch("nVert",&nVert_);
	evTree_->Branch("evWeight",&evWeight_);
	for(int i=0;i<storeNVert; i++){
		evTree_->Branch(Form("MVA%d",i),&MVA_[i]);
		evTree_->Branch(Form("dZ%d",i),&dZ_[i]);
		evTree_->Branch(Form("diphoM%d",i),&diphoM_[i]);
		evTree_->Branch(Form("diphoCosTheta%d",i),&diphoCosTheta_[i]);
		evTree_->Branch(Form("diphoDeltaPhi%d",i),&diphoDeltaPhi_[i]);
		evTree_->Branch(Form("diphoPt%d",i),&diphoPt_[i]);
	}


    if(PADEBUG)	cout << "InitRealMicroAnalysis END"<<endl;
}

// ----------------------------------------------------------------------------------------------------
void MicroAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
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
    // int diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,false, &smeared_pho_energy[0] );
    // Better do this without applying the ID
    int diphoton_id = 0;
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

    for (int vi=0;vi<l.vtx_std_n;vi++){
    	ptasym_ = vtxAna_.ptasym(vi);
    	ptbal_ = vtxAna_.ptbal(vi);
    	logsumpt2_ = vtxAna_.logsumpt2(vi);

    	TVector3 *curVtx = (TVector3*)l.vtx_std_xyz->At(vi);
    	dZToGen_ = (*curVtx-*genVtx).Mag();
    	dZtoClosest_ = closest_reco[vi];

    	isClosestToGen_ = (vi == closest_idx);

    	uTree_->Fill();
    }

    //per event tree
    vector<int> & rankedVtxs = (*l.vtx_std_ranked_list)[diphoton_id];
    vtxAna_.preselection(rankedVtxs);
    vtxAna_.evaluate(*tmvaReader_,tmvaMethod);
    /// for (int vi=0;vi<l.vtx_std_n;vi++){
    /// MVA_.assign(-10,storeNVert);
    dZTrue_ = ( *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[0]) - *genVtx).Z();
    for (size_t vi=0;vi<rankedVtxs.size();vi++){
    	if(vi>=storeNVert) break;
    	MVA_[vi] = vtxAna_.mva(rankedVtxs[vi]);
	dZ_[vi] = ( *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[vi]) - *(TVector3*)l.vtx_std_xyz->At(rankedVtxs[0])).Z();
    }
    evTree_->Fill();


}
