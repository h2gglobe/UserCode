#include "../interface/MicroAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
MicroAnalysis::MicroAnalysis()  : 
    name_("MicroAnalysis"),
    pho1_p4_(0), pho2_p4_(0), dipho_p4_(0)
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
	uFile_->Close();
}

// ----------------------------------------------------------------------------------------------------
void MicroAnalysis::Init(LoopAll& l) 
{
	StatAnalysis::Init(l);

	cout << "Initializing MicroAnalysis" << endl;
	if(PADEBUG) cout << "InitRealMicroAnalysis START"<<endl;
	//TODO open ttree
	uFile_ = TFile::Open(l.histFileName,"recreate");
	uTree_ = new TTree("utree","MicroAnalysis Tree");

	uTree_->Branch("dZtoGen",&dZtoGen_,"dZtoGen/F");
	uTree_->Branch("nVert",&nVert_,"nVert/I");
	uTree_->Branch("pho1_p4",&pho1_p4_,32000,0);
	uTree_->Branch("pho2_p4",&pho2_p4_,32000,0);
	uTree_->Branch("dipho1_p4",&dipho_p4_,32000,0);
	uTree_->Branch("ptasym",&ptasym_,"ptasym/F");
	uTree_->Branch("ptbal",&ptbal_,"ptbal/F");
	uTree_->Branch("logsumpt2",&logsumpt2_,"logsumpt2/F");
	uTree_->Branch("isClosestToGen",&isClosestToGen_,"isClosestToGen/O");

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
    int diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,false, &smeared_pho_energy[0] ); 
    // std::cerr << "Selected pair " << l.dipho_n << " " << diphoton_id << std::endl;
	if (diphoton_id <= -1 ) return;


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
    nVert_ = l.gv_n;
    pho1_p4_ = &lead_p4;
    pho2_p4_ = &sublead_p4;
    dipho_p4_ = &Higgs;

    TVector3 *genVtx = (TVector3*)l.gv_pos->At(0);
    int closest_idx = -1;
    float closest_dist = 1e6;
    for (int vi=0;vi<l.gv_n;vi++){
    	TVector3 *curVtx = (TVector3*)l.vtx_std_xyz->At(vi);
    	float dist = (*curVtx-*genVtx).Mag();
    	if (dist < closest_dist){
    		closest_idx = vi;
    		closest_dist = dist;
    	}
    }

    for (int vi=0;vi<l.gv_n;vi++){
    	ptasym_ = vtxAna_.ptasym(vi);
    	ptbal_ = vtxAna_.ptbal(vi);
    	logsumpt2_ = vtxAna_.logsumpt2(vi);

    	TVector3 *curVtx = (TVector3*)l.vtx_std_xyz->At(vi);

    	dZtoGen_ = (*curVtx-*genVtx).Mag();

    	isClosestToGen_ = (vi == closest_idx);

    	uTree_->Fill();
    }

}

