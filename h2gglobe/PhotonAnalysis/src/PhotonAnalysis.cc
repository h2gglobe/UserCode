#include "../interface/PhotonAnalysis.h"


#include "PhotonReducedInfo.h"
#include "Sorters.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#define PADEBUG 0

using namespace std;

// ----------------------------------------------------------------------------------------------------
PhotonAnalysis::PhotonAnalysis()  : 
    runStatAnalysis(false), doTriggerSelection(false),
    name_("PhotonAnalysis"),
    vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams),
    tmvaPerVtxMethod("BDTG"),
    tmvaPerVtxWeights(""),
    tmvaPerEvtMethod("evtBTG"),
    tmvaPerEvtWeights(""),
    energyCorrectionMethod("DaunceyAndKenzie"), energyCorrected(0), energyCorrectedError(0)
{
    addConversionToMva=true;
    mvaVertexSelection=false;
    useDefaultVertex=false;
    forcedRho = -1.;

    keepPP = true;
    keepPF = true;
    keepFF = true; 

    doSystematics = true;    

    zero_ = 0.; 
}

// ----------------------------------------------------------------------------------------------------
PhotonAnalysis::~PhotonAnalysis() 
{}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Term(LoopAll& l) 
{}

// ----------------------------------------------------------------------------------------------------
void readEnergyScaleOffsets(const std::string &fname, EnergySmearer::energySmearingParameters::eScaleVector &escaleOffsets, 
                            EnergySmearer::energySmearingParameters::phoCatVector &photonCategories, bool data=true
                            )
{
    // read in energy scale corrections to be applied in run ranges
    std::fstream in(fname.c_str());
    assert( in );
    char line[200];
    float EBHighR9, EBLowR9, EBm4HighR9, EBm4LowR9, EEHighR9, EELowR9; 
    char catname[200];
    float mineta, maxeta, minr9, maxr9, offset, err;
    int type; 
    int  first, last;
    do {
        in.getline( line, 200, '\n' );

        if( sscanf(line,"%d %d %f %f %f %f %f %f",&first, &last, &EBHighR9, &EBLowR9, &EBm4HighR9, &EBm4LowR9, &EEHighR9, &EELowR9) == 8 ) { 
            std::cerr << "Energy scale by run " <<  first<< " " <<  last<< " " <<  EBHighR9<< " " <<  EBLowR9 << " " <<  EBm4HighR9<< " " <<  EBm4LowR9<< " " <<  EEHighR9<< " " <<  EELowR9 << std::endl;
            
            assert( ! data );
            escaleOffsets.push_back(EnergyScaleOffset(first,last));
            escaleOffsets.back().scale_offset["EBHighR9"] = -1.*EBHighR9;
            escaleOffsets.back().scale_offset["EBLowR9"]  = -1.*EBLowR9;
            escaleOffsets.back().scale_offset["EBm4HighR9"] = -1.*EBm4HighR9;
            escaleOffsets.back().scale_offset["EBm4LowR9"]  = -1.*EBm4LowR9;
            escaleOffsets.back().scale_offset["EEHighR9"] = -1.*EEHighR9;
            escaleOffsets.back().scale_offset["EELowR9"]  = -1.*EELowR9;
            escaleOffsets.back().scale_offset_error["EBHighR9"] = 0.;
            escaleOffsets.back().scale_offset_error["EBLowR9"]  = 0.;
            escaleOffsets.back().scale_offset_error["EBm4HighR9"] = 0.;
            escaleOffsets.back().scale_offset_error["EBm4LowR9"]  = 0.;
            escaleOffsets.back().scale_offset_error["EEHighR9"] = 0.;
            escaleOffsets.back().scale_offset_error["EELowR9"]  = 0.;
        } else if( sscanf(line,"%s %d %f %f %f %f %d %d %f %f", &catname, &type, &mineta, &maxeta, &minr9, &maxr9, &first, &last, &offset, &err  ) == 10 ) { 
	    std::cerr << "Energy scale (or smering) by run " <<  catname << " " << type << " " << mineta << " " << maxeta << " " << minr9 << " " << maxr9 << " " << first << " " << last << " " << offset << " " << err << std::endl;
	    
	    assert( type>=0 && type<=2 );

            EnergySmearer::energySmearingParameters::eScaleVector::reverse_iterator escaleOffset = 
                find(escaleOffsets.rbegin(),escaleOffsets.rend(),std::make_pair(first,last));
            if( escaleOffset == escaleOffsets.rend() ) {
                std::cerr << "  adding new range range " << first << " " << last << std::endl;
                escaleOffsets.push_back(EnergyScaleOffset(first,last));
                escaleOffset = escaleOffsets.rbegin();
            }
            // chck if the category is already defined
            if( find(photonCategories.begin(), photonCategories.end(), std::string(catname) ) == photonCategories.end() ) {
                std::cerr << "  defining new category" << std::endl;
                photonCategories.push_back(PhotonCategory(mineta,maxeta,minr9,maxr9,(PhotonCategory::photon_type_t)type,catname));
            }
            // assign the scale offset and error for this category and this run range 
            escaleOffset->scale_offset[catname] = data ? -offset : offset;
            escaleOffset->scale_offset_error[catname] = err;
        }

    } while( in );
    
    in.close();
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::loadPuMap(const char * fname, TDirectory * dir, TH1 * target)
{
    std::fstream in(fname);
    assert( in );
    char line[200];
    int typid;
    char dname[200];
    do {
        in.getline( line, 200, '\n' );

        if( sscanf(line,"%d %s",&typid,dname) != 2 ) { continue; } 
        std::cerr << "Reading PU weights for sample " << typid << " from " << dname << std::endl;
        TDirectory * subdir = (TDirectory *)dir->Get(dname);
        assert( subdir != 0 );
        loadPuWeights(typid, subdir, target);
    } while ( in );
    in.close();
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::loadPuWeights(int typid, TDirectory * puFile, TH1 * target)
{
    cout<<"Reweighting events for pileup."<<endl;
    TH1 * hweigh = (TH1*) puFile->Get("weights");
    if( hweigh == 0 ) {
        hweigh = (TH1*) puFile->Get("NPUWeights");
    }
    if( hweigh != 0 ) { 
        cout<< " This is a pre-processed pileup reweighing file." <<endl;
        TH1 * gen_pu = (TH1*)puFile->Get("generated_pu");
        if( gen_pu == 0 ) {
            gen_pu = (TH1*)puFile->Get("NPUSource");
        }
	if( target != 0 ) {
	    hweigh->Reset("ICE");
	    for( int ii=1; ii<hweigh->GetNbinsX(); ++ii ) {
		hweigh->SetBinContent( ii, target->GetBinContent( target->FindBin( hweigh->GetBinCenter(ii) ) ) );
	    }
	    hweigh->Divide(hweigh, gen_pu, 1., 1./gen_pu->Integral() );
	} else { 
	    // Normalize weights such that the total cross section is unchanged
	    TH1 * eff = (TH1*)hweigh->Clone("eff");
	    eff->Multiply(gen_pu);
	    hweigh->Scale( gen_pu->Integral() / eff->Integral()  );
	    delete eff;
	}
        weights[typid].clear();
        for( int ii=1; ii<hweigh->GetNbinsX(); ++ii ) {
            weights[typid].push_back(hweigh->GetBinContent(ii)); 
        }
    } 
    std::cout << "pile-up weights: ["<<typid<<"]";
    std::copy(weights[typid].begin(), weights[typid].end(), std::ostream_iterator<double>(std::cout,","));
    std::cout << std::endl;
}

// ----------------------------------------------------------------------------------------------------
float PhotonAnalysis::getPuWeight(int n_pu, int sample_type, bool warnMe)
{
    if ( sample_type !=0 && puHist != "") {
        bool hasSpecificWeight = weights.find( sample_type ) != weights.end() ; 
        if( sample_type < 0 && !hasSpecificWeight && warnMe ) {
            std::cerr  << "WARNING no pu weights specific for sample " << sample_type << std::endl;
        }
        std::vector<double> & puweights = hasSpecificWeight ? weights[ sample_type ] : weights[0]; 
        if(n_pu<puweights.size()){
            return puweights[n_pu]; 
        }    
        else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
            cout <<"n_pu ("<< n_pu<<") too big ("<<puweights.size()<<") ["<< sample_type <<"], event will not be reweighted for pileup"<<endl;
        }
    }
    return 1.;
} 

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::applyGenLevelSmearings(double & genLevWeight, const TLorentzVector & gP4, int npu, int sample_type, BaseGenLevelSmearer * sys, float syst_shift)
{
    static int nwarnings=10;
    for(std::vector<BaseGenLevelSmearer*>::iterator si=genLevelSmearers_.begin(); si!=genLevelSmearers_.end(); si++){
	float genWeight=1;
	if( sys != 0 && *si == *sys ) { 
	    (*si)->smearEvent(genWeight, gP4, npu, sample_type, syst_shift );
	} else {
	    (*si)->smearEvent(genWeight, gP4, npu, sample_type, 0. );
	}
	if( genWeight < 0. ) {
	    if( syst_shift == 0. ) {
		std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
		assert(0);
	    } else { 
		if( nwarnings-- > 0 ) {
		    std::cout <<  "WARNING: negative during systematic scan in " << (*si)->name() << std::endl;
		}
		genWeight = 0.;
	    }
	}
	genLevWeight*=genWeight;
    }
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::applySinglePhotonSmearings(std::vector<float> & smeared_pho_energy, std::vector<float> & smeared_pho_r9, std::vector<float> & smeared_pho_weight,
						int cur_type, const LoopAll & l, const float * energyCorrected, const float * energyCorrectedError,
						BaseSmearer * sys, float syst_shift
						
    )
{
    static int nwarnings = 10;
    photonInfoCollection.clear();
    smeared_pho_energy.resize(l.pho_n,0.);
    smeared_pho_r9.resize(l.pho_n,0.);
    smeared_pho_weight.resize(l.pho_n,0.);
    for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
	
        std::vector<std::vector<bool> > p;
        PhotonReducedInfo phoInfo (
	    *((TVector3*)     l.sc_xyz->At(l.pho_scind[ipho])), 
	    ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), 
	    energyCorrected[ipho],
	    l.pho_isEB[ipho], l.pho_r9[ipho],
	    true, // WARNING  setting pass photon ID flag for all photons. This is safe as long as only selected photons are used
	    (energyCorrectedError!=0?energyCorrectedError[ipho]:0)
	    );
	
	int ieta, iphi;
	l.getIetaIPhi(ipho,ieta,iphi);
	phoInfo.addSmearingSeed( (unsigned int)l.sc_raw[l.pho_scind[ipho]] + abs(ieta) + abs(iphi) + l.run + l.event + l.lumis ); 
	phoInfo.setSphericalPhoton(l.CheckSphericalPhoton(ieta,iphi));

	// FIXME add seed to syst smearings
	
        float pweight = 1.;
        // smear MC. But apply energy corrections and scale adjustement to data 
        if( cur_type != 0 && doMCSmearing ) {
            for(std::vector<BaseSmearer *>::iterator si=photonSmearers_.begin(); si!= photonSmearers_.end(); ++si ) {
                float sweight = 1.;
		if( sys != 0 && *si == *sys ) {
		    // move the smearer under study by syst_shift
		    (*si)->smearPhoton(phoInfo,sweight,l.run,syst_shift);
		} else {
		    // for the other use the nominal points
		    (*si)->smearPhoton(phoInfo,sweight,l.run,0.);
		}
                if( sweight < 0. ) {
		    if( syst_shift == 0. ) {
			std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
			assert(0);
		    } else { 
			if( nwarnings-- > 0 ) {
			    std::cout <<  "WARNING: negative during systematic scan in " << (*si)->name() << std::endl;
			}
			sweight = 0.;
		    }
                }
                pweight *= sweight;
            }
        } else if( cur_type == 0 ) {
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
        photonInfoCollection.push_back(phoInfo);
    }
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::fillDiphoton(TLorentzVector & lead_p4, TLorentzVector & sublead_p4, TLorentzVector & Higgs,
				  float & lead_r9, float & sublead_r9, TVector3 *& vtx, const float * energy,
				  const LoopAll & l, int diphoton_id)
{
    lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], energy);
    sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], energy);
    lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
    sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
    Higgs = lead_p4 + sublead_p4;    
    vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::applyDiPhotonSmearings(TLorentzVector & Higgs, TVector3 & vtx, int category, int cur_type, const TVector3 & truevtx, 
					    float & evweight, float & idmva1, float & idmva2,
					    BaseDiPhotonSmearer * sys, float syst_shift)
{
    static int nwarnings=10;
    float pth = Higgs.Pt();
    for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
        float rewei=1.;
	if( sys != 0 && *si == *sys ) { 
	    (*si)->smearDiPhoton( Higgs, vtx, rewei, category, cur_type, truevtx, idmva1, idmva2, syst_shift );
	} else {
	    (*si)->smearDiPhoton( Higgs, vtx, rewei, category, cur_type, truevtx, idmva1, idmva2, 0. );
	}
        if( rewei < 0. ) {
	    if( syst_shift == 0. ) {
		std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
		assert(0);
	    } else { 
		if( nwarnings-- > 0 ) {
		    std::cout <<  "WARNING: negative during systematic scan in " << (*si)->name() << std::endl;
		}
		rewei = 0.;
	    }
        }
        evweight *= rewei;
    }
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Init(LoopAll& l) 
{
    if(PADEBUG) 
        cout << "InitRealPhotonAnalysis START"<<endl;

    if(energyCorrectionMethod=="DaunceyAndKenzie"){
        energyCorrected     = (l.pho_residCorrEnergy);
        energyCorrectedError= (l.pho_residCorrResn);
    }else if(energyCorrectionMethod=="Bendavid"){
        energyCorrected     = (l.pho_regr_energy);
        energyCorrectedError= (l.pho_regr_energyerr);
    }else if(energyCorrectionMethod=="BendavidOTF"){
        energyCorrected     = (l.pho_regr_energy_otf);
        energyCorrectedError= (l.pho_regr_energyerr_otf);

        //  }else if(energyCorrectionMethod=="PFRegression"){
    }else{
        assert(doEcorrectionSmear==false);
    }
    if (doEcorrectionSmear) std::cout << "using energy correction type: " << energyCorrectionMethod << std::endl;
    else                    std::cout << "NOT using energy correction (sbattogiu)"<< std::endl;

    if( vtxVarNames.empty() ) {
        vtxVarNames.push_back("ptbal"), vtxVarNames.push_back("ptasym"), vtxVarNames.push_back("logsumpt2");
    }
    
    /// // trigger

    // /cdaq/physics/Run2011/5e32/v4.2/HLT/V2
    triggerSelections.push_back(TriggerSelection(160404,161176));
    triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");

    // /cdaq/physics/Run2011/5e32/v6.1/HLT/V1
    triggerSelections.push_back(TriggerSelection(161216,165633));
    triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon20_R9Id_Photon18_R9Id_v");

    // /cdaq/physics/Run2011/1e33/v2.3/HLT/V1
    triggerSelections.push_back(TriggerSelection(165970,166967));
    triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_R9Id_v");

    // /cdaq/physics/Run2011/1e33/v2.3/HLT/V1
    triggerSelections.push_back(TriggerSelection(167039,173198));
    triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_R9Id_v");

    // /cdaq/physics/Run2011/1e33/v2.3/HLT/V1
    triggerSelections.push_back(TriggerSelection(165970,166967));
    triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_R9Id_v");

    // /cdaq/physics/Run2011/1e33/v2.3/HLT/V3
    triggerSelections.push_back(TriggerSelection(167039,173198));
    triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_R9Id_v");

    // /cdaq/physics/Run2011/3e33/v1.1/HLT/V1
    triggerSelections.push_back(TriggerSelection(173236,178380));
    triggerSelections.back().addpath("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_v");
    triggerSelections.back().addpath("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_CaloIdXL_IsoXL_v");
    triggerSelections.back().addpath("HLT_Photon26_R9Id_Photon18_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_R9Id_v");

    // /cdaq/physics/Run2011/5e33/v1.4/HLT/V3
    triggerSelections.push_back(TriggerSelection(178420,-1));
    triggerSelections.back().addpath("HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_Mass60_v");
    triggerSelections.back().addpath("HLT_Photon26_CaloIdXL_IsoXL_Photon18_R9IdT_Mass60_v");
    triggerSelections.back().addpath("HLT_Photon26_R9IdT_Photon18_CaloIdXL_IsoXL_Mass60_v");
    triggerSelections.back().addpath("HLT_Photon26_R9IdT_Photon18_R9IdT_Mass60_v");
    triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v");
    triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
    triggerSelections.back().addpath("HLT_Photon36_R9Id_Photon22_R9Id_v");

    // n-1 plots for VBF tag 2011 
    l.SetCutVariables("cut_VBFLeadJPt",       &myVBFLeadJPt);
    l.SetCutVariables("cut_VBFSubJPt",        &myVBFSubJPt);
    l.SetCutVariables("cut_VBF_Mjj",          &myVBF_Mjj);
    l.SetCutVariables("cut_VBF_dEta",         &myVBFdEta);
    l.SetCutVariables("cut_VBF_Zep",          &myVBFZep);
    l.SetCutVariables("cut_VBF_dPhi",         &myVBFdPhi);
    l.SetCutVariables("cut_VBF_Mgg0",         &myVBF_Mgg);
    l.SetCutVariables("cut_VBF_Mgg2",         &myVBF_Mgg);
    l.SetCutVariables("cut_VBF_Mgg4",         &myVBF_Mgg);
    l.SetCutVariables("cut_VBF_Mgg10",        &myVBF_Mgg);
    l.SetCutVariables("cut_VBF_Mgg4_100_180",        &myVBF_Mgg);
    l.SetCutVariables("cut_VBF_Mgg2_100_180",        &myVBF_Mgg);

    // n-1 plots for VH hadronic tag 2011
    l.SetCutVariables("cut_VHhadLeadJPt",      &myVHhadLeadJPt);
    l.SetCutVariables("cut_VHhadSubJPt",       &myVHhadSubJPt);
    l.SetCutVariables("cut_VHhad_Mjj",         &myVHhad_Mjj);
    l.SetCutVariables("cut_VHhad_dEta",        &myVHhaddEta);
    l.SetCutVariables("cut_VHhad_Zep",         &myVHhadZep);
    l.SetCutVariables("cut_VHhad_dPhi",        &myVHhaddPhi);
    l.SetCutVariables("cut_VHhad_Mgg0",        &myVHhad_Mgg);
    l.SetCutVariables("cut_VHhad_Mgg2",        &myVHhad_Mgg);
    l.SetCutVariables("cut_VHhad_Mgg4",        &myVHhad_Mgg);
    l.SetCutVariables("cut_VHhad_Mgg10",        &myVHhad_Mgg);
    l.SetCutVariables("cut_VHhad_Mgg2_100_160",        &myVHhad_Mgg);
    l.SetCutVariables("cut_VHhad_Mgg4_100_160",        &myVHhad_Mgg);

    // n-1 plot for ClassicCats
    l.SetCutVariables("cutnm1hir9EB_r9",             &sublead_r9);
    l.SetCutVariables("cutnm1hir9EB_isoOverEt",      &sublead_isoOverEt);
    l.SetCutVariables("cutnm1hir9EB_badisoOverEt",   &sublead_badisoOverEt);
    l.SetCutVariables("cutnm1hir9EB_trkisooet",      &sublead_trkisooet);
    l.SetCutVariables("cutnm1hir9EB_sieie",          &sublead_sieie);
    l.SetCutVariables("cutnm1hir9EB_drtotk",         &sublead_drtotk);
    l.SetCutVariables("cutnm1hir9EB_hovere",         &sublead_hovere);
    l.SetCutVariables("cutnm1hir9EB_Mgg",            &sublead_mgg);

    l.SetCutVariables("cutnm1lor9EB_r9",             &sublead_r9);
    l.SetCutVariables("cutnm1lor9EB_isoOverEt",      &sublead_isoOverEt);
    l.SetCutVariables("cutnm1lor9EB_badisoOverEt",   &sublead_badisoOverEt);
    l.SetCutVariables("cutnm1lor9EB_trkisooet",      &sublead_trkisooet);
    l.SetCutVariables("cutnm1lor9EB_sieie",          &sublead_sieie);
    l.SetCutVariables("cutnm1lor9EB_drtotk",         &sublead_drtotk);
    l.SetCutVariables("cutnm1lor9EB_hovere",         &sublead_hovere);
    l.SetCutVariables("cutnm1lor9EB_Mgg",            &sublead_mgg);

    l.SetCutVariables("cutnm1hir9EE_r9",             &sublead_r9);
    l.SetCutVariables("cutnm1hir9EE_isoOverEt",      &sublead_isoOverEt);
    l.SetCutVariables("cutnm1hir9EE_badisoOverEt",   &sublead_badisoOverEt);
    l.SetCutVariables("cutnm1hir9EE_trkisooet",      &sublead_trkisooet);
    l.SetCutVariables("cutnm1hir9EE_sieie",          &sublead_sieie);
    l.SetCutVariables("cutnm1hir9EE_drtotk",         &sublead_drtotk);
    l.SetCutVariables("cutnm1hir9EE_hovere",         &sublead_hovere);
    l.SetCutVariables("cutnm1hir9EE_Mgg",            &sublead_mgg);

    l.SetCutVariables("cutnm1lor9EE_r9",             &sublead_r9);
    l.SetCutVariables("cutnm1lor9EE_isoOverEt",      &sublead_isoOverEt);
    l.SetCutVariables("cutnm1lor9EE_badisoOverEt",   &sublead_badisoOverEt);
    l.SetCutVariables("cutnm1lor9EE_trkisooet",      &sublead_trkisooet);
    l.SetCutVariables("cutnm1lor9EE_sieie",          &sublead_sieie);
    l.SetCutVariables("cutnm1lor9EE_drtotk",         &sublead_drtotk);
    l.SetCutVariables("cutnm1lor9EE_hovere",         &sublead_hovere);
    l.SetCutVariables("cutnm1lor9EE_Mgg",            &sublead_mgg);


    // CiC initialization
    // FIXME should move this to GeneralFunctions
    l.runCiC = true;
    const int phoNCUTS = LoopAll::phoNCUTS;
    const int phoCiC6NCATEGORIES = LoopAll::phoCiC6NCATEGORIES;
    const int phoCiC4NCATEGORIES = LoopAll::phoCiC4NCATEGORIES;
    const int phoNCUTLEVELS = LoopAll::phoNCUTLEVELS;

    for(int iLevel=0; iLevel<phoNCUTLEVELS; ++iLevel) {
        float cic6_cuts_lead[phoNCUTS][phoCiC6NCATEGORIES];
        float cic6_cuts_sublead[phoNCUTS][phoCiC6NCATEGORIES];
        float cic4_cuts_lead[phoNCUTS][phoCiC4NCATEGORIES];
        float cic4_cuts_sublead[phoNCUTS][phoCiC4NCATEGORIES];

        // get the cut values for the current photon CIC level 
        l.SetPhotonCutsInCategories((LoopAll::phoCiCIDLevel)iLevel, &cic6_cuts_lead[0][0], &cic6_cuts_sublead[0][0], &cic4_cuts_lead[0][0], &cic4_cuts_sublead[0][0] );
        
        // rearrange the returned values to arrays with more meaningful names
        float * cic6_cuts_arrays_lead[phoNCUTS] = {
            &l.cic6_cut_lead_isosumoet[0][0], 
            &l.cic6_cut_lead_isosumoetbad[0][0], 
            &l.cic6_cut_lead_trkisooet[0][0], 
            &l.cic6_cut_lead_sieie[0][0],
            &l.cic6_cut_lead_hovere[0][0],
            &l.cic6_cut_lead_r9[0][0], 
            &l.cic6_cut_lead_drtotk_25_99[0][0], 
            &l.cic6_cut_lead_pixel[0][0] 
        };
        
        float * cic6_cuts_arrays_sublead[phoNCUTS] = {
            &l.cic6_cut_sublead_isosumoet[0][0], 
            &l.cic6_cut_sublead_isosumoetbad[0][0], 
            &l.cic6_cut_sublead_trkisooet[0][0], 
            &l.cic6_cut_sublead_sieie[0][0], 
            &l.cic6_cut_sublead_hovere[0][0], 
            &l.cic6_cut_sublead_r9[0][0],
            &l.cic6_cut_sublead_drtotk_25_99[0][0], 
            &l.cic6_cut_sublead_pixel[0][0]
        };

        float * cic4_cuts_arrays_lead[phoNCUTS] = {
            &l.cic4_cut_lead_isosumoet[0][0], 
            &l.cic4_cut_lead_isosumoetbad[0][0], 
            &l.cic4_cut_lead_trkisooet[0][0], 
            &l.cic4_cut_lead_sieie[0][0],
            &l.cic4_cut_lead_hovere[0][0], 
            &l.cic4_cut_lead_r9[0][0], 
            &l.cic4_cut_lead_drtotk_25_99[0][0], 
            &l.cic4_cut_lead_pixel[0][0] 
        };
        
        float * cic4_cuts_arrays_sublead[phoNCUTS] = {
            &l.cic4_cut_sublead_isosumoet[0][0], 
            &l.cic4_cut_sublead_isosumoetbad[0][0], 
            &l.cic4_cut_sublead_trkisooet[0][0], 
            &l.cic4_cut_sublead_sieie[0][0], 
            &l.cic4_cut_sublead_hovere[0][0], 
            &l.cic4_cut_sublead_r9[0][0],
            &l.cic4_cut_sublead_drtotk_25_99[0][0], 
            &l.cic4_cut_sublead_pixel[0][0]
        };

        for(int iCut=0; iCut<phoNCUTS; ++iCut) {
            for(int iCat=0; iCat<phoCiC6NCATEGORIES; ++iCat) {
                cic6_cuts_arrays_lead[iCut][iLevel*phoCiC6NCATEGORIES+iCat] = cic6_cuts_lead[iCut][iCat];
                cic6_cuts_arrays_sublead[iCut][iLevel*phoCiC6NCATEGORIES+iCat] = cic6_cuts_sublead[iCut][iCat];
            }
            for(int iCat=0; iCat<phoCiC4NCATEGORIES; ++iCat) {
                cic4_cuts_arrays_lead[iCut][iLevel*phoCiC4NCATEGORIES+iCat] = cic4_cuts_lead[iCut][iCat];
                cic4_cuts_arrays_sublead[iCut][iLevel*phoCiC4NCATEGORIES+iCat] = cic4_cuts_sublead[iCut][iCat];
            }
        }
    } // end of loop over all photon cut levels
    
    //--------------------

    if( tmvaPerVtxWeights != ""  ) {
        tmvaPerVtxVariables_.push_back("ptbal"), tmvaPerVtxVariables_.push_back("ptasym"), tmvaPerVtxVariables_.push_back("logsumpt2");
        if( addConversionToMva ) {
            tmvaPerVtxVariables_.push_back("limPullToConv");
            tmvaPerVtxVariables_.push_back("nConv");
        }
        tmvaPerVtxReader_ = new TMVA::Reader( "!Color:!Silent" );
        HggVertexAnalyzer::bookVariables( *tmvaPerVtxReader_, tmvaPerVtxVariables_ );
        tmvaPerVtxReader_->BookMVA( tmvaPerVtxMethod, tmvaPerVtxWeights );
    } else {
        tmvaPerVtxReader_ = 0;
    }
    if( tmvaPerEvtWeights != "" ) {
        tmvaPerEvtReader_ = new TMVA::Reader( "!Color:!Silent" );
        HggVertexAnalyzer::bookPerEventVariables( *tmvaPerEvtReader_ );
        tmvaPerEvtReader_->BookMVA( tmvaPerEvtMethod, tmvaPerEvtWeights );
    } else {
        tmvaPerEvtReader_ = 0;
    }
    assert( !mvaVertexSelection || tmvaPerVtxReader_ != 0 );

    eSmearDataPars.categoryType = "2CatR9_EBEBm4EE";
    eSmearDataPars.byRun = true;
    //eSmearDataPars.n_categories = 4; //GF
    eSmearDataPars.n_categories = 6; //GF
    std::cerr << "Reading energy scale offsets " << scale_offset_file << std::endl;
    readEnergyScaleOffsets(scale_offset_file, eSmearDataPars.scale_offset_byrun, eSmearDataPars.photon_categories);
    // if the scale offset file defines the categories set the category type to automatic
    if( ! eSmearDataPars.photon_categories.empty() ) {
        eSmearDataPars.categoryType = "Automagic";
        eSmearDataPars.n_categories = -1;
    }
    // E resolution smearing NOT applied to data 
    eSmearDataPars.smearing_sigma["EBHighR9"] = 0.;
    eSmearDataPars.smearing_sigma["EBLowR9"]  = 0.;
    eSmearDataPars.smearing_sigma["EBm4HighR9"] = 0.;
    eSmearDataPars.smearing_sigma["EBm4LowR9"]  = 0.;
    eSmearDataPars.smearing_sigma["EEHighR9"] = 0.;
    eSmearDataPars.smearing_sigma["EELowR9"]  = 0.;
    // E resolution systematics NOT applied to data 
    eSmearDataPars.smearing_sigma_error["EBHighR9"] = 0.;
    eSmearDataPars.smearing_sigma_error["EBLowR9"]  = 0.;
    eSmearDataPars.smearing_sigma_error["EBm4HighR9"] = 0.;
    eSmearDataPars.smearing_sigma_error["EBm4LowR9"]  = 0.;
    eSmearDataPars.smearing_sigma_error["EEHighR9"] = 0.;
    eSmearDataPars.smearing_sigma_error["EELowR9"]  = 0.;
    
    // energy scale corrections to Data
    eScaleDataSmearer = new EnergySmearer( eSmearDataPars );
    eScaleDataSmearer->name("E_scale_data");
    eScaleDataSmearer->doEnergy(true);
    eScaleDataSmearer->scaleOrSmear(true);
    
    if( scale_offset_error_file.empty() ) {
        //eSmearPars.categoryType = "2CatR9_EBEE"; //GF
        eSmearPars.categoryType = "2CatR9_EBEBm4EE";
        eSmearPars.byRun = false;
        //eSmearPars.n_categories = 4; //GF
        eSmearPars.n_categories = 6;
        // E scale is shifted for data, NOT for MC 
        eSmearPars.scale_offset["EBHighR9"] = 0.;
        eSmearPars.scale_offset["EBLowR9"]  = 0.;
        eSmearPars.scale_offset["EBm4HighR9"] = 0.;
        eSmearPars.scale_offset["EBm4LowR9"]  = 0.;
        eSmearPars.scale_offset["EEHighR9"] = 0.;
        eSmearPars.scale_offset["EELowR9"]  = 0.;
        // E scale systematics are applied to MC, NOT to data
        eSmearPars.scale_offset_error["EBHighR9"] = scale_offset_error_EBHighR9;
        eSmearPars.scale_offset_error["EBLowR9"]  = scale_offset_error_EBLowR9;
        eSmearPars.scale_offset_error["EBm4HighR9"] = scale_offset_error_EBHighR9;
        eSmearPars.scale_offset_error["EBm4LowR9"]  = scale_offset_error_EBLowR9;
        eSmearPars.scale_offset_error["EEHighR9"] = scale_offset_error_EEHighR9;
        eSmearPars.scale_offset_error["EELowR9"]  = scale_offset_error_EELowR9;
        // E resolution smearing applied to MC 
        eSmearPars.smearing_sigma["EBHighR9"] = smearing_sigma_EBHighR9;
        eSmearPars.smearing_sigma["EBLowR9"]  = smearing_sigma_EBLowR9;
        eSmearPars.smearing_sigma["EBm4HighR9"] = smearing_sigma_EBm4HighR9;
        eSmearPars.smearing_sigma["EBm4LowR9"]  = smearing_sigma_EBm4LowR9;
        eSmearPars.smearing_sigma["EEHighR9"] = smearing_sigma_EEHighR9;
        eSmearPars.smearing_sigma["EELowR9"]  = smearing_sigma_EELowR9;
        // E resolution systematics applied to MC 
        eSmearPars.smearing_sigma_error["EBHighR9"] = smearing_sigma_error_EBHighR9;
        eSmearPars.smearing_sigma_error["EBLowR9"]  = smearing_sigma_error_EBLowR9;
        eSmearPars.smearing_sigma_error["EBm4HighR9"] = smearing_sigma_error_EBm4HighR9;
        eSmearPars.smearing_sigma_error["EBm4LowR9"]  = smearing_sigma_error_EBm4LowR9;
        eSmearPars.smearing_sigma_error["EEHighR9"] = smearing_sigma_error_EEHighR9;
        eSmearPars.smearing_sigma_error["EELowR9"]  = smearing_sigma_error_EELowR9;
        // error on photon corrections set to a fraction of the correction itself; number below is tentative (GF: push it to .dat)  
        eSmearPars.corrRelErr  = 0.5;
    } else {
        // Read energy scale errors and energy smaerings from dat files
        assert( ! scale_offset_error_file.empty() && ! smearing_file.empty() );

        // Use the same format used for the run-dependent energy corrections
        EnergySmearer::energySmearingParameters::eScaleVector tmp_scale_offset, tmp_smearing;
        EnergySmearer::energySmearingParameters::phoCatVector tmp_scale_cat, tmp_smearing_cat;
        readEnergyScaleOffsets(scale_offset_error_file, tmp_scale_offset, tmp_scale_cat,false);
        readEnergyScaleOffsets(smearing_file, tmp_smearing, tmp_smearing_cat,false);

        // make sure that the scale correction and smearing info is as expected
        assert( tmp_scale_offset.size() == 1); assert( tmp_smearing.size() == 1 );
        assert( ! tmp_smearing_cat.empty() );
        /// assert( tmp_smearing_cat == tmp_scale_cat );

        // copy the read info to the smarer parameters
        eSmearPars.categoryType = "Automagic";
        eSmearPars.byRun = false;
        eSmearPars.n_categories = tmp_smearing_cat.size();
        eSmearPars.photon_categories = tmp_smearing_cat;
        
        eSmearPars.scale_offset = tmp_scale_offset[0].scale_offset;
        eSmearPars.scale_offset_error = tmp_scale_offset[0].scale_offset_error;
        
        eSmearPars.smearing_sigma = tmp_smearing[0].scale_offset;
        eSmearPars.smearing_sigma_error = tmp_smearing[0].scale_offset_error;
    }
    
    // energy scale systematics to MC
    eScaleSmearer = new EnergySmearer( eSmearPars );
    eScaleSmearer->name("E_scale");
    eScaleSmearer->doEnergy(true);
    eScaleSmearer->scaleOrSmear(true);

    if( doEcorrectionSmear ) {
        eCorrSmearer = new EnergySmearer( eSmearPars );
        eCorrSmearer->name("E_corr");
        // activating pho corrections to this instance of EnergySmearer, implies that it won't touch Escale and Eresolution
        eCorrSmearer->doCorrections(true); 
    }
    
    if (l.typerun == 2 || l.typerun == 1) {
    }
    // MassResolution 
    massResolutionCalculator = new MassResolution();    
    /* -------------------------------------------------------------------------------------------
       Pileup Reweighting
       https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
       ----------------------------------------------------------------------------------------------  */
    if (puHist != "") {
        if(PADEBUG) 
            cout << "Opening PU file"<<endl;
        TFile* puFile = TFile::Open( puHist );
        if (puFile) {
	    TH1 * target = 0;
	    
	    if( puTarget != "" ) {
		TFile * puTargetFile = TFile::Open( puTarget ); 
		assert( puTargetFile != 0 );
		target = (TH1*)puTargetFile->Get("pileup");
		target->Scale( 1. / target->Integral() );
	    }
	    
            if( puMap != "" ) {
                loadPuMap(puMap, puFile, target); 
            } else {
                loadPuWeights(0, puFile, target);
            }
            puFile->Close();
        }
        else {
            cout<<"Error opening " <<puHist<<" pileup reweighting histogram, using 1.0"<<endl; 
            weights[0].resize(50);
            for (unsigned int i=0; i<weights[0].size(); i++) weights[0][i] = 1.0;
        }
        if(PADEBUG) 
            cout << "Opening PU file END"<<endl;
    } 

    // Load up instances of PhotonFix for local coordinate calculations
    /*  
        PhotonFix::initialise("4_2",photonFixDat);  
        std::cout << "Regression corrections from -> " << regressionFile.c_str() << std::endl;
        fgbr = new TFile(regressionFile.c_str(),"READ");
        fReadereb = (GBRForest*)fgbr->Get("EBCorrection");
        fReaderebvariance = (GBRForest*)fgbr->Get("EBUncertainty");  
        fReaderee = (GBRForest*)fgbr->Get("EECorrection");
        fReadereevariance = (GBRForest*)fgbr->Get("EEUncertainty");      
        fgbr->Close();
    */
    
    // -------------------------------------------------------------------- 
    if(PADEBUG) 
        cout << "InitRealPhotonAnalysis END"<<endl;

    // FIXME book of additional variables

}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
    if(PADEBUG) 
        cout << "Analysis START"<<endl;
    pho_presel.clear();

    //remove process ID 18 from gamma+jet to avoid double counting with born+box
    if (l.itype[l.current]==3 && l.process_id==18) return;

    //apply pileup reweighting
    unsigned int n_pu = l.pu_n;
    float weight =1.;
    if (l.itype[l.current] !=0 && puHist != "") {
        std::vector<double> & puweights = weights.find( l.itype[l.current] ) != weights.end() ? weights[ l.itype[l.current] ] : weights[0]; 
        if(n_pu<puweights.size()){
            weight *= puweights[n_pu];
        }    
        else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
            cout<<"n_pu ("<< n_pu<<") too big ("<<puweights.size()<<"), event will not be reweighted for pileup"<<endl;
        }
    }

    if( pho_presel.empty() ) { 
        PreselectPhotons(l,jentry);
    }

    if( pho_acc.size() < 2 || pho_et[ pho_acc[0] ] < presel_scet1 ) return;

    int leadLevel=LoopAll::phoSUPERTIGHT, subLevel=LoopAll::phoSUPERTIGHT;

    //Fill histograms to use as denominator (pre-selection only) and numerator (selection applied)
    //for photon ID efficiency calculation.  To avoid ambiguities concerning vertex choice, use only 
    //events with one diphoton pair (close to 100% of signal events)
    if (l.dipho_n==1) {

        int ivtx = l.dipho_vtxind[0];
        int lead = l.dipho_leadind[0];
        int sublead = l.dipho_subleadind[0];

        TLorentzVector lead_p4 = l.get_pho_p4(lead,ivtx,&corrected_pho_energy[0]); 
        TLorentzVector sublead_p4 = l.get_pho_p4(sublead,ivtx,&corrected_pho_energy[0]); 
        float leadEta = ((TVector3 *)l.sc_xyz->At(l.pho_scind[lead]))->Eta();
        float subleadEta = ((TVector3 *)l.sc_xyz->At(l.pho_scind[sublead]))->Eta();

        //apply pre-selection
        bool passpresel = true;
        if(lead_p4.Pt() < 40. || sublead_p4.Pt() < 30. || 
           fabs(leadEta) > 2.5 || fabs(subleadEta) > 2.5 || 
           ( fabs(leadEta) > 1.4442 && fabs(leadEta) < 1.566 ) || ( fabs(subleadEta) > 1.4442 && fabs(subleadEta) < 1.566 ))
            passpresel = false;
        if (lead != sublead && passpresel) {

            int leadpho_category = l.PhotonCategory(lead, 2, 2);
            int subleadpho_category = l.PhotonCategory(sublead, 2, 2);

            //Fill eta and pt distributions after pre-selection only (efficiency denominator)
            l.FillHist("pho1_pt_presel",0,lead_p4.Pt(), weight);
            l.FillHist("pho2_pt_presel",0,sublead_p4.Pt(), weight);
            l.FillHist("pho1_eta_presel",0,leadEta, weight);
            l.FillHist("pho2_eta_presel",0,subleadEta, weight);

            l.FillHist("pho1_pt_presel",leadpho_category+1,lead_p4.Pt(), weight);
            l.FillHist("pho2_pt_presel",subleadpho_category+1,sublead_p4.Pt(), weight);
            l.FillHist("pho1_eta_presel",leadpho_category+1,leadEta, weight);
            l.FillHist("pho2_eta_presel",subleadpho_category+1,subleadEta, weight);

            //Apply single photon CiC selection and fill eta and pt distributions (efficiency numerator)
            std::vector<std::vector<bool> > ph_passcut;
            if( l.PhotonCiCSelectionLevel(lead, ivtx, ph_passcut, 4, 0, &corrected_pho_energy[0]) >=  (LoopAll::phoCiCIDLevel) leadLevel) {
                l.FillHist("pho1_pt_sel",0,lead_p4.Pt(), weight);
                l.FillHist("pho1_eta_sel",0,leadEta, weight);
                l.FillHist("pho1_pt_sel",leadpho_category+1,lead_p4.Pt(), weight);
                l.FillHist("pho1_eta_sel",leadpho_category+1,leadEta, weight);
            }
            if( l.PhotonCiCSelectionLevel(sublead, ivtx, ph_passcut, 4, 1, &corrected_pho_energy[0]) >=  (LoopAll::phoCiCIDLevel) subLevel ) {
                l.FillHist("pho2_pt_sel",0,sublead_p4.Pt(), weight);
                l.FillHist("pho2_eta_sel",0,subleadEta, weight);
                l.FillHist("pho2_pt_sel",subleadpho_category+1,sublead_p4.Pt(), weight);
                l.FillHist("pho2_eta_sel",subleadpho_category+1,subleadEta, weight);
            }
        }
    }

    //Apply diphoton CiC selection
    int dipho_id = l.DiphotonCiCSelection((LoopAll::phoCiCIDLevel) leadLevel, (LoopAll::phoCiCIDLevel) subLevel, 40., 30.0, 4, applyPtoverM, &corrected_pho_energy[0]);

    if (dipho_id > -1){
        std::pair<int,int> diphoton_index = std::make_pair(l.dipho_leadind[dipho_id],l.dipho_subleadind[dipho_id]);
        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[dipho_id], l.dipho_vtxind[dipho_id], &corrected_pho_energy[0]);
        TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[dipho_id], l.dipho_vtxind[dipho_id], &corrected_pho_energy[0]);
        TLorentzVector Higgs = lead_p4 + sublead_p4;    
      
        //// TLorentzVector *lead_p4 = (TLorentzVector*)l.pho_p4->At(diphoton_index.first);
        //// TLorentzVector *sublead_p4 = (TLorentzVector*)l.pho_p4->At(diphoton_index.second);
        //// TLorentzVector Higgs = *lead_p4 + *sublead_p4;

        // Calculate the Rest Frame Quantities:
        TLorentzVector *rest_lead_p4 = (TLorentzVector*) lead_p4.Clone();
        TLorentzVector *rest_sublead_p4 = (TLorentzVector*) sublead_p4.Clone();
        TVector3 higgsBoostVector = Higgs.BoostVector();
        rest_lead_p4->Boost(higgsBoostVector);        // lead photon in Higgs Decay Frame
        rest_sublead_p4->Boost(higgsBoostVector);    // sub-lead photon in Higgs Decay Frame

        // Calculate variables to be filled
        float decayAngle  = l.DeltaPhi(lead_p4.Phi(),sublead_p4.Phi());
        float helicityAngle = fabs(lead_p4.E()-sublead_p4.E())/Higgs.P();
        float maxeta = lead_p4.Eta();
        if (fabs(sublead_p4.Eta())>fabs(lead_p4.Eta())) maxeta = sublead_p4.Eta();
        float HT = lead_p4.Pt() + sublead_p4.Pt();
        float mH = Higgs.M();


        //Fill histograms according to diphoton or single photon category, as appropriate

        int dipho_category = l.DiphotonCategory(diphoton_index.first, diphoton_index.second, Higgs.Pt(), 2, 2, 2);
        int leadpho_category = l.PhotonCategory(diphoton_index.first, 2, 2);
        int subleadpho_category = l.PhotonCategory(diphoton_index.second, 2, 2);
      
    
        //Only fill histograms for QCD and GJet if bkgCat is not promptprompt
        //if ((l.itype[l.current]==1 || l.itype[l.current]==2 || l.itype[l.current]==3) && bkgCat==promptprompt) return;

        l.FillHist("pt",0, Higgs.Pt(), weight);
        l.FillHist("ptOverM",0, Higgs.Pt()/Higgs.M(), weight);
        l.FillHist("eta",0, Higgs.Eta(), weight);
        l.FillHist("decayAngle",0, decayAngle, weight);
        l.FillHist("helicityAngle",0, helicityAngle, weight);
        l.FillHist("pho1_pt",0,lead_p4.Pt(), weight);
        l.FillHist("pho2_pt",0,sublead_p4.Pt(), weight);
        l.FillHist("pho1_ptOverM",0,lead_p4.Pt()/Higgs.M(), weight);
        l.FillHist("pho2_ptOverM",0,sublead_p4.Pt()/Higgs.M(), weight);
        l.FillHist("pho1_eta",0,lead_p4.Eta(), weight);
        l.FillHist("pho2_eta",0,sublead_p4.Eta(), weight);
        l.FillHist("pho_r9",0,l.pho_r9[diphoton_index.first], weight);
        l.FillHist("pho_r9",0,l.pho_r9[diphoton_index.second], weight);
        l.FillHist("maxeta",0,maxeta, weight);
        l.FillHist("ht",0, HT, weight);
        l.FillHist2D("ht_vs_m",0, HT, mH, weight);

        l.FillHist("pt",dipho_category+1, Higgs.Pt(), weight);
        l.FillHist("ptOverM",dipho_category+1, Higgs.Pt()/Higgs.M(), weight);
        l.FillHist("eta",dipho_category+1, Higgs.Eta(), weight);
        l.FillHist("decayAngle",dipho_category+1, decayAngle, weight);
        l.FillHist("helicityAngle",dipho_category+1, helicityAngle, weight);
        l.FillHist("pho1_pt",dipho_category+1,lead_p4.Pt(), weight);
        l.FillHist("pho2_pt",dipho_category+1,sublead_p4.Pt(), weight);
        l.FillHist("pho1_ptOverM",dipho_category+1,lead_p4.Pt()/Higgs.M(), weight);
        l.FillHist("pho2_ptOverM",dipho_category+1,sublead_p4.Pt()/Higgs.M(), weight);
        l.FillHist("pho1_eta",dipho_category+1,lead_p4.Eta(), weight);
        l.FillHist("pho2_eta",dipho_category+1,sublead_p4.Eta(), weight);
        l.FillHist("pho_r9",dipho_category+1,l.pho_r9[diphoton_index.first], weight);
        l.FillHist("pho_r9",dipho_category+1,l.pho_r9[diphoton_index.second], weight);
        l.FillHist("maxeta",dipho_category+1,maxeta, weight);
        l.FillHist("ht",dipho_category+1, HT, weight);
        l.FillHist2D("ht_vs_m",dipho_category+1, mH, HT, weight);

        //Fill separately for low and high mass sidebands and for signal region

        if (mH>100 && mH<110) {
            l.FillHist("pt_mlow",0, Higgs.Pt(), weight);
            l.FillHist("ptOverM_mlow",0, Higgs.Pt()/Higgs.M(), weight);
            l.FillHist("eta_mlow",0, Higgs.Eta(), weight);
            l.FillHist("decayAngle_mlow",0, decayAngle, weight);
            l.FillHist("helicityAngle_mlow",0, helicityAngle, weight);
            l.FillHist("pho1_pt_mlow",0,lead_p4.Pt(), weight);
            l.FillHist("pho2_pt_mlow",0,sublead_p4.Pt(), weight);
            //l.FillHist("pho1_pt_mlow",0,lead_p4.Pt()*R_low, weight);
            //l.FillHist("pho2_pt_mlow",0,sublead_p4.Pt()*R_low, weight);
            l.FillHist("pho1_ptOverM_mlow",0,lead_p4.Pt()/Higgs.M(), weight);
            l.FillHist("pho2_ptOverM_mlow",0,sublead_p4.Pt()/Higgs.M(), weight);
            l.FillHist("pho1_eta_mlow",0,lead_p4.Eta(), weight);
            l.FillHist("pho2_eta_mlow",0,sublead_p4.Eta(), weight);
            l.FillHist("pho_r9_mlow",0,l.pho_r9[diphoton_index.first], weight);
            l.FillHist("pho_r9_mlow",0,l.pho_r9[diphoton_index.second], weight);
            l.FillHist("maxeta_mlow",0,maxeta, weight);
            l.FillHist("ht_mlow",0, HT, weight);

            l.FillHist("pt_mlow",dipho_category+1, Higgs.Pt(), weight);
            l.FillHist("ptOverM_mlow",dipho_category+1, Higgs.Pt()/Higgs.M(), weight);
            l.FillHist("eta_mlow",dipho_category+1, Higgs.Eta(), weight);
            l.FillHist("decayAngle_mlow",dipho_category+1, decayAngle, weight);
            l.FillHist("helicityAngle_mlow",dipho_category+1, helicityAngle, weight);
            l.FillHist("pho1_pt_mlow",dipho_category+1,lead_p4.Pt(), weight);
            l.FillHist("pho2_pt_mlow",dipho_category+1,sublead_p4.Pt(), weight);
            l.FillHist("pho1_ptOverM_mlow",dipho_category+1,lead_p4.Pt()/Higgs.M(), weight);
            l.FillHist("pho2_ptOverM_mlow",dipho_category+1,sublead_p4.Pt()/Higgs.M(), weight);
            l.FillHist("pho1_eta_mlow",dipho_category+1,lead_p4.Eta(), weight);
            l.FillHist("pho2_eta_mlow",dipho_category+1,sublead_p4.Eta(), weight);
            l.FillHist("pho_r9_mlow",dipho_category+1,l.pho_r9[diphoton_index.first], weight);
            l.FillHist("pho_r9_mlow",dipho_category+1,l.pho_r9[diphoton_index.second], weight);
            l.FillHist("maxeta_mlow",dipho_category+1,maxeta, weight);
            l.FillHist("ht_mlow",dipho_category+1, HT, weight);
        } else if (mH>126 && mH<145) {
            l.FillHist("pt_mhigh",0, Higgs.Pt(), weight);
            l.FillHist("ptOverM_mhigh",0, Higgs.Pt()/Higgs.M(), weight);
            l.FillHist("eta_mhigh",0, Higgs.Eta(), weight);
            l.FillHist("decayAngle_mhigh",0, decayAngle, weight);
            l.FillHist("helicityAngle_mhigh",0, helicityAngle, weight);
            l.FillHist("pho1_pt_mhigh",0,lead_p4.Pt(), weight);
            l.FillHist("pho2_pt_mhigh",0,sublead_p4.Pt(), weight);
            //l.FillHist("pho1_pt_mhigh",0,lead_p4.Pt()*R_high, weight);
            //l.FillHist("pho2_pt_mhigh",0,sublead_p4.Pt()*R_high, weight);
            l.FillHist("pho1_ptOverM_mhigh",0,lead_p4.Pt()/Higgs.M(), weight);
            l.FillHist("pho2_ptOverM_mhigh",0,sublead_p4.Pt()/Higgs.M(), weight);
            l.FillHist("pho1_eta_mhigh",0,lead_p4.Eta(), weight);
            l.FillHist("pho2_eta_mhigh",0,sublead_p4.Eta(), weight);
            l.FillHist("pho_r9_mhigh",0,l.pho_r9[diphoton_index.first], weight);
            l.FillHist("pho_r9_mhigh",0,l.pho_r9[diphoton_index.second], weight);
            l.FillHist("maxeta_mhigh",0,maxeta, weight);
            l.FillHist("ht_mhigh",0, HT, weight);

            l.FillHist("pt_mhigh",dipho_category+1, Higgs.Pt(), weight);
            l.FillHist("ptOverM_mhigh",dipho_category+1, Higgs.Pt()/Higgs.M(), weight);
            l.FillHist("eta_mhigh",dipho_category+1, Higgs.Eta(), weight);
            l.FillHist("decayAngle_mhigh",dipho_category+1, decayAngle, weight);
            l.FillHist("helicityAngle_mhigh",dipho_category+1, helicityAngle, weight);
            l.FillHist("pho1_pt_mhigh",dipho_category+1,lead_p4.Pt(), weight);
            l.FillHist("pho2_pt_mhigh",dipho_category+1,sublead_p4.Pt(), weight);
            l.FillHist("pho1_ptOverM_mhigh",dipho_category+1,lead_p4.Pt()/Higgs.M(), weight);
            l.FillHist("pho2_ptOverM_mhigh",dipho_category+1,sublead_p4.Pt()/Higgs.M(), weight);
            l.FillHist("pho1_eta_mhigh",dipho_category+1,lead_p4.Eta(), weight);
            l.FillHist("pho2_eta_mhigh",dipho_category+1,sublead_p4.Eta(), weight);
            l.FillHist("pho_r9_mhigh",dipho_category+1,l.pho_r9[diphoton_index.first], weight);
            l.FillHist("pho_r9_mhigh",dipho_category+1,l.pho_r9[diphoton_index.second], weight);
            l.FillHist("maxeta_mhigh",dipho_category+1,maxeta, weight);
            l.FillHist("ht_mhigh",dipho_category+1, HT, weight);
        } else if (mH>=110 && mH<=126) {
            l.FillHist("pt_msig",0, Higgs.Pt(), weight);
            l.FillHist("ptOverM_msig",0, Higgs.Pt()/Higgs.M(), weight);
            l.FillHist("eta_msig",0, Higgs.Eta(), weight);
            l.FillHist("decayAngle_msig",0, decayAngle, weight);
            l.FillHist("helicityAngle_msig",0, helicityAngle, weight);
            l.FillHist("pho1_pt_msig",0,lead_p4.Pt(), weight);
            l.FillHist("pho2_pt_msig",0,sublead_p4.Pt(), weight);
            l.FillHist("pho1_ptOverM_msig",0,lead_p4.Pt()/Higgs.M(), weight);
            l.FillHist("pho2_ptOverM_msig",0,sublead_p4.Pt()/Higgs.M(), weight);
            l.FillHist("pho1_eta_msig",0,lead_p4.Eta(), weight);
            l.FillHist("pho2_eta_msig",0,sublead_p4.Eta(), weight);
            l.FillHist("pho_r9_msig",0,l.pho_r9[diphoton_index.first], weight);
            l.FillHist("pho_r9_msig",0,l.pho_r9[diphoton_index.second], weight);
            l.FillHist("maxeta_msig",0,maxeta, weight);
            l.FillHist("ht_msig",0, HT, weight);

            l.FillHist("pt_msig",dipho_category+1, Higgs.Pt(), weight);
            l.FillHist("ptOverM_msig",dipho_category+1, Higgs.Pt()/Higgs.M(), weight);
            l.FillHist("eta_msig",dipho_category+1, Higgs.Eta(), weight);
            l.FillHist("decayAngle_msig",dipho_category+1, decayAngle, weight);
            l.FillHist("helicityAngle_msig",dipho_category+1, helicityAngle, weight);
            l.FillHist("pho1_pt_msig",dipho_category+1,lead_p4.Pt(), weight);
            l.FillHist("pho2_pt_msig",dipho_category+1,sublead_p4.Pt(), weight);
            l.FillHist("pho1_ptOverM_msig",dipho_category+1,lead_p4.Pt()/Higgs.M(), weight);
            l.FillHist("pho2_ptOverM_msig",dipho_category+1,sublead_p4.Pt()/Higgs.M(), weight);
            l.FillHist("pho1_eta_msig",dipho_category+1,lead_p4.Eta(), weight);
            l.FillHist("pho2_eta_msig",dipho_category+1,sublead_p4.Eta(), weight);
            l.FillHist("pho_r9_msig",dipho_category+1,l.pho_r9[diphoton_index.first], weight);
            l.FillHist("pho_r9_msig",dipho_category+1,l.pho_r9[diphoton_index.second], weight);
            l.FillHist("maxeta_msig",dipho_category+1,maxeta, weight);
            l.FillHist("ht_msig",dipho_category+1, HT, weight);
        }
    }

    if(PADEBUG) 
        cout<<"myFillHistRed END"<<endl;

}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
    vtxAna_.setBranchAdresses(t,"vtx_std_");
    vtxAna_.getBranches(t,"vtx_std_",s);
}


// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::PreselectPhotons(LoopAll& l, int jentry) 
{
    // Photon preselection
    pho_acc.clear();
    pho_presel.clear();
    pho_presel_lead.clear();
    pho_et.clear();
    l.pho_matchingConv->clear();

    // Nominal smearing
    corrected_pho_energy.clear(); corrected_pho_energy.resize(l.pho_n,0.); 
    int cur_type = l.itype[l.current];

    // EDIT - 4 Dec 2011 NWardle Latest Ntuple Production uses New and Correct Regression so no need to calculate on the FLY corrections
    // TEMPORARY FIX TO CALCULATE CORRECTED ENERGIES SINCE REGRESSION WAS NOT STORED IN NTUPLES 
    // The following Fills the arrays with the ON-THE-FLY calculations
    //GetRegressionCorrections(l);  // need to pass LoopAll
    // -------------------------------------------------------------------------------------------//
    for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
        std::vector<std::vector<bool> > p;
        PhotonReducedInfo phoInfo (
                                   *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])),
                                   ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(),
                                   energyCorrected[ipho],
                                   l.pho_isEB[ipho],
                                   l.pho_r9[ipho],
                                   false,
                                   (energyCorrectedError!=0?energyCorrectedError[ipho]:0)
                                   );
        float pweight = 1.;
        float sweight = 1.;
        float eta = fabs(((TVector3 *)l.sc_xyz->At(l.pho_scind[ipho]))->Eta());
        if( doEcorrectionSmear )  { 
            eCorrSmearer->smearPhoton(phoInfo,sweight,l.run,0.); 
        }
        if( cur_type == 0 ) {          // correct energy scale in data
            float ebefore = phoInfo.energy();
            eScaleDataSmearer->smearPhoton(phoInfo,sweight,l.run,0.);
            pweight *= sweight;
            /// std::cerr << std::setprecision(3) << "adjusting energy scale " << " run: "<< l.run << " " << phoInfo.caloPosition().Eta() << " " 
            ///       << phoInfo.r9()<< " before: " << ebefore;
            /// std::cerr << " after: " << phoInfo.energy() << " 1.-ratio " << 100.*(1. - phoInfo.energy()/ebefore) << "%"  << std::endl;
        }
        // apply mc-derived photon corrections, to data and MC alike
        corrected_pho_energy[ipho] = phoInfo.energy();
        l.pho_genmatched[ipho]=GenMatchedPhoton( l, ipho);
    }

    for(int ipho=0; ipho<l.pho_n; ++ipho) {

        // match all photons in the original tree with the conversions from the merged collection and save the indices
        int iConv  =l.matchPhotonToConversion(ipho);
        if ( iConv>=0 )
            (*l.pho_matchingConv).push_back(l.matchPhotonToConversion(ipho));
        else
            (*l.pho_matchingConv).push_back(-1);

        // TLorentzVector * p4 = (TLorentzVector *) l.pho_p4->At(ipho);
        TLorentzVector p4 = l.get_pho_p4(ipho,0,&corrected_pho_energy[0]);
        // float eta  = fabs(((TVector3 *) l.pho_calopos->At(ipho))->Eta());
        float eta = fabs(((TVector3 *)l.sc_xyz->At(l.pho_scind[ipho]))->Eta());
        // photon et wrt 0,0,0
        float et = p4.Pt();
        pho_et.push_back(et);
        /// std::cerr << " " << p4->Pt() << " " << et << " " << eta;
      
        if( p4.Pt() < presel_scet2 || (eta>1.4442 && eta<1.566) || eta>presel_maxeta ) { 
            /// std::cerr << std::endl;
            continue;  
        }
        /// std::cerr << "keeping " << ipho << std::endl;
        pho_acc.push_back(ipho);
      
        bool isEB = l.pho_isEB[ipho];
        float & ecaliso = isEB ? presel_ecaliso_eb : presel_ecaliso_ee;
        float & sieie = isEB ? presel_sieie_eb : presel_sieie_ee;
        if( l.pho_ecalsumetconedr03[ipho] >= ecaliso ||  l.pho_sieie[ipho] >= sieie || l.pho_hoe[ipho] >= presel_hoe ) {
            continue;
        }
          
        //FIXME trigger matching
        pho_presel.push_back(ipho);
    } 

    std::sort(pho_acc.begin(),pho_acc.end(),
              SimpleSorter<float,std::greater<float> >(&pho_et[0]));
    std::sort(pho_presel.begin(),pho_presel.end(),
              SimpleSorter<float,std::greater<float> >(&pho_et[0]));

    if( pho_presel.size() > 1 ) {
        for(size_t ipho=0; ipho<pho_presel.size()-1; ++ipho ) {
            /// assert( ((TLorentzVector *)l.pho_p4->At(pho_presel[ipho]))->Pt() >= ((TLorentzVector *)l.pho_p4->At(pho_presel[ipho+1]))->Pt() );
            assert( pho_et[pho_presel[ipho]] >= pho_et[pho_presel[ipho+1]] );
        }
    }
    if( pho_acc.size()>1 ) {
        for(size_t ipho=0; ipho<pho_acc.size()-1; ++ipho ) {
            /// assert( ((TLorentzVector *)l.pho_p4->At(pho_acc[ipho]))->Pt() >= ((TLorentzVector *)l.pho_p4->At(pho_acc[ipho+1]))->Pt() );
            assert( pho_et[pho_acc[ipho]] >= pho_et[pho_acc[ipho+1]] );
        }
    }
    
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::FillReductionVariables(LoopAll& l, int jentry) 
{
    if(PADEBUG) 
        cout<<"myFillReduceVar START"<<endl;
    
    PreselectPhotons(l,jentry);
        
    if(PADEBUG) 
        cout<<"myFillReduceVar END"<<endl;

}

// ----------------------------------------------------------------------------------------------------
bool PhotonAnalysis::SelectEventsReduction(LoopAll& l, int jentry) 
{

    if(PADEBUG)  cout << " ****************** SelectEventsReduction " << endl;
    // require at least two reconstructed photons to store the event
    if( pho_acc.size() < 2 || l.get_pho_p4( pho_acc[0], 0, &corrected_pho_energy[0] ).Pt() < presel_scet1 ) { return false; }
    
    
    vtxAna_.clear();
    l.vtx_std_ranked_list->clear();
    l.dipho_vtx_std_sel->clear();
    l.vtx_std_ranked_list->clear();
    l.vtx_std_evt_mva->clear();
    l.vtx_std_sel=0;
    float maxSumPt = 0.;
    l.dipho_n = 0;

    // fill ID variables
    if( forcedRho >= 0. ) {
        l.rho = forcedRho;
    }
    l.FillCICInputs();
    l.FillCIC();

    if(l.itype[l.current]<0) {
        bool foundHiggs=FindHiggsObjects(l);
        if(PADEBUG)  cout << " foundHiggs? "<<foundHiggs<<std::endl;
    } else {
        SetNullHiggs(l);
    }
  


    if( pho_presel.size() < 2 ) {
        l.vtx_std_ranked_list->push_back( std::vector<int>() );
        for(int ii=0;ii<l.vtx_std_n; ++ii) { l.vtx_std_ranked_list->back().push_back(ii); }
        l.vtx_std_sel = 0;
    } else {
        // fully combinatorial vertex selection
        std::vector<std::pair<int,int> > diphotons;
        for(size_t ip=0; ip<pho_presel.size(); ++ip) {
            for(size_t jp=ip+1; jp<pho_presel.size(); ++jp) {
                diphotons.push_back( std::make_pair( pho_presel[ip], pho_presel[jp] ) );
            }
        }
        for(size_t id=0; id<diphotons.size(); ++id ) {
            
            int ipho1 = diphotons[id].first;
            int ipho2 = diphotons[id].second;
            
            if(PADEBUG)        cout << " SelectEventsReduction going to fill photon info " << endl;
            PhotonInfo pho1=l.fillPhotonInfos(ipho1,vtxAlgoParams.useAllConversions,&corrected_pho_energy[0]);
            PhotonInfo pho2=l.fillPhotonInfos(ipho2,vtxAlgoParams.useAllConversions,&corrected_pho_energy[0]);
            if(PADEBUG) cout << " SelectEventsReduction done with fill photon info " << endl;
            
            l.vertexAnalysis(vtxAna_, pho1, pho2 );
            // make sure that vertex analysis indexes are in synch 
            assert( (int)id == vtxAna_.pairID(ipho1,ipho2) );
            
            l.vtx_std_ranked_list->push_back( l.vertexSelection(vtxAna_, vtxConv_, pho1, pho2, vtxVarNames, mvaVertexSelection, 
                                                                tmvaPerVtxReader_, tmvaPerVtxMethod) );
            if( tmvaPerEvtReader_ ) {
                float vtxEvtMva = vtxAna_.perEventMva( *tmvaPerEvtReader_, tmvaPerEvtMethod, l.vtx_std_ranked_list->back() );
                l.vtx_std_evt_mva->push_back(vtxEvtMva);
            }
            if( l.vtx_std_ranked_list->back().size() != 0 && ! useDefaultVertex ) {  
                l.dipho_vtx_std_sel->push_back( (l.vtx_std_ranked_list)->back()[0] );
            } else {
                l.dipho_vtx_std_sel->push_back(0);
                std::cerr << "NO VERTEX SELECTED " << l.event << " " << l.run << " " << diphotons[id].first << " " << diphotons[id].second << std::endl;
            }
            l.dipho_n = id+1;
            l.dipho_leadind[id] = diphotons[id].first;
            l.dipho_subleadind[id] = diphotons[id].second;
            l.dipho_vtxind[id] = l.dipho_vtx_std_sel->back();
            
            TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[id], l.dipho_vtxind[id], &corrected_pho_energy[0] );
            TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[id], l.dipho_vtxind[id], &corrected_pho_energy[0] );
            l.dipho_sumpt[id] = lead_p4.Pt() + sublead_p4.Pt();
            
            if( l.dipho_sumpt[id] > maxSumPt ) {
                l.vtx_std_sel = l.dipho_vtx_std_sel->back();
                maxSumPt = l.dipho_sumpt[id];
            }
        }
    }
    

    return true;
}

// ----------------------------------------------------------------------------------------------------

bool PhotonAnalysis::SkimEvents(LoopAll& l, int jentry)
{
    l.b_pho_n->GetEntry(jentry);
    if( l.pho_n < 2 ) {
        return false;
    }

    // do not run trigger selection on MC
    int filetype = l.itype[l.current];
    bool skipTrigger = !doTriggerSelection || filetype != 0 || triggerSelections.empty();
    if( ! skipTrigger ) {
        // get the trigger selection for this run 
        l.b_run->GetEntry(jentry);
        std::vector<TriggerSelection>::iterator isel = find(triggerSelections.begin(), triggerSelections.end(), l.run );
        if(isel == triggerSelections.end() ) {
            std::cerr << "No trigger selection for run " << l.run << "defined" << std::endl;
            return true;
        }

        // get the trigegr data
        l.b_hlt1_bit->GetEntry(jentry);
        l.b_hlt_path_names_HLT1->GetEntry(jentry);
        if( !  isel->pass(*(l.hlt_path_names_HLT1),*(l.hlt1_bit)) ) {
            /// std::cerr << "failed "  << std::endl;
            return false;
        }
        //// l.countersred[trigCounter_]++;
    }
    
    if( l.typerun == l.kReduce || l.typerun == l.kFillReduce ) {
        //// if( filetype == 2 ) { // photon+jet
        ////    l.b_process_id->GetEntry(jentry);
        ////    if( l.process_id == 18 ) {
        ////        return false;
        ////    }
        //// }
        
        if( filetype != 0 && ! (keepPP && keepPF && keepFF) ) {
            l.b_gp_n->GetEntry(jentry);
            l.b_gp_mother->GetEntry(jentry);
            l.b_gp_status->GetEntry(jentry);
            l.b_gp_pdgid->GetEntry(jentry);
            l.b_gp_p4->GetEntry(jentry);
            
            int np = 0;
            for(int ip=0;ip<l.gp_n;++ip) {
                if( l.gp_status[ip] != 1 || l.gp_pdgid[ip] != 22 ) { 
                    continue; 
                }
                TLorentzVector * p4 = (TLorentzVector*) l.gp_p4->At(ip);
                if( p4->Pt() < 20. || fabs(p4->Eta()) > 3. ) { continue; }
                int mother_id = abs( l.gp_pdgid[ l.gp_mother[ip] ] );
                if( mother_id <= 25 ) { ++np; }
                if( np >= 2 ) { break; }
            }
            if( np >= 2 && ! keepPP ) { return false; }
            if( np == 1 && ! keepPF ) { return false; }
            if( np == 0 && ! keepFF ) { return false; }
        }
    }

    return true;
}

// ----------------------------------------------------------------------------------------------------

bool PhotonAnalysis::SelectEvents(LoopAll& l, int jentry) 
{
    return true;
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
    vtxAna_.branches(outputTree,"vtx_std_");    

    l.pho_matchingConv = new  std::vector<int>();
    l.Branch_pho_matchingConv(outputTree);
    
    l.vtx_std_evt_mva = new std::vector<float>();
    l.vtx_std_ranked_list = new std::vector<std::vector<int> >();
    l.pho_tkiso_recvtx_030_002_0000_10_01 = new std::vector<std::vector<float> >();
    l.pho_cic6cutlevel_lead = new std::vector<std::vector<Short_t> >();
    l.pho_cic6passcuts_lead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_cic6cutlevel_sublead = new std::vector<std::vector<Short_t> >();
    l.pho_cic6passcuts_sublead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_cic4cutlevel_lead = new std::vector<std::vector<Short_t> >();
    l.pho_cic4passcuts_lead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.pho_cic4cutlevel_sublead = new std::vector<std::vector<Short_t> >();
    l.pho_cic4passcuts_sublead = new std::vector<std::vector<std::vector<UInt_t> > >();
    l.dipho_vtx_std_sel =  new std::vector<int>();

    l.Branch_vtx_std_evt_mva(outputTree);
    l.Branch_vtx_std_ranked_list(outputTree);
    l.Branch_vtx_std_sel(outputTree);
    l.Branch_pho_tkiso_recvtx_030_002_0000_10_01(outputTree);
    l.Branch_pho_tkiso_badvtx_040_002_0000_10_01(outputTree);
    l.Branch_pho_tkiso_badvtx_id(outputTree);
    l.Branch_pho_drtotk_25_99(outputTree);

    l.Branch_dipho_n(outputTree);
    l.Branch_dipho_leadind(outputTree);
    l.Branch_dipho_subleadind(outputTree);
    l.Branch_dipho_vtxind(outputTree);
    l.Branch_dipho_sumpt(outputTree);
    
    l.Branch_pho_cic6cutlevel_lead( outputTree );
    l.Branch_pho_cic6passcuts_lead( outputTree );
    l.Branch_pho_cic6cutlevel_sublead( outputTree );
    l.Branch_pho_cic6passcuts_sublead( outputTree );
    l.Branch_pho_cic4cutlevel_lead( outputTree );
    l.Branch_pho_cic4passcuts_lead( outputTree );
    l.Branch_pho_cic4cutlevel_sublead( outputTree );
    l.Branch_pho_cic4passcuts_sublead( outputTree );

    l.Branch_pho_genmatched(outputTree);
    l.Branch_pho_regr_energy_otf(outputTree);
    l.Branch_pho_regr_energyerr_otf(outputTree);
    
    l.Branch_gh_gen2reco1( outputTree );
    l.Branch_gh_gen2reco2( outputTree );
    l.Branch_gh_vbfq1_pdgid( outputTree );
    l.Branch_gh_vbfq2_pdgid( outputTree );
    l.Branch_gh_vh_pdgid( outputTree );
    l.Branch_gh_vh1_pdgid( outputTree );
    l.Branch_gh_vh2_pdgid( outputTree );
    l.Branch_gh_higgs_p4( outputTree );
    l.Branch_gh_pho1_p4( outputTree );
    l.Branch_gh_pho2_p4( outputTree );
    l.Branch_gh_vbfq1_p4( outputTree );
    l.Branch_gh_vbfq2_p4( outputTree );
    l.Branch_gh_vh1_p4( outputTree );
    l.Branch_gh_vh2_p4( outputTree );

}


// ---------------------------------------------------------------------------------------------------------------------------------------------
int PhotonAnalysis::DiphotonMVASelection(LoopAll &l, HggVertexAnalyzer & vtxAna, Float_t & diphoMVA,  
                                         Float_t minTightMVA, Float_t minLooseMVA, Float_t leadPtMin, 
                                         Float_t subleadPtMin, std::string type, int ncategories,
                                         bool applyPtoverM, float *pho_energy_array, bool split) {

    //rho=0;// CAUTION SETTING RHO TO 0 FOR 2010 DATA FILES (RHO ISN'T IN THESE FILES)
    int selected_lead_index = -1;
    int selected_sublead_index = -1;
    float selected_lead_pt = -1;
    float selected_sublead_pt = -1;
  
    std::vector<int> passingLoose_dipho;
    std::vector<int> passingTight_dipho;
    std::vector<float> passingLoose_sumpt;
    std::vector<float> passingTight_sumpt;
    std::map<int,float> passing_diphomva;
    //std::vector<float> passing_leadphomva;
    //std::vector<float> passing_subphomva;
    if(PADEBUG)  std::cout <<"type "<< type<<std::endl;
    for(int idipho = 0; idipho < l.dipho_n; ++idipho ) {
        int ivtx = l.dipho_vtxind[idipho];
        int lead = l.dipho_leadind[idipho];
        int sublead = l.dipho_subleadind[idipho];
      
        if( lead == sublead ) { continue; }

        TLorentzVector lead_p4 = l.get_pho_p4(lead,ivtx,pho_energy_array); 
        TLorentzVector sublead_p4 = l.get_pho_p4(sublead,ivtx,pho_energy_array); 
      
        float leadEta = fabs(((TVector3 *)l.sc_xyz->At(l.pho_scind[lead]))->Eta());
        float subleadEta = fabs(((TVector3 *)l.sc_xyz->At(l.pho_scind[sublead]))->Eta());
        float m_gamgam = (lead_p4+sublead_p4).M();
        float pt_gamgam = (lead_p4+sublead_p4).Pt();

        // FIXME call to (*vtx_std_evt_mva)[ivtx] crashes program
        //std::cout<<"(*vtx_std_evt_mva)["<<ivtx<<"] "<< (*vtx_std_evt_mva)[ivtx]<<std::endl; 
        float vtxProb = vtxAna.vertexProbability((*l.vtx_std_evt_mva)[ivtx]);
        //float vtxProb = vtxAna.vertexProbability(0.0);
        //float vtxProb = 0.3;

        if( leadEta > 2.5 || subleadEta > 2.5 || 
            ( leadEta > 1.4442 && leadEta < 1.566 ) ||
            ( subleadEta > 1.4442 && subleadEta < 1.566 ) ) { continue; }
     
        // VBF cuts smoothly on lead pt/M but on straight pt on sublead to save sig eff and avoid HLT turn-on
        if(split){
            if ( lead_p4.Pt()/m_gamgam < leadPtMin/120. || sublead_p4.Pt()< subleadPtMin ) { continue; }
        }else{
            if( applyPtoverM ) {
                if ( lead_p4.Pt()/m_gamgam < leadPtMin/120. || sublead_p4.Pt()/m_gamgam < subleadPtMin/120. ) { continue; }
            } else {
                if ( lead_p4.Pt() < leadPtMin || sublead_p4.Pt() < subleadPtMin ) { continue; }
            }
        }

        if(PADEBUG)  std::cout << "getting photon ID MVA" << std::endl;
        std::vector<std::vector<bool> > ph_passcut;
        float leadmva = l.photonIDMVA(lead, ivtx, lead_p4, type.c_str()); 
        float submva = l.photonIDMVA(sublead, ivtx, sublead_p4, type.c_str());
        if(PADEBUG)  std::cout << "lead  leadmva  sublead  submva  "<<lead<<"  "<<leadmva<<"  "<<sublead<<"  "<<submva << std::endl;
        if(PADEBUG)  std::cout << "lead_p4.Pt()  leadEta  sublead_p4.Pt()  subleadEta  "<<lead_p4.Pt()<<"  "<<leadEta<<"  "<<sublead_p4.Pt()<<"  "<<subleadEta << std::endl;
        if(PADEBUG)  std::cout << "mass "<<m_gamgam<<std::endl;
        if(leadmva>minLooseMVA && submva>minLooseMVA) {
            passingLoose_dipho.push_back(idipho);
            passingLoose_sumpt.push_back(lead_p4.Et()+sublead_p4.Et());
            if(leadmva>minTightMVA && submva>minTightMVA) {
                passingTight_dipho.push_back(idipho);
                passingTight_sumpt.push_back(lead_p4.Et()+sublead_p4.Et());
            }
        } else { continue; } 

        if(PADEBUG)  std::cout << "got photon ID MVA" << std::endl;
     
        if(PADEBUG)  std::cout << "getting di-photon MVA" << std::endl;

        massResolutionCalculator->Setup(l,&photonInfoCollection[lead],&photonInfoCollection[sublead],idipho,eSmearPars,nR9Categories,nEtaCategories);
        //massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,lead,sublead,idipho,pt_gamgam, m_gamgam,eSmearPars,nR9Categories,nEtaCategories);

        float sigmaMrv = massResolutionCalculator->massResolutionCorrVtx();
        float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
        float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
    
        // FIXME  These things will move to MassResolution eventually...
        if(type=="UCSD") {
            float leta=fabs( ((TLorentzVector*)l.pho_p4->At(lead))->Eta() );
            float seta=fabs( ((TLorentzVector*)l.pho_p4->At(sublead))->Eta() );

            //MARCO FIXED
            float leadErr = GetSmearSigma(leta,l.pho_r9[lead]);
            float subleadErr = GetSmearSigma(seta,l.pho_r9[sublead]);

            Float_t t_dmodz = l.getDmOverDz(lead, sublead, pho_energy_array);
            float z_gg = ((TVector3*)(l.vtx_std_xyz->At(ivtx)))->Z();
            // will be updated to this
            //sigmaMwv = fabs(t_dmodz)*(sqrt(pow(double(l.bs_sigmaZ), 2) + pow(double(z_gg), 2)))/m_gamgam;
            sigmaMwv = (t_dmodz)*(sqrt(pow(double(l.bs_sigmaZ), 2) + pow(double(z_gg), 2)));
      
            sigmaMeonly = 0.5*sqrt( pow(l.pho_regr_energyerr[lead]/l.pho_regr_energy[lead],2)+pow(leadErr,2)+pow(l.pho_regr_energyerr[sublead]/l.pho_regr_energy[sublead],2)+pow(subleadErr,2) )*m_gamgam;
        }

        passing_diphomva[idipho]=
            (l.diphotonMVA(lead, sublead, ivtx, 
                           vtxProb, lead_p4, sublead_p4, 
                           sigmaMrv, sigmaMwv, sigmaMeonly, type.c_str()));
        if(PADEBUG)  std::cout << "got di-photon MVA" << std::endl;
        //passing_leadphomva.push_back(leadmva);
        //passing_subphomva.push_back(submva);
    }
  
    if( passingLoose_dipho.empty() ) { return -1; }

    if( passingTight_dipho.empty() ) {
        std::sort(passingLoose_dipho.begin(),passingLoose_dipho.end(),
                  SimpleSorter<float,std::greater<double> >(&passingLoose_sumpt[0]));
        // set reference to diphoton MVA
        diphoMVA = passing_diphomva[passingLoose_dipho[0]];

        if(PADEBUG)  std::cout<<"exiting DiphotonMVAsel"<<std::endl;
        return passingLoose_dipho[0];
    } else {
        std::sort(passingTight_dipho.begin(),passingTight_dipho.end(),
                  SimpleSorter<float,std::greater<double> >(&passingTight_sumpt[0]));
        // set reference to diphoton MVA
        diphoMVA = passing_diphomva[passingTight_dipho[0]];

        if(PADEBUG)  std::cout<<"exiting DiphotonMVAsel"<<std::endl;
        return passingTight_dipho[0];
    }

}

// ---------------------------------------------------------------------------------------------------------------------------------------------
int PhotonAnalysis::DiphotonMVAEventClass(LoopAll &l, float diphoMVA, int nCat, std::string type, int EBEB){
    if(PADEBUG)  std::cout<<"DiphotonMVAEventClass 1"<<std::endl;
  
    int eventClass = -1;

    float class5threshMIT[4]  = { 0.05,  0.55, 0.72, 0.89 };
    float class6threshUCSD[6] = { -0.4, -0.0356,  0.3889, 0.592, 0.6669,  0.7583 };
    // first 2 for (ebee+eeee) and last 6 for ebeb
    float class8threshUCSD[8]={-0.7, -0.11, -0.4, -0.0356,  0.3889, 0.592, 0.6669,  0.7583 };
    float class2threshUCSD[2]={-0.2, 0.3 };

    if(PADEBUG)  std::cout<<"DiphotonMVAEventClass diphoMVA:  "<<diphoMVA<<std::endl;
    if(PADEBUG)  std::cout<<"DiphotonMVAEventClass nCat:  "<<nCat<<std::endl;
    if(PADEBUG)  std::cout<<"DiphotonMVAEventClass type:  "<<type<<std::endl;
    if(PADEBUG)  std::cout<<"DiphotonMVAEventClass EBEB:  "<<EBEB<<std::endl;
  
    if(nCat==4){
        if(PADEBUG)  std::cout<<"DiphotonMVAEventClass 3"<<std::endl;
        for(int ithresh=0; ithresh<nCat; ithresh++){
            if(PADEBUG)  std::cout <<"eventClass "<<eventClass
                                   <<"\tclass5threshMIT["<<ithresh<<"] "<<class5threshMIT[ithresh]<<std::endl;
            if(class5threshMIT[ithresh]<diphoMVA && type=="MIT") eventClass++;
        }
    } else if(nCat==6){
        for(int ithresh=0; ithresh<nCat; ithresh++){
            if(class6threshUCSD[ithresh]<diphoMVA && type=="UCSD") eventClass++;
        }
    } else if(nCat==8){
        if(type=="UCSD") {
            for(int ithresh=0; ithresh<nCat; ithresh++){
                if(class8threshUCSD[ithresh]<diphoMVA && EBEB==1 && ithresh>1) {
                    if(PADEBUG)  std::cout <<"passing EBEB eventClass "<<eventClass<<std::endl;
                    eventClass++;
                }
        
                if(class8threshUCSD[ithresh]<diphoMVA && EBEB!=1 && ithresh<2) {
                    if(PADEBUG)  std::cout <<"passing !EBEB eventClass "<<eventClass<<std::endl;
                    eventClass++;
                }
            }
            if(eventClass!=-1 && EBEB==1) eventClass=eventClass+2;
        }
    } else if(nCat==2) {
        for(int ithresh=0; ithresh<nCat; ithresh++){
            if(class2threshUCSD[ithresh]<diphoMVA && type=="UCSD") eventClass++;
        }
    }

 

    if(PADEBUG)  std::cout<<"eventClass "<<eventClass<<std::endl;
    return eventClass;
}

void PhotonAnalysis::ResetAnalysis(){
}

float PhotonAnalysis::GetSmearSigma(float eta, float r9, int epoch){
    // not using epoch for now
    // r9cat + nr9Cat*etaCat
    float smear_nov14[] = {0.99e-2, 1.00e-2, 1.57e-2, 2.17e-2, 3.30e-2, 3.26e-2, 3.78e-2, 3.31e-2}; 
    float sigma = -1;
    if(epoch==0) {
        int nr9Cat=2;
        int r9Cat=(int)(r9<0.94);
        int nEtaCat=4;
        int etaCat=(int)(eta>1.0) + (int)(eta>1.5) + (int)(eta>2.0);
        sigma=smear_nov14[(int)(r9Cat + nr9Cat*etaCat)];
    }

    return sigma;
}
    
    
void PhotonAnalysis::SetNullHiggs(LoopAll& l){
  
    l.gh_higgs_p4 = new TLorentzVector(0,0,0,0);
  
    l.gh_pho1_p4 = new TLorentzVector(0,0,0,0);
    l.gh_pho2_p4 = new TLorentzVector(0,0,0,0);

    l.gh_vbfq1_p4 = new TLorentzVector(0,0,0,0);
    l.gh_vbfq2_p4 = new TLorentzVector(0,0,0,0);
    l.gh_vbfq1_pdgid=-10000;
    l.gh_vbfq2_pdgid=-10000;
  
    l.gh_vh1_p4 = new TLorentzVector(0,0,0,0);
    l.gh_vh2_p4 = new TLorentzVector(0,0,0,0);
    l.gh_vh_pdgid=-10000;
    l.gh_vh1_pdgid=-10000;
    l.gh_vh2_pdgid=-10000;

}


bool PhotonAnalysis::FindHiggsObjects(LoopAll& l){
  
    int higgsind=-1;
    int mc1=-1;
    int mc2=-1;
    int i1=-1;
    int i2=-1;

    l.FindMCHiggsPhotons( higgsind,  mc1,  mc2,  i1,  i2 );

  
    if(higgsind!=-1) l.gh_higgs_p4 = new TLorentzVector(*((TLorentzVector*)((TLorentzVector*) l.gp_p4->At(higgsind))->Clone()));
    else l.gh_higgs_p4 = new TLorentzVector(0,0,0,0);
  
    l.gh_gen2reco1=i1;
    l.gh_gen2reco2=i2;
  
  
    if(mc1!=-1) l.gh_pho1_p4 = new TLorentzVector(*((TLorentzVector*)((TLorentzVector*) l.gp_p4->At(mc1))->Clone()));
    else l.gh_pho1_p4 = new TLorentzVector(0,0,0,0);

    if(mc2!=-1) l.gh_pho2_p4 = new TLorentzVector(*((TLorentzVector*)((TLorentzVector*) l.gp_p4->At(mc2))->Clone()));
    else l.gh_pho2_p4 = new TLorentzVector(0,0,0,0);

  

    int vbfq1=-100;
    int vbfq2=-100;

    int vh=-100;
    int vh1=-100;
    int vh2=-100;

    if(higgsind!=-1){
        l.FindMCVBF(higgsind,vbfq1,vbfq2);
        l.FindMCVH(higgsind,vh,vh1,vh2);

        //std::cout<<"higgsind vbfq1 vbfq2 vh gp_pdgid[vh] vh1 vh2  "<<higgsind<<"  "<<vbfq1<<"  "<<vbfq2<<"  "<<vh<<"  "<<l.gp_pdgid[vh]<<"  "<<vh1<<"  "<<vh2<<std::endl;
    }
  
  
    if(vh==-100) l.gh_vh_pdgid=-10000;
    else l.gh_vh_pdgid=l.gp_pdgid[vh];

    if(vh1==-100) l.gh_vh1_pdgid=-10000;
    else l.gh_vh1_pdgid=l.gp_pdgid[vh1];

    if(vh2==-100) l.gh_vh2_pdgid=-10000;
    else l.gh_vh2_pdgid=l.gp_pdgid[vh2];
  
    if(vh1==-100) l.gh_vh1_p4 = new TLorentzVector(0,0,0,0);
    else l.gh_vh1_p4 = new TLorentzVector(*((TLorentzVector*)((TLorentzVector*) l.gp_p4->At(vh1))->Clone()));

    if(vh2==-100) l.gh_vh2_p4 = new TLorentzVector(0,0,0,0);
    else l.gh_vh2_p4 = new TLorentzVector(*((TLorentzVector*)((TLorentzVector*) l.gp_p4->At(vh2))->Clone()));


  
    if(vbfq1==-100) l.gh_vbfq1_pdgid=-10000;
    else l.gh_vbfq1_pdgid=l.gp_pdgid[vbfq1];

    if(vbfq2==-100) l.gh_vbfq2_pdgid=-10000;
    else l.gh_vbfq2_pdgid=l.gp_pdgid[vbfq2];

    if(vbfq1==-100) l.gh_vbfq1_p4 = new TLorentzVector(0,0,0,0);
    else l.gh_vbfq1_p4 = new TLorentzVector(*((TLorentzVector*)((TLorentzVector*) l.gp_p4->At(vbfq1))->Clone()));

    if(vbfq2==-100) l.gh_vbfq2_p4 = new TLorentzVector(0,0,0,0);
    else l.gh_vbfq2_p4 = new TLorentzVector(*((TLorentzVector*)((TLorentzVector*) l.gp_p4->At(vbfq2))->Clone()));

    return (higgsind != -1);

}



Bool_t PhotonAnalysis::GenMatchedPhoton(LoopAll& l, int ipho){
    Bool_t is_prompt = false;
    TLorentzVector* phop4 = (TLorentzVector*) l.pho_p4->At(ipho);

    for(int ip=0;ip<l.gp_n;++ip) {
        if( l.gp_status[ip] != 1 || l.gp_pdgid[ip] != 22 ) {
            continue;
        }
        TLorentzVector * p4 = (TLorentzVector*) l.gp_p4->At(ip);
        if( p4->Pt() < 20. || fabs(p4->Eta()) > 3. ) { continue; }
        int mother_id = abs( l.gp_pdgid[ l.gp_mother[ip] ] );
        if( mother_id <= 25 ) {
            float dr = phop4->DeltaR(*p4);
            if (dr<0.2) {
                is_prompt = true;
                break;
            }
        }
    }
    return is_prompt;

}


bool PhotonAnalysis::ClassicCatsNm1Plots(LoopAll& l, int diphoton_nm1_id, float* smeared_pho_energy, float eventweight, float myweight){
    
    bool pass = false;
    if(diphoton_nm1_id==-1) return pass;

    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_nm1_id], l.dipho_vtxind[diphoton_nm1_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_nm1_id], l.dipho_vtxind[diphoton_nm1_id], &smeared_pho_energy[0]);
    TLorentzVector diphoton = lead_p4+sublead_p4;
    
    int photon_category   = l.PhotonCategory(l.dipho_subleadind[diphoton_nm1_id],2,2);
    sublead_r9            = l.pho_r9[l.dipho_subleadind[diphoton_nm1_id]];
    sublead_trkisooet     = (*l.pho_tkiso_recvtx_030_002_0000_10_01)[l.dipho_subleadind[diphoton_nm1_id]][l.dipho_vtxind[diphoton_nm1_id]]*50/sublead_p4.Et();
    sublead_isoOverEt     = (l.pho_hcalsumetconedr04[l.dipho_subleadind[diphoton_nm1_id]]
                             +  l.pho_ecalsumetconedr03[l.dipho_subleadind[diphoton_nm1_id]]
                             +  (*l.pho_tkiso_recvtx_030_002_0000_10_01)[l.dipho_subleadind[diphoton_nm1_id]][l.dipho_vtxind[diphoton_nm1_id]]
                             - 0.17*l.rho)*50/sublead_p4.Et();
    sublead_badisoOverEt  = (l.pho_hcalsumetconedr04[l.dipho_subleadind[diphoton_nm1_id]]
                             +  l.pho_ecalsumetconedr04[l.dipho_subleadind[diphoton_nm1_id]]
                             +  l.pho_tkiso_badvtx_040_002_0000_10_01[l.dipho_subleadind[diphoton_nm1_id]]
                             - 0.52*l.rho)*50/sublead_p4.Et();
    
    sublead_sieie         = l.pho_sieie[l.dipho_subleadind[diphoton_nm1_id]];
    sublead_drtotk        = l.pho_drtotk_25_99[l.dipho_subleadind[diphoton_nm1_id]];
    sublead_hovere        = l.pho_hoe[l.dipho_subleadind[diphoton_nm1_id]];
    sublead_mgg           = diphoton.M();
    
    int applyCutsType = 15 + photon_category;
    
    pass = l.ApplyCutsFill(0,applyCutsType, eventweight, myweight);

    return pass;
}


bool PhotonAnalysis::ElectronTag2011(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy, bool nm1, float eventweight, float myweight){
    bool tag = false;

    if(diphotonVHlep_id==-1) return tag;
    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
    
    int elInd = l.ElectronSelection(lead_p4, sublead_p4, l.dipho_vtxind[diphotonVHlep_id]);
    if(elInd!=-1) tag = true;
    
    return tag;
}



bool PhotonAnalysis::MuonTag2011(LoopAll& l, int diphotonVHlep_id, float* smeared_pho_energy, bool nm1, float eventweight, float myweight){
    bool tag = false;

    if(diphotonVHlep_id==-1) return tag;
        
    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHlep_id], l.dipho_vtxind[diphotonVHlep_id], &smeared_pho_energy[0]);
    
    int muonInd = l.MuonSelection(lead_p4, sublead_p4, l.dipho_vtxind[diphotonVHlep_id]);
    if(muonInd!=-1) tag = true;
    
    return tag;
}



bool PhotonAnalysis::VBFTag2011(LoopAll& l, int diphoton_id, float* smeared_pho_energy, bool nm1, float eventweight, float myweight){
    bool tag = false;

    if(diphoton_id==-1) return tag;
    float jet1ptcut =15.0;
    float jet2ptcut =15.0;
  
  
    TLorentzVector lead_p4    = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
  	TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
          
    std::pair<int, int> jets = l.Select2HighestPtJets(lead_p4, sublead_p4, jet1ptcut, jet2ptcut );
    if(jets.first==-1 or jets.second==-1) return tag;
    
    TLorentzVector diphoton = lead_p4+sublead_p4;
  
    TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.first);
    TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.second);
    TLorentzVector dijet = (*jet1) + (*jet2);
    
    myVBFLeadJPt= jet1->Pt();
    myVBFSubJPt = jet2->Pt();
    myVBF_Mjj   = dijet.M();
    myVBFdEta   = fabs(jet1->Eta() - jet2->Eta());
    myVBFZep    = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
    myVBFdPhi   = fabs(diphoton.DeltaPhi(dijet));
    myVBF_Mgg   = diphoton.M();
  
    if(nm1){
        tag = l.ApplyCutsFill(0,1, eventweight, myweight);
    } else {
        tag = l.ApplyCuts(0,1);
    }
  
    return tag;
}

bool PhotonAnalysis::VHhadronicTag2011(LoopAll& l, int diphotonVHhad_id, float* smeared_pho_energy, bool nm1, float eventweight, float myweight){
    bool tag = false;

    if(diphotonVHhad_id==-1) return tag;
    float jet1ptcut =15.0;
    float jet2ptcut =15.0;
  
  
    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphotonVHhad_id], l.dipho_vtxind[diphotonVHhad_id], &smeared_pho_energy[0]);
    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphotonVHhad_id], l.dipho_vtxind[diphotonVHhad_id], &smeared_pho_energy[0]);
    
    std::pair<int, int> jets = l.Select2HighestPtJets(lead_p4, sublead_p4, jet1ptcut, jet2ptcut );
    if(jets.first==-1 or jets.second==-1) return tag;
  
    TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.first);
    TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(jets.second);
    TLorentzVector dijet = (*jet1) + (*jet2);
    
    TLorentzVector diphoton = lead_p4+sublead_p4;
    
    myVHhadLeadJPt = jet1->Pt();
    myVHhadSubJPt = jet2->Pt();
    myVHhad_Mjj = dijet.M();
    myVHhaddEta = fabs(jet1->Eta() - jet2->Eta());
    myVHhadZep  = fabs(diphoton.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
    myVHhaddPhi = fabs(diphoton.DeltaPhi(dijet));
    myVHhad_Mgg =diphoton.M();

    
    if(nm1){
        tag = l.ApplyCutsFill(0,2, eventweight, myweight);
    } else {
        tag = l.ApplyCuts(0,2);
    }
  
    return tag;
}




// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
