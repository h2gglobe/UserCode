#ifndef __MICROANALYSIS__
#define __MICROANALYSIS__

#include "StatAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "RooContainer.h"
#include "RooGenFunction.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"

#include "TMVA/Reader.h"

#define speedOfLight 299792458.


// ------------------------------------------------------------------------------------
class MicroAnalysis : public StatAnalysis
{
public:
	
	MicroAnalysis();
	virtual ~MicroAnalysis();
	
	virtual const std::string & name() const { return name_; };
	
	// LoopAll analysis interface implementation
	void Init(LoopAll&);
	void Term(LoopAll&);
	
	virtual void Analysis(LoopAll&, Int_t);

	TString uFileName;
	std::string tmvaMethod;
	std::string tmvaWeights;
	TString roofitProbs;
	TString timePerfHist;
	float   timeOptimismFactor;
	Int_t   storeNVert;
	// bool    useMvaRanking,addConversionToMva;

protected:
	
	float getTimeResol(float absDeltaEta, bool iseb1, bool iseb2);
	float getDeltaTof(TVector3 &posLead, TVector3 &posSubLead, TVector3 &posVertex);
	float getExtraTravelTime(TVector3 &posSC, TVector3 &posVertex);

	std::string name_;
	TFile *uFile_;

	//Variables for the tree
	// All floats because that's what TMVA likes
	TTree *uTree_;
	UShort_t	nVert_;
	UShort_t	nPU_;
	TLorentzVector *pho1_;
	TLorentzVector *pho2_;
	TLorentzVector *dipho_;
	Float_t absDeltaEta_;
	Float_t ptasym_;
	Float_t ptbal_;
	Float_t logsumpt2_;
	Float_t tofCorrTdiff_;
	Float_t dZToGen_;
	Float_t dZtoClosest_;
	Float_t	evWeight_;
	Bool_t isClosestToGen_;
	Float_t nConv_;
	Float_t convCompat_;
	Float_t dZToConv_;
	Float_t pullToConv_;
	Float_t rprod_;
	Float_t convRprod_;
	std::vector<std::string> vtxVarNames_;
	std::vector<float> vtxVars_;
	
	// TMVA reader
	vector<string> tmvaVariables_, rankVariables_;
	TMVA::Reader *tmvaReader_;
	TTree *evTree_;
	Float_t dZTrue_, zRMS_, zTrue_, rTrue_, alphaTrue_;
	Int_t category_;
	Float_t mTrue_, mTrueVtx_;
	vector<float> alpha_;
	vector<float> MVA_;
	vector<float> rankprod_;
	vector<float> convRankprod_;
	vector<float> prob_;
	vector<float> corrProb_;
	vector<float> wrongProb_;
	vector<float> zRMSn_;
	vector<float> dZ_;
	vector<float> diphoM_;
	vector<float> diphoCosTheta_;
	vector<float> diphoCosDeltaPhi_;
	vector<float> diphoPt_;
	vector<float> dZToConvV_;
	vector<float> pullToConvV_;
	
	RooGenFunction * vtxProb_, * sigProb_, * bkgProb_;
	TFile * rooFile_;
	TFile * rooFileTime_;
	
	TH2F* dtVSdEtaEBEB_, * dtVSdEtaEBEE_, * dtVSdEtaEEEE_;
	std::vector<TH1F> dtEBEB_, dtEBEE_,  dtEEEE_;

};



#endif
