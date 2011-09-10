#ifndef __MICROANALYSIS__
#define __MICROANALYSIS__

#include "StatAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
#include "TTree.h"

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

	//Variables for the tree
	UShort_t	nVert_;
	UShort_t	nPU_;
	TLorentzVector *pho1_;
	TLorentzVector *pho2_;
	TLorentzVector *dipho_;
	Float_t ptasym_;
	Float_t ptbal_;
	Float_t logsumpt2_;
	Float_t dZToGen_;
	Float_t dZtoClosest_;
	Float_t	evWeight_;
	Bool_t isClosestToGen_;


	TString uFileName;

protected:
	std::string name_;
	TTree *uTree_;
	TFile *uFile_;
};

#endif
