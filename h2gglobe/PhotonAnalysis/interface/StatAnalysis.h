#ifndef __STATANALYSIS__
#define __STATANALYSIS__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
#include "PdfWeightSmearer.h"
#include "InterferenceSmearer.h"
#include <iostream>
#include <fstream>
#include "math.h"

// ------------------------------------------------------------------------------------
class StatAnalysis : public PhotonAnalysis 
{
 public:
    
    StatAnalysis();
    virtual ~StatAnalysis();
    
    virtual const std::string & name() const { return name_; };
    
    // LoopAll analysis interface implementation
    void Init(LoopAll&);
    void Term(LoopAll&);
    
    void GetBranches(TTree *, std::set<TBranch *>& );
    
    virtual bool SelectEvents(LoopAll&, int);
    virtual void ResetAnalysis();
    virtual bool Analysis(LoopAll&, Int_t);
    
    std::string efficiencyFile;
    std::string efficiencyFileMVA;
    std::string effPhotonCategoryType;
    int effPhotonNCat;

    EfficiencySmearer::efficiencySmearingParameters effSmearPars;
    DiPhoEfficiencySmearer::diPhoEfficiencySmearingParameters diPhoEffSmearPars;

    double GetDifferentialKfactor(double, int);

    void FillSignalLabelMap(LoopAll &l);
    std::string GetSignalLabel(int id, LoopAll &l) ;
    
    bool unblind;
    
    bool  doEscaleSyst, doEresolSyst, doPhotonIdEffSyst, doVtxEffSyst, doR9Syst, doTriggerEffSyst, doKFactorSyst;
    bool  doEscaleSmear, doEresolSmear, doPhotonIdEffSmear, doVtxEffSmear, doR9Smear, doTriggerEffSmear, 
	doKFactorSmear, doInterferenceSmear;
    float systRange;
    int   nSystSteps;   
    //int   nEtaCategories, nR9Categories, nPtCategories;
    std::vector<int> cicCutLevels;
    std::vector<int> bkgPolOrderByCat;
    bool doMcOptimization;
    bool fillOptTree;
    bool doFullMvaFinalTree;

    bool splitwzh;

    void fillOpTree(LoopAll& l, const TLorentzVector & lead_p4, const TLorentzVector & sublead_p4, float *smeared_pho_energy, Float_t vtxProb,
		     std::pair<int, int> diphoton_index, Int_t diphoton_id, Float_t phoid_mvaout_lead, Float_t phoid_mvaout_sublead,
		     Float_t weight, Float_t evweight, Float_t mass, Float_t sigmaMrv, Float_t sigmaMwv,
		     const TLorentzVector & Higgs, Float_t diphobdt_output, Int_t category, bool VBFevent, Float_t myVBF_Mjj, Float_t myVBFLeadJPt, 
		     Float_t myVBFSubJPt, Int_t nVBFDijetJetCategories, bool isSyst, std::string name1);

    void fillOpTree(LoopAll& l, const TLorentzVector & lead_p4, const TLorentzVector & sublead_p4, Float_t vtxProb,
		    std::pair<int, int> diphoton_index, Int_t diphoton_id, Float_t phoid_mvaout_lead, Float_t phoid_mvaout_sublead,
		    Float_t weight, Float_t mass, Float_t sigmaMrv, Float_t sigmaMwv,
		    const TLorentzVector & Higgs, Float_t diphobdt_output, Int_t category, bool VBFevent, Float_t myVBF_Mjj, Float_t myVBFLeadJPt, 
		    Float_t myVBFSubJPt, Int_t nVBFDijetJetCategories, bool isSyst, std::string name1);

    int nDataBins;  
    bool scaleClusterShapes, scaleR9Only;
    bool dumpAscii, dumpMcAscii;
    float phoidMvaCut;
    std::vector<double> zeePtBinLowEdge, zeePtWeight;
    std::vector<int> sigPointsToBook;
    
    std::string kfacHist;
    std::string pdfWeightHist;

    TH1D *thm110,*thm120,*thm130,*thm140;

 protected:
    // Factorized Analysis method + loop over systematics
    virtual bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, 
			      int & category, int & diphoton_id,
			      bool & isCorrectVertex,float &kinematic_bdtout,
			      bool isSyst=false, 
			      float syst_shift=0., bool skipSelection=false,
			      BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0); 

    virtual void FillRooContainer(LoopAll& l, int cur_type, float mass, float diphotonMVA, int category, float weight, 
				  bool isCorrectVertex, int diphoton_ind=-1);
    virtual void AccumulateSyst(int cur_type, float mass, float diphotonMVA, int category, float weight,
				std::vector<double> & mass_errors,
				std::vector<double> & mva_errors,
				std::vector<int>    & categories,
				std::vector<double> & weights);
    virtual void FillRooContainerSyst(LoopAll& l, const std::string & name,int cur_type,
				      std::vector<double> & mass_errors, std::vector<double> & mva_errors,
				      std::vector<int>    & categories, std::vector<double> & weights, int diphoton_id=-1);
    
    bool VHmuevent, VHelevent, VBFevent, VHhadevent,VHhadBtagevent, VHmetevent,TTHhadevent,TTHlepevent,tHqLeptonicevent;  //met at analysis step
    bool VHlep1event, VHlep2event;
    int VHelevent_cat;
    int VHmuevent_cat;
    int VHmetevent_cat;
    double genLevWeight; 
    float sigmaMrv, sigmaMwv;
    int vbfIjet1, vbfIjet2;

    void buildBkgModel(LoopAll& l, const std::string & postfix);

    std::vector<float> smeared_pho_energy;
    std::vector<float> smeared_pho_r9;
    std::vector<float> smeared_pho_weight;

    void  computeExclusiveCategory(LoopAll & l, int & category, std::pair<int,int> diphoton_index, float pt, float diphoBDT=1. );	
    int  categoryFromBoundaries(std::vector<float> & v, float val);
    int  categoryFromBoundaries2D(std::vector<float> & v1, std::vector<float> & v2, std::vector<float> & v3, float val1, float val2, float val3);

    void fillControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, const TLorentzVector & Higgs, 
			  float lead_r9, float sublead_r9, int diphoton_index, 
              int category, bool rightvtx, float evweight , TVector3* vtx, LoopAll &, 
              int muVtx=-1, int mu_ind=-1, int elVtx=-1, int el_ind=-1, float bdtoutput=-10);

    void fillSignalEfficiencyPlots(float weight, LoopAll & l );

    void rescaleClusterVariables(LoopAll &l);

    EfficiencySmearer *idEffSmearer, *r9Smearer;
    DiPhoEfficiencySmearer *vtxEffSmearer, *triggerEffSmearer;
    KFactorSmearer * kFactorSmearer;
    PdfWeightSmearer * pdfWeightSmearer;
    InterferenceSmearer * interferenceSmearer;
    
    std::string name_;
    std::map<int,std::string> signalLabels;
    float nevents, sumwei, sumaccept, sumsmear, sumev; 
    
    int nInclusiveCategories_;
    int nCategories_;
    int nPhotonCategories_;
    int diPhoCounter_;
    int nVBFCategories  ; 
    int nVHhadCategories; 
    int nVHhadBtagCategories; 
    int nTTHhadCategories; 
    int nTTHlepCategories; 
    int nVHlepCategories; 
    int nVHmetCategories; 
    int ntHqLeptonicCategories;
    
    // RooStuff
    RooContainer *rooContainer;

    ofstream eventListText;
    //vector<double> weights;
    TFile *kfacFile;
    
};

#endif


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
