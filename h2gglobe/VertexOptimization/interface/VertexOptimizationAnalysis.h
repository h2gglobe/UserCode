#ifndef __VertexOptimizationAnalysis__
#define __VertexOptimizationAnalysis__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis/interface/PhotonAnalysis.h"
#include "PhotonAnalysis/interface/StatAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
#include <iostream>
#include <fstream>
#include "math.h"

// ------------------------------------------------------------------------------------
class VertexOptimizationAnalysis : public StatAnalysis 
{
 public:
    
    VertexOptimizationAnalysis();
    virtual ~VertexOptimizationAnalysis();
    
    virtual const std::string & name() const { return name_; };
    
    // LoopAll analysis interface implementation
    void Init(LoopAll&);
    void Term(LoopAll&);
    
    void ReducedOutputTree(LoopAll &l, TTree * outputTree);

    virtual bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, int & category, int & diphoton_id,
			      bool & isCorrectVertex,
			      bool isSyst=false, 
			      float syst_shift=0., bool skipSelection=false,
			      BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0); 
    
    TString uFileName;
    
private:
    std::vector<std::string> vtxVarNames_;
    std::vector<float> vtxVars_;
    
    TFile * uFile_;
    TTree * uTree_;
    TH1 * hMinBiasSpecturm_;
    bool isClosestToGen_;
    int   nPU_, nVert_;
    float evWeight_;
    
};

#endif


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
