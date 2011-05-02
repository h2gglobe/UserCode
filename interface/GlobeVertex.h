#ifndef GLOBEVERTEX_H
#define GLOBEVERTEX_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"


#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <iostream>

class GlobeVertex {
 public:
  
  GlobeVertex(const edm::ParameterSet&, const char*);
  virtual ~GlobeVertex() {};

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  // variables
  Int_t vtx_n;

  TClonesArray *bs_xyz;
  Float_t bs_sigmaZ;
  Float_t bs_x0Error;
  Float_t bs_y0Error;
  Float_t bs_z0Error;
  Float_t bs_sigmaZ0Error;

  TClonesArray *vtx_xyz;
  TClonesArray *vtx_dxdydz;
  TClonesArray *vtx_vectorp3;
  Float_t vtx_scalarpt[MAX_VERTICES];
  Float_t vtx_x2dof[MAX_VERTICES];
  Float_t vtx_ndof[MAX_VERTICES];

  Int_t vtx_ntks[MAX_VERTICES];
  std::vector<std::vector<short> >* vtx_tkind;
  std::vector<std::vector<float> >* vtx_tkweight;

 private:
  const char* nome;
  GlobeCuts *gCUT;

  edm::InputTag vertexColl;
  edm::InputTag trackColl;
  edm::InputTag bsColl;
  int debug_level;
};

#endif
