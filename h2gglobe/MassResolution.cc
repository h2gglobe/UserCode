#include "MassResolution.h"

//----------------------------------------------------------//
// Project:   MassResolution
// Author:    Matt Kenzie (matthew.william.kenzie@cern.ch)
// Modified:  25/08/2011
// Admins:    Matth Kenzie (matthew.william.kenzie@cern.ch)
//---------------------------------------------------------//

/*
See MassResolution.h for instructions
*/

MassResolution::MassResolution(){}

void MassResolution::Setup(LoopAll &l, PhotonReducedInfo &leadInfo, PhotonReducedInfo &subleadInfo,int diphoton_index,EnergySmearer::energySmearingParameters eSmearPars, int nR9Categories, int nEtaCategories) 
//  leadPhoton(leadInfo),
//  subleadPhoton(subleadInfo)
{

  leadPhoton = leadInfo;	
  subleadPhoton = subleadInfo;	
  vertex = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_index]);
  vtx_dxdydz = (TVector3*)l.vtx_std_dxdydz->At(l.dipho_vtxind[diphoton_index]);

  lead_sc_pos    = leadPhoton.caloPosition();
  sublead_sc_pos = subleadPhoton.caloPosition();

  lead_Eres    = leadPhoton.corrEnergyErr();
  sublead_Eres = subleadPhoton.corrEnergyErr();

  lead_r9      = leadPhoton.r9();
  sublead_r9   = subleadPhoton.r9();

  TLorentzVector lead_p4=leadPhoton.p4(vertex->X(),vertex->Y(),vertex->Z());
  TLorentzVector sublead_p4=subleadPhoton.p4(vertex->X(),vertex->Y(),vertex->Z());

  higgsMass=(lead_p4+sublead_p4).M();

}
// return the mass resolution given correct vertex
double MassResolution::massResolutionCorrVtx(){

  TLorentzVector lead_p4=leadPhoton.p4(vertex->X(),vertex->Y(),vertex->Z());
  TLorentzVector sublead_p4=subleadPhoton.p4(vertex->X(),vertex->Y(),vertex->Z());

  double lead_E = lead_p4.E();
  double sublead_E = sublead_p4.E();
  double alpha = lead_p4.Angle(sublead_p4.Vect());
  double lead_sig = leadPhotonResolution();
  double sublead_sig = subleadPhotonResolution();
  double alpha_sig = angleResolutionCorrVtx();
  
  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

}
// return the mass resolution wrong vertex
double MassResolution::massResolutionWrongVtx(){
  
  TLorentzVector lead_p4=leadPhoton.p4(vertex->X(),vertex->Y(),vertex->Z());
  TLorentzVector sublead_p4=subleadPhoton.p4(vertex->X(),vertex->Y(),vertex->Z());

  double lead_E = lead_p4.E();
  double sublead_E = sublead_p4.E();
  double alpha = lead_p4.Angle(sublead_p4.Vect());
  double lead_sig = leadPhotonResolution();
  double sublead_sig = subleadPhotonResolution();
  double alpha_sig = angleResolutionWrongVtx();
  
  return 0.5*higgsMass*TMath::Sqrt(((lead_sig*lead_sig)/(lead_E*lead_E))+((sublead_sig*sublead_sig)/(sublead_E*sublead_E))+((alpha_sig*alpha_sig)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

}

// return energy contribution to mass resolution only
double MassResolution::massResolutionEonly() {

  TLorentzVector lead_p4=leadPhoton.p4(vertex->X(),vertex->Y(),vertex->Z());
  TLorentzVector sublead_p4=subleadPhoton.p4(vertex->X(),vertex->Y(),vertex->Z());

  double lead_E = lead_p4.E();
  double sublead_E = sublead_p4.E();
  double lead_sig = leadPhotonResolution();
  double sublead_sig = subleadPhotonResolution();

  return 0.5*higgsMass*TMath::Sqrt((lead_sig*lead_sig)/(lead_E*lead_E)+(sublead_sig*sublead_sig)/(sublead_E*sublead_E));
}


// return angle resolution given the vertex choice is correct
double MassResolution::angleResolutionCorrVtx() {
  return propagateDz(dzResolutionCorrVtx());
}

// return angle resolution given the vertex choice is wrong
double MassResolution::angleResolutionWrongVtx() {
  return propagateDz(dzResolutionWrongVtx());
}

// return angle resolution given a convolution of correct/wrong vertex as func of higgsPt
double MassResolution::angleResolution() {
  return propagateDz(dzResolution());
}
// return lead photon resolution without smearing
double MassResolution::leadPhotonResolutionNoSmear() {
  return lead_Eres;
}
// return sublead photon resolution without smearing
double MassResolution::subleadPhotonResolutionNoSmear() {
  return sublead_Eres;
}
// return lead photon resolution 
double MassResolution::leadPhotonResolution() {
  TLorentzVector lead_p4=leadPhoton.p4(vertex->X(),vertex->Y(),vertex->Z());
  bool sphericalLeadPhoton_=leadPhoton.isSphericalPhoton();
  return getPhotonResolution(lead_p4.E(),lead_Eres,lead_r9, lead_sc_pos.Eta(),lead_iDet,sphericalLeadPhoton_);
}
// return sublead photon resolution
double MassResolution::subleadPhotonResolution() {
  TLorentzVector sublead_p4=subleadPhoton.p4(vertex->X(),vertex->Y(),vertex->Z());
  bool sphericalSubleadPhoton_=subleadPhoton.isSphericalPhoton();
  return getPhotonResolution(sublead_p4.E(),sublead_Eres,sublead_r9,sublead_sc_pos.Eta(),sublead_iDet,sphericalSubleadPhoton_);
}

// Actually compute resolution given a photon
double MassResolution::getPhotonResolution(double photonEnergy, double photonResolution, double r9,  double scEta, bool iDet, bool ispherical) {


  // Get the photon-category sigma
  std::string myCategory="";
  if (_eSmearPars.categoryType=="Automagic") 
    {
	    EnergySmearer::energySmearingParameters::phoCatVectorConstIt vit = find(_eSmearPars.photon_categories.begin(), _eSmearPars.photon_categories.end(), 
										    std::make_pair( fabs((float)scEta), (float)r9 ) );
	    if( vit ==  _eSmearPars.photon_categories.end() ) {
		    //cout << "PhotonCheck " << r9 << " " << phoCat << " " << scEta << " " << iDet << " " << endl;
		    std::cerr << "Could not find energy scale correction for this photon " << (float)scEta << " " <<  (float)r9 << std::endl;
		    assert( 0 );
	    }
	    myCategory = vit->name;
    } 
  else if (_eSmearPars.categoryType=="2CatR9_EBEE")
    {
      if (iDet)
	myCategory+="EB";
      else
	myCategory+="EE";
      
      if (r9>=0.94)
	myCategory+="HighR9";
      else
	myCategory+="LowR9";
    }
  else if (_eSmearPars.categoryType=="2CatR9_EBEBm4EE")
    {
      if (iDet && fabs(scEta)      < 1.)
	myCategory+="EB";
      else if (iDet && fabs(scEta) > 1.)
	myCategory+="EBm4";
      else
	myCategory+="EE";
      
      if (r9>=0.94)
	myCategory+="HighR9";
      else
	myCategory+="LowR9";
    }
  else if (_eSmearPars.categoryType=="EBEE")
    {
      if (iDet)
	myCategory+="EB";
      else
	myCategory+="EE";
    }
  else
    {
      std::cout << "Unknown categorization. No category name is returned" << std::endl;
    }



    double categoryResolution = ispherical ? 0.0045*photonEnergy : _eSmearPars.smearing_sigma[myCategory]*photonEnergy;	
    return TMath::Sqrt(categoryResolution*categoryResolution + photonResolution*photonResolution);

}

//return dz resolution given correct vertex (used 10mm)
double MassResolution::dzResolutionCorrVtx() {
  return 0.1;
}
//return dz resolution given wrong vertex (using sqrt(2)*5.8cm)
double MassResolution::dzResolutionWrongVtx() {
  return TMath::Sqrt(2.)*5.8;
}
  
//return dz resolution from dz wrong and dz right (stored in TGraph as func of higgsPt)
double MassResolution::dzResolution() { 
  return dz;
}

// propagate error on z to error on angle
double MassResolution::propagateDz(double dz){

  TLorentzVector lead_p4=leadPhoton.p4(vertex->X(),vertex->Y(),vertex->Z());
  TLorentzVector sublead_p4=subleadPhoton.p4(vertex->X(),vertex->Y(),vertex->Z());

  double alpha = lead_p4.Angle(sublead_p4.Vect());
  if (alpha!= sublead_p4.Angle(lead_p4.Vect())) std::cout << "Error: Angle between photons not consistent" << std::endl;
  
  double x1 = lead_sc_pos.X()-vertex->X();
  double y1 = lead_sc_pos.Y()-vertex->Y();
  double z1 = lead_sc_pos.Z()-vertex->Z();
 
  double x2 = sublead_sc_pos.X()-vertex->X();
  double y2 = sublead_sc_pos.Y()-vertex->Y();
  double z2 = sublead_sc_pos.Z()-vertex->Z();
 
  double r1 = TMath::Sqrt(x1*x1+y1*y1+z1*z1);
  double r2 = TMath::Sqrt(x2*x2+y2*y2+z2*z2);

  double cos_term = TMath::Cos(lead_p4.Phi()-sublead_p4.Phi());
  double sech1 = SecH(lead_p4.Eta());
  double sech2 = SecH(sublead_p4.Eta());
  double tanh1 = TanH(lead_p4.Eta());
  double tanh2 = TanH(sublead_p4.Eta());

  double numerator1 = sech1*(sech1*tanh2-tanh1*sech2*cos_term);
  double numerator2 = sech2*(sech2*tanh1-tanh2*sech1*cos_term);
  double denominator = 1. - tanh1*tanh2-sech1*sech2*cos_term;

  double ResTerm = (-1.*dz/denominator)*(numerator1/r1 + numerator2/r2);

  double angleResolution = ResTerm*(1.-TMath::Cos(alpha))/TMath::Sin(alpha);

  return angleResolution;

}

// utility functions
double MassResolution::SecH(double x){
  return 1.0/TMath::CosH(x);
}

double MassResolution::TanH(double x){
  return TMath::TanH(x);
}

