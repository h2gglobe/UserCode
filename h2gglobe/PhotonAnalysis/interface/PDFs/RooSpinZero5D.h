/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


#include <vector>
#include <string>

using namespace std;


class RooSpinZero5D {
public:
  RooSpinZero5D();
  RooSpinZero5D(const char *name, const char *title,
	   vector<double>allVars, 
	   vector<double>allPars);
  RooSpinZero5D(const char *name, const char *title,
	   double _h1,
	   double _h2,
	   double _hs,
	   double _Phi,
	   double _Phi1,
	   double _mZZ,
	   double _fppVal,
	   double _fmmVal,
	   double _fpmVal,
	   double _fp0Val,
	   double _f0mVal,
	   double _phippVal,
	   double _phimmVal,
	   double _phipmVal,
	   double _phip0Val,
	   double _phi0mVal,
	   double _fz1Val,
	   double _fz2Val,
	   double _R1Val,
	   double _R2Val,
	   double _para2,
	   double _para4,
	   double _para6,
	   double _para8,
	   double _acca0,
	   double _acca1,
	   double _acca2,
	   double _acca4,
	   double _a2,
	   double _a4,
	   double _cutOff,
	   double _g,
	   double _b2,
	   double _b4,
	   double _N);

  RooSpinZero5D(const RooSpinZero5D& other, const char* name=0) ;
  //virtual TObject* clone(const char* newname) const { return new RooSpinZero5D(*this,newname); }
  inline virtual ~RooSpinZero5D() { }
  //Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  //Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  double evaluate() const;
  void SetAcceptanceParameters();
  void setVars(vector<double> allVars);
  void setVars(double newcostheta1,double newcostheta2,
	       double newcosthetas,double newphi,
	       double newphi1, double newmZZ);
  void setParams(vector<double> allPars);

protected:

  string name_;
  string title_;
  double h1 ;
  double h2 ;
  double hs ;
  double Phi ;
  double Phi1 ;
  double mZZ ; 
  double fppVal ;
  double fmmVal ;
  double fpmVal ;
  double fp0Val ;
  double f0mVal ;
  double phippVal ;
  double phimmVal ;
  double phipmVal ;
  double phip0Val ;
  double phi0mVal ;
  double fz1Val ;
  double fz2Val ;
  double R1Val ;
  double R2Val ;
  double para2 ;
  double para4 ;
  double para6 ;
  double para8 ;
  double acca0 ;
  double acca1 ;
  double acca2 ;
  double acca4 ;
  double a2 ;
  double a4 ;
  double cutOff ;
  double g ;
  double b2 ;
  double b4 ;
  double N  ;
  //----------------------
  double para2_0;
  double para2_1;
  double para2_2;
  double para4_0;
  double para4_1;
  double para4_2;
  double para6_0;
  double para6_1;
  double para6_2;
  double para8_0;
  double para8_1;
  double para8_2;
  double acca0_0;
  double acca0_1;
  double acca0_2;
  double acca1_0;
  double acca1_1;
  double acca1_2;
  double acca2_0;
  double acca2_1;
  double acca2_2;
  double acca4_0;
  double acca4_1;
  double acca4_2;
  double a2_0;
  double a2_1;
  double a2_2;
  double a4_0;
  double a4_1;
  double a4_2;
  double cutOff_0;
  double cutOff_1;
  double cutOff_2;
  double g_0;
  double g_1;
  double g_2;
  double b2_0;
  double b2_1;
  double b2_2;
  double b4_0;
  double b4_1;
  double b4_2;
  //----------------------


private:

  //ClassDef(RooSpinZero5D,0) // Your description goes here...
};
 

