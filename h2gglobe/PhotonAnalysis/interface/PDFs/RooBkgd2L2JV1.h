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

#ifndef ROO_BKGD2L2JV1
#define ROO_BKGD2L2JV1

#include <vector>
#include <string>

using namespace std;

class RooBkgd2L2JV1{
public:
  RooBkgd2L2JV1();
  RooBkgd2L2JV1(const char *name, const char *title, 
		double mycostheta1, double mycostheta2, 
		double mycosthetas, double myphi, 
		double myphi1, double mZZ, 
		vector<double>allPars);
  RooBkgd2L2JV1(const char *name, const char *title,
	     double _h1,	  
	     double _h2,	  
	     double _hs, 	  
	     double _phi,	  
	     double _phi1,	  
	     double _mZZ,	  
	     double _h1param0,
	     double _h1param2,
	     double _h1param4,
	     double _h2param0,
             double _h2param2,
             double _h2param4,
             double _h2g,	  
             double _h2cutOff,
	     double _hsparam2,
             double _hsparam4,
             double _hsacca2,
             double _hsacca4, 
	     double _pacca2, 
             double _pacca4, 
             double _pacca6, 
             double _pacca8, 
             double _p1acca2,
             double _p1acca4,
	     double _p1acca6,
             double _p1acca8 	     
	     );

  RooBkgd2L2JV1(const RooBkgd2L2JV1& other, const char* name=0) ;
  //virtual TObject* clone(const char* newname) const { return new RooBkgd2L2JV1(*this,newname); }
  inline virtual ~RooBkgd2L2JV1() { }
  double evaluate() const;
  void setVars(vector<double> allVars);
  void setVars(double newcostheta1, double newcostheta2, 
	       double newcosthetas, double newphi,
	       double newphi1, double mZZ);
  void setPars(vector<double> allPars);
  void SetParameters();

protected:

  string name_;
  string title_;
  double h1;     	  
  double h2;     	  
  double hs;     	  
  double phi;    	  
  double phi1;   	  
  double mZZ;    
  double h1param0;
  double h1param2;
  double h1param4;
  double h2param0;
  double h2param2;
  double h2param4;
  double h2g;    
  double h2cutOff;
  double hsparam2;
  double hsparam4;
  double hsacca2;
  double hsacca4;
  double pacca2; 
  double pacca4; 
  double pacca6; 
  double pacca8; 
  double p1acca2;
  double p1acca4;
  double p1acca6;
  double p1acca8;
  
  double h1param0_0;
  double h1param0_1;
  double h1param0_2;

  double h1param2_0;
  double h1param2_1;
  double h1param2_2;
  
  double h1param4_0;
  double h1param4_1;
  double h1param4_2;

  double h2param0_0;
  double h2param0_1;
  double h2param0_2;


  double h2param2_0;
  double h2param2_1;
  double h2param2_2;

  double h2param4_0;
  double h2param4_1;
  double h2param4_2;

  double h2g_0;    
  double h2g_1;    
  double h2g_2;    

  double h2cutOff_0;
  double h2cutOff_1;
  double h2cutOff_2;

  double hsparam2_0;
  double hsparam2_1;
  double hsparam2_2;

  double hsparam4_0;
  double hsparam4_1;
  double hsparam4_2;

  double hsacca2_0;
  double hsacca2_1;
  double hsacca2_2;

  double hsacca4_0;
  double hsacca4_1;
  double hsacca4_2;

  double pacca2_0; 
  double pacca2_1; 
  double pacca2_2; 

  double pacca4_0; 
  double pacca4_1; 
  double pacca4_2; 

  double pacca6_0; 
  double pacca6_1; 
  double pacca6_2; 

  double pacca8_0; 
  double pacca8_1; 
  double pacca8_2; 

  double p1acca2_0;
  double p1acca2_1;
  double p1acca2_2;

  double p1acca4_0;
  double p1acca4_1;
  double p1acca4_2;

  double p1acca6_0;
  double p1acca6_1;
  double p1acca6_2;

  double p1acca8_0;
  double p1acca8_1;
  double p1acca8_2;

private:

  //ClassDef(RooBkgd2L2JV1,0) // Your description goes here...
};
 
#endif
