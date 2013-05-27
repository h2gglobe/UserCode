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

#ifndef ROO_BPHIONE
#define ROO_BPHIONE

#include <vector>
#include <string>

 
class RooBphione {
public:
  RooBphione();
  RooBphione(const char *name, const char *title,double myphi1,vector<double>allPars);
  RooBphione(const char *name, const char *title,
    	     double _x,
             double _para2,
             double _para4,
             double _acca2,
             double _acca4,
             double _acca6,
	     double _acca8);
  RooBphione(const RooBphione& other, const char* name=0) ;
  //  virtual TObject* clone(const char* newname) const { return new RooBphione(*this,newname); }
  inline virtual ~RooBphione() { }
  double evaluate() const ;
  void setVars(double newphi1);
  void setPars(vector<double> allPars);

protected:
  string name_;
  string title_;
  double x ;//phistarone for background
  double para2 ;
  double para4 ;
  double acca2 ;
  double acca4 ;
  double acca6 ;
  double acca8 ;
  


private:

  //  ClassDef(RooBphione,0) // Your description goes here...
};
 
#endif
