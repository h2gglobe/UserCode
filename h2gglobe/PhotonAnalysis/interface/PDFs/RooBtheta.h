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

#ifndef ROO_BTHETA
#define ROO_BTHETA
#include <vector>
#include <string>

class RooBtheta{
public:
  RooBtheta();
  RooBtheta(const char *name, const char *title,double costheta,vector<double>allPars);
  RooBtheta(const char *name, const char *title,
 	    double _x,
            double _para2,
            double _para4,
	    double _acca2,
	    double _acca4);
  RooBtheta(const RooBtheta& other, const char* name=0) ;
  //  virtual TObject* clone(const char* newname) const { return new RooBtheta(*this,newname); }
  inline virtual ~RooBtheta() { }
  double evaluate() const ;
  void setVars(double newcostheta);
  void setPars(vector<double> allPars);

protected:
  string name_;
  string title_;
  double x ;
  double para2 ;
  double para4 ;
  double acca2 ;
  double acca4 ;
  


private:

};
 
#endif
