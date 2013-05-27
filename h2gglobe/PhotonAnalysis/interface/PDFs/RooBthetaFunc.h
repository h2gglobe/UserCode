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

#ifndef ROO_BTHETAFUNC
#define ROO_BTHETAFUNC

#include <vector>
#include <string>

class RooBthetaFunc{
public:
  RooBthetaFunc();
  RooBthetaFunc(const char *name, const char *title,vector<double> allVars,vector<double>allPars);
  RooBthetaFunc(const char *name, const char *title,
 	        double _x,
                double _y,
                double _acca0,
	        double _acca2,
	        double _acca4);
  RooBthetaFunc(const RooBthetaFunc& other, const char* name=0) ;
  // virtual TObject* clone(const char* newname) const { return new RooBthetaFunc(*this,newname); }
  inline virtual ~RooBthetaFunc() { }
  double evaluate() const ;
  void setVars(vector<double> allVars);
  void setPars(vector<double> allPars);

protected:
  string name_;
  string title_;
  double x ;
  double y ;
  double acca0 ;
  double acca2 ;
  double acca4 ;
  


private:

};
 
#endif
