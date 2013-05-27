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
 
class RooBkg2l2jHL {
public:
  RooBkg2l2jHL();
  RooBkg2l2jHL(const char *name, const char *title,double myx,vector<double>allPars);
  RooBkg2l2jHL(const char *name, const char *title,
 	        double _x,
                double _acca0,
	        double _acca2,
	        double _acca4,
		double _g,
		double _cutOff);

  RooBkg2l2jHL(const RooBkg2l2jHL& other, const char* name=0) ;
  //  virtual TObject* clone(const char* newname) const { return new RooBkg2l2jHL(*this,newname); }
  inline virtual ~RooBkg2l2jHL() { }
  double evaluate() const ;
  void setVars(double newx);
  void setPars(vector<double> allPars);

protected:
 string name_;
  string title_;
  double x ;
  double acca0 ;
  double acca2 ;
  double acca4 ;
  double g ;
  double cutOff ;


private:

};
 
