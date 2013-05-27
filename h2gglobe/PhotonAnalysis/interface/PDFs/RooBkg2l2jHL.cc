#include <Riostream.h>
#include <math.h>

#include "RooBkg2l2jHL.h"


RooBkg2l2jHL::RooBkg2l2jHL():
  name_("bkg2l2j1"),
  title_("bkg2l2jtitle"),
  x(0.0),
  acca0(0.0),
  acca2(0.0),
  acca4(0.0),
  g(0.0),
  cutOff(0.0)
{
}

RooBkg2l2jHL::RooBkg2l2jHL(const char *name, const char *title,double myx,vector<double>allPars):
  name_(name),
  title_(title)
{
  setVars( myx);
  setPars(allPars);
}


RooBkg2l2jHL::RooBkg2l2jHL(const char *name, const char *title, 
			      double _x,
			      double _acca0,
			      double _acca2,
			      double _acca4,
			      double _g,
			      double _cutOff) :
  name_(name),
   title_(title),                       
   x(_x),
   acca0(_acca0),
   acca2(_acca2),
   acca4(_acca4),
   g(_g),		 
   cutOff(_cutOff) 
 { 
 } 


 RooBkg2l2jHL::RooBkg2l2jHL(const RooBkg2l2jHL& other, const char* name) :  
   name_(name),
   x(other.x),
   acca0(other.acca0),
   acca2(other.acca2),
   acca4(other.acca4),
   g(other.g),		 
   cutOff(other.cutOff) 
 { 
 } 



 double RooBkg2l2jHL::evaluate() const 
 { 
   return (1.0+acca0*x*x)*(1.0+acca2*x*x+acca4*x*x*x*x)/(1 + exp((-cutOff + x)/g));
 } 

void RooBkg2l2jHL::setVars(double newx){
  x=newx;
}
void RooBkg2l2jHL::setPars(vector<double> allPars){

  if(int(allPars.size()) < 5){
    std::cout<<"ERROR in RooBkg2l2jHL::setPars! Size of input params less than minimum! "<<allPars.size()<<" < 6"<<std::endl;

  }
  acca0=allPars.at(0) ;
  acca2=allPars.at(1) ;
  acca4=allPars.at(2) ;
  g=allPars.at(3);
  cutOff=allPars.at(4);
}
