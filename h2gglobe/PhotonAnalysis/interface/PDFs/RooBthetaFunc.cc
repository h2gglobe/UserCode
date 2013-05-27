#include <Riostream.h>
#include <math.h>

#include "RooBthetaFunc.h"


// ClassImp(RooBthetaFunc) 

RooBthetaFunc::RooBthetaFunc():
  name_("bthetafunc"),
  title_("bthetafunctitle"),
  x(0.0),
  acca0(0.0),
  acca2(0.0),
  acca4(0.0)
{
}

RooBthetaFunc::RooBthetaFunc(const char *name, const char *title,vector<double>allVars,vector<double>allPars):
  name_(name),
  title_(title)
{
  setVars(allVars);
  setPars(allPars);
}


 RooBthetaFunc::RooBthetaFunc(const char *name, const char *title, 
			      double _x,
			      double _y,
			      double _acca0,
			      double _acca2,
			      double _acca4) :
   name_(name),
   title_(title),
   x(_x),
   y(_y),
   acca0(_acca0),
   acca2(_acca2),
   acca4(_acca4)
 { 
 } 


 RooBthetaFunc::RooBthetaFunc(const RooBthetaFunc& other, const char* name) :  
   name_(name),
   x(other.x),
   y(other.y),
   acca0(other.acca0),
   acca2(other.acca2),
   acca4(other.acca4)
 { 
 } 



 double RooBthetaFunc::evaluate() const 
 { 
   return (1.0+acca0*x*x)*(1.0+acca0*y*y)
         *(1.0+acca2*x*x+acca4*x*x*x*x)*(1.0+acca2*y*y+acca4*y*y*y*y);
 } 


void RooBthetaFunc::setVars(vector<double> allVars){
  x=allVars.at(0);
  y=allVars.at(1);
}

void RooBthetaFunc::setPars(vector<double> allPars){

  if(int(allPars.size()) <3 ){
    std::cout<<"ERROR in RooBphione::setPars! Size of input params less than minimum! "<<allPars.size()<<" < 4"<<std::endl;

  }
  acca0=allPars.at(0) ;
  acca2=allPars.at(1) ;
  acca4=allPars.at(2) ;

}
