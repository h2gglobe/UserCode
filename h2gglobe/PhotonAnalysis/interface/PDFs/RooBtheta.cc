#include <Riostream.h>
#include <math.h>

#include "RooBtheta.h"


RooBtheta::RooBtheta():
  name_("btheta"),
  title_("bthetatitle"),
  x(0.0),
  para2(0.0),
  para4(0.0),
  acca2(0.0),
  acca4(0.0)
{
}

RooBtheta::RooBtheta(const char *name, const char *title,double costheta,vector<double>allPars):
  name_(name),
  title_(title)
{
  setVars( costheta);
  setPars(allPars);
}

 RooBtheta::RooBtheta(const char *name, const char *title, 
		      double _x,
		      double _para2,
                      double _para4,
		      double _acca2,
		      double _acca4) :
   name_(name),
   title_(title),
   x(_x),
   para2(_para2),
   para4(_para4),
   acca2(_acca2),
   acca4(_acca4)
 { 
 } 


 RooBtheta::RooBtheta(const RooBtheta& other, const char* name) :  
   name_(name),
   x(other.x),
   para2(other.para2),
   para4(other.para4),
   acca2(other.acca2),
   acca4(other.acca4)
 { 
 } 



 double RooBtheta::evaluate() const 
 { 
   return (1.0+para2*x*x+para4*x*x*x*x)*(1.0+acca2*x*x+acca4*x*x*x*x);
 } 

void RooBtheta::setVars(double newcostheta){
  x=newcostheta;
}

void RooBtheta::setPars(vector<double> allPars){

  if(int(allPars.size()) <4 ){
    std::cout<<"ERROR in RooBphione::setPars! Size of input params less than minimum! "<<allPars.size()<<" < 4"<<std::endl;

  }
  para2=allPars.at(0) ;
  para4=allPars.at(1) ;
  acca2=allPars.at(2) ;
  acca4=allPars.at(3) ;

}
