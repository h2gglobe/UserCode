#include <Riostream.h>
#include <math.h>
#include "RooBphione.h"

// ClassImp(RooBphione) 
RooBphione::RooBphione():
  name_("bphi1"),
  title_("bphi1title"),
  x(0.0),
  para2(0.0),
  para4(0.0),
  acca2(0.0),
  acca4(0.0),
  acca6(0.0),
  acca8(0.0)
{
}

RooBphione::RooBphione(const char *name, const char *title,double myphi1,vector<double>allPars):
  name_(name),
  title_(title)
{
  setVars( myphi1);
  setPars(allPars);
}
 RooBphione::RooBphione(const char *name, const char *title, 
			double _x,
			double _para2,
			double _para4,
			double _acca2,
			double _acca4,
			double _acca6,
			double _acca8) :
   name_(name),
   title_(title),
   x(_x),
   para2(_para2),
   para4(_para4),
   acca2(_acca2),
   acca4(_acca4),
   acca6(_acca6),
   acca8(_acca8)
 { 
 } 


 RooBphione::RooBphione(const RooBphione& other, const char* name) :  
   name_(name),
   title_("dummytitle"),
   x(other.x),
   para2(other.para2),
   para4(other.para4),
   acca2(other.acca2),
   acca4(other.acca4),
   acca6(other.acca6),
   acca8(other.acca8)
 { 
 } 



 double RooBphione::evaluate() const 
 { 
   double accp = acca2*cos(x)+acca4*cos(2.0*x)+acca6*cos(4.0*x)+acca8;

   double heli = 1.0+para2*pow(x,2)+para4*pow(x,4);

   double heliCorr = heli*accp;

   if((heliCorr<=0.0)) heliCorr = 0.0;

   return heliCorr;

 } 

void RooBphione::setVars(double newphi1){
  x=newphi1;
}
void RooBphione::setPars(vector<double> allPars){

  if(int(allPars.size()) < 6){
    std::cout<<"ERROR in RooBphione::setPars! Size of input params less than minimum! "<<allPars.size()<<" < 6"<<std::endl;

  }
  para2=allPars.at(0) ;
  para4=allPars.at(1) ;
  acca2=allPars.at(2) ;
  acca4=allPars.at(3) ;
  acca6=allPars.at(4) ;
  acca8=allPars.at(5) ;

}
