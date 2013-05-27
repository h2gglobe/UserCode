#include <Riostream.h>
#include <math.h>
#include "RooBkgd2L2JV1.h"

using namespace std;

//ClassImp(RooBkgd2L2JV1) 

RooBkgd2L2JV1::RooBkgd2L2JV1():
  name_("background"),
  title_("backgroundtitle"),
  h1(0.0),
  h2(0.0),
  hs(0.0),
  phi(0.0),
  phi1(0.0),
  mZZ(0.0),
  h1param0(0.0),
  h1param2(0.0),
  h1param4(0.0),
  h2param0(0.0),
  h2param2(0.0),
  h2param4(0.0),
  h2g(0.0),
  h2cutOff(0.0),
  hsparam2(0.0),
  hsparam4(0.0),
  hsacca2(0.0),
  hsacca4(0.0),
  pacca2(0.0),
  pacca4(0.0),
  pacca6(0.0),
  pacca8(0.0),
  p1acca2(0.0),
  p1acca4(0.0),
  p1acca6(0.0),
  p1acca8(0.0)
{
}

RooBkgd2L2JV1::RooBkgd2L2JV1(const char *name, const char *title,
			     double mycostheta1, double mycostheta2,
			     double mycosthetas, double myphi1, 
			     double myphi, double mymZZ,
			     vector<double>allPars):
  name_(name),
  title_(title)
{
  setVars(mycostheta1,mycostheta2,mycosthetas,myphi,myphi1,mymZZ);
  setPars(allPars);
}

RooBkgd2L2JV1::RooBkgd2L2JV1(const char *name, const char *title, 
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
             double _p1acca8):                              
  name_(name),
  title_(title),
  h1(_h1),
  h2(_h2),
  hs(_hs),
  phi(_phi),
  phi1(_phi1),
  mZZ(_mZZ),
  h1param0(_h1param0),
  h1param2(_h1param2),
  h1param4(_h1param4),
  h2param0(_h2param0),
  h2param2(_h2param2),
  h2param4(_h2param4),
  h2g(_h2g),
  h2cutOff(_h2cutOff),
  hsparam2(_hsparam2),
  hsparam4(_hsparam4),
  hsacca2(_hsacca2),
  hsacca4(_hsacca4),
  pacca2(_pacca2),
  pacca4(_pacca4),
  pacca6(_pacca6),
  pacca8(_pacca8),
  p1acca2(_p1acca2),
  p1acca4(_p1acca4),
  p1acca6(_p1acca6),
  p1acca8(_p1acca8)
{ 
} 


 RooBkgd2L2JV1::RooBkgd2L2JV1(const RooBkgd2L2JV1& other, const char* name) :  
   name_(name),
   title_("dummytitle"),
   h1(other.h1),
   h2(other.h2),
   hs(other.hs),
   phi(other.phi),
   phi1(other.phi1),
   mZZ(other.mZZ),
   h1param0(other.h1param0),
   h1param2(other.h1param2),
   h1param4(other.h1param4),
   h2param0(other.h2param0),
   h2param2(other.h2param2),
   h2param4(other.h2param4),
   h2g(other.h2g),
   h2cutOff(other.h2cutOff),
   hsparam2(other.hsparam2),
   hsparam4(other.hsparam4),
   hsacca2(other.hsacca2),
   hsacca4(other.hsacca4),
   pacca2(other.pacca2),
   pacca4(other.pacca4),
   pacca6(other.pacca6),
   pacca8(other.pacca8),
   p1acca2(other.p1acca2),
   p1acca4(other.p1acca4),
   p1acca6(other.p1acca6),
   p1acca8(other.p1acca8)
{ 
} 



 double RooBkgd2L2JV1::evaluate() const 
 { 
   double h1Func = (1.0+h1param0*h1*h1)*(1.0+h1param2*h1*h1+h1param4*h1*h1*h1*h1);
   //   std::cout << "h1Func: " << h1Func << endl;
   double h2Func = (1.0+h2param0*h2*h2)*(1.0+h2param2*h2*h2+h2param4*h2*h2*h2*h2)
     /(1 + exp((-h2cutOff + h2)/h2g));
   //   std::cout <<"h2Func: " << h2Func <<endl;
   double hsFunc =  (1.0+hsparam2*hs*hs+hsparam4*hs*hs*hs*hs)*(1.0+hsacca2*hs*hs+hsacca4*hs*hs*hs*hs);
   //   std::cout <<"hsFunc: " << hsFunc <<endl;
   double phiFunc = (pacca2*cos(phi)+pacca4*cos(2.0*phi)+pacca6*cos(4.0*phi)+pacca8);
   //   std::cout <<"phiFunc: " << phiFunc <<endl;
   double phi1Func= (p1acca2*cos(phi1)+p1acca4*cos(2.0*phi1)+p1acca6*cos(4.0*phi1)+p1acca8);
   //   std::cout <<"phi1Func: " << phi1Func <<endl;

   double h1_Norm = 2*(1 + (h1param0+h1param2)/3 + (h1param4+h1param0*h1param2)/5
		       + h1param0*h1param4/7);
   //   std::cout <<"h1Norm: " << h1_Norm <<endl;
   double h2_Norm = (1 + (h2param0+h2param2)/3 + (h2param4+h2param0*h2param2)/5
		     + h2param0*h2param4/7);
   //   std::cout <<"h2Norm: " << h2_Norm<<endl;
   double hs_Norm = 2*(1 + (hsacca2+hsparam2)/3 + (hsacca4+hsparam2*hsacca2+hsparam4)/5
		       + (hsparam2*hsacca4+hsparam4*hsacca2)/7 + hsparam4*hsacca4/9);
   //   std::cout <<"hsNorm: " << hs_Norm <<endl;
   double phi_Norm = 2*3.1415*pacca8;
   //   std::cout <<"phiNorm: " << phi_Norm <<endl;
   double phi1_Norm= 2*3.1415*p1acca8;
   //   std::cout <<"phi1Norm: " << phi1_Norm <<endl;

   if(h2Func<0){ 
     //std::cout << "JUST SO YOU KNOW: h2FUNC was less than zero, specifically: " << h2Func << " at h2: " << h2 << std::endl;
     h2Func=.0001;
   }
   if(hsFunc<0){
     //std::cout << "JUST SO YOU KNOW: hsFUNC was less than zero, specifically: " << hsFunc << std::endl;
     hsFunc=.0001;
   }  
   if(h1Func<0){
     //std::cout << "JUST SO YOU KNOW: h1FUNC was less than zero, specifically: " << h1Func << std::endl;
     h1Func=.0001;
   }
   if(phiFunc<0){
     //std::cout << "JUST SO YOU KNOW: phiFUNC was less than zero, specifically: " << phiFunc << std::endl;
     phiFunc=.0001;
   }
   if(phi1Func<0){
     //std::cout << "JUST SO YOU KNOW: phi1FUNC was less than zero, specifically: " << phi1Func << std::endl;
     phi1Func=.0001;
   }   
   //std::cout << "Bkg Prod: " << h1Func*h2Func*hsFunc*phiFunc*phi1Func/(h1_Norm*h2_Norm*hs_Norm*phi_Norm*phi1_Norm) << endl;
   return(h1Func*h2Func*hsFunc*phiFunc*phi1Func/(h1_Norm*h2_Norm*hs_Norm*phi_Norm*phi1_Norm));
 } 


void RooBkgd2L2JV1::setVars(vector<double> allVars){
  if(int(allVars.size())<6){
    std::cout<<"ERROR from RooBkgd2L2JV1::setVars. Size of vector with measurables is less than 6: "<<allVars.size()<<std::endl;
  }

  h1=allVars.at(0);
  h2=allVars.at(1);
  hs=allVars.at(2);
  phi=allVars.at(3);
  phi1=allVars.at(4);
  mZZ=allVars.at(5);
}


void RooBkgd2L2JV1::setVars(double newcostheta1, double newcostheta2,
			    double newcosthetas, double newphi,
			    double newphi1, double newmZZ){
  h1=newcostheta1;
  h2=newcostheta2;
  hs=newcosthetas;
  phi=newphi;
  phi1=newphi1;
  mZZ=newmZZ;
}

void RooBkgd2L2JV1::setPars(vector<double> allPars){
  if(int(allPars.size())<20){
    std::cout<<"ERROR in RooBkgd2L2JV1::setPars! Size of input params less than minimum: " 
	     << allPars.size() << std::endl;
  }

  h1param0=allPars.at(0); 
  h1param2=allPars.at(1);
  h1param4=allPars.at(2);
  h2param0=allPars.at(3);
  h2param2=allPars.at(4);
  h2param4=allPars.at(5); 
  h2g=allPars.at(6);
  h2cutOff=allPars.at(7);
  hsparam2=allPars.at(8);
  hsparam4=allPars.at(9);
  hsacca2=allPars.at(10);
  hsacca4=allPars.at(11);
  pacca2=allPars.at(12);
  pacca4=allPars.at(13);   
  pacca6=allPars.at(14);
  pacca8=allPars.at(15);
  p1acca2=allPars.at(16);
  p1acca4=allPars.at(17);
  p1acca6=allPars.at(18);
  p1acca8=allPars.at(19);
	    
}  


void RooBkgd2L2JV1::SetParameters(){

  h1param0_0=-0.0852538;
  h1param0_1=0.002857;
  h1param0_2=0.;
  
  h1param2_0=-0.26096;
  h1param2_1=0.000886;
  h1param2_2=0.;
  
  h1param4_0=-0.272046;
  h1param4_1=-0.001431;
  h1param4_2=0.;

  h2param0_0=-.941379;
  h2param0_1=-.000114;
  h2param0_2=0.;

  h2param2_0=-2.54539;
  h2param2_1=0.01936;
  h2param2_2=0.;

  h2param4_0=.26186;
  h2param4_1=-0.01569;
  h2param4_2=0.;

  h2g_0=.210251;    
  h2g_1=-.0000388;    
  h2g_2=0.0;    

  h2cutOff_0=.68897;
  h2cutOff_1=.000250291;
  h2cutOff_2=0.;

  hsparam2_0=-0.11984;
  hsparam2_1=-.0001719;
  hsparam2_2=0.;

  hsparam4_0=-1.3136;
  hsparam4_1=.009102;
  hsparam4_2=0.;

  hsacca2_0=-.11984;
  hsacca2_1=-.00017198;
  hsacca2_2=0.;

  hsacca4_0=-1.3136;
  hsacca4_1=.009102;
  hsacca4_2=0.;

  pacca2_0=17.613; 
  pacca2_1=-0.06057; 
  pacca2_2=0.; 

  pacca4_0=-119.16; 
  pacca4_1=.235575; 
  pacca4_2=0.; 

  pacca6_0=-100.991; 
  pacca6_1=.252066; 
  pacca6_2=0.; 

  pacca8_0=57.7097; 
  pacca8_1=2.0291; 
  pacca8_2=0.; 

  p1acca2_0=249.457;
  p1acca2_1=-1.41561;
  p1acca2_2=.00193596;

  p1acca4_0=99.1997;
  p1acca4_1=-1.1647;
  p1acca4_2=.00133026;

  p1acca6_0=-383.76;
  p1acca6_1=1.80341;
  p1acca6_2=-.001999;

  p1acca8_0=10000.;  
  p1acca8_1=0.;
  p1acca8_2=0.;

  p1acca2 = p1acca2_0 + p1acca2_1*mZZ + p1acca2_2*mZZ*mZZ;
  p1acca4 = p1acca4_0 + p1acca4_1*mZZ + p1acca4_2*mZZ*mZZ;
  p1acca6 = p1acca6_0 + p1acca6_1*mZZ + p1acca6_2*mZZ*mZZ;
  p1acca8 = p1acca8_0 + p1acca8_1*mZZ + p1acca8_2*mZZ*mZZ;

  pacca2 = pacca2_0 + pacca2_1*mZZ + pacca2_2*mZZ*mZZ; 
  pacca4 = pacca4_0 + pacca4_1*mZZ + pacca4_2*mZZ*mZZ; 
  pacca6 = pacca6_0 + pacca6_1*mZZ + pacca6_2*mZZ*mZZ; 
  pacca8 = pacca8_0 + pacca8_1*mZZ + pacca8_2*mZZ*mZZ; 

  hsacca2 = hsacca2_0 + hsacca2_1*mZZ + hsacca2_2*mZZ*mZZ;
  hsacca4 = hsacca4_0 + hsacca4_1*mZZ + hsacca4_2*mZZ*mZZ;

  hsparam2 = hsparam2_0 + hsparam2_1*mZZ + hsparam2_2*mZZ*mZZ;
  hsparam4 = hsparam4_0 + hsparam4_1*mZZ + hsparam4_2*mZZ*mZZ;

  h2cutOff = h2cutOff_0 + h2cutOff_1*mZZ + h2cutOff_2*mZZ*mZZ;
  h2g = h2g_0 + h2g_1*mZZ + h2g_2*mZZ*mZZ;

  h2param0 = h2param0_0 + h2param0_1*mZZ + h2param0_2*mZZ*mZZ;
  h2param2 = h2param2_0 + h2param2_1*mZZ + h2param2_2*mZZ*mZZ;
  h2param4 = h2param4_0 + h2param4_1*mZZ + h2param4_2*mZZ*mZZ;

  h1param0 = h1param0_0 + h1param0_1*mZZ + h1param0_2*mZZ*mZZ;
  h1param2 = h1param2_0 + h1param2_1*mZZ + h1param2_2*mZZ*mZZ;
  h1param4 = h1param4_0 + h1param4_1*mZZ + h1param4_2*mZZ*mZZ;

}
