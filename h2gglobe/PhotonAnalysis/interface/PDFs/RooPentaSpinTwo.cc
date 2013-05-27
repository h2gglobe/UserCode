#include <Riostream.h>
#include <math.h>

#include "RooPentaSpinTwo.h"




RooPentaSpinTwo::RooPentaSpinTwo():
  name_("Penta1"),
   title_("Penta1Title"),
   h1(0.),
   h2(0.0),
   Phi(0.0),
   hs(0.0),
   Phi1(0.0),
   fppVal(0.0),
   fmmVal(0.0),
   fpmVal(0.0),
   fp0Val(0.0),
   f0mVal(0.0),
   phippVal(0.0),
   phimmVal(0.0),
   phipmVal(0.0),
   phip0Val(0.0),
   phi0mVal(0.0),
   fz1Val(0.0),
   fz2Val(0.0),
   R1Val(0.0),
   R2Val(0.0),
   para2(0.0),
   para4(0.0),
   para6(0.0),
   para8(0.0),
   acca0(0.0),
   acca1(0.0),
   acca2(0.0),
   acca4(0.0),
   a2(0.0),
   a4(0.0),
   cutOff(0.0),
   g(0.0),
   b2(0.0),
   b4(0.0),
   N(0.0)
{
}



RooPentaSpinTwo::RooPentaSpinTwo(const char *name, const char *title, vector<double> allVars, vector<double> allPars):
  name_(name),
  title_(title){

  setVars(allVars);
  setParams(allPars);
}

 RooPentaSpinTwo::RooPentaSpinTwo(const char *name, const char *title, 
				  double _h1,
				  double _h2,
				  double _Phi,
                                  double _hs,
				  double _Phi1,
				  double _fppVal,
				  double _fmmVal,
				  double _fpmVal,
				  double _fp0Val,
				  double _f0mVal,
				  double _phippVal,
				  double _phimmVal,
				  double _phipmVal,
				  double _phip0Val,
				  double _phi0mVal,
				  double _fz1Val,
				  double _fz2Val,
				  double _R1Val,
				  double _R2Val,
				  double _para2,
				  double _para4,
				  double _para6,
				  double _para8,
				  double _acca0,
				  double _acca1,
				  double _acca2,
				  double _acca4,
				  double _a2,
				  double _a4,
				  double _cutOff,
				  double _g,
                                  double _b2,
				  double _b4,
				  double _N):  
   name_(name),
   title_(title),
   h1(_h1),
   h2(_h2),
   Phi(_Phi),
   hs(_hs),
   Phi1(_Phi1),
   fppVal(_fppVal),
   fmmVal(_fmmVal),
   fpmVal(_fpmVal),
   fp0Val(_fp0Val),
   f0mVal(_f0mVal),
   phippVal(_phippVal),
   phimmVal(_phimmVal),
   phipmVal(_phipmVal),
   phip0Val(_phip0Val),
   phi0mVal(_phi0mVal),
   fz1Val(_fz1Val),
   fz2Val(_fz2Val),
   R1Val(_R1Val),
   R2Val(_R2Val),
   para2(_para2),
   para4(_para4),
   para6(_para6),
   para8(_para8),
   acca0(_acca0),
   acca1(_acca1),
   acca2(_acca2),
   acca4(_acca4),
   a2(_a2),
   a4(_a4),
   cutOff(_cutOff),
   g(_g),
   b2(_b2),
   b4(_b4),
   N(_N)
 { 
 } 


 RooPentaSpinTwo::RooPentaSpinTwo(const RooPentaSpinTwo& other, const char* name) :  
   name_(name),
   h1(other.h1),
   h2(other.h2),
   Phi(other.Phi),
   hs(other.hs),
   Phi1(other.Phi1),
   fppVal(other.fppVal),
   fmmVal(other.fmmVal),
   fpmVal(other.fpmVal),
   fp0Val(other.fp0Val),
   f0mVal(other.f0mVal),
   phippVal(other.phippVal),
   phimmVal(other.phimmVal),
   phipmVal(other.phipmVal),
   phip0Val(other.phip0Val),
   phi0mVal(other.phi0mVal),
   fz1Val(other.fz1Val),
   fz2Val(other.fz2Val),
   R1Val(other.R1Val),
   R2Val(other.R2Val),
   para2(other.para2),
   para4(other.para4),
   para6(other.para6),
   para8(other.para8),
   acca0(other.acca0),
   acca1(other.acca1),
   acca2(other.acca2),
   acca4(other.acca4),
   a2(other.a2),
   a4(other.a4),
   cutOff(other.cutOff),
   g(other.g),
   b2(other.b2),
   b4(other.b4),
   N(other.N)
 {
 } 

void RooPentaSpinTwo::setVars(vector<double> allVars){

  if(int(allVars.size())<5){
    std::cout<<"ERROR from RooPentaSpinTwo::setVars. Size of vector with measurables is less than 5: "<<allVars.size()<<std::endl;
  }

  h1=  allVars.at(0);
  h2=  allVars.at(1);
  Phi= allVars.at(2);
  hs=  allVars.at(3);
  Phi1=allVars.at(4);
}

void RooPentaSpinTwo::setVars(double newh1, double newh2, double newphi, double newhs, double newphi1 ){

  h1=  newh1;
  h2=  newh2;
  Phi= newphi;
  hs=  newhs;
  Phi1=newphi1;
}

void RooPentaSpinTwo::setParams(vector<double> allPars){

  if(int(allPars.size())<29){
    std::cout<<"ERROR from RooPentaSpinTwo::setParams. Size of vector with params is less than # of params: "<<allPars.size()<<" < 29"<<std::endl;
  }


  fppVal=allPars.at(0);
  fmmVal=allPars.at(1);
  fpmVal=allPars.at(2);
  fp0Val=allPars.at(3);
  f0mVal=allPars.at(4);
  phippVal=allPars.at(5);
  phimmVal=allPars.at(6);
  phipmVal=allPars.at(7);
  phip0Val=allPars.at(8);
  phi0mVal=allPars.at(9);
  fz1Val=allPars.at(10);
  fz2Val=allPars.at(11);
  R1Val=allPars.at(12);
  R2Val=allPars.at(13);
  para2=allPars.at(14);
  para4=allPars.at(15);
  para6=allPars.at(16);
  para8=allPars.at(17);
  acca0=allPars.at(18);
  acca1=allPars.at(19);
  acca2=allPars.at(20);
  acca4=allPars.at(21);
  a2=allPars.at(22);
  a4=allPars.at(23);
  cutOff=allPars.at(24);
  g=allPars.at(25);
  b2=allPars.at(26);
  b4=allPars.at(27);
  N=allPars.at(28);
  



}//end setParams

 double RooPentaSpinTwo::evaluate() const 
 { 
   double shs = sqrt(1-hs*hs);
   double sh1 = sqrt(1-h1*h1);
   double sh2 = sqrt(1-h2*h2);
   
   if ((1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal) < 0) return 1e-9;
   
   double term1Coeff = (2.-2.*fz1Val+fz2Val-6.*(2.-4.*fz1Val-fz2Val)*pow(hs,2)+3.*(6.-10.*fz1Val-5.*fz2Val)*pow(hs,4));
   double term1A = 4.*(1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal)*pow(sh1,2)*pow(sh2,2);
   double term1B = (fppVal+fmmVal)*((1.+h1*h1)*(1.+h2*h2)+4.*R1Val*R2Val*h1*h2);
   double term1C = -2.*(fppVal-fmmVal)*(R1Val*h1*(1.+h2*h2)+R2Val*h2*(1.+h1*h1));
   double term1D = 4.*sqrt(fppVal*(1-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*(R1Val-h1)*(R2Val-h2)*sh1*sh2*cos(Phi+phippVal);
   double term1E = 4.*sqrt(fmmVal*(1-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*(R1Val+h1)*(R2Val+h2)*sh1*sh2*cos(Phi-phimmVal);
   double term1F = 2.*sqrt(fppVal*fmmVal)*pow(sh1,2)*pow(sh2,2)*cos(2.*Phi+phippVal-phimmVal);
   double term1 = term1Coeff*(term1A+term1B+term1C+term1D+term1E+term1F);

   double term2Coeff = 8.*(fz1Val+fz2Val+3.*(2.-3.*fz1Val-2.*fz2Val)*pow(hs,2)-(6.-10.*fz1Val-5.*fz2Val)*pow(hs,4));
   double term2A = (fp0Val+f0mVal)*(1.-h1*h1*h2*h2) - (fp0Val-f0mVal)*(R1Val*h1*pow(sh2,2)+R2Val*h2*pow(sh1,2));
   double term2B = 2.*sqrt(fp0Val*f0mVal)*sh1*sh2*(R1Val*R2Val-h1*h2)*cos(Phi+phip0Val-phi0mVal);
   double term2 = term2Coeff*(term2A+term2B);
   
   double term3Coeff = -8.*(fz1Val-fz2Val+(6.-10.*fz1Val-5.*fz2Val)*pow(hs,2))*pow(shs,2)*sh1*sh2*cos(Phi + 2.*Phi1);
   double term3A = (fp0Val+f0mVal)*(R1Val*R2Val+h1*h2)-(fp0Val-f0mVal)*(R1Val*h2+R2Val*h1)+2.*sqrt(fp0Val*f0mVal)*sh1*sh2*cos(Phi+phip0Val-phi0mVal);
   double term3 = term3Coeff*term3A;
   
   double term4Coeff = 6.-2.*fz1Val-5.*fz2Val-6.*(2.-2.*fz1Val-3.*fz2Val)*pow(hs,2)+(6.-10.*fz1Val-5.*fz2Val)*pow(hs,4);
   double term4A = fpmVal*((1.+h1*h1)*(1.+h2*h2)-4.*R1Val*R2Val*h1*h2);
   double term4 = term4Coeff*term4A;

   double term5 = pow(shs,4)*(6.-10.*fz1Val-5.*fz2Val)*fpmVal*pow(sh1,2)*pow(sh2,2)*cos(2.*Phi+4.*Phi1);
   
   double term6Coeff = (-1.)*sqrt(6.)*(2.-2.*fz1Val-3.*fz2Val-(6.-10.*fz1Val-5.*fz2Val)*pow(hs,2))*pow(shs,2);
   double term6A = 2.*sqrt(fpmVal*(1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*sh1*sh2*((R1Val-h1)*(R2Val+h2)*cos(Phi+2.*Phi1-phipmVal)+(R1Val+h1)*(R2Val-h2)*cos(Phi+2.*Phi1+phipmVal));
   double term6B = sqrt(fpmVal*fmmVal)*(sh1*sh1*(1.+2.*R2Val*h2+h2*h2)*cos(2.*Phi1-phipmVal+phimmVal)+sh2*sh2*(1.+2.*R1Val*h1+h1*h1)*cos(2.*Phi+2.*Phi1+phipmVal-phimmVal));
   double term6C = sqrt(fppVal*fpmVal)*(sh1*sh1*(1.-2.*R2Val*h2+h2*h2)*cos(2.*Phi1+phipmVal-phippVal)+sh2*sh2*(1.-2.*R1Val*h1+h1*h1)*cos(2.*Phi+2.*Phi1-phipmVal+phippVal));
   double term6 = term6Coeff*(term6A+term6B+term6C);
   
   /// new mixing terms
   double term7Coeff = -4.*sqrt(3.)*(2.-4.*fz1Val-fz2Val-(6.-10.*fz1Val-5.*fz2Val)*pow(hs,2))*hs*shs;
   double term7A = sqrt(fmmVal*f0mVal)*(sh1*(R1Val+h1)*(1.+2.*R2Val*h2+h2*h2)*cos(Phi1+phimmVal-phi0mVal)-sh2*(R2Val+h2)*(1.+2.*R1Val*h1+h1*h1)*cos(Phi+Phi1-phimmVal+phi0mVal));
   double term7B = sqrt(fmmVal*fp0Val)*(sh1*sh1*sh2*(R2Val+h2)*cos(Phi-Phi1-phimmVal+phip0Val)-sh2*sh2*sh1*(R1Val+h1)*cos(2.*Phi+Phi1-phimmVal+phip0Val));
   double term7C = -sqrt(fppVal*f0mVal)*(sh1*sh1*sh2*(R2Val-h2)*cos(Phi-Phi1+phippVal-phi0mVal)-sh2*sh2*sh1*(R1Val-h1)*cos(2.*Phi+Phi1+phippVal-phi0mVal));
   double term7D = -sqrt(fppVal*fp0Val)*(sh1*(R1Val-h1)*(1.-2.*R2Val*h2+h2*h2)*cos(Phi1+phip0Val-phippVal)-sh2*(R2Val-h2)*(1.-2.*R1Val*h1+h1*h1)*cos(Phi+Phi1+phippVal-phip0Val));
   double term7E = -2.*sqrt(f0mVal*(1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*(sh1*sh2*sh2*(R1Val+h1)*cos(Phi1+phi0mVal)-sh2*sh1*sh1*(R2Val+h2)*cos(Phi+Phi1-phi0mVal));
   double term7F = 2.*sqrt(fp0Val*(1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*(sh1*sh2*sh2*(R1Val-h1)*cos(Phi1-phip0Val)-sh2*sh1*sh1*(R2Val-h2)*cos(Phi+Phi1+phip0Val));
   double term7 = term7Coeff*(term7A+term7B+term7C+term7D+term7E+term7F);

   double term8Coeff = 2.*sqrt(2.)*hs*shs*(6.-6.*fz1Val-9.*fz2Val-(6.-10.*fz1Val-5.*fz2Val)*pow(hs,2));
   double term8A = sqrt(fpmVal*f0mVal)*(sh1*(R1Val-h1)*(1.+2.*R2Val*h2+h2*h2)*cos(Phi1-phipmVal+phi0mVal)-sh2*(R2Val-h2)*(1.+2.*R1Val*h1+h1*h1)*cos(Phi+Phi1+phipmVal-phi0mVal));
   double term8B = sqrt(fpmVal*fp0Val)*(sh2*(R2Val+h2)*(1.-2.*R1Val*h1+h1*h1)*cos(Phi+Phi1-phipmVal+phip0Val)-sh1*(R1Val+h1)*(1.-2.*R2Val*h2+h2*h2)*cos(Phi1+phipmVal-phip0Val));
   double term8 = term8Coeff*(term8A+term8B);
   
   double term9Coeff = -2.*sqrt(2.)*hs*pow(shs,3)*(6.-10.*fz1Val-5.*fz2Val);
   double term9A = sqrt(fpmVal*f0mVal)*(sh1*sh1*sh2*(R2Val+h2)*cos(Phi+3.*Phi1-phipmVal+phi0mVal)-sh2*sh2*sh1*(R1Val+h1)*cos(2.*Phi+3.*Phi1+phipmVal-phi0mVal));
   double term9B = sqrt(fpmVal*fp0Val)*(sh1*sh2*sh2*(R1Val-h1)*cos(2.*Phi+3.*Phi1-phipmVal+phip0Val)-sh2*sh1*sh1*(R2Val-h2)*cos(Phi+3.*Phi1+phipmVal-phip0Val));
   double term9 = term9Coeff*(term9A+term9B);

   // signs of the interference terms are flipped!!!!
   if (true){
     term7 *= (-1.);
     term8 *= (-1.);
     term9 *= (-1.);
   }
   
   double sum = term1+term2+term3+term4+term5+term6+term7+term8+term9;

   double accp = (1.0+para2*pow(hs,2)+para4*pow(hs,4)+para6*pow(hs,6)+para8*pow(hs,8))
                  *(1+b2*pow(h1,2)+b4*pow(h1,4))
                  *(acca0+acca1*cos(Phi1)+acca2*cos(2.0*Phi1)+acca4*cos(4.0*Phi1))
                  *(1.0+a2*pow(h2,2)+a4*pow(h2,4))/(1.0+exp((h2-cutOff)/g));
   if(sum*accp<=0){sum=accp=.0001;}
   return sum*accp ;
 } 



/******************
int RooPentaSpinTwo::getAnalyticalIntegral(RooArgSet& allVars, 
                                           RooArgSet& analVars, const char* //rangeName
                                          ) const
{
  if (matchArgs(allVars,analVars,h1,Phi1,Phi)) return 1 ;
  return 0 ;
}

double RooPentaSpinTwo::analyticalIntegral(int code, const char* rangeName) const
{
  switch(code)
    {
    case 1:
      {
	double phi00Val=0;
	double f00;
	f00=1-fmmVal-fppVal-2*f0mVal-2*fp0Val-2*fpmVal;
	double term1 =  (-2*acca0*3.1415*(4*(35 + 7*b2 + 3*b4)*f00*3.1415*(-1 + pow(h2,2)) - 2*(35 + 14*b2 + 9*b4)*3.1415*(fmmVal + fppVal + 2*(fmmVal - fppVal)*R2Val*h2 + (fmmVal + fppVal)*pow(h2,2)))*(2 - 2*fz1Val + fz2Val + 6*(-2 + 4*fz1Val + fz2Val)*pow(hs,2) - 3*(-6 + 10*fz1Val + 5*fz2Val)*pow(hs,4)))/105;
	  double term2 = (8*acca0*3.1415*3.1415*(7*(5*b2 + 3*(5 + b4))*(f0mVal + fp0Val) + 2*(35 + 7*b2 + 3*b4)*(f0mVal - fp0Val)*R2Val*h2 - (35 + 21*b2 + 15*b4)*(f0mVal + fp0Val)*pow(h2,2))*(fz1Val + fz2Val - 3*(-2 + 3*fz1Val + 2*fz2Val)*pow(hs,2) + (-6 + 10*fz1Val + 5*fz2Val)*pow(hs,4)))/105;
	  double term3 = (-8*acca2*(35 + 7*b2 + 3*b4)*sqrt(f0mVal)*sqrt(fp0Val)*3.1415*3.1415*(1 - pow(h2,2))*(-1 + pow(hs,2))*(-fz1Val + fz2Val + (-6 + 10*fz1Val + 5*fz2Val)*pow(hs,2))*cos(phi0mVal - phip0Val))/105;
	  double term4 = (4*acca0*(35 + 14*b2 + 9*b4)*fpmVal*3.1415*3.1415*(1 + pow(h2,2))*(6 - 2*fz1Val - 5*fz2Val + 6*(-2 + 2*fz1Val + 3*fz2Val)*pow(hs,2) + (6 - 10*fz1Val - 5*fz2Val)*pow(hs,4)))/105;
	  double term5 = (sqrt(2/3)*acca2*(35 + 7*b2 + 3*b4)*sqrt(fpmVal)*3.1415*3.1415*(-1 + pow(hs,2))*(2 - 6*pow(hs,2) + fz2Val*(-3 + 5*pow(hs,2)) + 2*fz1Val*(-1 + 5*pow(hs,2)))*(sqrt(fmmVal)*(1 + 2*R2Val*h2 + pow(h2,2))*cos(phimmVal - phipmVal) + sqrt(fppVal)*(1 - 2*R2Val*h2 + pow(h2,2))*cos(phipmVal - phippVal)))/35;
	  double term6 = -(acca1*(8 + 2*b2 + b4)*3.1415*3.1415*3.1415*R1Val*hs*sqrt(3 - 3*pow(hs,2))*(-2 + 4*fz1Val + fz2Val + (6 - 10*fz1Val - 5*fz2Val)*pow(hs,2))*(2*sqrt(f00*f0mVal)*(-1 + pow(h2,2))*cos(phi00Val - phi0mVal) + sqrt(f0mVal)*sqrt(fmmVal)*(1 + 2*R2Val*h2 + pow(h2,2))*cos(phi0mVal - phimmVal) + 2*sqrt(f00*fp0Val)*cos(phi00Val - phip0Val) - 2*sqrt(f00*fp0Val)*pow(h2,2)*cos(phi00Val - phip0Val) - sqrt(fp0Val)*sqrt(fppVal)*cos(phip0Val - phippVal) + 2*sqrt(fp0Val)*sqrt(fppVal)*R2Val*h2*cos(phip0Val - phippVal) - sqrt(fp0Val)*sqrt(fppVal)*pow(h2,2)*cos(phip0Val - phippVal)))/16;
	  double term7 = -(acca1*(8 + 2*b2 + b4)*sqrt(fpmVal)*3.1415*3.1415*3.1415*R1Val*hs*sqrt(1 - pow(hs,2))*(6 - 6*fz1Val - 9*fz2Val + (-6 + 10*fz1Val + 5*fz2Val)*pow(hs,2))*(sqrt(f0mVal)*(1 + 2*R2Val*h2 + pow(h2,2))*cos(phi0mVal - phipmVal) + sqrt(fp0Val)*(-1 + 2*R2Val*h2 - pow(h2,2))*cos(phip0Val - phipmVal)))/(16*sqrt(2));

	  return N*(term1+term2+term3+term4+term5+term6+term7)*(1+para2*hs*hs+para4*hs*hs*hs*hs+para6*pow(hs,6)+para8*pow(hs,8))*(1+a2*h2*h2+a4*h2*h2*h2*h2)/(1 + exp((-cutOff + h2)/g));
      }
    }

assert(0) ;
return 0 ;
}




*******************/
