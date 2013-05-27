#include <vector>
#include <string>
using namespace std;
class RooPentaSpinTwo{
public:
  RooPentaSpinTwo();
  RooPentaSpinTwo(const char *name, const char *title, std::vector<double>allVars, std::vector<double>allPars);

  RooPentaSpinTwo(const char *name, const char *title,
	          double _h1,//costheta1
        	  double _h2,//costheta2
              	  double _Phi,//phi
                  double _hs,//costhetastar
                  double _Phi1,//phistar1
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
		  double _N);
  RooPentaSpinTwo(const RooPentaSpinTwo& other, const char* name=0) ;
  // virtual TObject* clone(const char* newname) const { return new RooPentaSpinTwo(*this,newname); }
  inline virtual ~RooPentaSpinTwo() { }
  //  int getAnalyticalIntegral(RooArgSet& allPars, RooArgSet& analVars, const char* rangeName=0) const ;
  //  double analyticalIntegral(int code, const char* rangeName=0) const ;
  double evaluate() const ;
  void setParams(vector<double> allPars);
  void setVars(vector<double>allVars);
  void setVars(double newh1, double newh2, double newphi, double newhs, double newphi1);

protected:

  string name_;
  string title_;
  double h1 ;
  double h2 ;
  double Phi ;
  double hs ;
  double Phi1 ;
  double fppVal ;
  double fmmVal ;
  double fpmVal ;
  double fp0Val ;
  double f0mVal ;
  double phippVal ;
  double phimmVal ;
  double phipmVal ;
  double phip0Val ;
  double phi0mVal ;
  double fz1Val ;
  double fz2Val ;
  double R1Val ;
  double R2Val ;
  double para2 ;
  double para4 ;
  double para6 ;
  double para8 ;
  double acca0 ;
  double acca1 ;
  double acca2 ;
  double acca4 ;
  double a2 ;
  double a4 ;
  double cutOff ;
  double g ;
  double b2 ;
  double b4 ;
  double N  ;


private:

};
 
