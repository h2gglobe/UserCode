#ifndef SAMPLECONTAINER
#define SAMPLECONTAINER

#include <string>

class SampleContainer {

 public:
  SampleContainer();
  ~SampleContainer();
  
  void computeWeight(float);
     
  float weight;
  float event_weight;
  int itype;
  int ind;
  int histoplotit;
  std::string filesshortnam;
  long long int ntot;
  int nred;
  float lumi; 
  float xsec;
  float kfactor; 
  float scale;
  float lumireal;
  
 private:

};

#endif
