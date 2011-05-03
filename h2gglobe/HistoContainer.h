#ifndef HISTOCONTAINER
#define HISTOCONTAINER

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <RooDataHist.h>
#include <RooRealVar.h>
#include <map>
#include <string>

class HistoContainer {

 public:
  HistoContainer();
  HistoContainer(int, std::string);
  ~HistoContainer();
    
  void Add(char *, int, int, float, float);
  void Add(char *, int, int, float, float, int, float, float);
  void Add(char *, int, int, float, float, float, float);

  void Fill(std::string, int, float);
  void Fill(std::string, int, float, float);
  
  void Fill2D(std::string, int, float, float);
  void Fill2D(std::string, int, float, float, float);
  
  void Save();
 
  int getHistVal();
  void setHistVal(int);
  void setHistNam(std::string);
  void setScale(float);
  std::string ModifiedName(char*, int);
  float total_scale;

 private:
  void writeRooDataHist(TH1F*);
  int histVal;
  std::string histNam;
  std::map<std::string, std::vector<TH1F> > h1;
  std::map<std::string, std::vector<TH2F> > h2;
  std::map<std::string, std::vector<TProfile> > hp;
};

#endif
