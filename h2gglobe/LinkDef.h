#ifndef ROOT_TREE_VECTOR_LINKDEF_H 
#define ROOT_TREE_VECTOR_LINKDEF_H 1

#ifdef __CINT__

//#pragma link off all classes;

#pragma link C++ class std::vector<std::vector<unsigned int> >+;
#pragma link C++ class std::vector<std::vector<unsigned short> >+;
#pragma link C++ class std::vector<std::vector<int> >+;
#pragma link C++ class std::vector<std::vector<short> >+;
#pragma link C++ class std::vector<std::vector<float> >+;
//#pragma link C++ class std::vector<unsigned short>+;
#pragma link C++ class std::vector<std::string>+;
#pragma link C++ class std::vector<std::vector<std::vector<unsigned int> > >+;
#pragma link C++ class std::pair<std::string, std::string>+;

// #pragma link C++ defined_in "CategoryOptimizer/interface/CategoryOptimizer.h";
#include "TmpLinkDef.h"

#endif

#endif // ROOT_TREE_VECTOR_LINKDEF_H
