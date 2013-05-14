#include <sys/stat.h>
#include <iostream> 
#include <string> 
#include <vector> 
#include <algorithm> 
#include <stdexcept> 

#include <cstdio> 
#include <cstdlib> 

#include <boost/regex.hpp>

int fexist(char*);

class ExecCommand { 
public: 
  ExecCommand( const std::string &cmd ) { 
    FILE *pfd = popen( cmd.c_str(), "r" ); 
    if ( pfd == 0 ) { 
      throw std::runtime_error( "Command or process could not be executed." ); 
    } 
    while (!feof(pfd)) { 
      char buf[ 1024 ] = {0}; 
      if (fgets(buf, sizeof(buf), pfd) > 0) { 
	m_results.push_back( std::string(buf) ); 
      } 
    } 
    pclose( pfd ); 
  } 

  const std::vector<std::string>& getResults() const { return m_results; } 
  int getTag() {
    int result = 0;
    std::string tag = "None";
    if (m_results.size() > 0) {
      const char* pattern = "[A-Za-z]V(\\d+)_(\\d+)_(\\d+)";
      boost::regex re(pattern);
      boost::cmatch matches;
      if (boost::regex_search(m_results[0].c_str(), matches, re)) {
	result = atoi((matches[1]+matches[2]+matches[3]).c_str());
      }
    }

    return result;
  }
  
private:
  std::vector<std::string> m_results; 
}; 
