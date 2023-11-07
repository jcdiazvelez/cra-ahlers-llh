#include <boost/format.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/foreach.hpp>
#include <sstream>
#include <fstream>  
#include <iomanip>  
#include <string>
#include <stdexcept>
#include <iter-lhreco-proj/pickle.h>


using boost::format;

void 
pickle::dump(const std::vector<double>& map, std::string filename, std::string delimit) 
{ 
        std::ofstream ofs; 
        ofs.open (filename.c_str(), std::ofstream::out | std::ofstream::app); 
        BOOST_FOREACH(const double& entry, map){ 
                ofs << std::setprecision(5) << entry; 
                ofs << delimit; 
        } 
        ofs.close(); 
}

void 
pickle::load(std::vector<double>& map, std::string filename)
{ 
        std::ifstream ifs; 
        ifs.open (filename.c_str(), std::ifstream::in); 
        std::string delim;
        double entry;
        while (ifs.good()) { 
                std::cin >> std::dec >> entry; 
                std::cin >> delim;
                map.push_back(entry); 
        } 
        ifs.close(); 
}
