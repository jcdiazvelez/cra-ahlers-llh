#ifndef PICKLE_INCLUDE
#define PICKLE_INCLUDE

namespace pickle { 
        void 
        dump(const std::vector<double>& map, std::string filename, std::string delimit="\n");

        void 
        load(std::vector<double>& map, std::string filename);
}


#endif // PICKLE_INCLUDE
