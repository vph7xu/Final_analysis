//...................................................................
//Utility.h - add some features to for global use
//


#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <vector>
#include <unordered_map>


class Utility {
public:
	    std::unordered_map<std::string,double>
        readSimpleKVGlobal(const std::string& file) const;
};


#endif //UTILITY_H
