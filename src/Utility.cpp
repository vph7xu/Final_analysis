#include "Utility.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <functional>
#include <utility>
#include <algorithm>

using std::string; using std::vector; using std::unordered_map;


unordered_map<string,double>Utility::readSimpleKVGlobal(const string& f) const
{
    unordered_map<string,double> m;
    std::ifstream in(f);
    if(!in){ std::cerr << "[PhysCalc] cannot open " << f << " "; return m; }

    string line;
    while(std::getline(in,line)){
        if(line.empty()) continue;
        std::stringstream ss(line);
        string key; ss >> key;
        if(key[0]=='#') continue;             // allow comment lines

        // optional '=' token — eat it if present
        if(ss.peek()=='='){ ss.get(); }
        else{
            // maybe spaces before '='
            char c; ss >> c; if(c!='=') ss.unget(); else { /* '=' already consumed */ }
        }
        double val; ss >> val;
        m[key] = val;
    }
    return m;
}
