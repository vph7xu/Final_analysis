#include "PhysicsAsymmetryMulti.h"

#include <TFile.h>
#include <TNamed.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <functional>  // std::function explicit instantiation
#include <utility>     // std::pair

using std::cerr;

// -----------------------------------------------------------------------------
void PhysicsAsymmetryMulti::addSource(const std::string& f,const std::string& k)
{
    sources_.push_back({f,k});
}

// -----------------------------------------------------------------------------
CorrEntry PhysicsAsymmetryMulti::loadPair(const std::string& file,const std::string& key)
{
    CorrEntry c;
    TFile f(file.c_str());
    if(f.IsOpen()){
        if(auto n = dynamic_cast<TNamed*>(f.Get(key.c_str()))){
            std::stringstream ss(n->GetTitle()); char comma;
            ss>>c.v>>comma>>c.e; return c; }
    }
    // fall‑back to TXT: expect two numbers on first line
    std::ifstream in(file); if(!in) return c;
    std::string line; std::getline(in,line); std::stringstream ss(line);
    ss>>c.v>>c.e; return c;
}

// -----------------------------------------------------------------------------
template<typename Expr>
bool PhysicsAsymmetryMulti::compute(const Expr& expr)
{
    // load every requested correction once
    for(auto& pr : sources_){ corr_[pr.second] = loadPair(pr.first, pr.second); }

    // let user‑supplied expression build A and σ from the map
    auto res = expr(corr_);
    A_   = res.first;
    dA_  = res.second;
    bool ok = !(std::isnan(A_)||std::isnan(dA_));
    if(!ok) cerr<<"[PhysAsymMulti] expression returned NaN\n";
    return ok;
}

// explicit instantiation of template method for linker
template bool PhysicsAsymmetryMulti::compute<std::function<std::pair<double,double>(const std::unordered_map<std::string, CorrEntry>&)>>(const std::function<std::pair<double,double>(const std::unordered_map<std::string, CorrEntry>&)> &);

