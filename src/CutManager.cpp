// -----------------------------------------------------------------------------
// CutManager.cpp — explicit low/high keys, supports double & int branches
// -----------------------------------------------------------------------------
#include "CutManager.h"

#include <unordered_map>
#include <vector>
#include <string>
#include <exception>

// -----------------------------------------------------------------------------
// Pointer maps (no auto _L/_H generation; keys listed explicitly below)
// -----------------------------------------------------------------------------
namespace {
    const std::unordered_map<std::string, double BranchVars::*> ptrDbl = {
        {"W2",        &BranchVars::W2},
        {"dx",        &BranchVars::dx},
        {"dy",        &BranchVars::dy},
        {"eHCAL",     &BranchVars::eHCAL},
        {"coin_time", &BranchVars::coin_time}
    };

    const std::unordered_map<std::string, int BranchVars::*> ptrInt = {
        {"runnum",   &BranchVars::runnum},
        {"helicity", &BranchVars::helicity}
    };

    struct CutEntry {
        const char* var;       // name as in maps above
        const char* lowKey;    // key in JSON/TXT for low  threshold
        const char* highKey;   // key in JSON/TXT for high threshold
    };

    // Here you define which variables carry which threshold keys
    const std::vector<CutEntry> cutTable = {
        {"W2",        "W2_L",        "W2_H"},
        {"dx",        "dx_L",        "dx_H"},
        {"dy",        "dy_L",        "dy_H"},
        {"eHCAL",     "eHCAL_L",     "eHCAL_H"},
        {"coin_time", "coin_time_L", "coin_time_H"},
        {"runnum",    "runnum_L",    "runnum_H"}
        // add helicity cut if desired
    };
}

// safeGet — fetch numeric value if exists and convertible to double
static bool safeGet(const CutConfig& cfg, const char* key, double& out)
{
    if (!cfg.contains(key)) return false;
    try { out = cfg[key]; }
    catch (const std::exception&) { return false; }
    return true;
}

// helper function
bool CutManager::inRange(double val, const char* kL, const char* kH) const
{
    double lo, hi;
    if (safeGet(cfg_, kL, lo) && val < lo) return false;
    if (safeGet(cfg_, kH, hi) && val > hi) return false;
    return true;
}

bool CutManager::passAll(const BranchVars& v) const
{
    // ---- hard veto: keep only helicity ±1 ----
    if (v.helicity != 1 && v.helicity != -1) return false;
    
    // ---- hard veto: global cuts ----
    if (v.ntrack<1) return false; //track 
    if (abs(v.vz)>0.27) return false; //vertex cut
    if (v.eHCAL<0.025) return false; //hcal energy cut
    if (v.ePS<0.2) return false; //preshower energy cut

    if (abs(((v.ePS+v.eSH)/v.trP)-1.0)>0.2) return false; //eoverp cut

    for (const auto& ce : cutTable) {
        // pick pointer map depending on variable type
        double val = 0.0;
        if (auto it = ptrDbl.find(ce.var); it != ptrDbl.end()) {
            val = v.*(it->second);
        } else if (auto it2 = ptrInt.find(ce.var); it2 != ptrInt.end()) {
            val = static_cast<double>(v.*(it2->second));
        } else {
            continue; // unknown variable name → skip
        }
        if (!inRange(val, ce.lowKey, ce.highKey)) return false;
    }
    return true;
}

