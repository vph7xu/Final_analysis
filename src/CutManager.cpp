#include "CutManager.h"

namespace {
    /** Pointer-to-member lookup so we can query BranchVars generically */
    const std::unordered_map<std::string, double BranchVars::*> ptr = {
        {"W2",        &BranchVars::W2},
        {"dx",        &BranchVars::dx},
        {"dy",        &BranchVars::dy},
        {"eHCAL",     &BranchVars::eHCAL},
        {"coin_time", &BranchVars::coin_time}
        //{"runnum",    &BranchVars::runnum}
    };
}

bool CutManager::pass(const BranchVars& v,
                      const char* varName,
                      const char* lowKey,
                      const char* highKey) const
{
    const double val = v.*ptr.at(varName);
    if (cfg_.contains(lowKey)  && val < cfg_[lowKey])  return false;
    if (cfg_.contains(highKey) && val > cfg_[highKey]) return false;
    return true;
}

