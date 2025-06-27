#include "BranchVars.h"
#include <TError.h>


std::unordered_map<std::string, void*> BranchVars::addrMap_() {
    return {
        {"W2",        &W2},
        {"dx",        &dx},
        {"dy",        &dy},
        {"eHCAL",     &eHCAL},
        {"coin_time", &coin_time},
        {"runnum",    &runnum},
        // --- extend when you add a new data member! ---
    };
}

void BranchVars::attach(TTree* t)
{
    for (auto& kv : addrMap_()) {
        if (t->GetBranch(kv.first.c_str()))
            t->SetBranchAddress(kv.first.c_str(), kv.second);
        else
            Warning("BranchVars::attach", "Missing branch %s — value will stay zero.", kv.first.c_str());
    }
}

