#include "BranchVars.h"
#include <TError.h>


std::unordered_map<std::string, void*> BranchVars::addrMap_() {
    return {
        {"ebeam",     &ebeam},
        {"trPx",      &trPx},
        {"trPy",      &trPy},
        {"trPz",      &trPz},
        {"etheta",    &etheta},
        {"W2",        &W2},
        {"dx",        &dx},
        {"dy",        &dy},
        {"eHCAL",     &eHCAL},
        {"coin_time", &coin_time},
        {"runnum",    &runnum},
        {"helicity",  &helicity},
        {"IHWP",      &IHWP},
        {"He3Pol",    &He3Pol},
        {"err_He3Pol",&err_He3Pol},
        {"ePS",       &ePS},
        {"vz",        &vz},
        {"eSH",       &eSH},
        {"trP",       &trP},
        {"grinch_track",        &grinch_track},
        {"grinch_clus_size",    &grinch_clus_size},
        {"datetime",            &datetime},
        {"Q2",        &Q2},
        {"pN_expect", &pN_expect},
        {"ntrack",    &ntrack},
        {"trP_sbs",   &trP_sbs},
        {"ntrack_sbs",&ntrack_sbs},
        {"vz_sbs",    &vz_sbs},
        {"theta_pq",  &theta_pq},
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

