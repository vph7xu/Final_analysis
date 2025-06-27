#ifndef BRANCH_VARS_H
#define BRANCH_VARS_H

#include <TTree.h>
#include <string>
#include <unordered_map>

/**  Thin POD-style holder for all leaves you ever read.
 *   Call attach(tree) once per file; after that just use .W2, .dx, …
 */
class BranchVars {
public:
    //  ==== public leaf members (extend as you wish) ====
    double W2        = 0;   ///< reconstructed W²
    double dx        = 0;   ///< SBS-HCAL horizontal residual
    double dy        = 0;
    double eHCAL     = 0;
    double coin_time = 0;
    int    runnum    = 0;
    // … add more here, no other file changes needed …

    /** Bind all data members to the branches of a TTree. */
    void attach(TTree* t);

private:
    /** Helper so we can iterate over the members generically. */
    std::unordered_map<std::string, void*> addrMap_();
};

#endif  // BRANCH_VARS_H

