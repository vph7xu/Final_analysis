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
    double ebeam     = 0.0;
    double trPx      = 0.0;
    double trPy      = 0.0;
    double trPz      = 0.0;
    double etheta    = 0.0;
    double W2        = 0.0;   ///< reconstructed W²
    double dx        = 0.0;   ///< SBS-HCAL horizontal residual
    double dy        = 0.0;
    double eHCAL     = 0.0;
    double coin_time = 0.0;
    int    runnum    = 0;
    int    helicity  = 0;
    int    IHWP      = 0;
    double He3Pol    = 0.0;
    double ePS       = 0.0;
    double vz        = 0.0;
    double eSH       = 0.0;
    double trP       = 0.0;
    double grinch_track     = 0.0;
    double grinch_clus_size = 0.0;
    double Q2        = 0.0;
    double pN_expect = 0.0;
    double ntrack    = 0.0;
    double trP_sbs   = 0.0;
    double ntrack_sbs       = 0.0;
    double vz_sbs           = 0.0;
    double theta_pq         = 0.0;
    TDatime *datetime = nullptr;
    double err_He3Pol= 0.0; 
    double thtgt = 0.0;
    double thetabend = 0.0;
    int    hcal_nclus = 0;
    double hcal_clus_e[100];
    double hcal_clus_x[100];
    double hcal_clus_y[100];
    double hcal_clus_atime[100];
    double hcal_clus_nblk[100];
    
    // … add more here, no other file changes needed …

    /** Bind all data members to the branches of a TTree. */
    void attach(TTree* t);

private:
    /** Helper so we can iterate over the members generically. */
    std::unordered_map<std::string, void*> addrMap_();
};

#endif  // BRANCH_VARS_H

