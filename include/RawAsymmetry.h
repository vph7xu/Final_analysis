// -----------------------------------------------------------------------------
// RawAsymmetry.h  —  overall spectrum + per‑run asymmetry module
// -----------------------------------------------------------------------------
// * Requires BranchVars to expose: runnum (int) and helicity (int ±1)
// * Uses CutManager for numeric cuts and optional RunQuality veto
// * Writes:
//     – ROOT file with overall Aexp histogram
//     – text file: run  N_plus  N_minus  A_raw  dA_stat
// -----------------------------------------------------------------------------
#ifndef RAW_ASYMMETRY_H
#define RAW_ASYMMETRY_H

#include "BranchVars.h"
#include "CutManager.h"
#include "RunQuality.h"

#include <TChain.h>
#include <TH1D.h>
#include <unordered_map>
#include <string>

class RawAsymmetry {
public:
    /**
     * Constructor.
     * @param cuts      Numeric cut manager (must out‑live this object)
     * @param rq        Optional run‑quality pointer (nullptr ⇒ veto disabled)
     * @param nbins     Histogram bins
     * @param xmin,xmax Histogram x‑range
     * @param outRoot   ROOT output file name
     * @param outTxt    Text output file name
     */
    RawAsymmetry(const CutManager& cuts,
                 const RunQuality* rq  = nullptr,
                 const char * kin = "GEN3_He3",
                 int    nbins   = 100,
                 double xmin    = -4,
                 double xmax    =  3,
                 const char* outRoot = "raw_asymmetry.root",
                 const char* outTxt  = "raw_asym_per_run.txt");

    /** Loop over entries (BranchVars must already be attached) */
    void process(TChain& ch, BranchVars& v);

private:
    const CutManager& cuts_;
    const RunQuality* rq_;   // nullable
    const char* kin_;

    TH1D        h_;          // overall Aexp spectrumcutsIdx
    std::string rootName_;
    std::string txtName_;

    // run → {N_plus, N_minus}
    std::unordered_map<int, std::pair<long long, long long>> counts_;
};

#endif // RAW_ASYMMETRY_H

