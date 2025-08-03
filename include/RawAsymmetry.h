// -----------------------------------------------------------------------------
// RawAsymmetry.h  —  overall spectrum + per‑run asymmetry module
// -----------------------------------------------------------------------------
// * Requires BranchVars to expose: runnum (int) and helicity (int ±1)
// * Uses CutManager for numeric cuts and optional RunQuality veto
// * Writes:
//     – ROOT file with overall Aexp histogram
//     – text file: run  N_plus  N_minus  A_raw  dA_stat
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// RawAsymmetry.h  — beam‑ and target‑polarisation aware asymmetry module
// -----------------------------------------------------------------------------
#ifndef RAW_ASYMMETRY_H
#define RAW_ASYMMETRY_H

#include "BranchVars.h"
#include "CutManager.h"
#include "AnalysisCuts.h"
#include "RunQuality.h"

#include <TChain.h>
#include <TH1D.h>
#include <TDatime.h>

#include <string>
#include <vector>
#include <unordered_map>

// Row of Beam_pol.csv
struct BeamRow {
    TDatime start, end;   // validity window
    double  val;          // polarisation
    double  err;          // ±σ
};

struct Counts { long long Np=0, Nm=0; };

struct PolAcc {                       // accumulate per run
    double sumBeam=0,  sumBeamErr2=0;
    double sumHe3 =0,  sumHe3Err2 =0;
    long   w      =0;                 // weight (events)
};

struct PolAccSys {
    // Event-weighted means (use w=1 or charge/live-time if desired)
    double sumWb  = 0.0, sumWPb = 0.0;  // beam: Σw, Σ(w*P)
    double sumWt  = 0.0, sumWPt = 0.0;  // target: Σw, Σ(w*P)

    // Blockwise totals for systematic combination
    // Beam block id is the CSV row index (interval)
    std::unordered_map<int, double>       beam_blk_sumWP;  // idx -> Σ(w*P)
    std::unordered_map<int, double>       beam_blk_err;    // idx -> abs error

    // Target block key is made from (P, err) rounded to 1e-4 to group identical calibrations
    std::unordered_map<long long, double> targ_blk_sumWP;  // key -> Σ(w*P)
    std::unordered_map<long long, double> targ_blk_err;    // key -> abs error
};

class RawAsymmetry {
public:
    RawAsymmetry(const CutManager&  cuts,
                 const AnalysisCuts& acuts,
                 const RunQuality*  rq   = nullptr,
                 const char*        kin  = "GEN3_He3",
                 int    nbins = 100,
                 double xmin  = -4,
                 double xmax  =  3,
                 const char* rootFile = "raw_asymmetry.root",
                 const char* txtFile  = "raw_asymmetry.txt");

    void process(TChain& ch, BranchVars& v);

private:
    // helpers
    std::vector<BeamRow> readBeamCSV(const std::string&) const;
    std::pair<double,double> beamAt(const std::vector<BeamRow>&, const TDatime&) const;

    const CutManager&  cuts_;
    const AnalysisCuts& c_;
    const RunQuality*  rq_;
    const char*      kin_;

    TH1D        h_;           // Aexp spectrum
    std::string rootF_, txtF_;

    std::unordered_map<int, Counts> counts_; // run → N±
    //std::unordered_map<int, PolAcc> pols_;   // run → pol sums //old implementation we don't use this anymore
    std::unordered_map<int, PolAccSys> pols_;
};

#endif // RAW_ASYMMETRY_H

