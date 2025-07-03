// -----------------------------------------------------------------------------
// AccidentalCorrection.h —  asymmetry vs coincidence-time module
// derive uniform bins from AnalysisCuts    (cuts, binWidth)
// using cuts.coin_ac_L / coin_ac_H         (defaults: 5 ns bin width)
// -----------------------------------------------------------------------------
#ifndef ACCIDENTAL_CORRECTION_H
#define ACCIDENTAL_CORRECTION_H

#include "BranchVars.h"
#include "AnalysisCuts.h"
#include "RunQuality.h"
#include <TChain.h>
#include <string>
#include <vector>
#include <functional>

class AccidentalCorrection {
public:

    // derive [lo,hi] from AnalysisCuts.coin_ac_L/H with fixed width
    AccidentalCorrection(const AnalysisCuts& cuts,
                         const RunQuality* rq  = nullptr,    			
                         const char* kin = "GEN3_He3",
                         double binWidth = 5.0,                  // ns
                         const char* rootFile = "asym_vs_cointime.root",
                         const char* csvFile  = "asym_vs_cointime.csv");

    void process(TChain& ch, BranchVars& v);

private:
    const AnalysisCuts& c_;
    const RunQuality* rq_;   // nullable
    std::vector<double> edges_;
    const char* kin_;

    std::string rootName_;
    std::string csvName_;

    std::vector<long long> cntP_, cntM_;

    void initBins();
};

#endif // ACCIDENTAL_CORRECTION_H

