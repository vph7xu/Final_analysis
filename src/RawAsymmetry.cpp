#include "RawAsymmetry.h"

#include <algorithm>  // std::max
#include <iostream>

RawAsymmetry::RawAsymmetry(const CutManager& cut,
                           int nbins, double xmin, double xmax,
                           const char* outFile)
    : cut_     (cut)
    , hA_      ("hA", "raw asym.; A_{exp}; counts", nbins, xmin, xmax)
    , outFile_ (outFile)
{}


void RawAsymmetry::process(TChain& ch, const BranchVars& v)
{
    const Long64_t nEntries = ch.GetEntries();
    const Long64_t step     = 100;//std::max<Long64_t>(1, nEntries / 50);   // 2 %

    std::cout << "\n[ RawAsym ] processing " << nEntries << " events …\n";

    for (Long64_t i = 0; i < nEntries; ++i) {
        ch.GetEntry(i);
        if (!cut_.passAll(v)) continue;

        const double Aexp = (v.dy != 0) ? v.dx / v.dy : 0.0;
        hA_.Fill(Aexp);

        if (i % step == 0 || i == nEntries - 1) {
            double frac = double(i + 1) / nEntries;
            int barw = 42, pos = barw * frac;
            std::cout << '\r' << '[';
            for (int j = 0; j < barw; ++j)
                std::cout << (j < pos ? '=' : (j == pos ? '>' : ' '));
            std::cout << "] " << int(frac * 100) << " %" << std::flush;
        }
    }
    std::cout << "\nDone.\n";

    TFile fout(outFile_.c_str(), "RECREATE");
    hA_.Write();
    fout.Close();
}

