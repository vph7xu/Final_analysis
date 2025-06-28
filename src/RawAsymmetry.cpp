#include "RawAsymmetry.h"

#include <TFile.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>

RawAsymmetry::RawAsymmetry(const CutManager& cuts,
                           const RunQuality* rq,
                           int nbins, double xmin, double xmax,
                           const char* outRoot,
                           const char* outTxt)
    : cuts_   (cuts),
      rq_     (rq),
      h_      ("hA", "raw asym.; A_{exp}; counts", nbins, xmin, xmax),
      rootName_(outRoot),
      txtName_ (outTxt) {}

void RawAsymmetry::process(TChain& ch, BranchVars& v)
{
    const Long64_t nEntries = ch.GetEntries();
    const Long64_t step     = 100;//std::max<Long64_t>(1, nEntries / 50);

    std::cout << "\n[ RawAsymmetry ] processing " << nEntries << " events …\n";

    for (Long64_t i = 0; i < nEntries; ++i) {
        ch.GetEntry(i);

        // numeric cuts
        if (!cuts_.passAll(v)) continue;
        // run-quality veto (if provided)
        if (rq_ && (!rq_->helicityOK(v.runnum) || !rq_->mollerOK(v.runnum))) continue;

        // Accumulate N+ / N- per run (assumes BranchVars.helicity exists)
        auto &pair = counts_[v.runnum];
        if (v.helicity > 0) ++pair.first;      // N+
        else if (v.helicity < 0) ++pair.second; // N-

        // overall spectrum (optional)
        //const double Aexp = (v.dy != 0) ? v.dx / v.dy : 0.0;
        h_.Fill(v.dx);

        // progress bar
        if (i % step == 0 || i == nEntries - 1) {
            double frac = double(i + 1) / nEntries;
            int barw = 42, pos = static_cast<int>(barw * frac);
            std::cout << '\r' << '[';
            for (int j = 0; j < barw; ++j)
                std::cout << (j < pos ? '=' : (j == pos ? '>' : ' '));
            std::cout << "] " << static_cast<int>(frac * 100) << " %" << std::flush;
        }
    }
    std::cout << "\nDone.\n";

    // ---- save histogram ----
    {
        TFile fout(rootName_.c_str(), "RECREATE");
        h_.Write();
        fout.Close();
    }

    // ---- compute & write per-run asymmetry & statistical error ----
    std::ofstream txt(txtName_);
    txt << "#run  N_plus  N_minus  A_raw  dA_stat\n";
    for (const auto& kv : counts_) {
        int run   = kv.first;
        long long Np = kv.second.first;
        long long Nm = kv.second.second;
        long long sum = Np + Nm;
        double A = (sum > 0) ? double(Np - Nm) / sum : 0.0;
        double dA = (sum > 0) ? 2.0 * std::sqrt(double(Np) * Nm) / (sum * sum) : 0.0; // 2*sqrt(N+ N-)/(N+ + N-)^2
        txt << run << " " << Np << " " << Nm << " " << A << " " << dA << "\n";
    }
    txt.close();

    std::cout << "Per-run asymmetries written to " << txtName_ << "\n";
}
