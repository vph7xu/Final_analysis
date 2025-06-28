#include "PlotDXDY.h"

#include <TFile.h>
#include <iostream>
#include <cmath>

PlotDXDY::PlotDXDY(const std::string& cfgFile, const char* outRoot)
    : h2_("hDXDY", "dy vs dx; dx; dy", 120, -0.06, 0.06, 120, -0.06, 0.06),
      outFile_(outRoot)
{
    if (!cfg_.load(cfgFile))
        std::cerr << "[PlotDXDY] Warning: could not read cuts from " << cfgFile << ". Using full range.\n";
}

bool PlotDXDY::passesDXDY(const BranchVars& v) const
{
    auto in = [&](double val, const char* loK, const char* hiK){
        double lo, hi;
        bool hasLo = cfg_.contains(loK);
        bool hasHi = cfg_.contains(hiK);
        if (hasLo && val < cfg_[loK]) return false;
        if (hasHi && val > cfg_[hiK]) return false;
        return true;
    };

    if (!in(v.dx, "dx_L", "dx_H")) return false;
    if (!in(v.dy, "dy_L", "dy_H")) return false;

    if (cfg_.contains("helicity")) {
        int wanted = static_cast<int>(cfg_["helicity"]);
        if (wanted != 999 && v.helicity != wanted) return false;
    }
    return true;
}

void PlotDXDY::process(TChain& ch, BranchVars& v)
{
    Long64_t n = ch.GetEntries();
    std::cout << "[PlotDXDY] looping over " << n << " events\n";
    for (Long64_t i=0;i<n;++i){
        ch.GetEntry(i);
        if (!passesDXDY(v)) continue;
        h2_.Fill(v.dx, v.dy);
    }
    TFile f(outFile_.c_str(), "RECREATE");
    h2_.Write();
    f.Close();
    std::cout << "[PlotDXDY] histogram written to " << outFile_ << "\n";
}

