// -----------------------------------------------------------------------------
// PlotDXDY.h — 2‑D residual (dx,dy) plotting module with custom cut set
// -----------------------------------------------------------------------------
// * Reads its own CutConfig (can be DIFFERENT from the one in CutManager)
// * Example JSON keys used:
//       "dx_L" "dx_H"  (double)
//       "dy_L" "dy_H"  (double)
//       "helicity"     (int  ±1  or 999 for inclusive)
// * Fills TH2D and writes to ROOT.
// -----------------------------------------------------------------------------
#ifndef PLOT_DXDY_H
#define PLOT_DXDY_H

#include "BranchVars.h"
#include "CutConfig.h"
#include <TChain.h>
#include <TH2D.h>
#include <string>

class PlotDXDY {
public:
    /**
     * @param cfgFile   JSON/YAML/TXT defining dx/dy cuts specifically for this plot
     * @param outRoot   File to write the TH2D into (defaults plot_dxdy.root)
     */
    explicit PlotDXDY(const std::string& cfgFile,
                      const char* outRoot = "plot_dxdy.root");

    void process(TChain& ch, BranchVars& v);

private:
    CutConfig  cfg_;
    TH2D       h2_;
    std::string outFile_;

    bool passesDXDY(const BranchVars& v) const;
};

#endif // PLOT_DXDY_H
