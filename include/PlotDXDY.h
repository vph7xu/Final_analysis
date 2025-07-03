// -----------------------------------------------------------------------------
// PlotDXDY.h — 2D (dx,dy) residual plot module
//   * Consumes AnalysisCuts (single-source-of-truth) instead of its own JSON
//   * Cuts applied: dx_L/H, dy_L/H, helicity (±1 or inclusive 999)
//   * Produces TH2D saved to a ROOT file
// -----------------------------------------------------------------------------
#ifndef PLOT_DXDY_H
#define PLOT_DXDY_H

#include "BranchVars.h"
#include "AnalysisCuts.h"
#include "RunQuality.h"

#include <TChain.h>
#include <TH2D.h>
#include <string>

class PlotDXDY {
public:
    PlotDXDY(const AnalysisCuts& cuts,
             const RunQuality* rq  = nullptr,
             const char* kin = "GEN3_He3",
             const char* outRoot = "plot_dxdy.root"
             );

    void process(TChain& ch, BranchVars& v);

private:
    const AnalysisCuts& c_;
    const RunQuality* rq_;   // nullable
    TH1D        hvz_;
    TH1D        hePS_;
    TH1D        heHCAL_;
    TH1D        hW2_;
    TH1D        hcointime_;
    TH2D        hDXDY_;
    TH2D        hDXW2_;
    TH2D        hDYW2_;
    TH1D        hDX_;
    TH1D        hDY_;
    std::string outFile_;
    const char* kin_;

    bool passes(const BranchVars& v) const;
};

#endif // PLOT_DXDY_H

