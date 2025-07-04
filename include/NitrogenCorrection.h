// -----------------------------------------------------------------------------
// PionCorrection.h — pion‑region raw asymmetry & correction factor
// -----------------------------------------------------------------------------
// * Uses BranchVars fields: eHCAL (cluster energy), helicity, runnum
// * Region is defined by AnalysisCuts.pion_L / pion_H.
// * An external selector (lambda or CutManager predicate) performs all
//   pre‑selection (e.g. passAll, run‑quality veto, kinematic cuts).
// * Produces:
//     – Optionally prints to stdout
// -----------------------------------------------------------------------------
#ifndef NITROGEN_CORRECTION_H
#define NITROGEN_CORRECTION_H

#include "BranchVars.h"
#include "BranchVarsSim.h"
#include "AnalysisCuts.h"
#include "RunQuality.h"
#include <TChain.h>
#include <TH1D.h>
#include <TF1.h>
#include <string>
#include <functional>

class NitrogenCorrection {
public:
    /**
     * @param cuts       AnalysisCuts with pion_L / pion_H defined (GeV)
     * @param selector   bool predicate to accept events before pion window
     * @param rootFile   output ROOT file (default pion_correction.root)
     */
    NitrogenCorrection(const AnalysisCuts& cuts,
                   const char* kin = "GEN3_He3",
                   const char* rootFile = "nitrogen_correction.root");

    void process(TChain& ch_QE_sim, TChain& ch_N2_sim, BranchVarsSim& vQE, BranchVarsSim& vN2);

private:
    const AnalysisCuts& c_;
    const char* kin_;
    
};

#endif // NITROGEN_CORRECTION_H

