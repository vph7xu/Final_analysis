// -----------------------------------------------------------------------------
// PionCorrection.h — pion‑region raw asymmetry & correction factor
// -----------------------------------------------------------------------------
// * Uses BranchVars fields: eHCAL (cluster energy), helicity, runnum
// * Region is defined by AnalysisCuts.pion_L / pion_H.
// * An external selector (lambda or CutManager predicate) performs all
//   pre‑selection (e.g. passAll, run‑quality veto, kinematic cuts).
// * Produces:
//     – ROOT file storing a TNamed("pion_asym","value,error")
//     – Optionally prints to stdout
// -----------------------------------------------------------------------------
#ifndef PION_CORRECTION_H
#define PION_CORRECTION_H

#include "BranchVars.h"
#include "BranchVarsSim.h"
#include "AnalysisCuts.h"
#include "RunQuality.h"
#include <TChain.h>
#include <TH1D.h>
#include <TF1.h>
#include <string>
#include <functional>

class PionCorrection {
public:
    /**
     * @param cuts       AnalysisCuts with pion_L / pion_H defined (GeV)
     * @param selector   bool predicate to accept events before pion window
     * @param rootFile   output ROOT file (default pion_correction.root)
     */
    PionCorrection(const AnalysisCuts& cuts,
                   const RunQuality* rq  = nullptr,
                   const char* kin = "GEN3_He3",
                   const char* rootFile = "pion_correction.root");

    void process(TChain& ch, TChain& ch_QE_sim, TChain& ch_pim_sim, BranchVars& v, BranchVarsSim& vQE, BranchVarsSim& vPim);

    double asym()  const { return asym_; }
    double error() const { return err_;  }

private:
    const AnalysisCuts& c_;
    const RunQuality* rq_;
    const char* kin_;

    std::string outFile_;

    long long Np_ = 0;
    long long Nm_ = 0;

    double asym_ = 0.0;
    double err_  = 0.0;
    
    TH1D* performFit(TH1D* h_data, TH1D* h_pion, TH1D* h_QE,
                 double &pion_weight, double &qe_weight);
    
};

#endif // PION_CORRECTION_H

