#ifndef RAW_ASYMMETRY_H
#define RAW_ASYMMETRY_H

#include "BranchVars.h"
#include "CutManager.h"

#include <TChain.h>
#include <TH1D.h>
#include <TFile.h>

class RawAsymmetry {
public:
    /** Configure histogram and (optionally) the output file name */
    RawAsymmetry(const CutManager& cut,
                 int    nbins  = 120,
                 double xmin   = -1.2,
                 double xmax   =  1.2,
                 const char* outFile = "raw_asymmetry.root");

    /** Run over a chain that already has BranchVars attached */
    void process(TChain& ch, const BranchVars& v);

private:
    const CutManager& cut_;
    TH1D              hA_;
    std::string       outFile_;
};

#endif

