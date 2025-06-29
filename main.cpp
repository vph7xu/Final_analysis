#include "BranchVars.h"
#include "CutConfig.h"
#include "CutManager.h"
#include "RunQuality.h"
#include "RawAsymmetry.h"
#include "AnalysisCuts.h"
#include "PlotDXDY.h"

#include <TChain.h>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <root-file[s]> <cuts.json> [run_quality.json|csv]\n";
        return 1;
    }

    const bool haveRQ = (argc >= 4);
    const int  cutsIdx = argc - (haveRQ ? 2 : 1);   // position of cuts.json
    const char* cutsFile = argv[cutsIdx];
    const char* rqFile   = haveRQ ? argv[argc - 1] : nullptr;

    /* ---------- 1. Build the chain ---------- */
    TChain ch("Tout");               // adjust tree name if needed
    for (int i = 1; i < cutsIdx; ++i)
        ch.Add(argv[i]);

    /* ---------- 2. Attach branches ---------- */
    BranchVars v;
    v.attach(&ch);

    /* ---------- 3. Numeric cuts ---------- */
    CutConfig cfg;
    if (!cfg.load(cutsFile)) {
        std::cerr << "Cannot read cuts file " << cutsFile << '\n';
        return 2;
    }
    CutManager cuts(cfg);

    /* ---------- 4. Optional run-quality ---------- */
    RunQuality rq; const RunQuality* rqPtr = nullptr;
    if (haveRQ) {
        auto isJson = [](const char* name){
            std::string_view sv(name);
            return sv.size() >= 5 && sv.substr(sv.size() - 5) == ".json";
        };
        bool ok = isJson(rqFile) ? rq.loadJSON(rqFile)
                                 : rq.loadCSV (rqFile);
        if (!ok) {
            std::cerr << "[warning] Could not load run-quality file " << rqFile
                      << "; disabling run-quality veto.\n";
        } else {
            rqPtr = &rq;
        }
    }

    /* ---------- 5. Run RawAsymmetry module ---------- */
    RawAsymmetry mod(cuts, rqPtr);   // default histogram settings inside class
    mod.process(ch, v);

    AnalysisCuts icuts(cutsFile);   // load once

    PlotDXDY dxdy(icuts);                  // uses dx/dy/helicity from cuts
    dxdy.process(ch, v);


    return 0;
}


