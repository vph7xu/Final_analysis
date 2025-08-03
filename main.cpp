#include "BranchVars.h"
#include "BranchVarsSim.h"
#include "CutConfig.h"
#include "CutManager.h"
#include "RunQuality.h"
#include "RawAsymmetry.h"
#include "AnalysisCuts.h"
#include "PlotDXDY.h"
#include "AccidentalCorrection.h"
#include "PionCorrection.h"
#include "InelasticCorrection.h"
#include "NitrogenCorrection.h"
#include "PhysicsAsymmetryMulti.h"
#include "PhysicsAsymmetryCalc.h"
#include "AvgKinematics.h"
#include "GenExtraction.h"

#include <TChain.h>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <data-rootfile> <simQE-rootfile> <simpim-rootfile> <siminelastic-rootfile> <simN2-rootfile> <cuts.json> [run_quality.json|csv] Kinematic_point\n";
        return 1;
    }

    const bool haveRQ = (argc >= 9);
    const int  cutsIdx = argc - (haveRQ ? 3 : 2);   // position of cuts.json
    const char* cutsFile = argv[cutsIdx];
    const char* kin = argv[argc - (haveRQ ? 1 : 0)];
    const char* rqFile   = haveRQ ? argv[argc - 2] : nullptr;

    /* ---------- 1. Build the chain ---------- */
    TChain ch("Tout");               // adjust tree name if needed
    //for (int i = 1; i < cutsIdx; ++i)
    ch.Add(argv[1]); // data file

    TChain chsimQE("Tout");
    chsimQE.Add(argv[2]);

    TChain chsimPim("Tout");
    chsimPim.Add(argv[3]);

    TChain chsiminelastic("Tout");
    chsiminelastic.Add(argv[4]);

    TChain chsimN2("Tout");
    chsimN2.Add(argv[5]);

    /*TChain ch_sim_QE("Tout");
    ch_sim_QE.Add(argv[2]);

    TChain ch_sim_pim("Tout");
    ch_sim_pim.Add(argv[3]);

    TChain ch_sim_inelastic("Tout");
    ch_sim_inelastic.Add(argv[4]);
    */

    /* ---------- 2. Attach branches ---------- */
    BranchVars v;
    v.attach(&ch);

    BranchVarsSim vsimQE;
    vsimQE.attach(&chsimQE);

    BranchVarsSim vsimPim;
    vsimPim.attach(&chsimPim);

    BranchVarsSim vsiminelastic;
    vsiminelastic.attach(&chsiminelastic);

    BranchVarsSim vsimN2;
    vsimN2.attach(&chsimN2);

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

    AnalysisCuts icuts(cutsFile);   // load once

    /* ---------- 5. Run RawAsymmetry module ---------- */
    //RawAsymmetry mod(cuts, icuts, rqPtr, kin);   // default histogram settings inside class
    //mod.process(ch, v);

    //PlotDXDY dxdy(icuts, rqPtr, kin);                  // uses dx/dy/helicity from cuts
    //dxdy.process(ch, v);

    //AccidentalCorrection AccidentalCorrection(icuts, rqPtr, kin);
    //AccidentalCorrection.process(ch,v);

    //PionCorrection PionCorrection(icuts, rqPtr,kin);
    //PionCorrection.process(ch, chsimQE, chsimPim, v, vsimQE, vsimPim);

    //NitrogenCorrection NitrogenCorrection(icuts,kin);
    //NitrogenCorrection.process(chsimQE,chsimN2,vsimQE,vsimN2);

    //InelasticCorrection InelasticCorrection(icuts,rqPtr, kin);
    //InelasticCorrection.process(ch, chsimQE, chsiminelastic,v,vsimQE,vsiminelastic);

    PhysicsAsymmetryCalc phys(kin,1);
    phys.run();

    AvgKinematics avgkin(icuts,rqPtr,kin);
    avgkin.process(ch,v);

    GenExtraction gen(kin,1);
    gen.process();

    return 0;
}