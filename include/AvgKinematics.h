// -----------------------------------------------------------------------------
// AvgKinematics.h — averages of Q², ε, Px, Pz using unified AnalysisCuts JSON
// -----------------------------------------------------------------------------
#ifndef AVG_KINEMATICS_H
#define AVG_KINEMATICS_H

#include "BranchVars.h"
#include "AnalysisCuts.h"
#include "RunQuality.h"
#include <TChain.h>
#include <TVector3.h>
#include <string>

class AvgKinematics {
public:
    AvgKinematics(const AnalysisCuts& cuts,
                  const RunQuality*   rq       = nullptr,
                  const std::string&  kinTag   = "GEN3_He3"
                  );

    bool process(TChain& ch, BranchVars& v);

private:
    bool good(const BranchVars& v) const;
    struct Ang { double h=0, v=0; };
    Ang fieldCSV() const;                   // read DB/Field_Meas.csv

    std::string      kin_;
    const AnalysisCuts& cuts_;
    const RunQuality*   rq_{};
    TVector3         tgtPolDir_;
};

#endif
