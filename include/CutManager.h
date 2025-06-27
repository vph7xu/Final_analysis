#ifndef CUT_MANAGER_H
#define CUT_MANAGER_H

#include "BranchVars.h"
#include "CutConfig.h"

/** Evaluates all enabled cuts for a given event */
class CutManager {
public:
    explicit CutManager(const CutConfig& c) : cfg_(c) {}

    /** true → event survives */
    bool passAll(const BranchVars& v) const
    {
        return pass(v, "W2",        "W2_L",     "W2_H")
            && pass(v, "dx",        "dx_L",     "dx_H")
            && pass(v, "dy",        "dy_L",     "dy_H")
            && pass(v, "eHCAL",     "eHCAL_L",  "eHCAL_H")
            && pass(v, "coin_time", "coinL",    "coinH");
    }

private:
    const CutConfig& cfg_;

    /** Generic range test:  cfg[keyLow] < var < cfg[keyHigh]   (if the keys exist) */
    bool pass(const BranchVars& v,
              const char* varName,
              const char* lowKey,
              const char* highKey) const;
};

#endif

