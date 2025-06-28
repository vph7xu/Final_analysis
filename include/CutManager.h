// -----------------------------------------------------------------------------
// CutManager.h — numeric cut evaluator
//   * Supports double and int branches (maps are in .cpp)
//   * Always vetoes events whose BranchVars.helicity is not ±1.
// -----------------------------------------------------------------------------
#ifndef CUT_MANAGER_H
#define CUT_MANAGER_H

#include "BranchVars.h"
#include "CutConfig.h"

class CutManager {
public:
    explicit CutManager(const CutConfig& cfg) : cfg_(cfg) {}

    /**
     * Return true if the event passes all enabled cuts and helicity is ±1.
     */
    bool passAll(const BranchVars& v) const;

private:
    const CutConfig& cfg_;

    // helper: check val vs. cfg_[lowKey / highKey]
    bool inRange(double val, const char* lowKey, const char* highKey) const;
};

#endif // CUT_MANAGER_H
