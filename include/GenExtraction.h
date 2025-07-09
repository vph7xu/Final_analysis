// -----------------------------------------------------------------------------
// GenExtraction.h  —  solve for λ = GE/GM from physics‐asym summary & kinematic
// averages.  Uses plain key = value TXT files; no dependency on parse.h.
// -----------------------------------------------------------------------------
#ifndef GEN_EXTRACTION_H
#define GEN_EXTRACTION_H

#include <string>

class GenExtraction {
public:
    /**
     * @param kinTag   base tag (e.g. "GEN3"); physics file is
     *                  txt/<kin>_He3_physics_neutron_asymmetry_summary_eHCAL_cut_<flag>.txt
     * @param cutFlag  true → uses files with suffix "cut_1", false → "cut_0"
     */
    GenExtraction(const std::string& kinTag, bool cutFlag);

    /// read the two TXT inputs, compute λ and its errors, write results file
    bool process();

private:
    /// mini helper: fetch a single numeric value from a key=value TXT file
    static double readKV(const std::string& file, const std::string& key);

    std::string kin_;      // e.g. "GEN3"
    std::string flag_;     // "0" or "1"
};

#endif // GEN_EXTRACTION_H

