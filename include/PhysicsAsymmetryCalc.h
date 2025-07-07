// -----------------------------------------------------------------------------
// PhysicsAsymmetryCalc.h  —  wrap PhysicsAsymmetryMulti for one‑shot use
// -----------------------------------------------------------------------------
// * Constructor takes the kinematic tag (e.g. "GEN3_He3") and the eHCAL‑cut
//   flag (bool -> 0/1) so it knows the standard file locations used in your
//   workflow:
//       DB/corrections/<kin>_corrections_eHCAL_cut_<flag>.txt
//       txt/<kin>_raw_neutron_asymmetry_results_eHCAL_cut_<flag>.txt
//       txt/<kin>_average_polarization_results_eHCAL_cut_<flag>.txt
//       … plus all per‑correction ROOT/TXT outputs already produced.
//
// * run() loads every file, evaluates the hard‑wired formula (identical to
//   the standalone script you posted), and writes the final physics asymmetry
//   to       txt/<kin>_physics_neutron_asymmetry_summary_eHCAL_cut_<flag>.txt
//
// * After run() you may query value() / error() programmatically.
// -----------------------------------------------------------------------------
#ifndef PHYSICS_ASYM_CALC_H
#define PHYSICS_ASYM_CALC_H

#include "PhysicsAsymmetryMulti.h"
#include <string>
#include <vector>
#include <unordered_map>

// one row from the run-by-run raw-asymmetry table
struct RunRow {
    int    run;
    long   Np, Nm;
    double Araw, dAraw;   // raw asymmetry ±stat
    double beam, dBeam;   // beam polarisation ±σ
    double tgt , dTgt;    // ³He polarisation ±σ
};

class PhysicsAsymmetryCalc {
public:
    PhysicsAsymmetryCalc(const std::string& kinTag, bool eHCALcut);
    bool run();                      // returns true on success

    double value() const { return A_; }
    double error() const { return dA_; }

private:

    std::vector<RunRow> readRawTable(const std::string& file) const;
    std::unordered_map<std::string,double>
        readSimpleKV(const std::string& file) const;

    std::string kin_; std::string flagStr_;
    PhysicsAsymmetryMulti multi_;
    double A_=0, dA_=0;
};

#endif // PHYSICS_ASYM_CALC_H

