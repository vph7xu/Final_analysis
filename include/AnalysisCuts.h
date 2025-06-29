/***************************************************************
 *  AnalysisCuts.h   — simple POD with direct members
 *  ------------------------------------------------------------
 *  • Call AnalysisCuts cuts("GEN3.json");
 *  • Then read   cuts.W2_L , cuts.coin_H , cuts.helicity , …
 *  • Unknown / missing keys stay at their default (0.0)
 **************************************************************/
#ifndef ANALYSIS_CUTS_H
#define ANALYSIS_CUTS_H

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <string>

struct AnalysisCuts
{
    /* ---- public data members (add more as you need) ---- */
    double W2_L      = 0, W2_H      = 0;
    double dx_L      = 0, dx_H      = 0;
    double dy_L      = 0, dy_H      = 0;
    double eHCAL_L   = 0, eHCAL_H   = 0;
    double coin_L    = 0, coin_H    = 0;
    double coin_ac_L = 0, coin_ac_H = 0;
    double runnum_L  = 0, runnum_H  = 0;
    int    helicity  = 999;                 // 999 ⇒ accept ±1

    /* ---- ctor: load once from the master JSON ---- */
    explicit AnalysisCuts(const std::string& file)
    {
        std::ifstream f(file);
        if (!f) { std::cerr << "[AnalysisCuts] cannot open " << file << '\n'; return; }

        nlohmann::json j;  f >> j;

        auto fetch = [&](auto& ref, const char* key){
            if (!j.contains(key)) return;
            const auto& v = j[key];
            try {
                if      (v.is_number_float())   ref = v.get<double>();
                else if (v.is_number_integer()) ref = static_cast<double>(v.get<int64_t>());
                else if (v.is_string())         ref = std::stod(v.get<std::string>());
            } catch (...) {/* ignore malformed */}
        };

        fetch(W2_L,      "W2_L");      fetch(W2_H,      "W2_H");
        fetch(dx_L,      "dx_L");      fetch(dx_H,      "dx_H");
        fetch(dy_L,      "dy_L");      fetch(dy_H,      "dy_H");
        fetch(eHCAL_L,   "eHCAL_L");   fetch(eHCAL_H,   "eHCAL_H");
        fetch(coin_L,    "coin_L");    fetch(coin_H,    "coin_H");
        fetch(coin_ac_L, "coin_ac_L"); fetch(coin_ac_H, "coin_ac_H");
        fetch(runnum_L,  "runnum_L");  fetch(runnum_H,  "runnum_H");

        if (j.contains("helicity"))
            helicity = j["helicity"].get<int>();
    }
};

#endif  // ANALYSIS_CUTS_H

