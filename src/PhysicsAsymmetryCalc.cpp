#include "PhysicsAsymmetryCalc.h"
#include "parse.h"   // simple key=value reader

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <functional>
#include <utility>
#include <algorithm>

using std::string; using std::vector; using std::unordered_map;

// -----------------------------------------------------------------------------
PhysicsAsymmetryCalc::PhysicsAsymmetryCalc(const string& kin, bool cut)
    : kin_(kin), flagStr_(cut?"1":"0") {}

// ---------------------------- helper readers ---------------------------------
vector<RunRow> PhysicsAsymmetryCalc::readRawTable(const string& file) const
{
    vector<RunRow> rows; std::ifstream in(file);
    if(!in){ std::cerr << "[PhysCalc] missing " << file << '\n'; return rows; }

    string head; std::getline(in, head);          // skip header
    string line;
    while(std::getline(in, line)){
        std::stringstream ss(line); RunRow r;
        ss >> r.run >> r.Np >> r.Nm >> r.Araw >> r.dAraw
           >> r.beam >> r.dBeam >> r.tgt >> r.dTgt;
        rows.push_back(r);
    }
    return rows;
}


unordered_map<string,double>PhysicsAsymmetryCalc::readSimpleKV(const string& f) const
{
    unordered_map<string,double> m;
    std::ifstream in(f);
    if(!in){ std::cerr << "[PhysCalc] cannot open " << f << " "; return m; }

    string line;
    while(std::getline(in,line)){
        if(line.empty()) continue;
        std::stringstream ss(line);
        string key; ss >> key;
        if(key[0]=='#') continue;             // allow comment lines

        // optional '=' token — eat it if present
        if(ss.peek()=='='){ ss.get(); }
        else{
            // maybe spaces before '='
            char c; ss >> c; if(c!='=') ss.unget(); else { /* '=' already consumed */ }
        }
        double val; ss >> val;
        m[key] = val;
    }
    return m;
} 

// ---------------------------------- run() ------------------------------------
bool PhysicsAsymmetryCalc::run()
{
    // ------- files -----------------------------------------------------------
    string rawFile = "corrections/raw_asymmetry_" + kin_ + ".txt";
    auto   corrFile = [&](const char* stem){
        return string("corrections/") + stem + "Correction_" + kin_ + ".txt"; };

    const string accFile = corrFile("Accidental");
    const string inelFile= corrFile("Inelastic");
    const string nitFile = corrFile("Nitrogen");
    const string pionFile= corrFile("Pion");

    // -------- per‑run raw table ---------------------------------------------
    auto runs = readRawTable(rawFile);
    if (runs.empty()) return false;

    // -------- read TXT correction files directly -------------------------------
    const auto accMap  = readSimpleKV(accFile);
    const auto pionMap = readSimpleKV(pionFile);
    const auto inelMap = readSimpleKV(inelFile);
    const auto nitMap  = readSimpleKV(nitFile);

    auto V=[&](const auto& M,const char* k){ auto it=M.find(k); return it!=M.end()? it->second : 0.0; };

    // constant neutron polarisation
    const double Pn = 0.96;

    // central values
    double Aacc = V(accMap ,"A_acc");
    double Api  = V(pionMap,"A_pi");
    double Ain  = V(inelMap,"A_in");
    double Ap   = V(inelMap,"A_p");
    double Afsi = 0.0;            // not provided yet

    double facc = V(accMap ,"f_acc");
    double fpi  = V(pionMap,"f_pi");
    double fin  = V(inelMap,"f_in");
    double fp   = V(inelMap,"f_p");
    double ffsi = 0.0;
    double fN2  = V(nitMap ,"f_N2");

    // errors
    double errAacc = V(accMap ,"err_A_acc");
    double errApi  = V(pionMap,"err_A_pi");
    double errAin  = V(inelMap,"err_A_in");
    double errAp   = V(inelMap,"err_A_p");
    double errAfsi = 0.0;

    double errfacc = V(accMap ,"err_f_acc");
    double errfpi  = V(pionMap,"err_f_pi");
    double errfin  = V(inelMap,"err_f_in");
    double errfp   = V(inelMap,"err_f_p");
    double errffsi = 0.0;
    double errfN2  = V(nitMap ,"err_f_N2");

    const double fnDil = 1 - facc - fN2 - fpi - fin - fp - ffsi;

    

    // -------------- avg pol placeholders (adjust as needed) -----------------
    double avgHe3   = 0.45, dAvgHe3  = 0.05;
    double avgBeam  = 0.85, dAvgBeam = 0.05;
    double dPn      = 0.05;                      // uncertainty on Pn

    /*    // ---------- avg polarisation constants (FIX THIS) ----------------------------------
    //auto avgPol = readSimpleKV(avgPolFile);
    double avgHe3   = 0.45;//avgPol["avg_He3pol"];
    double dAvgHe3  = 0.05;//avgPol["err_avg_He3pol"];
    double avgBeam  = 0.85;//avgPol["avg_beampol"];
    double dAvgBeam = 0.05;//avgPol["err_avg_beampol"];
    double dPn      = 0.05;//avgPol["err_avg_Pn"];   // err on Pn
    */


    // -------------- per‑run output -----------------------------------------
    std::ofstream per("txt/physics_neutron_asymmetry_results_per_run_"+kin_+".txt");
    per << "#run A_phys dA_stat\n";

    double num=0, den=0; const double eps=1e-10;

    for(const auto& r: runs){
        const double denom = r.tgt*0.01 * r.beam*0.01 * Pn * fnDil;
        if(std::fabs(denom)<eps) continue;

        const double Aphys = (r.Araw - facc*Aacc - fpi*Api - fin*Ain - fp*Ap - ffsi*Afsi)/denom;
        const double dA    = r.dAraw / denom;

        per << r.run << ' ' << Aphys << ' ' << dA << '\n';

        if(dA>eps){ num += Aphys/(dA*dA); den += 1.0/(dA*dA); }
    }
    per.close();

    if(den==0){ std::cerr << "[PhysCalc] no valid runs\n"; return false; }
    A_ = num/den; dA_ = 1.0/std::sqrt(den);

    // -------------- systematic error ---------------------------------------
    double errSys = std::sqrt(
          ( facc*facc*errAacc*errAacc + fpi*fpi*errApi*errApi +
            fin*fin*errAin*errAin   + fp*fp*errAp*errAp + ffsi*ffsi*errAfsi*errAfsi )/(Pn*Pn*fnDil*fnDil)
        + std::pow(A_*errfN2/fnDil,2)
        + std::pow((Pn*A_-Aacc)*errfacc/(Pn*fnDil),2)
        + std::pow((Pn*A_-Api )*errfpi /(Pn*fnDil),2)
        + std::pow((Pn*A_-Ain )*errfin /(Pn*fnDil),2)
        + std::pow((Pn*A_-Ap  )*errfp  /(Pn*fnDil),2)
        + std::pow((Pn*A_-Afsi)*errffsi/(Pn*fnDil),2)
        + A_*A_*( dAvgHe3*dAvgHe3 + (dPn/Pn)*(dPn/Pn) + (dAvgBeam/avgBeam)*(dAvgBeam/avgBeam) ) );

    // -------------- summary -------------------------------------------------
    std::ofstream sum("txt/physics_neutron_asymmetry_summary_"+kin_+".txt");
    sum << "A_phys_avg " << A_ << ' ' << dA_ << '\n'
        << "err_sys     " << errSys << '\n';

    std::cout << "[PhysCalc] A_phys_avg = " << A_ << " ± " << dA_ << " (stat) ± " << errSys << " (sys)\n";
    std::cout << "Ain="<<Ain<<" Ap="<<Ap<<" Api="<<Api<<" Aacc="<<Aacc<<std::endl;


        /* --- quick diagnostic printout ------------------------------------------ */
    auto dumpMap = [](const std::string& title,
                      const std::unordered_map<std::string,double>& m)
    {
        std::cout << "\n=== " << title << " contents ===\n";
        for (const auto& kv : m)
            std::cout << kv.first << " = " << kv.second << '\n';
    };

    dumpMap("Accidental", accMap);
    dumpMap("Pion      ", pionMap);
    dumpMap("Inelastic ", inelMap);
    dumpMap("Nitrogen  ", nitMap);
    /* ------------------------------------------------------------------------ */


    return true;
}
