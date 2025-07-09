#include "GenExtraction.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

using std::string;

namespace {
string pathPhys(const string& kin,const string& fl)
{ return "txt/physics_neutron_asymmetry_summary_"+kin+".txt"; }
string pathKin (const string& kin,const string& fl)
{ return "txt/average_kinematic_values_"+kin+".txt"; }
string pathOut (const string& kin,const string& fl)
{ return "txt/results_eHCAL_cut_"+kin+".txt"; }
}

// -----------------------------------------------------------------------------
GenExtraction::GenExtraction(const string& tag,bool cut)
    : kin_(tag), flag_(cut?"1":"0") {}

// -----------------------------------------------------------------------------
double GenExtraction::readKV(const string& file,const string& key)
{
    std::ifstream in(file);
    if(!in){ std::cerr << "[GenExtraction] cannot open " << file << '\n'; return NAN; }
    string line;
    while(std::getline(in,line)){
        if(line.empty() || line[0]=='#') continue;
        std::stringstream ss(line);
        string k; ss >> k;
        if(k!=key) continue;
        if(ss.peek()=='=') ss.get();          // optional '='
        else{
            char c; ss >> c; if(c!='=') ss.unget();
        }
        double v{0}; ss >> v; return v;
    }
    std::cerr << "[GenExtraction] key "<<key<<" not found in "<<file<<"\n";
    return NAN;
}

// -----------------------------------------------------------------------------
bool GenExtraction::process()
{
    const string fPhys = pathPhys(kin_,flag_);
    const string fKin  = pathKin (kin_,flag_);

    // physics‑asym summary
    const double Aphy   = readKV(fPhys,"Aphys");
    const double dAst   = readKV(fPhys,"err_Aphys_stat");
    const double dAsy   = readKV(fPhys,"err_Aphys_sys");

    // averages
    const double Q2     = readKV(fKin,"Q2_avg");
    const double tau    = readKV(fKin,"tau_avg");
    const double eps    = readKV(fKin,"epsilon_avg");
    const double Px     = readKV(fKin,"Px_avg");
    const double Pz     = readKV(fKin,"Pz_avg");

    if(std::isnan(Aphy)||std::isnan(Q2)) return false;

    // algebra (same as original standalone script)
    const double A = (eps/tau)*Aphy;
    const double B = std::sqrt(2*eps*(1-eps)/tau)*Px;
    const double C = Aphy + std::sqrt(1-eps*eps)*Pz;
    const double disc = B*B - 4*A*C;
    if(disc<=0){ std::cerr << "[GenExtraction] negative discriminant\n"; return false; }
    const double lambda = (-B + std::sqrt(disc))/(2*A);

    const double g  = (C/(A*std::sqrt(disc)) + lambda/A)*(eps/tau);
    const double h  = 1.0/disc;
    const double dLamStat = std::sqrt(g*g + h) * dAst;
    const double dLamSys  = std::sqrt(g*g + h) * dAsy;
    const double dLamTot  = std::hypot(dLamStat,dLamSys);

    // write
    std::ofstream out(pathOut(kin_,flag_));
    out << "lambda = " << lambda << '\n'
        << "sigma_lambda_stat = " << dLamStat << '\n'
        << "sigma_lambda_sys  = " << dLamSys  << '\n'
        << "sigma_lambda      = " << dLamTot  << '\n';

    std::cout << "[GenExtraction] λ = "<<lambda<<" ± "<<dLamStat
              << " (stat) ± "<<dLamSys<<" (sys) → "<<pathOut(kin_,flag_)<<"\n";
    return true;
}
