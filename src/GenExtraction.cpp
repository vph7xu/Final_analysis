#include "GenExtraction.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

using std::string;

namespace {
string pathPhys(const string& kin,const string& fl)
{ return "txt/"+kin+"/physics_neutron_asymmetry_summary_"+kin+".txt"; }
string pathKin (const string& kin,const string& fl)
{ return "txt/"+kin+"/average_kinematic_values_"+kin+".txt"; }
string pathAcc (const string& kin,const string& fl)
{ return "corrections/"+kin+"/AccidentalCorrection_"+kin+".txt"; }
string pathPion(const string& kin,const string& fl)
{ return "corrections/"+kin+"/PionCorrection_"+kin+".txt"; }
string pathInel(const string& kin,const string& fl)
{ return "corrections/"+kin+"/InelasticCorrection_"+kin+".txt"; }
string pathN2  (const string& kin,const string& fl)
{ return "corrections/"+kin+"/NitrogenCorrection_"+kin+".txt"; }
string pathOut (const string& kin,const string& fl)
{ return "txt/"+kin+"/results_eHCAL_cut_"+kin+".txt"; }
string pathSummary (const string& kin,const string& fl)
{ return "txt/"+kin+"/summary_"+kin+".txt"; }

// --- Neutron form-factor lookup (tabulated parameterization) ---------------
// Expected columns (whitespace separated, '#' comments allowed):
//   Q2  GEn/GD  dGEn/GD  dGEn_Par/GD  GMn/muGD  dGMn/muGD  dGMn_Par/muGD
// We only need GMn/muGD (and its uncertainty). Here "muGD" means (mu_n * G_D).
struct NeutronLookupTable {
    std::vector<double> q2;         // GeV^2
    std::vector<double> gmn_muGD;   // GMn/(mu_n * G_D)
    std::vector<double> dgmn_muGD;  // uncertainty on GMn/(mu_n * G_D)
};

inline double dipoleGD(double Q2){
    // Standard dipole: G_D(Q^2) = (1 + Q^2/0.71)^{-2}
    const double a = 1.0 + Q2/0.71;
    return 1.0/(a*a);
}

static bool readNeutronLookup(const string& file, NeutronLookupTable& tab){
    std::ifstream in(file);
    if(!in) return false;

    tab = NeutronLookupTable{};
    string line;
    while(std::getline(in,line)){
        if(line.empty()) continue;
        const auto firstNon = line.find_first_not_of(" \t\r\n");
        if(firstNon==string::npos) continue;
        if(line[firstNon]=='#') continue;

        std::stringstream ss(line);
        double Q2=0, GEn_GD=0, dGEn_GD=0, dGEnPar_GD=0, GMn_muGD=0, dGMn_muGD=0, dGMnPar_muGD=0;
        if(!(ss >> Q2 >> GEn_GD >> dGEn_GD >> dGEnPar_GD >> GMn_muGD >> dGMn_muGD >> dGMnPar_muGD)) continue;
        tab.q2.push_back(Q2);
        tab.gmn_muGD.push_back(GMn_muGD);
        // Prefer the parameterization uncertainty (last column) if sensible; otherwise use dGMn/muGD.
        const double s = (dGMnPar_muGD>0 ? dGMnPar_muGD : dGMn_muGD);
        tab.dgmn_muGD.push_back(std::abs(s));
    }
    return !tab.q2.empty();
}

static bool interpolateGMn_overMu(const NeutronLookupTable& tab, double Q2,
                                 double& GMn_overMu, double& dGMn_overMu){
    // Interpolate GMn/(mu_n*G_D) at Q2, then multiply by G_D(Q2).
    // Output is GMn/mu_n (i.e., divided by mu_n), so caller can apply mun.
    if(tab.q2.empty()) return false;
    const double GD = dipoleGD(Q2);

    // Clamp to the table range.
    if(Q2 <= tab.q2.front()){
        GMn_overMu  = tab.gmn_muGD.front()  * GD;
        dGMn_overMu = tab.dgmn_muGD.front() * GD;
        return true;
    }
    if(Q2 >= tab.q2.back()){
        GMn_overMu  = tab.gmn_muGD.back()  * GD;
        dGMn_overMu = tab.dgmn_muGD.back() * GD;
        return true;
    }

    auto it = std::upper_bound(tab.q2.begin(), tab.q2.end(), Q2);
    const size_t i1 = std::distance(tab.q2.begin(), it);
    const size_t i0 = i1 - 1;

    const double x0 = tab.q2[i0], x1 = tab.q2[i1];
    const double t  = (Q2 - x0)/(x1 - x0);

    const double y0 = tab.gmn_muGD[i0],  y1 = tab.gmn_muGD[i1];
    const double s0 = tab.dgmn_muGD[i0], s1 = tab.dgmn_muGD[i1];

    const double y  = y0 + t*(y1 - y0);
    const double s  = s0 + t*(s1 - s0);

    GMn_overMu  = y * GD;
    dGMn_overMu = std::abs(s) * GD;
    return true;
}
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

    static const double mun = -1.9130427; // +/- 5E-7 muN

    const string fPhys = pathPhys(kin_,flag_);
    const string fKin  = pathKin (kin_,flag_);
    const string fAcc  = pathAcc (kin_,flag_);
    const string fPion = pathPion(kin_,flag_);
    const string fInel = pathInel(kin_,flag_);
    const string fNitrogen   = pathN2  (kin_,flag_);

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

    // corrections
    const double facc   = readKV(fAcc,"f_acc");
    const double fpi    = readKV(fPion,"f_pi");
    const double fin    = readKV(fInel,"f_in");
    const double fp     = readKV(fInel,"f_p");
    const double fN2    = readKV(fNitrogen,"f_N2");
    const double ffsi   = 0.0287;  // not finalized

    const double errfacc = readKV(fAcc,"err_f_acc");
    const double errfpi  = readKV(fPion,"err_f_pi");
    const double errfin  = readKV(fInel,"err_f_in");
    const double errfp   = readKV(fInel,"err_f_p");
    const double errfN2  = readKV(fNitrogen,"err_f_N2");
    const double errffsi = 0.0026; // not finalized

    const double Aacc = readKV(fAcc,"A_acc");
    const double Api  = readKV(fPion,"A_pi");
    const double Ain  = readKV(fInel,"A_in");
    const double Ap   = readKV(fInel,"A_p");
    const double Afsi = 0.0003;    // not finalized

    const double errAacc = readKV(fAcc,"err_A_acc");
    const double errApi  = readKV(fPion,"err_A_pi");
    const double errAin  = readKV(fInel,"err_A_in");
    const double errAp   = readKV(fInel,"err_A_p");
    const double errAfsi = 0.0005; // not finalized


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
    const double dLamStat = std::sqrt(pow(g + std::sqrt(h),2)) * dAst;
    const double dLamSys  = std::sqrt(pow(g + std::sqrt(h),2)) * dAsy;
    const double dLamTot  = std::hypot(dLamStat,dLamSys);

    // --- NEW: read GMn(Q2) from neutron_lookup.dat and compute GEn ----------
    // The table provides GMn/(mu_n * G_D). We interpolate and convert to GMn.
    NeutronLookupTable tab;
    bool okTab = readNeutronLookup("World_data/neutron_lookup.dat", tab);
    if(!okTab) okTab = readNeutronLookup("neutron_look.dat", tab);

    double GMn_overMu = NAN, dGMn_overMu = NAN; // GMn/mu_n
    if(okTab) okTab = interpolateGMn_overMu(tab, Q2, GMn_overMu, dGMn_overMu);
    if(!okTab){
        std::cerr << "[GenExtraction] WARNING: could not read neutron_lookup.dat (or neutron_look.dat). "
                     "Skipping GEn extraction from GMn parameterization.\n";
    }

    const double GD = dipoleGD(Q2);
    const double GMn_abs  = okTab ? (mun * GMn_overMu) : NAN;
    const double dGMn_abs = okTab ? (std::abs(mun) * dGMn_overMu) : NAN;

    // GEn = (GEn/GMn) * GMn = lambda * GMn
    const double GEn_abs   = okTab ? (lambda * GMn_abs) : NAN;
    // Propagate uncertainties (assume GMn parameterization uncertainty is independent of lambda)
    const double dGEn_stat = okTab ? (std::abs(GMn_abs) * dLamStat) : NAN;
    const double dGEn_sys  = okTab ? (std::abs(GMn_abs) * dLamSys)  : NAN;
    const double dGEn_GMn  = okTab ? (std::abs(lambda) * dGMn_abs)  : NAN;
    const double dGEn_tot  = okTab ? std::hypot(std::hypot(dGEn_stat, dGEn_sys), dGEn_GMn) : NAN;

    // write
    std::ofstream out(pathOut(kin_,flag_));
    out << "lambda = " << lambda << '\n'
        << "sigma_lambda_stat = " << dLamStat << '\n'
        << "sigma_lambda_sys  = " << dLamSys  << '\n'
        << "sigma_lambda      = " << dLamTot  << '\n';

    if(okTab){
        out << "GD = " << GD << '\n'
            << "GMn = " << GMn_abs << '\n'
            << "sigma_GMn_param = " << dGMn_abs << '\n'
            << "GEn = " << GEn_abs << '\n'
            << "sigma_GEn_stat = " << dGEn_stat << '\n'
            << "sigma_GEn_sys  = " << dGEn_sys  << '\n'
            << "sigma_GEn_GMnParam = " << dGEn_GMn << '\n'
            << "sigma_GEn      = " << dGEn_tot << '\n';
    }

    std::ofstream sum(pathSummary(kin_,flag_));
    sum << "Q2 = " << Q2 << '\n'
        << "tau = " << tau << '\n'
        << "epsilon  = " << eps << '\n'
        << "Px      = " << Px << '\n'
        << "Pz      = " << Pz << '\n'
        << "Aacc      = " <<Aacc<<"+-" << errAacc << '\n'
        << "Api      = " <<Api  << "+-" << errApi << '\n'
        << "Ain      = " <<Ain  << "+-" << errAin << '\n'
        << "Ap       = " <<Ap   << "+-" << errAp  << '\n'
        << "Afsi     = " <<Afsi << "+-" << errAfsi << '\n'
        << "facc     = " <<facc << "+-" << errfacc << '\n'
        << "fpi      = " <<fpi  << "+-" << errfpi  << '\n'
        << "fin      = " <<fin  << "+-" << errfin  << '\n'
        << "fp       = " <<fp   << "+-" << errfp   << '\n'
        << "fN2      = " <<fN2  << "+-" << errfN2  << '\n'
        << "ffsi     = " <<ffsi << "+-" << errffsi << '\n'
        << "lambda = " << lambda << "+-" << dLamStat << "(stat) +-"<< dLamSys <<"(sys)"<< '\n'
        << "mu_n*G_E/G_M = " << lambda * mun << "+-" << dLamStat * mun <<"(stat) +-"<< dLamSys * mun <<"(sys)"<< '\n';

    if(okTab){
        sum << "GD = " << GD << '\n'
            << "GMn(param) = " << GMn_abs << "+-" << dGMn_abs << " (param)" << '\n'
            << "GEn = " << GEn_abs
            << "+-" << dGEn_stat << " (stat)"
            << " +-" << dGEn_sys  << " (sys)"
            << " +-" << dGEn_GMn  << " (GMn param)"
            << " => " << dGEn_tot << " (total)" << '\n';
    }

    std::cout << "[GenExtraction] λ = "<<lambda<<" ± "<<dLamStat
              << " (stat) ± "<<dLamSys<<" (sys) → "<<pathOut(kin_,flag_)<<"\n";
    return true;
}
