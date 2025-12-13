#include "GenExtraction.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

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

    // write
    std::ofstream out(pathOut(kin_,flag_));
    out << "lambda = " << lambda << '\n'
        << "sigma_lambda_stat = " << dLamStat << '\n'
        << "sigma_lambda_sys  = " << dLamSys  << '\n'
        << "sigma_lambda      = " << dLamTot  << '\n';

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

    std::cout << "[GenExtraction] λ = "<<lambda<<" ± "<<dLamStat
              << " (stat) ± "<<dLamSys<<" (sys) → "<<pathOut(kin_,flag_)<<"\n";
    return true;
}
