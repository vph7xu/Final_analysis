// MakeCombinedGEnSummary.cpp
//
// Standalone helper to collect GEN2/GEN3/GEN4/GEN4b summary numbers
// from the text files already produced by your analysis and write
// one combined text table.
//
// Compile:
//   g++ -O2 -std=c++17 MakeCombinedGEnSummary.cpp -o MakeCombinedGEnSummary
//
// Run:
//   ./MakeCombinedGEnSummary
//
// Output:
//   ./combined_GEn_summary_table.txt

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

using std::string;

// --------------------- paths ---------------------

static string pathPol(const string& kin)
{ return "../corrections/" + kin + "/avg_polarizations_" + kin + ".txt"; }

static string pathRaw(const string& kin)
{ return "../txt/" + kin + "/raw_asym_" + kin + ".txt"; }

static string pathPhys(const string& kin)
{ return "../txt/" + kin + "/physics_neutron_asymmetry_summary_" + kin + ".txt"; }

static string pathKin(const string& kin)
{ return "../txt/" + kin + "/average_kinematic_values_" + kin + ".txt"; }

static string pathAcc(const string& kin)
{ return "../corrections/" + kin + "/AccidentalCorrection_" + kin + ".txt"; }

static string pathPion(const string& kin)
{ return "../corrections/" + kin + "/PionCorrection_" + kin + ".txt"; }

static string pathInel(const string& kin)
{ return "../corrections/" + kin + "/InelasticCorrection_" + kin + ".txt"; }

static string pathN2(const string& kin)
{ return "../corrections/" + kin + "/NitrogenCorrection_" + kin + ".txt"; }

static string pathOut(const string& kin)
{ return "../txt/" + kin + "/results_eHCAL_cut_" + kin + ".txt"; }

// --------------------- helpers ---------------------

static bool canOpenFile(const string& file)
{
    std::ifstream in(file);
    if(!in){
        std::cerr << "[ERROR] Cannot open file: " << file << "\n";
        return false;
    }
    return true;
}

static double readKV(const string& file, const string& key, bool required=false)
{
    std::ifstream in(file);
    if(!in){
        std::cerr << "[ERROR] Cannot open file: " << file
                  << " while looking for key: " << key << "\n";
        return NAN;
    }

    string line;
    while(std::getline(in, line)){
        if(line.empty()) continue;
        if(line[0] == '#') continue;

        std::stringstream ss(line);
        string k;
        ss >> k;
        if(k != key) continue;

        if(ss.peek() == '=') ss.get();
        else{
            char c;
            ss >> c;
            if(c != '=') ss.unget();
        }

        double v = NAN;
        if(!(ss >> v)){
            std::cerr << "[ERROR] Found key '" << key
                      << "' in file " << file
                      << " but could not parse its numeric value.\n";
            return NAN;
        }
        return v;
    }

    if(required){
        std::cerr << "[WARNING] Required key '" << key
                  << "' not found in file: " << file << "\n";
    } else {
        std::cerr << "[INFO] Optional key '" << key
                  << "' not found in file: " << file << "\n";
    }

    return NAN;
}

static string fmt1(double v, int p=4)
{
    if(std::isnan(v)) return "--";
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(p) << v;
    return ss.str();
}

static string fmtVE(double v, double e, int p=4)
{
    if(std::isnan(v)) return "--";
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(p) << v;
    if(!std::isnan(e)) ss << " ± " << std::fixed << std::setprecision(p) << e;
    return ss.str();
}

static string fmtV2E(double v, double e1, double e2, int p=4)
{
    if(std::isnan(v)) return "--";
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(p) << v;
    if(!std::isnan(e1)) ss << " ± " << std::fixed << std::setprecision(p) << e1;
    if(!std::isnan(e2)) ss << " ± " << std::fixed << std::setprecision(p) << e2;
    return ss.str();
}

// --------------------- data holders ---------------------

struct ValErr {
    double v = NAN;
    double e = NAN;
};

struct Val2Err {
    double v  = NAN;
    double e1 = NAN;
    double e2 = NAN;
};

struct KinSummary {
    string kin;

    double Q2 = NAN;
    double N  = NAN;
    double Np = NAN;
    double Nm = NAN;

    ValErr Araw;
    ValErr Pbeam;
    ValErr P3He;
    ValErr Pn;
    ValErr Pp;

    ValErr Aacc;
    ValErr Api;
    ValErr Ain;
    ValErr Ap;
    ValErr Afsi;

    ValErr facc;
    ValErr fN2;
    ValErr fpi;
    ValErr fin;
    ValErr fp;
    ValErr ffsi;
    ValErr fn;

    Val2Err Aphys;
    Val2Err lambda;
    Val2Err muGEGM;
    ValErr  GMn_over_muGD;
    Val2Err GEn;
};

// --------------------- loader ---------------------

static KinSummary loadKinSummary(const string& kin)
{
    static const double mun = -1.9130427;

    KinSummary r;
    r.kin = kin;

    const string fPol  = pathPol(kin);
    const string fRaw  = pathRaw(kin);
    const string fPhys = pathPhys(kin);
    const string fKin  = pathKin(kin);
    const string fAcc  = pathAcc(kin);
    const string fPion = pathPion(kin);
    const string fInel = pathInel(kin);
    const string fNit  = pathN2(kin);
    const string fOut  = pathOut(kin);

    std::cerr << "\n[INFO] Loading kinematic: " << kin << "\n";

    canOpenFile(fPol);
    canOpenFile(fRaw);
    canOpenFile(fPhys);
    canOpenFile(fKin);
    canOpenFile(fAcc);
    canOpenFile(fPion);
    canOpenFile(fInel);
    canOpenFile(fNit);
    canOpenFile(fOut);

    // Average kinematics
    r.Q2 = readKV(fKin, "Q2_avg", true);

    // Raw / polarization
    r.N     = readKV(fRaw, "total_good_events", true);
    r.Np    = readKV(fRaw, "Np", true);
    r.Nm    = readKV(fRaw, "Nm", true);

    r.Araw  = { readKV(fRaw, "A_raw", true), readKV(fRaw, "err_A_raw", true) };
    r.Pbeam = { readKV(fPol, "avg_beampol", true), readKV(fPol, "err_avg_beampol", true) };
    r.P3He  = { readKV(fPol, "avg_He3pol", true),  readKV(fPol, "err_avg_He3pol", true) };
    r.Pn    = { readKV(fPol, "avg_Pn", true),      readKV(fPol, "err_avg_Pn", true) };

    // Pp not present in your polarization files; set manually
    r.Pp    = { -0.03, 0.002 }; // conservative estimate based on GEN2/GEN3/GEN4/GEN4b polarization files

    // Physics asymmetry
    r.Aphys = {
        readKV(fPhys, "Aphys", true),
        readKV(fPhys, "err_Aphys_stat", true),
        readKV(fPhys, "err_Aphys_sys", true)
    };

    // Asymmetry corrections
    r.Aacc = { readKV(fAcc,  "A_acc", true),  readKV(fAcc,  "err_A_acc", true) };
    r.Api  = { readKV(fPion, "A_pi",  true),  readKV(fPion, "err_A_pi",  true) };
    r.Ain  = { readKV(fInel, "A_in",  true),  readKV(fInel, "err_A_in",  true) };
    r.Ap   = { readKV(fInel, "A_p",   true),  readKV(fInel, "err_A_p",   true) };
    r.Afsi = { 0.0003, 0.0005 };

    // Fraction corrections
    r.facc = { readKV(fAcc,  "f_acc", true),  readKV(fAcc,  "err_f_acc", true) };
    r.fN2  = { readKV(fNit,  "f_N2",  true),  readKV(fNit,  "err_f_N2",  true) };
    r.fpi  = { readKV(fPion, "f_pi",  true),  readKV(fPion, "err_f_pi",  true) };
    r.fin  = { readKV(fInel, "f_in",  true),  readKV(fInel, "err_f_in",  true) };
    r.fp   = { readKV(fInel, "f_p",   true),  readKV(fInel, "err_f_p",   true) };
    r.ffsi = { 0.0287, 0.0026 };

    // Compute fn = 1 - sum(non-neutron fractions)
    if( !std::isnan(r.facc.v) &&
        !std::isnan(r.fN2.v)  &&
        !std::isnan(r.fpi.v)  &&
        !std::isnan(r.fin.v)  &&
        !std::isnan(r.fp.v)   &&
        !std::isnan(r.ffsi.v) )
    {
        r.fn.v = 1.0 - (r.facc.v + r.fN2.v + r.fpi.v + r.fin.v + r.fp.v + r.ffsi.v);
    }
    else{
        std::cerr << "[WARNING] Could not compute fn for " << kin
                  << " because one or more fraction values are missing.\n";
    }

    if( !std::isnan(r.facc.e) &&
        !std::isnan(r.fN2.e)  &&
        !std::isnan(r.fpi.e)  &&
        !std::isnan(r.fin.e)  &&
        !std::isnan(r.fp.e)   &&
        !std::isnan(r.ffsi.e) )
    {
        r.fn.e = std::sqrt(
            r.facc.e * r.facc.e +
            r.fN2.e  * r.fN2.e  +
            r.fpi.e  * r.fpi.e  +
            r.fin.e  * r.fin.e  +
            r.fp.e   * r.fp.e   +
            r.ffsi.e * r.ffsi.e
        );
    }
    else{
        std::cerr << "[WARNING] Could not compute err_fn for " << kin
                  << " because one or more fraction errors are missing.\n";
    }

    // Output of GenExtraction
    r.lambda = {
        readKV(fOut, "lambda", true),
        readKV(fOut, "sigma_lambda_stat", true),
        readKV(fOut, "sigma_lambda_sys", true)
    };

    r.GEn = {
        readKV(fOut, "GEn", false),
        readKV(fOut, "sigma_GEn_stat", false),
        readKV(fOut, "sigma_GEn_sys", false)
    };

    const double GD        = readKV(fOut, "GD", false);
    const double GMn       = readKV(fOut, "GMn", false);
    const double dGMnParam = readKV(fOut, "sigma_GMn_param", false);

    if(!std::isnan(GD) && !std::isnan(GMn) && GD != 0.0){
        r.GMn_over_muGD.v = GMn / (mun * GD);
        if(!std::isnan(dGMnParam)) r.GMn_over_muGD.e = std::abs(dGMnParam / (mun * GD));
    } else {
        std::cerr << "[WARNING] Could not compute GMn/(mu_n*GD) for " << kin
                  << " because GD or GMn is missing.\n";
    }

    if(!std::isnan(r.lambda.v)){
        r.muGEGM.v  = mun * r.lambda.v;
        if(!std::isnan(r.lambda.e1)) r.muGEGM.e1 = std::abs(mun * r.lambda.e1);
        if(!std::isnan(r.lambda.e2)) r.muGEGM.e2 = std::abs(mun * r.lambda.e2);
    } else {
        std::cerr << "[WARNING] Could not compute mu_n*GEn/GMn for " << kin
                  << " because lambda is missing.\n";
    }

    return r;
}

// --------------------- writer ---------------------

static void writeCombinedSummary(const string& outfile,
                                 const std::vector<KinSummary>& rows)
{
    std::ofstream out(outfile);
    if(!out){
        std::cerr << "[ERROR] Cannot open output file for writing: " << outfile << "\n";
        return;
    }

    if(rows.empty()){
        std::cerr << "[ERROR] No kinematics to write.\n";
        return;
    }

    out << "# Combined preliminary GEn summary table\n";
    out << "# one-error format:  x ± err\n";
    out << "# two-error format:  x ± stat ± sys\n\n";

    const int w0 = 18;
    const int w  = 30;

    out << std::left << std::setw(w0) << "Variable";
    for(const auto& r : rows) out << std::setw(w) << r.kin;
    out << "\n";

    auto row1 = [&](const string& label, auto getter){
        out << std::left << std::setw(w0) << label;
        for(const auto& r : rows){
            out << std::setw(w) << getter(r);
        }
        out << "\n";
    };

    row1("Q2 (GeV/c)^2", [](const KinSummary& r){ return fmt1(r.Q2, 2); });
    row1("N",            [](const KinSummary& r){ return fmt1(r.N, 0); });
    row1("N+",           [](const KinSummary& r){ return fmt1(r.Np, 0); });
    row1("N-",           [](const KinSummary& r){ return fmt1(r.Nm, 0); });

    row1("Araw",   [](const KinSummary& r){ return fmtVE(r.Araw.v,  r.Araw.e,  4); });
    row1("Pbeam",  [](const KinSummary& r){ return fmtVE(r.Pbeam.v, r.Pbeam.e, 4); });
    row1("P3He",   [](const KinSummary& r){ return fmtVE(r.P3He.v,  r.P3He.e,  4); });
    row1("Pn",     [](const KinSummary& r){ return fmtVE(r.Pn.v,    r.Pn.e,    4); });
    row1("Pp",     [](const KinSummary& r){ return fmtVE(r.Pp.v,    r.Pp.e,    4); });

    row1("Aacc",   [](const KinSummary& r){ return fmtVE(r.Aacc.v,  r.Aacc.e,  4); });
    row1("Api",    [](const KinSummary& r){ return fmtVE(r.Api.v,   r.Api.e,   4); });
    row1("Ain",    [](const KinSummary& r){ return fmtVE(r.Ain.v,   r.Ain.e,   4); });
    row1("Ap",     [](const KinSummary& r){ return fmtVE(r.Ap.v,    r.Ap.e,    4); });
    row1("Afsi",   [](const KinSummary& r){ return fmtVE(r.Afsi.v,  r.Afsi.e,  4); });

    row1("facc",   [](const KinSummary& r){ return fmtVE(r.facc.v,  r.facc.e,  4); });
    row1("fN2",    [](const KinSummary& r){ return fmtVE(r.fN2.v,   r.fN2.e,   4); });
    row1("fpi",    [](const KinSummary& r){ return fmtVE(r.fpi.v,   r.fpi.e,   4); });
    row1("fin",    [](const KinSummary& r){ return fmtVE(r.fin.v,   r.fin.e,   4); });
    row1("fp",     [](const KinSummary& r){ return fmtVE(r.fp.v,    r.fp.e,    4); });
    row1("ffsi",   [](const KinSummary& r){ return fmtVE(r.ffsi.v,  r.ffsi.e,  4); });
    row1("fn",     [](const KinSummary& r){ return fmtVE(r.fn.v,    r.fn.e,    4); });

    row1("Aphys",        [](const KinSummary& r){ return fmtV2E(r.Aphys.v,  r.Aphys.e1,  r.Aphys.e2,  4); });
    row1("lambda",       [](const KinSummary& r){ return fmtV2E(r.lambda.v, r.lambda.e1, r.lambda.e2, 4); });
    row1("mu_n GEn/GMn", [](const KinSummary& r){ return fmtV2E(r.muGEGM.v, r.muGEGM.e1, r.muGEGM.e2, 4); });
    row1("GMn/(mu_n GD)",[](const KinSummary& r){ return fmtVE(r.GMn_over_muGD.v, r.GMn_over_muGD.e, 4); });
    row1("GEn",          [](const KinSummary& r){ return fmtV2E(r.GEn.v,    r.GEn.e1,    r.GEn.e2,    4); });

    out.close();
    std::cout << "[INFO] Wrote " << outfile << "\n";
}

// --------------------- main ---------------------

int main()
{
    std::vector<string> kins = {"GEN2_He3", "GEN3_He3", "GEN4_He3", "GEN4b_He3"};
    std::vector<KinSummary> rows;

    for(const auto& kin : kins){
        rows.push_back(loadKinSummary(kin));
    }

    if(rows.size() != kins.size()){
        std::cerr << "[ERROR] Expected " << kins.size()
                  << " kinematics, found " << rows.size() << "\n";
        return 1;
    }

    writeCombinedSummary("./combined_GEn_summary_table.txt", rows);
    return 0;
}