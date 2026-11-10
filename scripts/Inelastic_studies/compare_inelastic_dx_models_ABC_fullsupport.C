// Standalone GEn-II 1D dx template/background comparison in ROOT C++.
// Run from Final_analysis repository root:
//   root -l -b -q 'compare_inelastic_dx_models_ABC_fullsupport.C+()'
//   root -l -b -q 'compare_inelastic_dx_models_ABC_fullsupport.C+("", "", "", "cuts/GEN3.json", "DB/run_quality.json", "out_plots/inelastic_dx_background_models_GEN3_ABC_fullsupport", 100, true, false, "2,3,4,5,6")'
// Empty path arguments keep the defaults defined below.
// This full-support variant fills extended MC templates before fitting a
// central dx window. A fitted shift no longer creates an internal hard cutoff.
// Requirements: ROOT with Minuit2 (ROOT 6.x), C++17. No external JSON library.

#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dxbg {

static const char* kData =
    "/media/vimukthi/Sata/Research_JLab/outfiles/pass3/intimehighest/"
    "QE_data_GEN3_sbs100p_nucleon_np_model2_sbstrackingon_intimehighest.root";
static const char* kQE =
    "/media/vimukthi/Sata/Research_JLab/outfiles/SIMC/noNRad_correct_raster/"
    "QE_sim_simc_He3eN_GEN3_sbs100p_nucleon_np_model2_elastic_simc.root";
static const char* kInelastic =
    "/media/vimukthi/Sata/Research_JLab/outfiles/inelastic_simulation_2025/"
    "QE_sim_GEN3_sbs100p_nucleon_np_model2_inelastic.root";

struct Config {
    double vzmax = 0.27;
    double eHCAL_L, EoverP_L, EoverP_H, W2_L, W2_H;
    double dy_L, dy_H, dx_L, dx_H, coin_L, coin_H, runnum_L, runnum_H;
};

std::string readText(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Could not open " + path);
    std::ostringstream s;
    s << f.rdbuf();
    return s.str();
}

double jsonNumber(const std::string& source, const std::string& name) {
    const std::regex pattern("\\\"" + name + "\\\"\\s*:\\s*([-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?)");
    std::smatch m;
    if (!std::regex_search(source, m, pattern))
        throw std::runtime_error("Missing numeric cut '" + name + "' in JSON");
    return std::stod(m[1].str());
}

Config readCuts(const std::string& path) {
    const auto j = readText(path);
    Config c;
    c.eHCAL_L = jsonNumber(j, "eHCAL_L");
    c.EoverP_L = jsonNumber(j, "EoverP_L");
    c.EoverP_H = jsonNumber(j, "EoverP_H");
    c.W2_L = jsonNumber(j, "W2_L");
    c.W2_H = jsonNumber(j, "W2_H");
    c.dy_L = jsonNumber(j, "dy_L");
    c.dy_H = jsonNumber(j, "dy_H");
    c.dx_L = jsonNumber(j, "dx_L");
    c.dx_H = jsonNumber(j, "dx_H");
    c.coin_L = jsonNumber(j, "coin_L");
    c.coin_H = jsonNumber(j, "coin_H");
    c.runnum_L = jsonNumber(j, "runnum_L");
    c.runnum_H = jsonNumber(j, "runnum_H");
    return c;
}

std::set<int> readGoodRuns(const std::string& path) {
    // Expected run_quality.json format:
    // { "1234": {"helicity":true,"moller":true}, ... }
    // Also accepts numeric 1/0 instead of true/false.
    const auto j = readText(path);
    const std::regex oneRun("\\\"([0-9]+)\\\"\\s*:\\s*\\{([^{}]*)\\}");
    const std::regex hel("\\\"helicity\\\"\\s*:\\s*(true|false|1|0)", std::regex_constants::icase);
    const std::regex mol("\\\"moller\\\"\\s*:\\s*(true|false|1|0)", std::regex_constants::icase);
    std::set<int> good;
    for (auto it = std::sregex_iterator(j.begin(), j.end(), oneRun);
         it != std::sregex_iterator(); ++it) {
        const std::string body = (*it)[2].str();
        std::smatch hm, mm;
        if (!std::regex_search(body, hm, hel) || !std::regex_search(body, mm, mol)) continue;
        auto isGood = [](const std::string& s) { return s == "1" || s == "true" || s == "TRUE" || s == "True"; };
        if (isGood(hm[1].str()) && isGood(mm[1].str())) good.insert(std::stoi((*it)[1].str()));
    }
    if (good.empty())
        throw std::runtime_error("No good runs parsed from " + path + ". Check run-quality JSON structure, or set applyRunQuality=false.");
    return good;
}

struct Input {
    std::vector<double> data, proton, neutron, inel;
    int nBins = 0;
    double xmin = 0.0, xmax = 0.0;
};

std::unique_ptr<TFile> openFile(const std::string& path) {
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) throw std::runtime_error("Cannot open ROOT file " + path);
    return f;
}

TTree* findTree(TFile* f) {
    for (const char* name : {"Tout", "T"}) {
        if (auto* t = dynamic_cast<TTree*>(f->Get(name))) return t;
    }
    throw std::runtime_error(std::string("File ") + f->GetName() + " has neither Tout nor T tree");
}

using Leaves = std::map<std::string, TLeaf*>;
Leaves enableLeaves(TTree* tree, const std::vector<std::string>& names) {
    Leaves leaves;
    tree->SetBranchStatus("*", 0);
    for (const auto& name : names) {
        TLeaf* leaf = tree->GetLeaf(name.c_str());
        if (!leaf) throw std::runtime_error(std::string(tree->GetName()) + " missing branch/leaf '" + name + "'");
        // SetBranchStatus on the leaf's parent branch name, which may differ
        // from the leaf name for older TTrees.
        tree->SetBranchStatus(leaf->GetBranch()->GetName(), 1);
        leaves[name] = leaf;
    }
    return leaves;
}

double get(const Leaves& leaves, const char* key) {
    return leaves.at(key)->GetValue(0);
}

void fillInput(Input& input, const std::string& path, const Config& c,
               const std::set<int>& goodRuns, bool useQuality,
               int which, bool matchedDy, double padding,
               Long64_t maxEntries, double preshiftP, double preshiftN, double preshiftI) {
    auto file = openFile(path);
    TTree* tree = findTree(file.get());
    std::vector<std::string> names = {"dx", "dy", "W2", "vz", "eHCAL", "ePS", "eSH", "trP"};
    if (which == 0) names.insert(names.end(), {"runnum", "ntrack", "coin_time", "helicity"});
    if (which != 0) names.push_back("weight");
    if (which == 1) names.push_back("fnucl");
    const Leaves leaves = enableLeaves(tree, names);
    const Long64_t total = tree->GetEntries();
    const Long64_t n = maxEntries > 0 ? std::min(total, maxEntries) : total;
    const double dyPad = (which == 0 || matchedDy) ? padding : 0.0;
    const double width = (input.xmax - input.xmin) / input.nBins;
    long long accepted = 0;
    std::cout << "[read] " << (which == 0 ? "data" : which == 1 ? "QE" : "inel")
              << " tree=" << tree->GetName() << " entries=" << n << " from " << path << std::endl;
    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);
        const double dx = get(leaves, "dx"), dy = get(leaves, "dy");
        const double trp = get(leaves, "trP");
        if (!(trp > 0) || !std::isfinite(dx) || !std::isfinite(dy)) continue;
        const double ep = (get(leaves, "ePS") + get(leaves, "eSH")) / trp;
        const double w2 = get(leaves, "W2");
        if (!std::isfinite(ep) || !std::isfinite(w2) ||
            std::abs(get(leaves, "vz")) > c.vzmax ||
            get(leaves, "eHCAL") < c.eHCAL_L ||
            get(leaves, "ePS") < 0.2 || ep < c.EoverP_L || ep > c.EoverP_H ||
            w2 <= c.W2_L || w2 >= c.W2_H ||
            dy <= c.dy_L - dyPad || dy >= c.dy_H + dyPad) continue;
        if (which == 0) {
            const int run = static_cast<int>(std::llround(get(leaves, "runnum")));
            if (run < c.runnum_L || run > c.runnum_H ||
                get(leaves, "ntrack") < 1 ||
                std::abs(get(leaves, "helicity")) != 1 ||
                get(leaves, "coin_time") < c.coin_L || get(leaves, "coin_time") > c.coin_H ||
                (useQuality && !goodRuns.count(run))) continue;
        }
        const double weight = which == 0 ? 1.0 : get(leaves, "weight");
        if (!std::isfinite(weight)) continue;
        std::vector<double>* dst = &input.data;
        double xval = dx;
        if (which == 1) {
            const int fnucl = static_cast<int>(std::llround(get(leaves, "fnucl")));
            if (fnucl == 1) { dst = &input.proton; xval += preshiftP; }
            else if (fnucl == 0) { dst = &input.neutron; xval += preshiftN; }
            else continue;
        } else if (which == 2) { dst = &input.inel; xval += preshiftI; }
        const int bin = static_cast<int>(std::floor((xval - input.xmin) / width));
        if (bin < 0 || bin >= input.nBins) continue;
        (*dst)[bin] += weight;
        ++accepted;
        if (i != 0 && i % 5000000 == 0) std::cout << "[read] " << i << "/" << n << std::endl;
    }
    std::cout << "[read] accepted in histogram: " << accepted << "/" << n << std::endl;
}

double total(const std::vector<double>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0);
}
void checkInput(const Input& in) {
    for (const auto& x : std::vector<std::pair<std::string, const std::vector<double>*>>{
         {"data", &in.data}, {"proton", &in.proton}, {"neutron", &in.neutron}, {"inelastic", &in.inel}}) {
        if (!std::isfinite(total(*x.second)) || total(*x.second) <= 0)
            throw std::runtime_error("Empty or invalid " + x.first + " histogram. Check cuts and paths.");
        for (double z : *x.second) if (z < -1e-9 * total(*x.second))
            throw std::runtime_error("Negative bins in " + x.first + " (weighted MC). Nonnegative model is inappropriate.");
    }
}

// Keep the FULL (guarded) MC template histogram available while evaluating
// the shifted shape ONLY in the requested fit bins. In particular, the bins
// immediately below the plotted xmin must exist before a positive dx shift.
// The shifted template is renormalized over fit bins, so A, A*B, A*C remain
// the component yields within the exact common fit window.
struct Shape {
    std::vector<double> bins, cdf;
    double sourceMin = 0.0, fitMin = 0.0, width = 0.0;
    int fitBins = 0;
    Shape(const std::vector<double>& raw, double sourceX0, double step,
          double fitX0, int nFit)
        : bins(raw), sourceMin(sourceX0), fitMin(fitX0),
          width(step), fitBins(nFit) {
        const double sum = total(bins);
        if (!(sum > 0) || fitBins <= 0 || !(width > 0))
            throw std::runtime_error("Cannot build shape from empty/invalid template");
        cdf.reserve(bins.size() + 1);
        cdf.push_back(0.0);
        for (double& z : bins) {
            z = std::max(0.0, z) / sum;
            cdf.push_back(cdf.back() + z);
        }
    }
    double integralTo(double x) const {
        if (x <= sourceMin) return 0.0;
        if (x >= sourceMin + width * bins.size()) return cdf.back();
        const double v = (x - sourceMin) / width;
        const int ix = std::min((int)bins.size()-1, std::max(0, (int)v));
        return cdf[ix] + (v - ix) * bins[ix];
    }
    std::vector<double> shifted(double delta) const {
        std::vector<double> out(fitBins, 0.0);
        for (int i = 0; i < fitBins; ++i) {
            // Positive delta means P_fit(x)=P_source(x-delta).
            const double left = fitMin + width * i - delta;
            out[i] = std::max(0.0, integralTo(left + width) - integralTo(left));
        }
        const double s = total(out);
        if (!(s > 0)) throw std::runtime_error("Shifted template has no support in fit window");
        for (double& z : out) z /= s;
        return out;
    }
};

// The nonnegative degree-n Bernstein basis columns are separately normalized to
// unit sums. Shape parameters mix these columns into one unit-normalized I_poly.
// The full model is ALWAYS A * (P + B * N + C * I), with A in events/bin
// times a unit-normalized shape and B,C dimensionless yield ratios.
double nChooseK(int n, int k) {
    double ans = 1;
    for (int j = 1; j <= k; ++j) ans *= double(n-j+1)/j;
    return ans;
}
std::vector<std::vector<double>> bernsteinBasis(int degree, int nb, double xmin, double xmax) {
    std::vector<std::vector<double>> basis(degree+1, std::vector<double>(nb));
    for (int k = 0; k <= degree; ++k) {
        for (int i = 0; i < nb; ++i) {
            const double t = (i + 0.5) / nb; // equal-width bins
            basis[k][i] = nChooseK(degree,k)*std::pow(t,k)*std::pow(1-t,degree-k);
        }
        const double s = total(basis[k]);
        for (auto& v : basis[k]) v /= s;
    }
    (void)xmin; (void)xmax;
    return basis;
}

struct Result {
    std::string name;
    bool ok = false;
    int npar = 0, ndf = 0;
    double pearson = 0, deviance = 0, reduced = 0, pvalue = 0, aic = 0;
    double nWindow = 0, bgWindow = 0, bgFraction = 0;
    std::vector<double> parameters, proton, neutron, background, model, residual, pull;
};

struct FitProblem {
    const std::vector<double>& data;
    const Shape& p, n, inel;
    int degree; // -1 = simulated inelastic shape, >=0 = Bernstein polynomial
    bool shift;
    double shiftLimit;
    std::vector<std::vector<double>> basis;
    double nData;

    FitProblem(const std::vector<double>& d, const Shape& a, const Shape& b, const Shape& c,
               int degree_, bool shift_, double limit, double xmin, double xmax)
        : data(d), p(a), n(b), inel(c), degree(degree_), shift(shift_), shiftLimit(limit),
          basis(degree_ < 0 ? std::vector<std::vector<double>>{} :
               bernsteinBasis(degree_, (int)d.size(), xmin, xmax)),
          nData(total(d)) {}

    // Parameter order:
    //  0: A [counts]   1: B = Nn/Np   2: C = NI/Np
    //  3 ... 3+(degree-1): polynomial log-shape-ratios (polynomial only)
    //  following parameters: dx_p, dx_n, and dx_i (simulation only)
    // Fix the final Bernstein weight to 1 before normalizing. This avoids
    // degeneracy between the overall polynomial normalization C and its shape.
    int nShapePars() const { return degree < 0 ? 0 : degree; }
    int shiftOffset() const { return 3 + nShapePars(); }
    int nPars() const { return shiftOffset() + (shift ? (degree < 0 ? 3 : 2) : 0); }

    std::vector<double> backgroundShape(const double* pars, double dxI) const {
        if (degree < 0) return inel.shifted(dxI);

        // I_poly(x) = sum_k w_k b_{k,degree}(x), sum_i I_poly(x_i)=1.
        // Softmax-like weights: w_k = exp(z_k)/(1+sum_j exp(z_j)),
        // with the last unnormalized weight fixed to one.
        const int nb = static_cast<int>(data.size());
        std::vector<double> out(nb, 0.0);
        std::vector<double> weights(degree+1, 1.0);
        for (int k = 0; k < degree; ++k) weights[k] = std::exp(pars[3+k]);
        const double sumW = total(weights);
        for (int k = 0; k <= degree; ++k) {
            const double w = weights[k] / sumW;
            for (int i = 0; i < nb; ++i) out[i] += w * basis[k][i];
        }
        return out;
    }

    void evaluate(const double* pars, std::vector<double>& vp,
                  std::vector<double>& vn, std::vector<double>& vb) const {
        const double A = pars[0];
        const double B = pars[1];
        const double C = pars[2];
        const int iShift = shiftOffset();
        const double dp = shift ? pars[iShift] : 0.0;
        const double dn = shift ? pars[iShift+1] : 0.0;
        const double di = (shift && degree < 0) ? pars[iShift+2] : 0.0;
        vp = p.shifted(dp);
        vn = n.shifted(dn);
        vb = backgroundShape(pars, di);

        // EXACTLY: M_i = A * [P_i + B*N_i + C*I_i], where each template
        // has unit sum over the fit range (including after shifts).
        for (size_t i = 0; i < data.size(); ++i) {
            vp[i] *= A;
            vn[i] *= A * B;
            vb[i] *= A * C;
        }
    }

    double poissonDeviance(const double* pars) const {
        std::vector<double> ap, an, ab;
        evaluate(pars, ap, an, ab);
        double out = 0.0;
        for (size_t i = 0; i < data.size(); ++i) {
            const double mu = std::max(1e-12, ap[i] + an[i] + ab[i]);
            const double obs = data[i];
            out += 2.0 * (mu - obs + (obs > 0 ? obs * std::log(obs/mu) : 0));
        }
        return out;
    }
};

Result fitOne(const std::string& label, const FitProblem& pr,
              double xmin, double xmax, double neutronLo, double neutronHi) {
    const int np = pr.nPars(), nb = (int)pr.data.size();
    const int nsh = pr.nShapePars(), iShift = pr.shiftOffset();
    Result best;
    best.name = label;
    double bestObjective = std::numeric_limits<double>::infinity();
    const std::vector<double> starts = pr.shift ? std::vector<double>{0.0, -0.025, 0.025} : std::vector<double>{0.0};
    for (double start : starts) {
        ROOT::Math::Functor fcn([&](const double* pars) { return pr.poissonDeviance(pars); }, np);
        std::unique_ptr<ROOT::Math::Minimizer> mini(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
        if (!mini) throw std::runtime_error("Minuit2 unavailable: check ROOT installation");
        mini->SetFunction(fcn);
        mini->SetMaxFunctionCalls(250000);
        mini->SetMaxIterations(10000);
        mini->SetTolerance(1e-5);
        mini->SetPrintLevel(0);
        // F(x) = A * (P(x-dp) + B*N(x-dn) + C*I(x-di)).
        // A, A*B, A*C are fitted component event yields over fit range.
        // The total expected count is A*(1+B+C).
        const double A0 = pr.nData / (1.0 + 0.20 + 0.30);
        mini->SetLimitedVariable(0, "A_proton_yield", A0,
                                 std::max(0.5, A0*0.002), 0.0, 3.0*pr.nData);
        mini->SetLimitedVariable(1, "B_N_over_P", 0.20, 0.005, 0.0, 20.0);
        mini->SetLimitedVariable(2, "C_I_over_P", 0.30, 0.005, 0.0, 20.0);
        for (int k = 0; k < nsh; ++k) {
            mini->SetLimitedVariable(3+k, "poly_log_weight_"+std::to_string(k),
                                     0.0, 0.05, -9.0, 9.0);
        }
        if (pr.shift) {
            const double s = std::min(pr.shiftLimit, std::abs(start));
            mini->SetLimitedVariable(iShift, "shift_p", std::copysign(s, start),
                                     0.001, -pr.shiftLimit, pr.shiftLimit);
            mini->SetLimitedVariable(iShift+1, "shift_n", -std::copysign(s, start),
                                     0.001, -pr.shiftLimit, pr.shiftLimit);
            if (pr.degree < 0)
                mini->SetLimitedVariable(iShift+2, "shift_inel", 0.0,
                                         0.001, -pr.shiftLimit, pr.shiftLimit);
        }
        const bool success = mini->Minimize();
        const double value = mini->MinValue();
        // Prefer converged fits; otherwise keep the smallest attained deviance.
        if (!std::isfinite(value)) continue;
        if ((success && !best.ok) || (success == best.ok && value < bestObjective)) {
            bestObjective = value;
            best.ok = success;
            best.parameters.assign(mini->X(), mini->X()+np);
        }
    }
    if (best.parameters.size() != static_cast<size_t>(np))
        throw std::runtime_error("Fitting failed for " + label + ": no finite objective result");
    best.npar = np;
    best.ndf = nb - np;
    if (best.ndf <= 0) throw std::runtime_error("Too few bins for " + label);
    pr.evaluate(best.parameters.data(), best.proton, best.neutron, best.background);
    best.model.resize(nb);
    best.pull.resize(nb);
    best.residual.resize(nb);
    const double width = (xmax - xmin) / nb;
    for (int i = 0; i < nb; ++i) {
        const double y = pr.data[i];
        const double mu = std::max(1e-12, best.proton[i] + best.neutron[i] + best.background[i]);
        best.model[i] = mu;
        best.residual[i] = y - mu;
        best.pull[i] = (y - mu)/std::sqrt(mu);
        best.pearson += (y-mu)*(y-mu)/mu;
        best.deviance += 2.0 * (mu - y + (y > 0 ? y*std::log(y/mu) : 0.0));
        const double x = xmin + (i + 0.5) * width;
        if (x >= neutronLo && x <= neutronHi) {
            best.nWindow += y;
            best.bgWindow += best.background[i];
        }
    }
    best.reduced = best.pearson / best.ndf;
    best.pvalue = TMath::Prob(best.pearson, best.ndf);
    best.aic = best.deviance + 2.0*best.npar;
    best.bgFraction = best.nWindow > 0 ? best.bgWindow/best.nWindow : std::numeric_limits<double>::quiet_NaN();
    return best;
}

std::unique_ptr<TH1D> makeHist(const std::string& name, const std::vector<double>& v,
                               double xmin, double xmax) {
    auto h = std::make_unique<TH1D>(name.c_str(), "", (int)v.size(), xmin, xmax);
    h->SetDirectory(nullptr);
    for (size_t i = 0; i < v.size(); ++i) h->SetBinContent((int)i+1, v[i]);
    return h;
}

void saveCanvas(TCanvas& canvas, const std::string& prefix) {
    canvas.SaveAs((prefix + ".pdf").c_str());
    canvas.SaveAs((prefix + ".png").c_str());
}

void drawOne(const Result& r, const std::vector<double>& data,
             double xmin, double xmax, double nlo, double nhi, const std::string& outdir,
             TFile* outputRoot) {
    const std::string tag = r.name;
    auto hd = makeHist("data_" + tag, data, xmin, xmax);
    auto hp = makeHist("proton_" + tag, r.proton, xmin, xmax);
    auto hn = makeHist("neutron_" + tag, r.neutron, xmin, xmax);
    auto hb = makeHist("background_" + tag, r.background, xmin, xmax);
    auto hm = makeHist("model_" + tag, r.model, xmin, xmax);
    auto hr = makeHist("residual_" + tag, r.residual, xmin, xmax);
    auto hpl = makeHist("pull_" + tag, r.pull, xmin, xmax);
    for (int ib = 1; ib <= hd->GetNbinsX(); ++ib)
        hd->SetBinError(ib, std::sqrt(std::max(0.0, data[ib-1])));

    outputRoot->cd();
    for (TH1D* h : {hd.get(), hp.get(), hn.get(), hb.get(), hm.get(), hr.get(), hpl.get()}) h->Write();

    TCanvas canvas(("c_" + tag).c_str(), ("dx fit " + tag).c_str(), 1150, 960);
    TPad top(("top_"+tag).c_str(), "", 0, 0.35, 1, 1);
    TPad middle(("middle_"+tag).c_str(), "", 0, 0.18, 1, 0.35);
    TPad bottom(("bottom_"+tag).c_str(), "", 0, 0, 1, 0.18);
    for (TPad* pad : {&top, &middle, &bottom}) {
        pad->SetLeftMargin(0.105); pad->SetRightMargin(0.035);
        pad->SetTicks(1,1);
        pad->Draw();
    }
    top.SetBottomMargin(0.02);
    middle.SetTopMargin(0.02); middle.SetBottomMargin(0.06);
    bottom.SetTopMargin(0.02); bottom.SetBottomMargin(0.30);

    top.cd();
    hd->SetTitle(Form("%s  |  Pearson #chi^{2}/ndf = %.1f/%d = %.3f;dx (m);Counts / bin",
                      tag.c_str(), r.pearson, r.ndf, r.reduced));
    hd->SetMarkerStyle(20); hd->SetMarkerSize(0.48); hd->SetLineColor(kBlack);
    hd->GetXaxis()->SetLabelSize(0);
    hd->GetYaxis()->SetTitleSize(0.045); hd->GetYaxis()->SetLabelSize(0.035);
    const double ymax = std::max(*std::max_element(data.begin(), data.end()),
                                 *std::max_element(r.model.begin(), r.model.end()));
    hd->SetMinimum(0); hd->SetMaximum(1.32*ymax);
    hd->Draw("E1");
    hp->SetLineColor(kMagenta+1); hp->SetLineWidth(2); hp->Draw("HIST SAME");
    hn->SetLineColor(kBlue+1); hn->SetLineWidth(2); hn->Draw("HIST SAME");
    hb->SetLineColor(kCyan+2); hb->SetLineWidth(2); hb->Draw("HIST SAME");
    hm->SetLineColor(kGreen+2); hm->SetLineWidth(3); hm->Draw("HIST SAME");
    hd->Draw("E1 SAME");
    TLegend legend(0.67, 0.60, 0.94, 0.88);
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.030);
    legend.AddEntry(hd.get(), "Data", "lep");
    legend.AddEntry(hm.get(), "Total model", "l");
    legend.AddEntry(hp.get(), "QE proton (SIMC)", "l");
    legend.AddEntry(hn.get(), "QE neutron (SIMC)", "l");
    legend.AddEntry(hb.get(), "Inelastic background", "l");
    legend.Draw();
    TLine lineL(nlo, 0, nlo, 0.95*ymax), lineH(nhi, 0, nhi, 0.95*ymax);
    for (TLine* line : {&lineL,&lineH}) {
        line->SetLineColor(kOrange+6); line->SetLineStyle(7); line->Draw();
    }

    middle.cd();
    hr->SetTitle(";dx (m);Data - Fit");
    hr->SetLineColor(kBlack); hr->SetLineWidth(2);
    hr->GetXaxis()->SetLabelSize(0);
    hr->GetYaxis()->SetTitleSize(0.12); hr->GetYaxis()->SetLabelSize(0.10);
    hr->GetYaxis()->SetTitleOffset(0.35); hr->GetYaxis()->SetNdivisions(504);
    hr->Draw("HIST");
    TLine zero1(xmin,0,xmax,0); zero1.SetLineStyle(2); zero1.Draw();

    bottom.cd();
    hpl->SetTitle(";dx (m);Pull");
    hpl->SetLineColor(kGray+2); hpl->SetMarkerStyle(20); hpl->SetMarkerSize(0.35);
    hpl->GetXaxis()->SetTitleSize(0.13); hpl->GetXaxis()->SetLabelSize(0.12);
    hpl->GetYaxis()->SetTitleSize(0.14); hpl->GetYaxis()->SetLabelSize(0.12);
    hpl->GetYaxis()->SetTitleOffset(0.3); hpl->GetYaxis()->SetNdivisions(505);
    hpl->SetMinimum(-6); hpl->SetMaximum(6);
    hpl->Draw("P");
    TLine z0(xmin,0,xmax,0), z3(xmin,3,xmax,3), zm3(xmin,-3,xmax,-3);
    for (TLine* l : {&z0,&z3,&zm3}) { l->SetLineStyle(2); l->Draw(); }
    canvas.cd();
    std::string taglower = tag;
    std::transform(taglower.begin(), taglower.end(), taglower.begin(), [](unsigned char c){return (char)std::tolower(c);});
    saveCanvas(canvas, outdir + "/fit_" + taglower);
}

void drawOverlays(const std::vector<Result>& results, const std::vector<double>& data,
                  double xmin, double xmax, const std::string& outdir) {
    std::vector<int> colors = {kGreen+2, kBlue+1, kMagenta+1, kOrange+7,
                               kCyan+2, kRed+1, kViolet+1, kGray+2};
    auto hd = makeHist("data_all", data, xmin, xmax);
    for (int i=1;i<=hd->GetNbinsX();++i) hd->SetBinError(i,std::sqrt(std::max(0.,data[i-1])));
    TCanvas cv("c_all_models", "Model overlay", 1100, 760);
    TPad top("p_all_top", "", 0, 0.26, 1, 1);
    TPad bot("p_all_bot", "", 0, 0, 1, 0.26);
    top.SetBottomMargin(0.025); top.SetLeftMargin(0.10); top.Draw();
    bot.SetTopMargin(0.015); bot.SetBottomMargin(0.26); bot.SetLeftMargin(0.10); bot.Draw();
    std::vector<std::unique_ptr<TH1D>> models, pulls;
    top.cd();
    hd->SetTitle("GEN3 dx fits: inelastic-background alternatives;dx (m);Counts / bin");
    hd->SetMaximum(1.3 * (*std::max_element(data.begin(), data.end())));
    hd->SetMarkerStyle(20); hd->SetMarkerSize(.4);
    hd->GetXaxis()->SetLabelSize(0);
    hd->Draw("E1");
    TLegend leg(0.65, 0.55, 0.94, 0.88);
    leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.029);
    leg.AddEntry(hd.get(), "Data", "lep");
    for (size_t k = 0; k < results.size(); ++k) {
        const auto& r=results[k];
        models.push_back(makeHist("overlay_"+r.name, r.model, xmin, xmax));
        models.back()->SetLineColor(colors[k%colors.size()]); models.back()->SetLineWidth(2);
        models.back()->Draw("HIST SAME");
        leg.AddEntry(models.back().get(),Form("%s: #chi^{2}_{red}=%.3f",r.name.c_str(),r.reduced),"l");
    }
    hd->Draw("E1 SAME"); leg.Draw();
    bot.cd();
    auto hAxes = makeHist("pull_axes", std::vector<double>(data.size(),0.0),xmin,xmax);
    hAxes->SetTitle(";dx (m);Pull"); hAxes->SetMinimum(-6); hAxes->SetMaximum(6);
    hAxes->GetXaxis()->SetTitleSize(.12); hAxes->GetXaxis()->SetLabelSize(.11);
    hAxes->GetYaxis()->SetTitleSize(.10); hAxes->GetYaxis()->SetLabelSize(.10);
    hAxes->GetYaxis()->SetTitleOffset(.47);
    hAxes->Draw("HIST");
    TLine z(xmin,0,xmax,0); z.SetLineStyle(2); z.Draw();
    for (size_t k=0; k<results.size(); ++k) {
        pulls.push_back(makeHist("overlay_pull_"+results[k].name, results[k].pull,xmin,xmax));
        pulls.back()->SetLineColor(colors[k%colors.size()]); pulls.back()->SetLineWidth(2);
        pulls.back()->Draw("HIST SAME");
    }
    cv.cd(); saveCanvas(cv,outdir+"/all_models_overlay");

    TCanvas cbg("c_backgrounds", "Background components", 1100,650);
    cbg.SetLeftMargin(.11); cbg.SetRightMargin(.04);
    TLegend lbg(.68,.60,.92,.88); lbg.SetBorderSize(0); lbg.SetFillStyle(0);
    std::vector<std::unique_ptr<TH1D>> backgrounds;
    double maxbg=0;
    for (auto& r:results) maxbg=std::max(maxbg,*std::max_element(r.background.begin(),r.background.end()));
    for (size_t k=0;k<results.size();++k) {
        backgrounds.push_back(makeHist("bg_overlay_"+results[k].name,results[k].background,xmin,xmax));
        auto* h=backgrounds.back().get();
        h->SetLineColor(colors[k%colors.size()]); h->SetLineWidth(2);
        h->SetTitle("Inelastic-background components;dx (m);Expected background counts / bin");
        h->SetMaximum(maxbg*1.25);
        h->Draw(k==0?"HIST":"HIST SAME");
        lbg.AddEntry(h,results[k].name.c_str(),"l");
    }
    lbg.Draw(); saveCanvas(cbg,outdir+"/inelastic_background_overlay");

    TCanvas cChi("c_chi2", "Reduced chi-square", 950,630);
    cChi.SetBottomMargin(.16);
    TH1D hChi("hReducedChi2", "Comparison of inelastic-background models;Background model;Pearson #chi^{2}/ndf",
              (int)results.size(),0,(int)results.size());
    hChi.SetStats(0); hChi.SetFillColor(kAzure-9); hChi.SetLineWidth(2);
    double ymax=0;
    for (size_t k=0;k<results.size();++k) {
        hChi.SetBinContent((int)k+1,results[k].reduced);
        hChi.GetXaxis()->SetBinLabel((int)k+1,results[k].name.c_str());
        ymax=std::max(ymax,results[k].reduced);
    }
    hChi.SetMinimum(0); hChi.SetMaximum(1.3*ymax);
    hChi.SetMarkerSize(1.5); hChi.Draw("HIST TEXT0");
    TLine one(0,1,(double)results.size(),1); one.SetLineStyle(2); one.SetLineColor(kRed+1); one.Draw();
    saveCanvas(cChi,outdir+"/reduced_chi2_comparison");
}

} // namespace dxbg

// ROOT entry point. Empty filename args use built-in GEN3 file paths.
// matchedDy=false reproduces your original 1D selections (data padded by 0.05m,
// simulation tight dy), while true applies the padded cut to both.
// fixShifts=true is useful for isolating background shape differences.
void compare_inelastic_dx_models_ABC_fullsupport(
    const char* dataFile = "",
    const char* qeFile = "",
    const char* inelFile = "",
    const char* cutsFile = "cuts/GEN3.json",
    const char* runQualityFile = "DB/run_quality.json",
    const char* outdir = "out_plots/inelastic_dx_background_models_GEN3_ABC_fullsupport",
    int bins = 200,
    bool matchedDy = false,
    bool fixShifts = false,
    const char* degreeList = "2,3,4,5",
    bool applyRunQuality = true,
    Long64_t maxEntriesPerFile = -1,
    double xmin = -2.2,
    double xmax = 0.8,
    double shiftLimit = 0.10,
    double fitxmin = 0.0,
    double fitxmax = 0.0,
    double templateGuard = 0.45) {
    using namespace dxbg;
    try {
        gROOT->SetBatch(kTRUE);
        gStyle->SetOptStat(0);
        if (bins < 16 || !(xmin < xmax)) throw std::runtime_error("Invalid bin count or dx range");
        if (!fixShifts && shiftLimit <= 0) throw std::runtime_error("shiftLimit must be > 0 unless fixShifts=true");
        const auto cuts = readCuts(cutsFile);
        const auto quality = applyRunQuality ? readGoodRuns(runQualityFile) : std::set<int>{};
        std::cout << "[config] good runs: " << quality.size() << "; matched dy: " << matchedDy
                  << "; shifts fixed: " << fixShifts << std::endl;
        gSystem->mkdir(outdir, true);
        // Fill data and templates over an EXTENDED dx support region with
        // exactly the same bin width as the displayed/fitted histogram.
        // This is important: cropping the MC before shifting forces its CDF
        // to zero outside xmin, causing the unphysical hard step in P(x-dp).
        const double step = (xmax-xmin)/bins;
        const double guardRequested = std::max(templateGuard, shiftLimit + 3.0*step);
        const int guardBins = (int)std::ceil(guardRequested/step);
        Input source;
        source.nBins = bins + 2*guardBins;
        source.xmin = xmin - guardBins*step;
        source.xmax = xmax + guardBins*step;
        source.data.assign(source.nBins,0.0);
        source.proton.assign(source.nBins,0.0);
        source.neutron.assign(source.nBins,0.0);
        source.inel.assign(source.nBins,0.0);
        std::cout << "[support] MC dx interval [" << source.xmin << ", "
                  << source.xmax << "] m, fit/display interval [" << xmin
                  << ", " << xmax << "] m; guard=" << guardBins*step
                  << " m each side; bin width=" << step << " m" << std::endl;
        fillInput(source, (dataFile && *dataFile) ? dataFile : kData, cuts, quality,
                  applyRunQuality, 0, matchedDy,0.05,maxEntriesPerFile,0,0,0);
        fillInput(source, (qeFile && *qeFile) ? qeFile : kQE, cuts, quality,
                  applyRunQuality, 1, matchedDy,0.05,maxEntriesPerFile,0,0,0);
        fillInput(source, (inelFile && *inelFile) ? inelFile : kInelastic, cuts, quality,
                  applyRunQuality, 2, matchedDy,0.05,maxEntriesPerFile,0,0,0);
        checkInput(source);
        std::cout << std::fixed << std::setprecision(2)
                  << "[hist] extended support totals: D=" << total(source.data)
                  << " P=" << total(source.proton) << " N=" << total(source.neutron)
                  << " I=" << total(source.inel) << std::endl;

        // Select data fit bins FROM the extended histogram, while leaving
        // all MC bins in source.proton/neutron/inel untouched.
        int first = guardBins, last = guardBins + bins - 1;
        if (fitxmax > fitxmin) {
            while (first <= last && source.xmin + (first+0.5)*step < fitxmin) ++first;
            while (last >= first && source.xmin + (last+0.5)*step > fitxmax) --last;
        }
        if (last-first+1 < 16)
            throw std::runtime_error("Too few bins inside requested fit interval");
        Input in;
        in.nBins = last-first+1;
        in.xmin = source.xmin + first*step;
        in.xmax = source.xmin + (last+1)*step;
        in.data.assign(source.data.begin()+first, source.data.begin()+last+1);
        if (total(in.data) <= 0)
            throw std::runtime_error("Empty data histogram inside requested fit window");
        std::cout << "[fit] dx interval [" << in.xmin << ", " << in.xmax
                  << "] m; data counts=" << total(in.data)
                  << "; bins=" << in.nBins << std::endl;

        // Use the full extended MC source as the input to shifted CDFs.
        // Resulting P, N, I are all renormalized to unit area in FIT bins.
        Shape sp(source.proton,source.xmin,step,in.xmin,in.nBins);
        Shape sn(source.neutron,source.xmin,step,in.xmin,in.nBins);
        Shape si(source.inel,source.xmin,step,in.xmin,in.nBins);
        std::set<int> degrees;
        std::stringstream degreeStream(degreeList ? degreeList : "2,3,4,5");
        std::string token;
        while (std::getline(degreeStream,token,',')) {
            if (token.empty()) continue;
            const int degree=std::stoi(token);
            if (degree<0 || degree>12) throw std::runtime_error("Polynomial degree must be 0..12");
            degrees.insert(degree);
        }
        if (degrees.empty()) throw std::runtime_error("No polynomial degrees supplied");
        std::vector<Result> results;
        results.reserve(1+degrees.size());
        std::vector<std::pair<std::string,int>> models={{"Simulation",-1}};
        for (int d:degrees) models.emplace_back("Pol"+std::to_string(d),d);
        for (const auto& model:models) {
            FitProblem problem(in.data,sp,sn,si,model.second,!fixShifts,shiftLimit,in.xmin,in.xmax);
            Result fit=fitOne(model.first,problem,in.xmin,in.xmax,cuts.dx_L,cuts.dx_H);
            std::cout << "[fit] " << std::setw(11) << fit.name
                      << " ok=" << fit.ok << " A=" << fit.parameters[0]
                      << " B=" << fit.parameters[1] << " C=" << fit.parameters[2]
                      << " chi2/ndf=" << fit.pearson << "/"
                      << fit.ndf << "=" << fit.reduced << " AIC=" << fit.aic
                      << " bg fraction in neutron dx=" << fit.bgFraction << std::endl;
            if (!fit.ok) std::cerr << "[warning] Minuit did not fully converge: " << fit.name << std::endl;
            results.push_back(std::move(fit));
        }
        TFile output((std::string(outdir)+"/fit_histograms.root").c_str(),"RECREATE");
        if (output.IsZombie()) throw std::runtime_error("Cannot create output ROOT file");
        for (auto& r:results) drawOne(r,in.data,in.xmin,in.xmax,cuts.dx_L,cuts.dx_H,outdir,&output);
        drawOverlays(results,in.data,in.xmin,in.xmax,outdir);
        output.Close();
        std::ofstream summary(std::string(outdir)+"/fit_summary.csv");
        summary << "model,converged,pearson_chi2,ndf,reduced_chi2,p_value,deviance,AIC,n_parameters,"
                   "background_full_fit,background_in_neutron_dx,data_in_neutron_dx,background_fraction_neutron_dx,"
                   "QEproton_full_fit,QEneutron_full_fit,A,B,C\n";
        summary<<std::setprecision(12);
        for (auto& r:results) summary << r.name << "," << r.ok << ","
            << r.pearson << "," << r.ndf << "," << r.reduced << "," << r.pvalue << ","
            << r.deviance << "," << r.aic << "," << r.npar << ","
            << total(r.background) << "," << r.bgWindow << "," << r.nWindow << ","
            << r.bgFraction << "," << total(r.proton) << "," << total(r.neutron)
            << "," << r.parameters[0] << "," << r.parameters[1] << "," << r.parameters[2] << "\n";
        std::ofstream pars(std::string(outdir)+"/fit_parameters.csv");
        pars << "model,parameter_index,name,value\n";
        for (const auto& r:results) {
            const int degree = r.name == "Simulation" ? -1 : std::stoi(r.name.substr(3));
            const int nsh = degree < 0 ? 0 : degree;
            for (size_t k=0; k<r.parameters.size(); ++k) {
                std::string pname;
                if (k == 0) pname = "A_proton_yield";
                else if (k == 1) pname = "B_neutron_over_proton";
                else if (k == 2) pname = "C_inelastic_over_proton";
                else if ((int)k < 3+nsh) pname = "log_bernstein_weight_" + std::to_string(k-3);
                else if ((int)k == 3+nsh) pname = "shift_proton_m";
                else if ((int)k == 4+nsh) pname = "shift_neutron_m";
                else pname = "shift_inelastic_m";
                pars << r.name << "," << k << "," << pname << ","
                     << std::setprecision(12) << r.parameters[k] << "\n";
            }
        }
        std::ofstream configText(std::string(outdir)+"/run_configuration.txt");
        configText << "Fit: F(x) = A * (P(x) + B*N(x) + C*I(x))\n"
                   << "A is the fitted proton yield (counts); B=Nn/Np; C=NI/Np.\n"
                   << "Polynomial I is a normalized Bernstein shape with degree independent log weights.\n"
                   << "Data ROOT: " << ((dataFile && *dataFile)?dataFile:kData) << "\n"
                   << "QE ROOT: " << ((qeFile && *qeFile)?qeFile:kQE) << "\n"
                   << "Inelastic ROOT: " << ((inelFile && *inelFile)?inelFile:kInelastic) << "\n"
                   << "Cuts: " << cutsFile << "\nRun quality: " << runQualityFile << "\n"
                   << "Run quality enabled: " << applyRunQuality << "\n"
                   << "Matched dy selection: " << matchedDy << "\n"
                   << "Shifts fixed: " << fixShifts << "\nShift limit (m): " << shiftLimit << "\n"
                   << "Fit x minimum: " << in.xmin << "\nFit x maximum: " << in.xmax << "\n"
                   << "Fit bins: " << in.nBins
                   << "\nTemplate support x minimum: " << source.xmin
                   << "\nTemplate support x maximum: " << source.xmax
                   << "\nTemplate guard requested (m): " << templateGuard
                   << "\nPolynomial degrees: " << degreeList << "\n"
                   << "Neutron dx window: [" << cuts.dx_L << ", " << cuts.dx_H << "]\n";
        std::cout << "[done] results in " << outdir << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << std::endl;
        throw;
    }
}
