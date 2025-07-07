// -----------------------------------------------------------------------------
// AccidentalAsymmetry.cpp  —  implementation
// -----------------------------------------------------------------------------
#include "AccidentalCorrection.h"

#include <TFile.h>
#include <TGraphErrors.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>

// ===== constructors ==========================================================

// derive from AnalysisCuts.coin_ac_L/H using fixed width
AccidentalCorrection::AccidentalCorrection(const AnalysisCuts& cuts,
                                           const RunQuality* rq,
                                           const char* kin,
                                           double binW,
                                           const char* root,
                                           const char* csv)
    : c_(cuts), kin_(kin), rq_(rq), rootName_(root), csvName_(csv)
{
    if (binW <= 0)
        throw std::runtime_error("PlotCoinTimeAsym: binWidth must be >0");
    double lo = c_.coin_ac_L-10;
    double hi = c_.coin_ac_H+10;
    if (hi <= lo)
        throw std::runtime_error("PlotCoinTimeAsym: coin_ac_L >= coin_ac_H");
    for (double x = lo; x <= hi + 1e-9; x += binW)
        edges_.push_back(x);
    if (edges_.back() < hi - 1e-9) edges_.push_back(hi);
    initBins();
}

// ===== private helpers =======================================================

void AccidentalCorrection::initBins()
{
    if (edges_.size() < 2)
        throw std::runtime_error("PlotCoinTimeAsym: need >=2 edges");
    cntP_.assign(edges_.size()-1, 0);
    cntM_.assign(edges_.size()-1, 0);
}

// ===== main loop ============================================================

void AccidentalCorrection::process(TChain& ch, BranchVars& v)
{
    Long64_t n = ch.GetEntries();
    const Long64_t step     = 100;
    double N_plus_acc = 0.0;
    double N_minus_acc = 0.0;
    double acc_events = 0.0;
    double QE_events = 0.0;

    std::cout << "[AccidentalCorrection] looping over " << n << " events\n";

    for (Long64_t i=0;i<n;++i){
        ch.GetEntry(i);

        if (rq_ && (!rq_->helicityOK(v.runnum) || !rq_->mollerOK(v.runnum))) continue;

        if (v.ntrack<1 || v.ePS<0.2 || abs(v.vz)>0.27 || v.eHCAL<c_.eHCAL_L || abs((v.ePS+v.eSH)/(v.trP)-1)>0.2||
            (c_.W2_L>v.W2 || v.W2>c_.W2_H) || (c_.dx_L>v.dx || v.dx>c_.dx_H) || (c_.dy_L>v.dy || v.dy>c_.dy_H) || 
            abs(v.helicity)!=1) continue;

        if ((c_.coin_ac_L<v.coin_time && v.coin_time<c_.coin_L-20) || (c_.coin_H+20<v.coin_time && v.coin_time<c_.coin_ac_H)){
            if (-1*v.helicity*v.IHWP*c_.Pkin_L == 1) N_plus_acc++;
            if (-1*v.helicity*v.IHWP*c_.Pkin_L == -1) N_minus_acc++;
        }

        if (c_.coin_L<v.coin_time && v.coin_time<c_.coin_H) QE_events++;
        if (c_.coin_L+35<v.coin_time && v.coin_time<c_.coin_H+35) acc_events++;

        auto it = std::upper_bound(edges_.begin(), edges_.end(), v.coin_time);
        if (it == edges_.begin() || it == edges_.end()) continue;
        size_t idx = static_cast<size_t>(it - edges_.begin() - 1);
        if (-1*v.helicity*v.IHWP*c_.Pkin_L == 1) ++cntP_[idx];
        else if (-1*v.helicity*v.IHWP*c_.Pkin_L == -1) ++cntM_[idx];

        // progress bar
        if (i % step == 0 || i == n - 1) {
            double frac = double(i + 1) / n;
            int barw = 42, pos = static_cast<int>(barw * frac);
            std::cout << '\r' << '[';
            for (int j = 0; j < barw; ++j)
                std::cout << (j < pos ? '=' : (j == pos ? '>' : ' '));
            std::cout << "] " << static_cast<int>(frac * 100) << " %" << std::flush;
        }
    }

    size_t nb = edges_.size()-1;
    std::vector<double> x(nb), y(nb), ex(nb), ey(nb);
    for (size_t i=0;i<nb;++i){
        double Np = static_cast<double>(cntP_[i]);
        double Nm = static_cast<double>(cntM_[i]);
        double sum = Np + Nm;
        double A   = (sum>0)? (Np-Nm)/sum : 0.0;
        double dA  = (sum>0)? 2.0*std::sqrt(Np*Nm)/(sum*sum*sum) : 0.0;
        x[i]  = 0.5*(edges_[i]+edges_[i+1]);
        ex[i] = 0.5*(edges_[i+1]-edges_[i]);
        y[i]  = A;
        ey[i] = dA;
    }



    TGraphErrors g(nb, x.data(), y.data(), nullptr, ey.data());
    g.SetName("gAsymVsCoinTime");
    g.SetTitle("Raw Asymmetry vs Coincidence Time; coin time; A_{raw}");

    TFile fout(Form("rootfiles/accidentals_%s.root",kin_)/*rootName_.c_str()*/, "RECREATE");
    g.Write();
    fout.Close();

    //if (!csvName_.empty()) {
        std::ofstream csv(Form("csv/accidentals_%s.csv",kin_));
        if (csv){
            csv << "#binCenter binHalfWidth A_raw dA\n";
            for (size_t i=0;i<nb;++i)
                csv << x[i] << "," << ex[i] << "," << y[i] << "," << ey[i] << "\n";
        }
    //}

    double f_acc = acc_events / QE_events;
    double A_acc = (N_plus_acc - N_minus_acc) / (N_plus_acc + N_minus_acc);

    double err_f_acc = std::sqrt(acc_events) / QE_events; // Poisson approx
    double err_A_acc = std::sqrt( (4.0 * N_plus_acc * N_minus_acc)
                                  / std::pow(N_plus_acc + N_minus_acc,3) );

    // Write accidental results to file
    
    std::ofstream txt(Form("corrections/AccidentalCorrection_%s.txt",kin_));
    txt << "N_plus_all      = " << N_plus_acc       << "\n";
    txt << "N_minus_all     = " << N_minus_acc      << "\n";
    txt << "A_acc           = " << A_acc            << "\n";
    txt << "err_A_acc       = " << err_A_acc        << "\n";
    txt << "accidental_events = " << acc_events << "\n";
    txt << "QE_events       = " << QE_events        << "\n";
    txt << "f_acc           = " << f_acc            << "\n";
    txt << "err_f_acc       = " << err_f_acc        << "\n";

    std::cout << "[AccidentalCorrection] graph written to " << rootName_ << "\n";

}

