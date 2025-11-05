// -----------------------------------------------------------------------------
// AccidentalAsymmetry.cpp  —  implementation
// -----------------------------------------------------------------------------
#include "AccidentalCorrection.h"

#include <TFile.h>
#include <TGraphErrors.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <TCanvas.h>
#include <TLine.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TSystem.h>   // for gSystem->mkdir
#include <TH2D.h>      // for TH2D
#include <TBox.h>      // for TBox (avoid incomplete type)

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

// Add this helper (e.g., above AccidentalCorrection::process)
static TGraphErrors* CalculateAsymmetry(std::vector<TH1D*>& helicityH,
                                        double lower_edge,
                                        double bin_width,
                                        const char* kin,
                                        bool flag_eHCAL_cut)
{
    const int nBins = static_cast<int>(helicityH.size());
    std::vector<double> xc(nBins), y(nBins), ey(nBins);

    gSystem->mkdir("txt", /*recursive*/true);
    std::ofstream outfile(Form("txt/%s_asymmetry_graph_cointime_binned_eHCAL_cut_%d.txt",
                               kin, int(flag_eHCAL_cut)));
    outfile << "cointime_bin_low_edge,Asymmetry_percent,Error_percent\n";

    for (int i=0;i<nBins;++i){
        TH1D* h = helicityH[i];
        const int bin_plus  = h->FindBin( 1.0);
        const int bin_minus = h->FindBin(-1.0);
        const double Np = h->GetBinContent(bin_plus);
        const double Nm = h->GetBinContent(bin_minus);
        const double S  = Np + Nm;

        double A = 0.0, sA = 0.0;
        if (S > 0.0){
            A  = (Np - Nm)/S;
            sA = std::sqrt((1-A*A)/S);//std::sqrt(4.0*Np*Nm)/std::pow(S,1.5); // correct σ_A
        }

        xc[i] = lower_edge + i*bin_width + 0.5*bin_width;
        y[i]  = 100.0*A;
        ey[i] = 100.0*sA;

        outfile << std::fixed << std::setprecision(6)
                << (lower_edge + i*bin_width) << "," << y[i] << "," << ey[i] << "\n";
    }
    outfile.close();

    auto* g = new TGraphErrors(nBins, xc.data(), y.data(), nullptr, ey.data());
    g->SetTitle("Helicity Asymmetry vs Cointime; Coincidence time (ns); Asymmetry (%)");
    g->SetMarkerStyle(21);
    g->SetMarkerColor(kBlue+1);
    g->SetLineColor(kBlue+1);
    g->SetLineWidth(3);
    return g;
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

    TH1D *h_cointime = new TH1D("h_cointime","Coincidence time (ns) ; Coincidence time (ns)" , 200,c_.coin_ac_L-20,c_.coin_ac_H+20);

    // Inside AccidentalCorrection::process(...) — add these NEW declarations after you create h_cointime:
    TH2D* h_dxdy_cut_W2_coin = new TH2D("h_dxdy_cut_W2_coin","dx-dy with W2+coin cuts;dy (m);dx (m)",200,-4,4, 200,-4,4);

    // cointime binning for asymmetry graph
    const double binMin   = c_.coin_ac_L - 20.0;
    const double binMax   = c_.coin_ac_H + 20.0;
    const double binWidth = 10.0;
    const int    nBins    = static_cast<int>((binMax - binMin)/binWidth) + 1;

    std::vector<TH1D*> Helicity_histograms;
    Helicity_histograms.reserve(nBins);
    for (int i=0;i<nBins;++i){
        const double lo = binMin + i*binWidth;
        const double hi = lo + binWidth;
        TH1D* h = new TH1D(Form("hel_bin_%d",i),
                           Form("Helicity for bin [%.1f,%.1f)",lo,hi),
                           5, -2.5, 2.5);
        Helicity_histograms.push_back(h);
    }

    // offset window edges (for boxes/lines)
    const double coin_time_offset_L = c_.coin_L + 35.0;
    const double coin_time_offset_H = c_.coin_H + 35.0;

    for (Long64_t i=0;i<n;++i){
        ch.GetEntry(i);

        if (rq_ && (!rq_->helicityOK(v.runnum) || !rq_->mollerOK(v.runnum))) continue;

        if (v.ntrack<1 || v.ePS<0.2 || abs(v.vz)>0.27 || v.eHCAL<c_.eHCAL_L || abs((v.ePS+v.eSH)/(v.trP)-1)>0.2||
            (c_.W2_L>v.W2 || v.W2>c_.W2_H) || (c_.dx_L>v.dx || v.dx>c_.dx_H) || (c_.dy_L>v.dy || v.dy>c_.dy_H) || 
            abs(v.helicity)!=1) continue;

        if ((pow((v.dy-0.0)/0.4,2)+pow((v.dx-0.0)/0.4,2))>1) continue;

        if ((c_.coin_ac_L<v.coin_time && v.coin_time<c_.coin_L-20) || (c_.coin_H+20<v.coin_time && v.coin_time<c_.coin_ac_H)){
            if (-1*v.helicity*v.IHWP*c_.Pkin_L == 1) N_plus_acc++;
            if (-1*v.helicity*v.IHWP*c_.Pkin_L == -1) N_minus_acc++;
        }

        if (c_.coin_L<v.coin_time && v.coin_time<c_.coin_H) QE_events++;
        if (c_.coin_L+35<v.coin_time && v.coin_time<c_.coin_H+35) acc_events++;

        h_cointime->Fill(v.coin_time);

        auto it = std::upper_bound(edges_.begin(), edges_.end(), v.coin_time);
        if (it == edges_.begin() || it == edges_.end()) continue;
        size_t idx = static_cast<size_t>(it - edges_.begin() - 1);
        if (-1*v.helicity*v.IHWP*c_.Pkin_L == 1) ++cntP_[idx];
        else if (-1*v.helicity*v.IHWP*c_.Pkin_L == -1) ++cntM_[idx];


        // per-bin helicity fill
        const int binIndex = static_cast<int>((v.coin_time - binMin)/binWidth);
        if (binIndex >= 0 && binIndex < nBins){
            const int sign = -1 * v.helicity * v.IHWP * c_.Pkin_L;
            Helicity_histograms[binIndex]->Fill(sign);
        }

        // dx-dy map within QE coin window
        if ( (c_.coin_L < v.coin_time) && (v.coin_time < c_.coin_H) ){
            h_dxdy_cut_W2_coin->Fill(v.dy, v.dx);
        }

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

    TFile fout(Form("rootfiles/%s/accidentals_%s.root",kin_,kin_)/*rootName_.c_str()*/, "RECREATE");
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

    double err_f_acc = std::sqrt(f_acc*(1-f_acc)/QE_events);//std::sqrt(acc_events) / QE_events; // Poisson approx
    double err_A_acc = std::sqrt((1-A_acc*A_acc)/(N_plus_acc+N_minus_acc));//std::sqrt( (4.0 * N_plus_acc * N_minus_acc)
                                 // / std::pow(N_plus_acc + N_minus_acc,3) );

    // Write accidental results to file
    
    std::ofstream txt(Form("corrections/%s/AccidentalCorrection_%s.txt",kin_,kin_));
    txt << "N_plus_all      = " << N_plus_acc       << "\n";
    txt << "N_minus_all     = " << N_minus_acc      << "\n";
    txt << "A_acc           = " << A_acc            << "\n";
    txt << "err_A_acc       = " << err_A_acc        << "\n";
    txt << "accidental_events = " << acc_events << "\n";
    txt << "QE_events       = " << QE_events        << "\n";
    txt << "f_acc           = " << f_acc            << "\n";
    txt << "err_f_acc       = " << err_f_acc        << "\n";

    std::cout << "[AccidentalCorrection] graph written to " << rootName_ << "\n";

    // Utility to draw a vertical dashed line spanning the full y-range of pad’s frame
    auto drawVLine_dashed = [](double x, double y1, double y2, int color=kGray+2) {
        TLine *L = new TLine(x, y1, x, y2);
        L->SetLineStyle(2);   // dashed
        L->SetLineWidth(3);
        L->SetLineColor(color);
        L->Draw("same");
    };

    auto drawVLine = [](double x, double y1, double y2, int color=kGray+2) {
        TLine *L = new TLine(x, y1, x, y2);
        L->SetLineStyle(1);   
        L->SetLineWidth(3);
        L->SetLineColor(color);
        L->Draw("same");
    };

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c","c",2400,1500);
    //TCanvas *c1 = new TCanvas("c1","c1",2400,1500);

    c->Divide(2,2);
    h_cointime->SetLineColor(kBlack);
    h_cointime->SetLineWidth(3);
    h_cointime->Draw();

    // Setup bins for cointime
    // const double binMin   = c_.coin_ac_L-30;
    // const double binMax   = c_.coin_ac_H+20;
    // const double binWidth = 10.0;
    // const int nBins = static_cast<int>((binMax - binMin) / binWidth) + 1;

    // span in y the visible histogram
    double yLow = 0.;                 // or -10. if you prefer matching your Box y1
    double yHigh = h_cointime->GetMaximum()*1.05;

    // interior bin edges
    for (int i=1;i<nBins;++i) {
        double xEdge = binMin + i*binWidth;
        drawVLine_dashed(xEdge, yLow, yHigh);
    }

    // optional “special” borders (QE and accidental window limits)
    drawVLine(c_.coin_L,      yLow, yHigh, kGreen);
    drawVLine(c_.coin_H,      yLow, yHigh, kGreen);
    drawVLine(c_.coin_L+35, yLow, yHigh, kRed);
    drawVLine(c_.coin_H+35, yLow, yHigh, kRed);


    //c1->Divide(2,2);

    c->Print(Form("images/%s/accidentals_plots_%s.pdf",kin_,kin_));
    c->Print(Form("images/%s/accidentals_plots_%s.png",kin_,kin_));

    // After writing root/csv/txt outputs — add this NEW canvas block:

    gStyle->SetOptStat(0);

    // helpers
    // auto drawVLine_dashed = [](double x, double y1, double y2, int color=kGray+2) {
    //     TLine *L = new TLine(x, y1, x, y2);
    //     L->SetLineStyle(2);
    //     L->SetLineWidth(2);
    //     L->SetLineColor(color);
    //     L->Draw("same");
    // };
    // auto drawVLine = [](double x, double y1, double y2, int color=kGray+2) {
    //     TLine *L = new TLine(x, y1, x, y2);
    //     L->SetLineStyle(1);
    //     L->SetLineWidth(3);
    //     L->SetLineColor(color);
    //     L->Draw("same");
    // };

    // boxes
    const double y1 = -10.0;
    const double y2 = h_cointime->GetMaximum()*1.05;

    TBox* box_offset = new TBox(coin_time_offset_L, y1, coin_time_offset_H, y2);
    box_offset->SetFillColorAlpha(6, 0.30);
    box_offset->SetLineColor(6);

    TBox* box_anti1 = new TBox(c_.coin_ac_L, y1, c_.coin_L-10, y2);
    box_anti1->SetFillColorAlpha(kRed, 0.30);
    box_anti1->SetLineColor(kRed);

    TBox* box_anti2 = new TBox(c_.coin_H+10, y1, c_.coin_ac_H, y2);
    box_anti2->SetFillColorAlpha(kRed, 0.30);
    box_anti2->SetLineColor(kRed);

    TBox* box_dxdy = new TBox(c_.dy_L, c_.dx_L, c_.dy_H, c_.dx_H);
    box_dxdy->SetFillStyle(0);
    box_dxdy->SetLineColor(kRed+1);
    box_dxdy->SetLineWidth(3);

    // canvas
    TCanvas *ccoin = new TCanvas("ccoin","ccoin", 3600, 3000);
    TCanvas *ccoin_1 = new TCanvas("ccoin_1","ccoin_1", 3600, 3000);
    ccoin->Divide(2,2);
    ccoin_1->Divide(1,1);
    // pad 1
    ccoin->cd(1);
    h_cointime->SetLineColor(kBlack);
    h_cointime->SetLineWidth(3);
    h_cointime->GetXaxis()->SetTitle("coincidence time (ns)");
    h_cointime->Draw("hist");
    box_offset->Draw("same");
    for (int i=1;i<nBins;++i) drawVLine_dashed(binMin + i*binWidth, 0.0, y2);
    drawVLine(c_.coin_L,          0.0, y2, kGreen+2);
    drawVLine(c_.coin_H,          0.0, y2, kGreen+2);
    drawVLine(coin_time_offset_L, 0.0, y2, kAzure+2);
    drawVLine(coin_time_offset_H, 0.0, y2, kAzure+2);

    // pad 2 with sub-pads
    ccoin->cd(2);
    ccoin_1->cd(1);
    TPad* p2_top = new TPad("p2_top","p2_top", 0.0, 0.30, 1.0, 1.0);
    p2_top->SetBottomMargin(0.02);
    p2_top->Draw();
    p2_top->cd();

    h_cointime->GetXaxis()->SetTitle("");
    h_cointime->Draw("hist");
    box_anti1->Draw("same");
    box_anti2->Draw("same");
    box_offset->Draw("same");
    for (int i=1;i<nBins;++i) drawVLine_dashed(binMin + i*binWidth, 0.0, y2);
    drawVLine(c_.coin_L,          0.0, y2, kGreen+1);
    drawVLine(c_.coin_H,          0.0, y2, kGreen+1);
    drawVLine(coin_time_offset_L, 0.0, y2, kRed+2);
    drawVLine(coin_time_offset_H, 0.0, y2, kRed+2);
    p2_top->Update();

    ccoin->cd(2);
    ccoin_1->cd(1);
    TPad* p2_bottom = new TPad("p2_bottom","p2_bottom", 0.0, 0.0, 1.0, 0.30);
    p2_bottom->SetTopMargin(0.02);
    p2_bottom->SetBottomMargin(0.25);
    p2_bottom->Draw();
    p2_bottom->cd();

    TGraphErrors* gAsymPct = CalculateAsymmetry(Helicity_histograms,
                                                binMin, binWidth, kin_, /*flag_eHCAL_cut=*/true);
    gAsymPct->Draw("AP");
    gAsymPct->GetXaxis()->SetRangeUser(binMin, binMax);
    gAsymPct->SetMinimum(-25);
    gAsymPct->SetMaximum( 25);
    gAsymPct->GetXaxis()->SetTitleSize(0.07);
    gAsymPct->GetXaxis()->SetLabelSize(0.06);
    gAsymPct->GetXaxis()->SetTitleOffset(0.95);
    gAsymPct->GetYaxis()->SetTitleSize(0.07);
    gAsymPct->GetYaxis()->SetLabelSize(0.06);
    gAsymPct->GetYaxis()->SetTitleOffset(0.8);
    gAsymPct->GetXaxis()->SetTitle("coincidence time (ns)");
    p2_bottom->Update();

    // pad 3
    ccoin->cd(3);
    h_dxdy_cut_W2_coin->SetContour(50);
    h_dxdy_cut_W2_coin->Draw("COLZ");
    box_dxdy->Draw("same");

    // save
    gSystem->mkdir("plots", true);
    //ccoin->SaveAs(Form("plots/%s_cointime_plots_for_accidentals_eHCAL_cut_%d.pdf", kin_, 1));
    ccoin->SaveAs(Form("images/%s/asymmetry_vs_cointime_%s_1.png", kin_, kin_));
    ccoin_1->SaveAs(Form("images/%s/asymmetry_vs_cointime_%s.png", kin_, kin_));
    //ccoin->SaveAs(Form("plots/%s_cointime_plots_for_accidentals_eHCAL_cut_%d.jpg", kin_, 1));

}

