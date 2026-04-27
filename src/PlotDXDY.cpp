#include "PlotDXDY.h"

#include <TFile.h>
#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TF1.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TMath.h>

PlotDXDY::PlotDXDY(const AnalysisCuts& cuts, const RunQuality* rq, const char* kin, const char* outRoot)
    : c_(cuts),
      rq_(rq),
      hvz_("hvz","Vertex z ; vz(m)",100,-0.5,0.5),
      hePS_("hePS","Pre-shower energy ; Pre-shower energy(GeV/c)",100,0.005,3),
      heoverP_("heoverP","E/p BBCal ; E/p",100,0.5,1.5),
      heHCAL_("heHCAL","HCal energy; HCal energy(GeV/c)",100,0.001,1),
      hW2_("hW2","W^{2};W^{2}((GeV/c)^{2})",100,-3,5),
      hcointime_("hcointime","Coincidence time; Coincidence time(ns)",100,c_.coin_L-30,c_.coin_H+30),
      hDXDY_("hDXDY", "dx vs dy; dy (m); dx (m)", 100, -4, 3, 100, -4, 3),
      hDXW2_("hDXW2","W^{2} vs dx; dx(m); W^{2} ((GeV/c)^{2})",100,-4,3,100,-3,5),
      hDYW2_("hDYW2","W^{2} vs dy; dy(m); W^{2} ((GeV/c)^{2})",100,-4,3,100,-3,5),
      hDX_("hDX", "dx ; dx(m)", 100, -4, 3),
      hDY_("hDY", "dy ; dy(m)", 100, -4, 3),
      hDXcointime_("hDXcointime_","cointime vs dx ; dx(m); cointime(ns)", 100,-4,3,50,c_.coin_L-10,c_.coin_H+10),
      hDYcointime_("hDYcointime_","cointime vs dy ; dy(m); cointime(ns)", 100,-4,3,50,c_.coin_L-10,c_.coin_H+10),
      hePSeoverP_("hePSeoverP_","E/p vs Pre-shower energy; Pre-shower energy(GeV/c) ; E/p ", 100,0.0,2,100,0.6,1.4),
      kin_(kin),
      outFile_(outRoot) {}

bool PlotDXDY::passes(const BranchVars& v) const
{
    // dx window
    /*if (c_.dx_L != 0 || c_.dx_H != 0) {
        if (v.dx < c_.dx_L || v.dx > c_.dx_H) return false;
    }
    // dy window
    if (c_.dy_L != 0 || c_.dy_H != 0) {
        if (v.dy < c_.dy_L || v.dy > c_.dy_H) return false;
    }
    // helicity selection (999 → accept ±1)
    if (c_.helicity != 999 && v.helicity != c_.helicity) return false;
    */

    if (v.ntrack<1) return false;
    if (abs(v.vz)>0.27) return false;
    if (v.ePS<0.2) return false;
    if (v.eHCAL<0.025) return false;
    if (abs((v.ePS+v.eSH)/v.trP -1)>0.2) return false;

    if (v.W2<c_.W2_L || v.W2>c_.W2_H) return false;
    if (v.coin_time<c_.coin_L || v.coin_time>c_.coin_H) return false;
    if (v.eHCAL<c_.eHCAL_L) return false;

    return true;
}

void PlotDXDY::process(TChain& ch, BranchVars& v)
{
    Long64_t n = ch.GetEntries();
    const Long64_t step     = 100;//std::max<Long64_t>(1, nEntries / 50);
    std::cout << "[PlotDXDY] looping over " << n << " events\n";
    for (Long64_t i = 0; i < n; ++i) {
        ch.GetEntry(i);

        if (rq_ && (!rq_->helicityOK(v.runnum) || !rq_->mollerOK(v.runnum))) continue;

        //if ((v.ntrack_sbs>0) && abs(v.vz_sbs)<0.27) continue; //sbs veto

        if (v.ntrack>0 && v.ePS>0.2 && v.eHCAL>0.025){
            hvz_.Fill(v.vz);
        }

        if (v.ntrack>0 && abs(v.vz)<0.27 && v.eHCAL>c_.eHCAL_L && (c_.coin_L<v.coin_time && v.coin_time<c_.coin_H)) {
            hePS_.Fill(v.ePS);
            hePSeoverP_.Fill(v.ePS, (v.ePS+v.eSH)/v.trP);
            heoverP_.Fill((v.ePS+v.eSH)/v.trP);

        //     if((c_.W2_L<v.W2 && v.W2<c_.W2_H) && (c_.dx_L<v.dx && v.dx<c_.dx_H) && (c_.dy_L<v.dy && v.dy<c_.dy_H) && ((pow((v.dy-c_.dy_c)/c_.dy_r,2)+pow((v.dx-c_.dx_c)/c_.dx_r,2))<1) 
        // && (c_.coin_L<v.coin_time && v.coin_time<c_.coin_H)){
        //         hePSeoverP_.Fill(v.ePS, (v.ePS+v.eSH)/v.trP);
        //         //heoverP_.Fill((v.ePS+v.eSH)/v.trP);
        //     }

        //     if((c_.W2_L<v.W2 && v.W2<c_.W2_H) && (c_.dx_L<v.dx && v.dx<c_.dx_H) && (c_.dy_L<v.dy && v.dy<c_.dy_H) && (((pow((v.dy-c_.dy_c)/c_.dy_r,2)+pow((v.dx-c_.dx_c)/c_.dx_r,2))<1)||(((pow((v.dy-c_.dy_P_c)/c_.dy_P_r,2)+pow((v.dx-c_.dx_P_c)/c_.dx_P_r,2))<1))) 
        // && (c_.coin_L<v.coin_time && v.coin_time<c_.coin_H)){
        //         //hePSeoverP_.Fill(v.ePS, (v.ePS+v.eSH)/v.trP);
        //         heoverP_.Fill((v.ePS+v.eSH)/v.trP);
        //     }
        }

        if (v.ntrack>0 && v.ePS>0.2 && abs(v.vz)<0.27){
            heHCAL_.Fill(v.eHCAL);
        }

        if (v.ntrack>0 && v.ePS>0.2 && abs(v.vz)<0.27 && v.eHCAL>c_.eHCAL_L && abs((v.ePS+v.eSH)/v.trP - 1)<0.2){
            if(c_.coin_L<v.coin_time && v.coin_time<c_.coin_H){
                hW2_.Fill(v.W2);
            }

            if((c_.W2_L<v.W2 && v.W2<c_.W2_H) && (c_.dx_L<v.dx && v.dx<c_.dx_H) && (c_.dy_L<v.dy && v.dy<c_.dy_H) && ((pow((v.dy-c_.dy_c)/c_.dy_r,2)+pow((v.dx-c_.dx_c)/c_.dx_r,2))<1)){
                hcointime_.Fill(v.coin_time);
            }

            if( (c_.coin_L<v.coin_time && v.coin_time<c_.coin_H) && (c_.dy_L<v.dy && v.dy<c_.dy_H) ){
                hDXW2_.Fill(v.dx,v.W2);
            }

            if( (c_.coin_L<v.coin_time && v.coin_time<c_.coin_H) && (c_.dx_L<v.dx && v.dx<c_.dx_H) ){
                hDYW2_.Fill(v.dy,v.W2);
            }
            if( (c_.W2_L<v.W2 && v.W2<c_.W2_H) && (c_.coin_L<v.coin_time && v.coin_time<c_.coin_H) && (c_.dx_L<v.dx && v.dx<c_.dx_H) ){
                hDY_.Fill(v.dy);
            }
            if((c_.W2_L<v.W2 && v.W2<c_.W2_H) && (c_.dy_L<v.dy && v.dy<c_.dy_H)){
                hDXcointime_.Fill(v.dx,v.coin_time);
            }
            if((c_.W2_L<v.W2 && v.W2<c_.W2_H) && (c_.dx_L<v.dx && v.dx<c_.dx_H)){
                hDYcointime_.Fill(v.dy,v.coin_time);
            }
        }

        if (!passes(v)) continue;
        hDXDY_.Fill(v.dy, v.dx);

        if(!passes(v) || v.dy<c_.dy_L || v.dy>c_.dy_H) continue;
        hDX_.Fill(v.dx);


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


    TFile f(Form("rootfiles/%s/plots_%s.root",kin_,kin_) /*outFile_.c_str()*/, "RECREATE");
    hvz_.Write();
    hePS_.Write();
    heHCAL_.Write();
    hW2_.Write();
    hcointime_.Write();
    hDXDY_.Write();
    hDXW2_.Write();
    hDYW2_.Write();
    hDX_.Write();
    hDY_.Write();
    f.Close();
    std::cout << "[PlotDXDY] histogram written to " << outFile_ << "\n";

    // One-plot-per-canvas version
    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    // gSystem->mkdir(Form("images/%s", kin_), true);

    // -------------------- basic styling --------------------
    auto Style1D = [](TH1* h){
    h->SetLineColor(kBlack);
    h->SetLineWidth(8);
    h->SetFillColorAlpha(kViolet+2, 0.3);
    };

    Style1D(&hvz_);
    Style1D(&hePS_);
    Style1D(&heoverP_);
    Style1D(&heHCAL_);
    Style1D(&hW2_);
    Style1D(&hcointime_);
    Style1D(&hDX_);
    Style1D(&hDY_);

    // Helper for 1D plots with cut lines
    // auto Save1DWithLines = [&](TH1& h, const char* cname, const char* outname,
    //                         std::vector<double> xlines)
    // {
    // TCanvas *cplot = new TCanvas(cname, cname, 3200, 2400);
    // cplot->cd();
    // h.Draw();

    // std::vector<TLine*> lines;
    // for(double x : xlines){
    //     TLine *ln = new TLine(x, h.GetMinimum(), x, h.GetMaximum());
    //     ln->SetLineColor(kRed);
    //     ln->SetLineWidth(8);
    //     ln->Draw("same");
    //     lines.push_back(ln);
    // }

    // cplot->SaveAs(Form("images/%s/%s.png", kin_, outname));
    // delete cplot;
    // };

    auto Save1DWithLines = [&](TH1& h, const char* cname, const char* outname,
                            std::vector<double> xlines,
                            bool doGausFit = false,
                            double fitMin = 0.0,
                            double fitMax = 0.0)
    {
        TCanvas *cplot = new TCanvas(cname, cname, 3200, 2400);
        cplot->cd();

        h.Draw();

        // draw cut lines
        std::vector<TLine*> lines;
        for(double x : xlines){
            TLine *ln = new TLine(x, 0.0, x, h.GetMaximum());
            ln->SetLineColor(kRed);
            ln->SetLineWidth(8);
            ln->Draw("same");
            lines.push_back(ln);
        }

        // optional Gaussian fit
        if(doGausFit && h.GetEntries() > 10){
            // if no range provided, use mean ± 2*RMS as a reasonable default
            if(fitMin >= fitMax){
                fitMin = h.GetMean() - 2.0*h.GetRMS();
                fitMax = h.GetMean() + 2.0*h.GetRMS();
            }

            TF1 *fg = new TF1(Form("fg_%s", cname), "gaus", fitMin, fitMax);
            fg->SetLineColor(kBlue+1);
            fg->SetLineWidth(6);

            // sensible initial parameters
            fg->SetParameters(h.GetMaximum(), h.GetMean(), h.GetRMS());

            h.Fit(fg, "RQ");  // R = use range, Q = quiet

            double mean    = fg->GetParameter(1);
            double sigma   = fg->GetParameter(2);
            double emean   = fg->GetParError(1);
            double esigma  = fg->GetParError(2);

            fg->Draw("same");

            TPaveText *pt = new TPaveText(0.60, 0.72, 0.88, 0.88, "NDC");
            pt->SetFillColor(0);
            pt->SetBorderSize(1);
            pt->SetTextAlign(12);
            pt->SetTextSize(0.035);
            pt->AddText(Form("Mean  = %.4f #pm %.4f", mean, emean));
            pt->AddText(Form("#sigma = %.4f #pm %.4f", sigma, esigma));
            pt->Draw("same");
        }

        cplot->SaveAs(Form("images/%s/%s.png", kin_, outname));
        delete cplot;
    };

    // Helper for 2D plots
    auto Save2D = [&](TH2& h, const char* cname, const char* outname, bool logz=false)
    {
    TCanvas *cplot = new TCanvas(cname, cname, 3200, 2400);
    cplot->cd();
    if(logz) gPad->SetLogz(1);
    h.Draw("COLZ");
    cplot->SaveAs(Form("images/%s/%s.png", kin_, outname));
    delete cplot;
    };

    // Helper for 2D dx-dy with ellipse
    auto SaveDXDYWithEllipse = [&](TH2& h, const char* cname, const char* outname)
    {
    TCanvas *cplot = new TCanvas(cname, cname, 3200, 2400);
    cplot->cd();
    gPad->SetLogz(1);
    h.Draw("COLZ");

    TEllipse *ell_dxdy = new TEllipse(c_.dy_c, c_.dx_c, c_.dy_r, c_.dx_r);
    ell_dxdy->SetFillStyle(0);
    ell_dxdy->SetLineColor(kBlack);
    ell_dxdy->SetLineWidth(3);
    ell_dxdy->SetLineStyle(2);
    ell_dxdy->Draw("same");

    cplot->SaveAs(Form("images/%s/%s.png", kin_, outname));
    delete cplot;
    };

    // -------------------- 1D plots --------------------
    Save1DWithLines(hvz_,       "cvz",       "vertex_z",         {-0.27, 0.27});
    Save1DWithLines(hePS_,      "cps",       "preshower_energy", {c_.pion_L});
    Save1DWithLines(heoverP_, "ceoverp", "eoverp_bbcal",{c_.EoverP_L, c_.EoverP_H});
    
    Save1DWithLines(heHCAL_,    "chcal",     "hcal_energy",      {c_.eHCAL_L});
    //Save1DWithLines(heoverP_,   "ceoverp",   "eoverp_bbcal",     {c_.EoverP_L, c_.EoverP_H});

    Save1DWithLines(hcointime_, "ccointime", "coincidence_time", {c_.coin_L, c_.coin_H});
    Save1DWithLines(hW2_,       "cw2",       "W2",               {c_.W2_L, c_.W2_H});
    Save1DWithLines(hDX_,       "cdx",       "dx",               {c_.dx_L, c_.dx_H});
    Save1DWithLines(hDY_,       "cdy",       "dy",               {c_.dy_L, c_.dy_H});

    // -------------------- 2D plots --------------------
    SaveDXDYWithEllipse(hDXDY_,     "cdxdy",        "dx_vs_dy");
    Save2D(hDXW2_,                  "cdxw2",        "W2_vs_dx", true);
    Save2D(hDYW2_,                  "cdyw2",        "W2_vs_dy", true);

    Save2D(hDXcointime_,            "cdxcointime",  "cointime_vs_dx", true);
    Save2D(hDYcointime_,            "cdycointime",  "cointime_vs_dy", true);
    Save2D(hePSeoverP_,             "ceps_eop",     "eoverp_vs_preshower", true);



    // TCanvas *c = new TCanvas("c","c",2400,1500);
    // TCanvas *c1 = new TCanvas("c1","c1",2400,1500);
    // TCanvas *c2 = new TCanvas("c2","c2",2400,1500); 
    // TCanvas *c3 = new TCanvas("c3","c3",2400,1500); 

    // gStyle->SetPalette(kRainBow);
    // gStyle->SetOptStat(0);

    // hvz_.SetLineColor(kBlack);
    // hePS_.SetLineColor(kBlack);
    // heoverP_.SetLineColor(kBlack);
    // heHCAL_.SetLineColor(kBlack);
    // hW2_.SetLineColor(kBlack);
    // hcointime_.SetLineColor(kBlack);
    // hDX_.SetLineColor(kBlack);
    // hDY_.SetLineColor(kBlack);


    // hvz_.SetLineWidth(3);
    // hePS_.SetLineWidth(3);
    // heoverP_.SetLineWidth(3);
    // heHCAL_.SetLineWidth(3);
    // hW2_.SetLineWidth(3);
    // hcointime_.SetLineWidth(3);
    // hDX_.SetLineWidth(3);
    // hDY_.SetLineWidth(3);

    // hvz_.SetFillColor(kViolet-9);
    // hePS_.SetFillColor(kViolet-9);
    // heoverP_.SetFillColor(kViolet-9);
    // heHCAL_.SetFillColor(kViolet-9);
    // hW2_.SetFillColor(kViolet-9);
    // hcointime_.SetFillColor(kViolet-9);
    // hDX_.SetFillColor(kViolet-9);
    // hDY_.SetFillColor(kViolet-9);

    // c->Divide(2,2);
    // c->cd(1);
    // TLine *l_vz = new TLine(-0.27, hvz_.GetMinimum(), -0.27, hvz_.GetMaximum());
    // TLine *h_vz = new TLine(0.27, hvz_.GetMinimum(), 0.27, hvz_.GetMaximum());
    // l_vz->SetLineColor(kRed);
    // h_vz->SetLineColor(kRed);
    // l_vz->SetLineWidth(3);
    // h_vz->SetLineWidth(3);
    // hvz_.Draw();
    // l_vz->Draw("same");
    // h_vz->Draw("same");
    
    // c->cd(2);
    // TLine *l_ps = new TLine(c_.pion_L, hePS_.GetMinimum(), c_.pion_L, hePS_.GetMaximum());
    // //TLine *h_ps = new TLine(0.2, hePS_.GetMinimum(), c_.ps_H, hePS_.GetMaximum());
    // l_ps->SetLineColor(kRed);
    // //h_ps->SetLineColor(kRed);
    // l_ps->SetLineWidth(3);
    // //h_ps->SetLineWidth(3);
    // hePS_.Draw();
    // l_ps->Draw("same");
    // //h_ps->Draw("same");

    // c->cd(3);
    // TLine *eHCAL_L = new TLine(c_.eHCAL_L, heHCAL_.GetMinimum(), c_.eHCAL_L, heHCAL_.GetMaximum());
    // eHCAL_L->SetLineColor(kRed);
    // eHCAL_L->SetLineWidth(3);
    // heHCAL_.Draw();
    // eHCAL_L->Draw("same");


    // c->cd(4);
    // TLine *eoverp_L = new TLine(c_.EoverP_L, 0, c_.EoverP_L, heoverP_.GetMaximum());
    // TLine *eoverp_H = new TLine(c_.EoverP_H, 0, c_.EoverP_H, heoverP_.GetMaximum());
    // eoverp_L->SetLineColor(kRed);
    // eoverp_H->SetLineColor(kRed);
    // eoverp_L->SetLineWidth(3);
    // eoverp_H->SetLineWidth(3);
    // heoverP_.Draw();
    // eoverp_L->Draw("same");
    // eoverp_H->Draw("same");



    // c2->Divide(2,2);
    // c2->cd(1);
    // TLine *coin_time_L = new TLine(c_.coin_L, hcointime_.GetMinimum(), c_.coin_L, hcointime_.GetMaximum());
    // TLine *coin_time_H = new TLine(c_.coin_H, hcointime_.GetMinimum(), c_.coin_H, hcointime_.GetMaximum());
    // coin_time_L->SetLineColor(kRed);
    // coin_time_H->SetLineColor(kRed);
    // coin_time_L->SetLineWidth(3);
    // coin_time_H->SetLineWidth(3);
    // hcointime_.Draw();
    // coin_time_L->Draw("same");
    // coin_time_H->Draw("same");

    // c2->cd(2);
    // TLine *W2_L = new TLine(c_.W2_L, hW2_.GetMinimum(), c_.W2_L, hW2_.GetMaximum());
    // TLine *W2_H = new TLine(c_.W2_H, hW2_.GetMinimum(), c_.W2_H, hW2_.GetMaximum());
    // W2_L->SetLineColor(kRed);
    // W2_H->SetLineColor(kRed);
    // W2_L->SetLineWidth(3);
    // W2_H->SetLineWidth(3);
    // hW2_.Draw();
    // W2_L->Draw("same");
    // W2_H->Draw("same");
    
    // c2->cd(3);
    // TLine *dx_L = new TLine(c_.dx_L, hDX_.GetMinimum(), c_.dx_L, hDX_.GetMaximum());
    // TLine *dx_H = new TLine(c_.dx_H, hDX_.GetMinimum(), c_.dx_H, hDX_.GetMaximum());
    // dx_L->SetLineColor(kRed);
    // dx_H->SetLineColor(kRed);
    // dx_L->SetLineWidth(3);
    // dx_H->SetLineWidth(3);
    // hDX_.Draw();
    // dx_L->Draw("same");
    // dx_H->Draw("same");

    // c2->cd(4);
    // TLine *dy_L = new TLine(c_.dy_L, hDY_.GetMinimum(), c_.dy_L, hDY_.GetMaximum());
    // TLine *dy_H = new TLine(c_.dy_H, hDY_.GetMinimum(), c_.dy_H, hDY_.GetMaximum());
    // dy_L->SetLineColor(kRed);
    // dy_H->SetLineColor(kRed);
    // dy_L->SetLineWidth(3);
    // dy_H->SetLineWidth(3);
    // hDY_.Draw();
    // dy_L->Draw("same");
    // dy_H->Draw("same");

    // // c1->SetLogz(1);

    // c1->Divide(2,2);
    // c1->cd(1);
    // gPad->SetLogz(1);
    // hDXDY_.Draw("COLZ");
    // TEllipse *ell_dxdy = new TEllipse(c_.dy_c, c_.dx_c, c_.dy_r, c_.dx_r);
    // ell_dxdy->SetFillStyle(0);     // transparent
    // ell_dxdy->SetLineColor(kBlack);
    // ell_dxdy->SetLineWidth(3);
    // ell_dxdy->SetLineStyle(2);     // dashed, optional
    // ell_dxdy->Draw("same");

    // c1->cd(2);
    // gPad->SetLogz(1);
    // hDXW2_.Draw("COLZ");
    // c1->cd(3);
    // gPad->SetLogz(1);
    // hDYW2_.Draw("COLZ");

    // // c3->SetLogz(1);

    // c3->Divide(2,2);
    // c3->cd(1);
    // gPad->SetLogz(1);
    // hDXcointime_.Draw("COLZ");
    // c3->cd(2);
    // gPad->SetLogz(1);
    // hDYcointime_.Draw("COLZ");
    // c3->cd(3);
    // gPad->SetLogz(1);
    // hePSeoverP_.Draw("COLZ");



    // c->Print(Form("images/%s/plots_%s.pdf(",kin_,kin_));
    // c1->Print(Form("images/%s/plots_%s.pdf",kin_,kin_));
    // c2->Print(Form("images/%s/plots_%s.pdf",kin_,kin_));
    // c3->Print(Form("images/%s/plots_%s.pdf)",kin_,kin_));
    
    // c->Print(Form("images/%s/plots_%s_1.png",kin_,kin_));
    // c1->Print(Form("images/%s/plots_%s_2.png",kin_,kin_));
    // c2->Print(Form("images/%s/plots_%s_3.png",kin_,kin_));
    // c3->Print(Form("images/%s/plots_%s_4.png",kin_,kin_));
}
