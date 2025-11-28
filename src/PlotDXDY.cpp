#include "PlotDXDY.h"

#include <TFile.h>
#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>

PlotDXDY::PlotDXDY(const AnalysisCuts& cuts, const RunQuality* rq, const char* kin, const char* outRoot)
    : c_(cuts),
      rq_(rq),
      hvz_("hvz","vz ; vz(m)",100,-0.5,0.5),
      hePS_("hePS","ePS ; ePS(GeV)",100,0.005,3),
      heHCAL_("heHCAL","eHCAL ; eHCAL(GeV)",100,0.001,3),
      hW2_("hW2","W^{2}(GeV^{2})",100,-3,5),
      hcointime_("hcointime","coincidence time; coin time(ns)",100,c_.coin_L-30,c_.coin_H+30),
      hDXDY_("hDXDY", "dy vs dx; dx (m); dy (m)", 100, -4, 3, 100, -4, 3),
      hDXW2_("hDXW2","W^{2} vs dx; dx(m); W^{2} (GeV^{2})",100,-4,3,100,-3,5),
      hDYW2_("hDYW2","W^{2} vs dy; dy(m); W^{2} (GeV^{2})",100,-4,3,100,-3,5),
      hDX_("hDX", "dx ; dx(m)", 100, -4, 3),
      hDY_("hDY", "dy ; dy(m)", 100, -4, 3),
      hDXcointime_("hDXcointime_","cointime vs dx ; dx(m); cointime(ns)", 100,-4,3,50,c_.coin_L-10,c_.coin_H+10),
      hDYcointime_("hDYcointime_","cointime vs dy ; dy(m); cointime(ns)", 100,-4,3,50,c_.coin_L-10,c_.coin_H+10),
      hePSeoverP_("hePSeoverP_","eoverp vs ePS; ePS(GeV) ; eoverp ", 100,0.0,2,100,0.6,1.4),
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

        if (v.ntrack>0 && v.ePS>0.2 && v.eHCAL>0.025){
            hvz_.Fill(v.vz);
        }

        if (v.ntrack>0 && abs(v.vz)<0.27 && v.eHCAL>c_.eHCAL_L && (c_.coin_L<v.coin_time && v.coin_time<c_.coin_H)) {
            hePS_.Fill(v.ePS);
            hePSeoverP_.Fill(v.ePS, (v.ePS+v.eSH)/v.trP);
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
        hDXDY_.Fill(v.dx, v.dy);

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


    TCanvas *c = new TCanvas("c","c",2400,1500);
    TCanvas *c1 = new TCanvas("c1","c1",2400,1500);
    TCanvas *c2 = new TCanvas("c2","c2",2400,1500); 
    TCanvas *c3 = new TCanvas("c3","c3",2400,1500); 

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    c->Divide(2,2);
    c->cd(1);
    hvz_.Draw();
    c->cd(2);
    hePS_.Draw();
    c->cd(3);
    heHCAL_.Draw();


    c1->Divide(2,2);
    c1->cd(1);
    hDXDY_.Draw("COLZ");
    c1->cd(2);
    hDXW2_.Draw("COLZ");
    c1->cd(3);
    hDYW2_.Draw("COLZ");

    c2->Divide(2,2);
    c2->cd(1);
    hcointime_.Draw();
    c2->cd(2);
    hW2_.Draw();
    c2->cd(3);
    hDX_.Draw();
    c2->cd(4);
    hDY_.Draw();

    c3->Divide(2,2);
    c3->cd(1);
    hDXcointime_.Draw("COLZ");
    c3->cd(2);
    hDYcointime_.Draw("COLZ");
    c3->cd(3);
    hePSeoverP_.Draw("COLZ");



    c->Print(Form("images/%s/plots_%s.pdf(",kin_,kin_));
    c1->Print(Form("images/%s/plots_%s.pdf",kin_,kin_));
    c2->Print(Form("images/%s/plots_%s.pdf",kin_,kin_));
    c3->Print(Form("images/%s/plots_%s.pdf)",kin_,kin_));
}
