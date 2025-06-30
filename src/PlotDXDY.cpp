#include "PlotDXDY.h"

#include <TFile.h>
#include <iostream>

PlotDXDY::PlotDXDY(const AnalysisCuts& cuts, const char* kin, const char* outRoot)
    : c_(cuts),
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

        if (v.ntrack>0 && v.ePS>0.2 && v.eHCAL>0.025){
            hvz_.Fill(v.vz);
        }

        if (v.ntrack>0 && abs(v.vz)<0.27 && v.eHCAL>0.025){
            hePS_.Fill(v.ePS);
        }

        if (v.ntrack>0 && v.ePS>0.2 && abs(v.vz)<0.27){
            heHCAL_.Fill(v.eHCAL);
        }

        if (v.ntrack>0 && v.ePS>0.2 && abs(v.vz)<0.27 && v.eHCAL>c_.eHCAL_L){
            if(c_.coin_L<v.coin_time && v.coin_time<c_.coin_H){
                hW2_.Fill(v.W2);
            }

            if((c_.W2_L<v.W2 && v.W2<c_.W2_H) && (c_.dx_L<v.dx && v.dx<c_.dx_H) && (c_.dy_L<v.dy && v.dy<c_.dy_H) ){
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

    TFile f(Form("plots_%s.root",kin_) /*outFile_.c_str()*/, "RECREATE");
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
}
