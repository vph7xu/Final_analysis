#include "RawAsymmetry.h"

#include <TFile.h>
#include <TH1D.h> 
#include <TCanvas.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <unordered_map>

// ---------- helpers ----------------------------------------------------------
static TDatime parseDT(const std::string& s)
{
    int y,m,d,h,min,sec;
    std::sscanf(s.c_str(), "%d-%d-%d %d:%d:%d", &y,&m,&d,&h,&min,&sec);
    return TDatime(y,m,d,h,min,sec);
}

// Identify beam interval and return its index for block grouping
struct BeamHit {
    int idx;       // CSV row index (block id)
    double val;    // polarization value
    double err;    // ABS systematic error (same units as val)
};

static BeamHit beamAtHit(const std::vector<BeamRow>& tbl, const TDatime& t)
{
    for (size_t i = 0; i < tbl.size(); ++i) {
        const auto& r = tbl[i];
        if (t >= r.start && t <= r.end) return BeamHit{ (int)i, r.val, r.err };
    }
    return BeamHit{ -1, -1.0, -1.0 };
}

// Make a stable key to group identical target (value, error) blocks
static inline long long makeKey(double val, double err)
{
    long long v = llround(val * 1e4); // round to 1e-4
    long long e = llround(err * 1e4);
    return (v << 32) ^ (e & 0xffffffffLL);
}

// ---------- constructor ------------------------------------------------------
RawAsymmetry::RawAsymmetry(const CutManager& cuts, const AnalysisCuts& acuts,
                           const RunQuality* rq, const char* kin,
                           int nb,double xlo,double xhi,
                           const char* root,const char* txt)
    : cuts_(cuts), c_(acuts), rq_(rq), kin_(kin),
      h_("hA","A_{exp}; A; events",nb,xlo,xhi),
      rootF_(root), txtF_(txt) {}

// ---------- CSV I/O ----------------------------------------------------------
std::vector<BeamRow> RawAsymmetry::readBeamCSV(const std::string& csv) const
{
    std::vector<BeamRow> v; std::ifstream f(csv);
    if(!f){ std::cerr<<"[RawAsym] cannot open "<<csv<<"\n"; return v; }
    std::string line; std::getline(f,line); // skip header
    while(std::getline(f,line)){
        std::stringstream ss(line); std::string a,b; double val,err;
        std::getline(ss,a,','); std::getline(ss,b,','); ss>>val; ss.ignore(1,','); ss>>err;
        v.push_back({parseDT(a),parseDT(b),val,err});
    }
    return v;
}

// Retained for compatibility (not used in the new accumulation)
std::pair<double,double> RawAsymmetry::beamAt(const std::vector<BeamRow>& tbl,
                                              const TDatime& t) const
{
    for(const auto& r:tbl) if(t>=r.start && t<=r.end) return {r.val,r.err};
    return {-1,-1};
}

// ---------- process ----------------------------------------------------------
void RawAsymmetry::process(TChain& ch, BranchVars& v)
{
    auto beamTbl = readBeamCSV("DB/Beam_pol.csv");
    if(beamTbl.empty()){ std::cerr<<"Beam_pol.csv missing → abort\n"; return; }

    Long64_t n = ch.GetEntries(); Long64_t step=100;
    std::cout<<"[RawAsym] processing "<<n<<" events…\n";

    double count_protons_plus =0.0; 
    double count_protons_minus =0.0;

    double count_neutrons_plus =0.0; 
    double count_neutrons_minus =0.0;

    double total_good_events =0.0;

    std::cout<<"eHCAL cut : "<<c_.eHCAL_L <<std::endl;
    std::cout<<"coin time cut : "<<c_.coin_L <<" to "<<c_.coin_H <<std::endl;
    std::cout<<"W2 cut : "<<c_.W2_L <<" to "<<c_.W2_H <<std::endl;
    std::cout<<"dx cut : "<<c_.dx_L <<" to "<<c_.dx_H <<std::endl;
    std::cout<<"dy cut : "<<c_.dy_L <<" to "<<c_.dy_H <<std::endl;
    std::cout<<"dx_P cut : "<<c_.dx_P_L <<" to "<<c_.dx_P_H <<std::endl;
    std::cout<<"dy_P cut : "<<c_.dy_P_L <<" to "<<c_.dy_P_H <<std::endl;
    std::cout<<"helicity cut : "<<c_.helicity <<std::endl;
    std::cout<<"dx cut center : "<<c_.dx_c <<" +/- "<<c_.dx_r <<std::endl;
    std::cout<<"dy cut center : "<<c_.dy_c <<" +/- "<<c_.dy_r <<std::endl;
    std::cout<<"dx_P cut center : "<<c_.dx_P_c <<" +/- "<<c_.dx_P_r <<std::endl;
    std::cout<<"dy_P cut center : "<<c_.dy_P_c <<" +/- "<<c_.dy_P_r <<std::endl;


    // --- histograms for basic cut checks ---
    TH1D h_ePS     ("h_ePS",      "ePS after basic cuts; ePS (GeV); events",           100, 0.0, 2.0);
    TH1D h_vz      ("h_vz",       "vz after basic cuts; vz (m); events",               120,-0.30,0.30);
    TH1D h_eHCAL   ("h_eHCAL",    "eHCAL after basic cuts; eHCAL (GeV); events",       100, 0.0, 5.0);
    TH1D h_EoverP  ("h_EoverP",   "(ePS+eSH)/p after basic cuts; E/p; events",         100, 0.0, 2.0);
    TH1D h_hel     ("h_hel",      "helicity after basic cuts; helicity; events",        5, -2.5,2.5);
    TH1D h_W2      ("h_W2",       "W^{2} after basic cuts; W^{2} (GeV^{2}); events",   100, 0.0,10.0);
    TH1D h_coin    ("h_coin",     "coin time after basic cuts; coin time (ns); events",100,100.0,200.0);

    TH1D h_dx_n    ("h_dx_n",     "dx (neutron window); dx (m); events",               100,-5.0,5.0);
    TH1D h_dy_n    ("h_dy_n",     "dy (neutron window); dy (m); events",               100,-5.0,5.0);
    TH1D h_dx_p    ("h_dx_p",     "dx (proton window); dx (m); events",                100,-5.0,5.0);
    TH1D h_dy_p    ("h_dy_p",     "dy (proton window); dy (m); events",                100,-5.0,5.0);
    // -------------------------------------------------------------------

    for(Long64_t i=0;i<n;++i){ ch.GetEntry(i);
        //if(!cuts_.passAll(v)) continue;
        if(rq_ && (!rq_->helicityOK(v.runnum)||!rq_->mollerOK(v.runnum))) continue;

        double EoverP = (v.ePS + v.eSH) / v.trP;

        double EoverPmone = EoverP -1.0;

        if (v.ntrack>0 && v.ePS>0.2 && abs(v.vz)<0.27 && v.eHCAL>c_.eHCAL_L && EoverPmone<0.2 && EoverPmone>-0.2 && abs(v.helicity)==1){
            // I don't know what's wrong with the E/p cut above, so I replaced it with this, absolute value did not work for some reason.
            // Probably its better to add E/p as a global cut anyway, just to be safe 
            //++total_good_events;

            // ---------- fill "before W2/coin/dx/dy" distributions ----------
            double EoverP = (v.ePS + v.eSH) / v.trP;
            h_ePS.Fill(v.ePS);
            h_vz.Fill(v.vz);
            h_eHCAL.Fill(v.eHCAL);
            h_EoverP.Fill(EoverP);
            h_hel.Fill(v.helicity);



            if((c_.W2_L<v.W2 && v.W2<c_.W2_H) && (c_.coin_L<v.coin_time && v.coin_time<c_.coin_H) 
                && (c_.dx_P_L<v.dx && v.dx<c_.dx_P_H) && (c_.dy_P_L<v.dy && v.dy<c_.dy_P_H)&& 
                ((pow((v.dy-c_.dy_P_c)/c_.dy_P_r,2)+pow((v.dx-c_.dx_P_c)/c_.dx_P_r,2))<1)){

                int helCorr_p = -1*v.helicity*v.IHWP*c_.Pkin_L;
                if(helCorr_p==1) ++count_protons_plus; else if(helCorr_p==-1) ++count_protons_minus;
                // dx,dy in proton window
                h_dx_p.Fill(v.dx);
                h_dy_p.Fill(v.dy);

            }

            else if( (c_.W2_L<v.W2 && v.W2<c_.W2_H) && (c_.coin_L<v.coin_time && v.coin_time<c_.coin_H) 
                && (c_.dx_L<v.dx && v.dx<c_.dx_H) && (c_.dy_L<v.dy && v.dy<c_.dy_H) 
                && ((pow((v.dy-c_.dy_c)/c_.dy_r,2)+pow((v.dx-c_.dx_c)/c_.dx_r,2))<1)){

                ++total_good_events;
                auto& cnt = counts_[v.runnum];
                int helCorr = -1*v.helicity*v.IHWP*c_.Pkin_L;
                if(helCorr==1) ++cnt.Np; else if(helCorr==-1) ++cnt.Nm;
                if(helCorr==1) ++count_neutrons_plus; else if(helCorr==-1) ++count_neutrons_minus;

                // NOTE: Axis title says A but you're filling dx; keep if intentional.
                h_.Fill(v.dx);

                h_W2.Fill(v.W2);
                h_coin.Fill(v.coin_time);

                // dx,dy in neutron window
                h_dx_n.Fill(v.dx);
                h_dy_n.Fill(v.dy);

                // --- Treat polarization errors as SYSTEMATIC (no 1/sigma^2 event weights) ---
                auto& pa = pols_[v.runnum];
                const double w_evt = 1.0; // replace with charge or live-time if desired

                // Beam: get CSV interval and accumulate event-weighted mean and block totals
                BeamHit bh = beamAtHit(beamTbl, *v.datetime);
                if (bh.idx >= 0 && bh.val > 0) {
                    pa.sumWb  += w_evt;
                    pa.sumWPb += w_evt * bh.val;
                    pa.beam_blk_sumWP[bh.idx] += w_evt * bh.val;
                    // store ABSOLUTE systematic error for that block (same units as val)
                    pa.beam_blk_err.emplace(bh.idx, bh.err);
                }

                // Target: group by (He3Pol, err_He3Pol). err_He3Pol is ABS systematic per block.
                if (v.He3Pol > 0 && v.err_He3Pol > 0) {
                    pa.sumWt  += w_evt;
                    pa.sumWPt += w_evt * v.He3Pol;
                    long long k = makeKey(v.He3Pol, v.err_He3Pol);
                    pa.targ_blk_sumWP[k] += w_evt * v.He3Pol;
                    pa.targ_blk_err.emplace(k, v.err_He3Pol);
                }

                if(i%step==0||i==n-1){ double f=double(i+1)/n; int bar=42,pos=int(bar*f);
                    std::cout<<'\r'<<'['; for(int j=0;j<bar;++j) std::cout<<(j<pos?'=':(j==pos?'>':' '));
                    std::cout<<"] "<<int(f*100)<<" %"<<std::flush; 
                    //std::cout<<"err_He3Pol = " << v.err_He3Pol<<std::endl;
                }
            }
        }
    }
    std::cout<<"\nDone.\n";

    // histogram
    { TFile f(Form("rootfiles/%s/raw_asymmetry_%s.root",kin_,kin_),"RECREATE"); h_.Write(); }

    // table
    std::ofstream out(Form("corrections/%s/raw_asymmetry_%s.txt",kin_,kin_));
    out<<"#run N+ N- A_raw dA_raw beamPol dBeam targetPol dTgt\n";

    double total_events = 0.0;
    double total_difference = 0.0;
    double total_Np = 0.0;
    double total_Nm = 0.0;

    double A_raw_verify_numerator = 0.0;
    double A_raw_verify_denominator = 0.0;

    std::vector<int> runs; runs.reserve(counts_.size());
    for (const auto& kv : counts_) runs.push_back(kv.first);
    std::sort(runs.begin(), runs.end());

    // Global accumulators for means and for systematic combination
    double SWb  = 0.0, SWPb  = 0.0; // beam: Σw, Σ(wP) across runs
    double SWt  = 0.0, SWPt  = 0.0; // target

    std::unordered_map<int, double>       beam_blk_sumWP_glob;
    std::unordered_map<int, double>       beam_blk_err_glob;   // abs error per beam block id
    std::unordered_map<long long, double> targ_blk_sumWP_glob;
    std::unordered_map<long long, double> targ_blk_err_glob;   // abs error per target block key

    for (int run : runs) {
        const auto& c = counts_.at(run);
        const auto& p = pols_.at(run);

        // Per-run raw asymmetry and its (stat) error
        double Np = static_cast<double>(c.Np);
        double Nm = static_cast<double>(c.Nm);
        double N  = Np + Nm;
        double A  = 0.0, dA = 0.0;
        if (N > 0.0) {
            A  = (Np - Nm) / N;
            dA = std::sqrt(std::max(0.0, (1.0 - A*A)/N));
            total_events     += N;
            total_difference += (Np - Nm);
            total_Np         += Np;
            total_Nm         += Nm;
            if (dA>0){
                A_raw_verify_numerator += (A/(dA*dA));
                A_raw_verify_denominator += (1/(dA*dA));
            }
        }

        // Per-run event-weighted polarization means
        double b  = (p.sumWb > 0) ? (p.sumWPb / p.sumWb) : -1.0;
        double t  = (p.sumWt > 0) ? (p.sumWPt / p.sumWt) : -1.0;

        // Per-run SYSTEMATIC errors (absolute), blockwise combination
        double db = -1.0, dt = -1.0;
        if (b > 0.0 && p.sumWPb > 0.0) {
            double S = 0.0;
            for (const auto& kvb : p.beam_blk_sumWP) {
                int    idx          = kvb.first;
                double sumWP_block  = kvb.second;
                double mu           = sumWP_block / p.sumWPb;  // contribution fraction
                double delta_frac   = p.beam_blk_err.at(idx) / b; // fractional sys of that block at Pbar scale
                S += (delta_frac) * (delta_frac); //removed mu may be should add back
            }
            db = b * std::sqrt(S); // absolute
        }
        if (t > 0.0 && p.sumWPt > 0.0) {
            double S = 0.0;
            for (const auto& kvt : p.targ_blk_sumWP) {
                long long key         = kvt.first;
                if(p.targ_blk_err.at(key)>0){
                    
                    double    sumWP_block = kvt.second;
                    double mu           = sumWP_block / p.sumWPt;
                    double delta_frac   = p.targ_blk_err.at(key) / t; // fractional at Pbar scale
                    S += (delta_frac ) * (delta_frac ); //mu removed for a check
                }
            }
            dt = t * std::sqrt(S/p.sumWt); // absolute // not sure revisit
        }

        out << run << " " << c.Np << " " << c.Nm << " "
            << A   << " " << dA   << " "
            << b   << " " << db   << " "
            << t   << " " << dt   << "\n";

        // Accumulate global means
        SWb  += p.sumWb;   SWPb += p.sumWPb;
        SWt  += p.sumWt;   SWPt += p.sumWPt;

        // Accumulate global blockwise sums for systematic combination
        for (const auto& kvb : p.beam_blk_sumWP) {
            int idx = kvb.first; double s = kvb.second;
            beam_blk_sumWP_glob[idx] += s;
            // store one representative abs error per block id
            if (!beam_blk_err_glob.count(idx)) beam_blk_err_glob[idx] = p.beam_blk_err.at(idx);
        }
        for (const auto& kvt : p.targ_blk_sumWP) {
            long long k = kvt.first; double s = kvt.second;
            targ_blk_sumWP_glob[k] += s;
            if (!targ_blk_err_glob.count(k)) targ_blk_err_glob[k] = p.targ_blk_err.at(k);
        }
    }

    out.close();
    std::cout<<"[RawAsym] table → "<<txtF_<<"\n";

    // Global raw asymmetry (stat only)
    double A_raw    = (total_events > 0.0) ? (total_difference/total_events) : 0.0;
    double err_Araw = (total_events > 0.0) ? std::sqrt(std::max(0.0, (1.0 - A_raw*A_raw)/total_events)) : 0.0;

    //Global raw proton asymmetry
    double total_protons = count_protons_plus + count_protons_minus;
    double diff_protons = count_protons_plus - count_protons_minus;

    double A_raw_p =  (total_protons>0.0) ? (diff_protons/total_protons) : 0.0;
    double err_Araw_p = (total_protons>0.0) ? std::sqrt(std::max(0.0,(1.0 - A_raw_p*A_raw_p)/total_protons)) : 0.0;

    // Global polarization means (event-weighted across runs)
    double total_avg_beam_polarization   = (SWb>0.0) ? (SWPb/SWb) : -1.0;
    double total_avg_target_polarization = (SWt>0.0) ? (SWPt/SWt) : -1.0;

    // Global SYSTEMATIC errors (absolute), blockwise combination across all runs
    double err_total_avg_beam_polarization   = -1.0;
    double err_total_avg_target_polarization = -1.0;

    double A_raw_verify = A_raw_verify_numerator/A_raw_verify_denominator;


    if (total_avg_beam_polarization > 0.0 && SWPb > 0.0) {
        double S = 0.0;
        for (const auto& kvb : beam_blk_sumWP_glob) {
            int    idx          = kvb.first;
            double sumWP_block  = kvb.second;
            double mu           = sumWP_block / SWPb;
            double delta_frac   = beam_blk_err_glob.at(idx) / total_avg_beam_polarization;
            S += (delta_frac) * (delta_frac);  //removed mu may be should add back
        }
        err_total_avg_beam_polarization = total_avg_beam_polarization * std::sqrt(S);
    }
    if (total_avg_target_polarization > 0.0 && SWPt > 0.0) {
        double S = 0.0;
        for (const auto& kvt : targ_blk_sumWP_glob) {
            long long key        = kvt.first;
            double    sumWP_block= kvt.second;
            double mu           = sumWP_block / SWPt;
            double delta_frac   = targ_blk_err_glob.at(key) / total_avg_target_polarization;
            S += (delta_frac) * (delta_frac ); //removed mu may be should add back
        }
        err_total_avg_target_polarization = total_avg_target_polarization * std::sqrt(S/SWt);// not sure revisit
    }

    // Write global polarization summary (convert to fraction if your downstream expects it)
    std::ofstream outpol(Form("corrections/%s/avg_polarizations_%s.txt",kin_,kin_));
    outpol<<"avg_beampol = "<<total_avg_beam_polarization*0.01<<"\n";
    outpol<<"err_avg_beampol = "<<err_total_avg_beam_polarization*0.01<<"\n";
    outpol<<"avg_He3pol = "<<total_avg_target_polarization*0.01<<"\n";
    outpol<<"err_avg_He3pol = "<<err_total_avg_target_polarization*0.01<<"\n";
    outpol<<"avg_Pn = "<<0.96<<"\n";
    outpol<<"err_avg_Pn = "<<0.005<<"\n";
    outpol.close();

    // Raw asymmetry summary
    std::ofstream rawAsym(Form("txt/%s/raw_asym_%s.txt",kin_,kin_));
    rawAsym<<"total_good_events = "<<total_good_events<<"\n";
    rawAsym<<"Np = "<<total_Np<<"\n";
    rawAsym<<"Nm = "<<total_Nm<<"\n";
    rawAsym<<"A_raw = "<<A_raw<<"\n";
    rawAsym<<"err_A_raw = "<<err_Araw<<"\n";
    rawAsym<<"Np_p = "<<count_protons_plus<<"\n";
    rawAsym<<"Np_m = "<<count_protons_minus<<"\n";
    rawAsym<<"A_raw_p = "<<A_raw_p<<"\n";
    rawAsym<<"err_Araw_p = "<<err_Araw_p<<"\n";
    rawAsym<<"Np_verify = "<<count_neutrons_plus<<"\n";
    rawAsym<<"Nm_verify = "<<count_neutrons_minus<<"\n";
    rawAsym<<"A_raw_verify = "<<A_raw_verify<<"\n";
    rawAsym<<"err_A_raw_verify = "<<std::sqrt(1/A_raw_verify_denominator)<<"\n";
    rawAsym.close();

    std::string pdfName = Form("plots/cuts_test_raw_asymmetry_%s.pdf", kin_);

    TCanvas c("c_rawasym", "RawAsymmetry QA", 1600, 1200);
    c.Divide(2,2);

    // ---------------- Page 1: basic spectrometer quantities --------------
    c.cd(1); h_ePS.Draw();
    c.cd(2); h_vz.Draw();
    c.cd(3); h_eHCAL.Draw();
    c.cd(4); h_EoverP.Draw();

    // First page: open the PDF with "("
    c.Print((pdfName + "(").c_str());

    // ---------------- Page 2: helicity / W2 / coin / (optional) dx_n ----
    c.Clear(); 
    c.Divide(2,2);

    c.cd(1); h_hel.Draw();
    c.cd(2); h_W2.Draw();
    c.cd(3); h_coin.Draw();
    c.cd(4); h_dx_n.Draw();   // neutron dx distribution

    // Middle page: normal print
    c.Print(pdfName.c_str());

    // ---------------- Page 3: remaining dx/dy distributions --------------
    c.Clear();
    c.Divide(2,2);

    c.cd(1); h_dy_n.Draw();   // neutron dy
    c.cd(2); h_dx_p.Draw();   // proton dx
    c.cd(3); h_dy_p.Draw();   // proton dy
    c.cd(4); h_.Draw();       // your original h_ (dx filled in neutron window)

    // Last page: close the PDF with ")"
    c.Print((pdfName + ")").c_str());

}
