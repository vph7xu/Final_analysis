#include "RawAsymmetry.h"

#include <TFile.h>
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

    for(Long64_t i=0;i<n;++i){ ch.GetEntry(i);
        if(!cuts_.passAll(v)) continue;
        if(rq_ && (!rq_->helicityOK(v.runnum)||!rq_->mollerOK(v.runnum))) continue;

        if (v.ntrack>0 && v.ePS>0.2 && std::abs(v.vz)<0.27 && v.eHCAL>c_.eHCAL_L){
            if( (c_.W2_L<v.W2 && v.W2<c_.W2_H) && (c_.coin_L<v.coin_time && v.coin_time<c_.coin_H) 
                && (c_.dx_L<v.dx && v.dx<c_.dx_H) && (c_.dy_L<v.dy && v.dy<c_.dy_H)){

                auto& cnt = counts_[v.runnum];
                int helCorr = -1*v.helicity*v.IHWP*c_.Pkin_L;
                if(helCorr==1) ++cnt.Np; else if(helCorr==-1) ++cnt.Nm;

                // NOTE: Axis title says A but you're filling dx; keep if intentional.
                h_.Fill(v.dx);

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
                    std::cout<<"] "<<int(f*100)<<" %"<<std::flush; }
            }
        }
    }
    std::cout<<"\nDone.\n";

    // histogram
    { TFile f(Form("rootfiles/raw_asymmetry_%s.root",kin_),"RECREATE"); h_.Write(); }

    // table
    std::ofstream out(Form("corrections/raw_asymmetry_%s.txt",kin_));
    out<<"#run N+ N- A_raw dA_raw beamPol dBeam targetPol dTgt\n";

    double total_events = 0.0;
    double total_difference = 0.0;
    double total_Np = 0.0;
    double total_Nm = 0.0;

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
                S += (delta_frac * mu) * (delta_frac * mu);
            }
            db = b * std::sqrt(S); // absolute
        }
        if (t > 0.0 && p.sumWPt > 0.0) {
            double S = 0.0;
            for (const auto& kvt : p.targ_blk_sumWP) {
                long long key         = kvt.first;
                double    sumWP_block = kvt.second;
                double mu           = sumWP_block / p.sumWPt;
                double delta_frac   = p.targ_blk_err.at(key) / t; // fractional at Pbar scale
                S += (delta_frac * mu) * (delta_frac * mu);
            }
            dt = t * std::sqrt(S); // absolute
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

    // Global polarization means (event-weighted across runs)
    double total_avg_beam_polarization   = (SWb>0.0) ? (SWPb/SWb) : -1.0;
    double total_avg_target_polarization = (SWt>0.0) ? (SWPt/SWt) : -1.0;

    // Global SYSTEMATIC errors (absolute), blockwise combination across all runs
    double err_total_avg_beam_polarization   = -1.0;
    double err_total_avg_target_polarization = -1.0;

    if (total_avg_beam_polarization > 0.0 && SWPb > 0.0) {
        double S = 0.0;
        for (const auto& kvb : beam_blk_sumWP_glob) {
            int    idx          = kvb.first;
            double sumWP_block  = kvb.second;
            double mu           = sumWP_block / SWPb;
            double delta_frac   = beam_blk_err_glob.at(idx) / total_avg_beam_polarization;
            S += (delta_frac * mu) * (delta_frac * mu);
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
            S += (delta_frac * mu) * (delta_frac * mu);
        }
        err_total_avg_target_polarization = total_avg_target_polarization * std::sqrt(S);
    }

    // Write global polarization summary (convert to fraction if your downstream expects it)
    std::ofstream outpol(Form("corrections/avg_polarizations_%s.txt",kin_));
    outpol<<"avg_beampol = "<<total_avg_beam_polarization*0.01<<"\n";
    outpol<<"err_avg_beampol = "<<err_total_avg_beam_polarization*0.01<<"\n";
    outpol<<"avg_He3pol = "<<total_avg_target_polarization*0.01<<"\n";
    outpol<<"err_avg_He3pol = "<<err_total_avg_target_polarization*0.01<<"\n";
    outpol<<"avg_Pn = "<<0.96<<"\n";
    outpol<<"err_avg_Pn = "<<0.005<<"\n";
    outpol.close();

    // Raw asymmetry summary
    std::ofstream rawAsym(Form("txt/raw_asym_%s.txt",kin_));
    rawAsym<<"Np = "<<total_Np<<"\n";
    rawAsym<<"Nm = "<<total_Nm<<"\n";
    rawAsym<<"A_raw = "<<A_raw<<"\n";
    rawAsym<<"err_A_raw = "<<err_Araw<<"\n";
    rawAsym.close();
}
