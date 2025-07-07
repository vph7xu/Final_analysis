#include "RawAsymmetry.h"

#include <TFile.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

// ---------- helper -----------------------------------------------------------
static TDatime parseDT(const std::string& s)
{
    int y,m,d,h,min,sec;
    std::sscanf(s.c_str(), "%d-%d-%d %d:%d:%d", &y,&m,&d,&h,&min,&sec);
    return TDatime(y,m,d,h,min,sec);
}

// ---------- constructor ------------------------------------------------------
RawAsymmetry::RawAsymmetry(const CutManager& cuts, const AnalysisCuts& a,
                           const RunQuality* rq, const char* kin,
                           int nb,double xlo,double xhi,
                           const char* root,const char* txt)
    : cuts_(cuts), c_(a), rq_(rq), kin_(kin),
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

        auto& cnt = counts_[v.runnum];
        int helCorr = -1*v.helicity*v.IHWP*c_.Pkin_L;
        if(helCorr==1) ++cnt.Np; else if(helCorr==-1) ++cnt.Nm;

        h_.Fill(v.dx);

        auto bp = beamAt(beamTbl,*v.datetime);
        auto& pa = pols_[v.runnum];
        if(bp.first>0){ pa.sumBeam += bp.first; pa.sumBeamErr2 += bp.second*bp.second; }
        pa.sumHe3 += v.He3Pol; pa.sumHe3Err2 += 0.05*0.05; // FIXME const err
        ++pa.w;

        if(i%step==0||i==n-1){ double f=double(i+1)/n; int bar=42,pos=int(bar*f);
            std::cout<<'\r'<<'['; for(int j=0;j<bar;++j) std::cout<<(j<pos?'=':(j==pos?'>':' '));
            std::cout<<"] "<<int(f*100)<<" %"<<std::flush; }
    }
    std::cout<<"\nDone.\n";

    // histogram
    { TFile f(Form("rootfiles/raw_asymmetry_%s.root",kin_),"RECREATE"); h_.Write(); }

    // table
    std::ofstream out(Form("corrections/raw_asymmetry_%s.txt",kin_));
    out<<"#run N+ N- A_raw dA_raw beamPol dBeam targetPol dTgt\n";

    std::vector<int> runs; runs.reserve(counts_.size());
    for (const auto& kv : counts_) runs.push_back(kv.first);
    std::sort(runs.begin(), runs.end());

    for (int run : runs) {
        const auto& c = counts_.at(run);
        const auto  p = pols_.at(run);

        long long sum = c.Np + c.Nm;
        double A  = 0.0, dA = 0.0;
        if (sum) {
            A  = double(c.Np - c.Nm) / sum;
            dA = 2.0 * std::sqrt(double(c.Np) * c.Nm) / (sum * sum);
        }

        double b  = p.w ? p.sumBeam / p.w                : -1;
        double db = p.w ? std::sqrt(p.sumBeamErr2) / p.w : -1;
        double t  = p.w ? p.sumHe3  / p.w                : -1;
        double dt = p.w ? std::sqrt(p.sumHe3Err2) / p.w  : -1;

        out << run << " " << c.Np << " " << c.Nm << " "
            << A   << " " << dA   << " "
            << b   << " " << db   << " "
            << t   << " " << dt   << "\n";
    }


    // table
    /*std::ofstream out(Form("raw_asymmetry_%s.txt",kin_));
    out<<"#run N+ N- A_raw dA_raw beamPol dBeam targetPol dTgt\n";

    for(auto& [run,c]:counts_){ auto p=pols_[run]; long long sum=c.Np+c.Nm;
        double A=0,dA=0; if(sum){ A=double(c.Np-c.Nm)/sum; dA=2*std::sqrt(double(c.Np)*c.Nm)/(sum*sum);}        
        double b  = p.w? p.sumBeam/p.w : -1;
        double db = p.w? std::sqrt(p.sumBeamErr2)/p.w : -1;
        double t  = p.w? p.sumHe3/p.w  : -1;
        double dt = p.w? std::sqrt(p.sumHe3Err2)/p.w : -1;
        out<<run<<" "<<c.Np<<" "<<c.Nm<<" "<<A<<" "<<dA<<" "<<b<<" "<<db<<" "<<t<<" "<<dt<<"\n";
    }*/
    std::cout<<"[RawAsym] table → "<<txtF_<<"\n";
}
