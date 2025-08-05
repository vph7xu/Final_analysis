#include "AvgKinematics.h"
//#include "parse.h"
//#include "models.h"
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

using std::string;

// ---------------- helper ----------------------------------------------------
AvgKinematics::Ang AvgKinematics::fieldCSV() const
{
    std::ifstream f("DB/Field_Meas.csv"); 
    string line; Ang a; 
    if(!f) return a;
    std::getline(f,line);
    while(std::getline(f,line)){
        std::stringstream ss(line);

        string tag,h,v;

        std::getline(ss,tag,',');
        std::getline(ss,h,','); 
        std::getline(ss,v,',');

        tag.erase(0, tag.find_first_not_of(" \t"));
        tag.erase(tag.find_last_not_of(" \t") + 1);

        if(tag==kin_){ 
            a.h=std::stod(h); 
            a.v=std::stod(v); 
            break; 
        }
    }

    std::cout<<"horizontal : "<<a.h <<" vertical : "<<a.v <<std::endl;

    return a;
}

// ---------------------------------------------------------------------------
AvgKinematics::AvgKinematics(const AnalysisCuts& c,const RunQuality* rq,const string& k)
    : cuts_(c), rq_(rq) ,kin_(k)
{
    auto ang=fieldCSV();
    tgtPolDir_.SetMagThetaPhi(1.0, ang.h*TMath::Pi()/180.0,
                                   ang.v*TMath::Pi()/180.0);
}

// ---------------------------------------------------------------------------
bool AvgKinematics::good(const BranchVars& v) const
{
    if (v.ntrack<1) return false;
    if (abs(v.vz)>0.27) return false;
    if (v.ePS<0.2) return false;
    if (v.eHCAL<0.025) return false;
    if (abs((v.ePS+v.eSH)/v.trP -1)>0.2) return false;

    if(cuts_.coin_L> v.coin_time || v.coin_time>cuts_.coin_H) return false;
    if(cuts_.W2_L  > v.W2       || v.W2      >cuts_.W2_H  )  return false;
    if(cuts_.dx_L  > v.dx       || v.dx      >cuts_.dx_H  )  return false;
    if(cuts_.dy_L  > v.dy       || v.dy      >cuts_.dy_H  )  return false;
    if(v.eHCAL < cuts_.eHCAL_L)                             return false;
    if(std::abs(v.helicity)!=1)                              return false;

    if(rq_ && (!rq_->helicityOK(v.runnum) || !rq_->mollerOK(v.runnum))) return false;
    return true;
}

// ---------------------------------------------------------------------------
bool AvgKinematics::process(TChain& ch, BranchVars& v)
{

    long long N=0; double Q2=0,tau=0,eps=0,Px=0,Pz=0; const double mN=0.938;

    Long64_t nentries = ch.GetEntries();
    const Long64_t step     = 100;

    std::cout<<"\n";
    std::cout << "[AvgKinematics] looping over " << nentries << " events\n";

    for(Long64_t i=0;i<ch.GetEntries();++i){ 
        ch.GetEntry(i); 
        if(!good(v)) continue;

        TLorentzVector Pe(0,0,v.ebeam,v.ebeam),PeP(v.trPx,v.trPy,v.trPz,v.trP);
        
        TLorentzVector q=Pe-PeP; TVector3 qv=q.Vect().Unit();
        
        TVector3 n=Pe.Vect().Cross(PeP.Vect()).Unit();
        
        double tauEv=v.Q2/(4*mN*mN);
        
        double epsEv=1.0/(1+2*(1+tauEv)*std::tan(v.etheta/2)*std::tan(v.etheta/2));
        
        double PxEv=n.Dot(qv.Cross(tgtPolDir_));
        
        double PzEv=qv.Dot(tgtPolDir_);
        
        ++N; 

        auto upd=[N](double& a,double val){ a+=(val-a)/N; };
        
        upd(Q2,v.Q2); upd(tau,tauEv); upd(eps,epsEv); upd(Px,PxEv); upd(Pz,PzEv); 

        // progress bar
        if (i % step == 0 || i == nentries - 1) {
            double frac = double(i + 1) / nentries;
            int barw = 42, pos = static_cast<int>(barw * frac);
            std::cout << '\r' << '[';
            for (int j = 0; j < barw; ++j)
                std::cout << (j < pos ? '=' : (j == pos ? '>' : ' '));
            std::cout << "] " << static_cast<int>(frac * 100) << " %" << std::flush;
        }

    }

    string out="txt/"+kin_+"/average_kinematic_values_"+kin_+".txt";
    
    std::ofstream txt(out);
    
    txt<<"Q2_avg = "<<Q2<<"\n"<<"tau_avg = "<<tau<<"\n"<<"epsilon_avg = "<<eps<<"\n"
       <<"Px_avg = "<<Px<<"\n"<<"Pz_avg = "<<Pz<<"\n";
    
    std::cout<<"[AvgKin] wrote "<<out<<" ("<<N<<" events)\n";
    
    return true;

}
