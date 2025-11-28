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
    if (abs(v.vz)>=0.27) return false;
    if (v.ePS<=0.2) return false;
    if (v.eHCAL<=cuts_.eHCAL_L) return false;
    if (abs((v.ePS+v.eSH)/v.trP -1)>=0.2) return false;

    if(cuts_.coin_L >= v.coin_time || v.coin_time >=cuts_.coin_H) return false;
    if(cuts_.W2_L  >= v.W2       || v.W2      >=cuts_.W2_H  )  return false;
    if(cuts_.dx_L  >= v.dx       || v.dx      >=cuts_.dx_H  )  return false;
    if(cuts_.dy_L  >= v.dy       || v.dy      >=cuts_.dy_H  )  return false;
    if(((pow((v.dy-cuts_.dy_c)/cuts_.dy_r,2)+pow((v.dx-cuts_.dx_c)/cuts_.dx_r,2))>=1))   return false;
    if(std::abs(v.helicity)!=1) return false;

    if(rq_ && (!rq_->helicityOK(v.runnum) || !rq_->mollerOK(v.runnum))) return false;
    return true;
}

// ---------------------------------------------------------------------------
bool AvgKinematics::process(TChain& ch, BranchVars& v)
{
    int total_good_events = 0; //redundant but whatever 

    long long N=0; double Q2=0,tau=0,eps=0,Px=0,Pz=0; const double mN=0.938;

    const int nexp = 6;
	double T_avg[nexp] = {0};
	double T_1_Q2 = 0;
	double T_1_tau = 0;
	double T_1_epsilon = 0;
	double T_1_Px = 0;
	double T_1_Pz = 0;

    Long64_t nentries = ch.GetEntries();
    const Long64_t step     = 100;

    std::cout<<"\n";
    std::cout << "[AvgKinematics] looping over " << nentries << " events\n";

    std::cout<<"dx cuts : "<<cuts_.dx_c<<" +/- "<<cuts_.dx_r<<std::endl;
    std::cout<<"dy cuts : "<<cuts_.dy_c<<" +/- "<<cuts_.dy_r<<std::endl;
    std::cout<<"coin time cuts : "<<cuts_.coin_L<<" to "<<cuts_.coin_H<<std::endl;
    std::cout<<"W2 cuts : "<<cuts_.W2_L<<" to "<<cuts_.W2_H<<std::endl;
    //std::cout<<"dx_P cuts : "<<cuts_.dx_P_c<<" +/- "<<cuts_.dx_P_r<<std::endl;
    //std::cout<<"dy_P cuts : "<<cuts_.dy_P_c<<" +/- "<<cuts_.dy_P_r<<std::endl;
    std::cout<<"helicity cut : "<<cuts_.helicity<<std::endl;
    std::cout<<"eHCAL cut : "<<cuts_.eHCAL_L<<std::endl;


    for(Long64_t i=0;i<ch.GetEntries();++i){ 
        ch.GetEntry(i); 


        //if(!good(v)) continue;

        if(rq_ && (!rq_->helicityOK(v.runnum)||!rq_->mollerOK(v.runnum))) continue;


        if (v.ntrack>0 && v.ePS>0.2 && abs(v.vz)<0.27 && v.eHCAL>cuts_.eHCAL_L && abs(((v.ePS+v.eSH)/v.trP)-1.0)<0.2 && abs(v.helicity)==1){
            
        if( (cuts_.W2_L<v.W2 && v.W2<cuts_.W2_H) && (cuts_.coin_L<v.coin_time && v.coin_time<cuts_.coin_H) 
                && (cuts_.dx_L<v.dx && v.dx<cuts_.dx_H) && (cuts_.dy_L<v.dy && v.dy<cuts_.dy_H) 
                && ((pow((v.dy-cuts_.dy_c)/cuts_.dy_r,2)+pow((v.dx-cuts_.dx_c)/cuts_.dx_r,2))<1)){

        ++total_good_events;

        TLorentzVector Pe(0,0,v.ebeam,v.ebeam),PeP(v.trPx,v.trPy,v.trPz,v.trP);
        
        TLorentzVector q=Pe-PeP; TVector3 qv=q.Vect().Unit();
        
        TVector3 n=Pe.Vect().Cross(PeP.Vect()).Unit();
        
        double tauEv=v.Q2/(4*mN*mN);
        double epsEv=1.0/(1+2*(1+tauEv)*std::tan(v.etheta/2)*std::tan(v.etheta/2));
        
        double PxEv=n.Dot(qv.Cross(tgtPolDir_));
        double PzEv=qv.Dot(tgtPolDir_);
        
		double B = -2 * sqrt(tauEv*(1+tauEv)) * std::tan(v.etheta/2) * PxEv;
		double C = -2 * tauEv * sqrt(1+tauEv+ pow((1+tauEv)*std::tan(v.etheta/2),2)) * std::tan(v.etheta/2) * PzEv;
		double D = tauEv/epsEv;

		double T_0 = C/D;
		double T_1 = B/D;
		double T_2 = -1*C / (D*D);
	  	double T_3 = -1*B / (D*D);
	  	double T_4 = C / (D*D*D);
	  	double T_5 = B / (D*D*D);

        ++N; 

        T_avg[0] += (T_0 - T_avg[0]) / N;
		T_avg[1] += (T_1 - T_avg[1]) / N;
		T_avg[2] += (T_2 - T_avg[2]) / N;
		T_avg[3] += (T_3 - T_avg[3]) / N;
		T_avg[4] += (T_4 - T_avg[4]) / N;
		T_avg[5] += (T_5 - T_avg[5]) / N;

        T_1_Q2 = T_1*v.Q2;
        T_1_tau = T_1*tauEv;
        T_1_epsilon = T_1*epsEv;
        T_1_Px = T_1*PxEv;
        T_1_Pz = T_1*PzEv; 

        

        auto upd=[N](double& a,double val){ a+=(val-a)/N; };
        
        //upd(Q2,v.Q2); upd(tau,tauEv); upd(eps,epsEv); upd(Px,PxEv); upd(Pz,PzEv);
        upd(Q2,T_1_Q2); upd(tau,T_1_tau); upd(eps,T_1_epsilon); upd(Px,T_1_Px); upd(Pz,T_1_Pz); 

        }
        }

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
    
    txt<<"Q2_avg = "<<Q2/T_avg[1]<<"\n"<<"tau_avg = "<<tau/T_avg[1]<<"\n"<<"epsilon_avg = "<<eps/T_avg[1]<<"\n"
       <<"Px_avg = "<<Px/T_avg[1]<<"\n"<<"Pz_avg = "<<Pz/T_avg[1]<<"\n"<<"total_good_events = "<<total_good_events<<"\n";
    
    std::cout<<"[AvgKin] wrote "<<out<<" ("<<N<<" events)\n";
    
    return true;

}
