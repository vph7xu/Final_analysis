#include "InelasticCorrection.h"
#include "Utility.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TNamed.h>
#include <TLegend.h>
#include <iostream>
#include <cmath>

using std::string; using std::vector; using std::unordered_map;

InelasticCorrection::InelasticCorrection(const AnalysisCuts& cuts,
                                         const RunQuality* rq,
                                         const char* kin,
                                         const char* rootFile)
    : c_(cuts), rq_(rq), kin_(kin), outFile_(rootFile) {}

// --- template shape fit -------------------------------------------------------
TH1D* InelasticCorrection::performFit(TH1D* hD, TH1D* hInel, TH1D* hQE_proton, TH1D* hQE_neutron,
                                      double& par0, double& par1, double& par2)
{
    TH1D *d=(TH1D*)hD->Clone("hD");    d->Scale(1./d->Integral());
    TH1D *inel=(TH1D*)hInel->Clone("hInel"); inel->Scale(1./inel->Integral());
    TH1D *qe_p=(TH1D*)hQE_proton->Clone("hQE_proton");     qe_p->Scale(1./qe_p->Integral());    
    TH1D *qe_n=(TH1D*)hQE_neutron->Clone("hQE_neutron");     qe_n->Scale(1./qe_n->Integral());

    TF1* f = new TF1("fdx",[&](double*x,double*par){
        int bin1 = qe_p->FindBin(x[0]);
        int bin2 = qe_n->FindBin(x[0]);
        int bin3 = inel->FindBin(x[0]);
        return par[0]*(qe_p->GetBinContent(bin1) + par[1]*qe_n->GetBinContent(bin2) + par[2]*inel->GetBinContent(bin3));
    }, d->GetXaxis()->GetXmin(), d->GetXaxis()->GetXmax(), 3);

    f->SetParameters(0.5,0.5); f->SetParLimits(0,0,1); f->SetParLimits(1,0,1);

    d->Fit(f,"RQ");

    par0 = f->GetParameter(0);
    par1 = f->GetParameter(1);
    par2 = f->GetParameter(2);

    TH1D* comb=(TH1D*)inel->Clone("hComb"); comb->Reset();
    for(int i=1;i<=comb->GetNbinsX();++i)
        comb->SetBinContent(i, par0*(qe_p->GetBinContent(i) + par1*qe_n->GetBinContent(i) + par2*inel->GetBinContent(i)));
    return comb;
}

// --- template shape fit -------------------------------------------------------
TH1D* InelasticCorrection::performFitW2(TH1D* hD, TH1D* hInel, TH1D* hQE_proton, TH1D* hQE_neutron,
                                      double& par0, double& par1, double Rnop)
{

    TH1D *d=(TH1D*)hD->Clone("hD");    d->Scale(1./d->Integral());
    TH1D *inel=(TH1D*)hInel->Clone("hInel"); inel->Scale(1./inel->Integral());
    TH1D *qe_p=(TH1D*)hQE_proton->Clone("hQE_proton");     qe_p->Scale(1./qe_p->Integral());    
    TH1D *qe_n=(TH1D*)hQE_neutron->Clone("hQE_neutron");     qe_n->Scale(1./qe_n->Integral());

    TF1* f = new TF1("fW2",[&](double*x,double*par){
        int bin1 = qe_p->FindBin(x[0]);
        int bin2 = qe_n->FindBin(x[0]);
        int bin3 = inel->FindBin(x[0]);
        return par[0]*(qe_p->GetBinContent(bin1) + Rnop*qe_n->GetBinContent(bin2) + par[1]*inel->GetBinContent(bin3));
    }, d->GetXaxis()->GetXmin(), d->GetXaxis()->GetXmax(), 2);

    //f->SetParameters(0.5,0.5); f->SetParLimits(0,0,1); f->SetParLimits(1,0,1);

    f->SetParameters(0.5, 0.5);
    f->SetParLimits(0, 0, 10000);   
    f->SetParLimits(1, 0, 10000);   
    //f->SetParLimits(2, 0, 1);   

    d->Fit(f,"RQ");

    par0 = f->GetParameter(0);
    par1 = f->GetParameter(1);
    //par2 = f->GetParameter(2);

    TH1D* comb=(TH1D*)inel->Clone("hComb"); comb->Reset();
    for(int i=1;i<=comb->GetNbinsX();++i)
        comb->SetBinContent(i, par0*(qe_p->GetBinContent(i) + Rnop*qe_n->GetBinContent(i) + par1*inel->GetBinContent(i)));
    return comb;

}

TH1D* InelasticCorrection::performFitW2_1(TH1D* hD, TH1D* hInel, TH1D* hQEp, TH1D* hQEn,
                                        double& alpha, double Rnop)
{
    TH1D *d=(TH1D*)hD->Clone("hDW2");      d->Scale(1.0/d->Integral());
    TH1D *pi=(TH1D*)hQEp->Clone("qepW2");  pi->Scale(1.0/pi->Integral());
    TH1D *ni=(TH1D*)hQEn->Clone("qenW2");  ni->Scale(1.0/ni->Integral());
    TH1D *ii=(TH1D*)hInel->Clone("inelW2");ii->Scale(1.0/ii->Integral());

    // unique name avoids collisions
    TF1* f = new TF1(Form("fW2_%p", d), [&](double* x, double* par){
        int bp = pi->FindBin(x[0]), bn = ni->FindBin(x[0]), bi = ii->FindBin(x[0]);
        double num = pi->GetBinContent(bp) + Rnop*ni->GetBinContent(bn) + par[0]*ii->GetBinContent(bi);
        return num / (1.0 + Rnop + par[0]);  // normalized mixture PDF
    }, d->GetXaxis()->GetXmin(), d->GetXaxis()->GetXmax(), 1);

    f->SetParameters(0,0.5);        // alpha start
    f->SetParLimits(0, 0.0, 100); // alpha >= 0
    d->Fit(f, "RQ");

    alpha = f->GetParameter(0);

    // Build combined PDF (unit-area) to draw
    TH1D* comb=(TH1D*)ii->Clone("hCombW2"); comb->Reset();
    for (int i=1;i<=comb->GetNbinsX();++i){
        double num = pi->GetBinContent(i) + Rnop*ni->GetBinContent(i) + alpha*ii->GetBinContent(i);
        comb->SetBinContent(i, num/(1.0+Rnop+alpha));
    }
    return comb;
}



// --- main driver --------------------------------------------------------------
void InelasticCorrection::process(TChain& ch, TChain& ch_QE, TChain& ch_inel,
                                  BranchVars& v, BranchVarsSim& vQE, BranchVarsSim& vInel)
{
    //if(c_.inel_W2_L==0 && c_.inel_W2_H==0){
    //    std::cerr << "[InelCorrection] inel_W2_L/H not set\n"; return; }

    //trying to readsome files

    Utility utility;

    auto dx_shifted_QE = [&](const BranchVarsSim& s){
        double d = s.dx;
        if (std::strcmp(kin_, "GEN3_He3") == 0 && s.fnucl==0) d -= 0.1;      // n
        else if (std::strcmp(kin_, "GEN4_He3") == 0 && s.fnucl==1) d -= 0.02;     // p
        else if (std::strcmp(kin_, "GEN4_He3")  == 0 && s.fnucl==0) d -= 0.07;
        else if (std::strcmp(kin_, "GEN4b_He3")== 0 && s.fnucl==0) d -= 0.025;    // n
        return d;
    };

    auto dx_shifted_Inel = [&](const BranchVarsSim& s){
        double d = s.dx;
        if (std::strcmp(kin_, "GEN3_He3") == 0)  d += 0.4;
        else if (std::strcmp(kin_, "GEN4_He3")==0)  d += 0.6;
        else if (std::strcmp(kin_, "GEN4b_He3")==0) d += 0.25;
        return d;
    };

    auto   corrFile = [&](const char* stem){
        return string("corrections/")+kin_+"/" + stem + "Correction_" + kin_ + ".txt"; };

    const string accFile = corrFile("Accidental");
    const string nitFile = corrFile("Nitrogen");
    const string pionFile= corrFile("Pion");

    // -------- read TXT correction files directly -------------------------------
    const auto accMap  = utility.readSimpleKVGlobal(accFile);
    const auto pionMap = utility.readSimpleKVGlobal(pionFile);
    const auto nitMap  = utility.readSimpleKVGlobal(nitFile);

    auto V=[&](const auto& M,const char* k){ auto it=M.find(k); return it!=M.end()? it->second : 0.0; };


    // central values
    double Aacc = V(accMap ,"A_acc");
    double Api  = V(pionMap,"A_pi");           // not provided yet

    double facc = V(accMap ,"f_acc");
    double fpi  = V(pionMap,"f_pi");
    double fN2  = V(nitMap ,"f_N2");

    // errors
    double errAacc = V(accMap ,"err_A_acc");
    double errApi  = V(pionMap,"err_A_pi");

    double errfacc = V(accMap ,"err_f_acc");
    double errfpi  = V(pionMap,"err_f_pi");
    double errfN2  = V(nitMap ,"err_f_N2");

    // histograms for fit (tight selection)
    TH1D hData ("hData",  "dx data; dx (m)", 100, -4.0, 3.0);
    TH1D hData_pos ("hData_pos",  "dx data (positive helicity); dx (m)", 100, -4.0, 3.0);
    TH1D hData_neg ("hData_neg",  "dx data (negative helicity); dx (m)", 100, -4.0, 3.0);
    TH1D hQE_proton("hQE_proton", "dx QE sim protons",          100, -4.0, 3.0);
    TH1D hQE_neutron("hQE_neutron", "dx QE sim neutrons",          100, -4.0, 3.0);
    TH1D hInelastic ("hInelastic",  "dx inelastic sim",        100, -4.0, 3.0);
    TH1D hInelastic_proton ("hInelastic_proton",  "dx inelastic sim protons",        100, -4.0, 3.0);
    TH1D hInelastic_neutron ("hInelastic_neutron",  "dx inelastic sim neutrons",        100, -4.0, 3.0);

    TH1D hData_W2("hData_W2",  "W2 data; W^{2} (GeV^{2})", 50, -1, 2);
    TH1D hQE_proton_W2("hQE_proton_W2", "W^{2} QE sim protons",          50, -1, 2);
    TH1D hQE_neutron_W2("hQE_neutron_W2", "W^{2} QE sim neutrons",          50, -1, 2);
    TH1D hInelastic_W2("hInelastic_W2",  "W^{2} inelastic sim",        50, -1, 2);

    // --- loop data
    Long64_t n=ch.GetEntries();
    std::cout << "\n"<<"[InelasticCorrection] looping over data " << n << " events\n";
    const Long64_t step     = 100;
    for(Long64_t i=0;i<n;++i){ 
        
        ch.GetEntry(i);
        
        if(rq_ && (!rq_->helicityOK(v.runnum)||!rq_->mollerOK(v.runnum))) continue;
   
        if( v.runnum<c_.runnum_L || v.runnum>c_.runnum_H ||
            v.ntrack<1 || abs(v.vz)>0.27 || v.eHCAL<c_.eHCAL_L || abs((v.ePS+v.eSH)/(v.trP)-1)>0.2 || v.ePS<0.2 ||
            (c_.coin_L>v.coin_time || v.coin_time>c_.coin_H) || 
            abs(v.helicity)!=1) continue;
        
        if(((pow((v.dy-0.0)/0.2,2)+pow((v.dx-0.0)/0.2,2))<=1)){
            hData_W2.Fill(v.W2);
        }

        if( v.runnum<c_.runnum_L || v.runnum>c_.runnum_H ||
            v.ntrack<1 || abs(v.vz)>0.27 || v.eHCAL<c_.eHCAL_L || abs((v.ePS+v.eSH)/(v.trP)-1)>0.2 || v.ePS<0.2 ||
            (c_.coin_L>v.coin_time || v.coin_time>c_.coin_H) || (c_.W2_L>v.W2 || v.W2>c_.W2_H) || (c_.dy_L>v.dy || v.dy>c_.dy_H) || 
            abs(v.helicity)!=1) continue; // no dx cut since we are fitting it
        
        hData.Fill(v.dx);
        
        if (v.helicity == 1){
            hData_pos.Fill(v.dx);
        }
        else if(v.helicity == -1){
            hData_neg.Fill(v.dx);
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
    
    // QE sim
    Long64_t nentries_QE=ch_QE.GetEntries();
    std::cout << "\n"<< "[InelasticCorrection] looping over QE sim " << nentries_QE << " events\n";
    for(Long64_t i=0;i<ch_QE.GetEntries();++i){ 

        ch_QE.GetEntry(i);

        if(/*vQE.ntrack<1 ||*/ abs(vQE.vz)>0.27 || vQE.eHCAL<c_.eHCAL_L || abs((vQE.ePS+vQE.eSH)/(vQE.trP)-1)>0.2 || vQE.ePS<0.2) continue;

        double dxq = dx_shifted_QE(vQE);

        if (((pow((vQE.dy-0.0)/0.2,2)+pow((dxq-0.0)/0.2,2))<=1)) {
            if (vQE.fnucl == 0) hQE_neutron_W2.Fill(vQE.W2, vQE.weight);  // see #3
            if (vQE.fnucl == 1) hQE_proton_W2.Fill(vQE.W2, vQE.weight);
        }


        if(/*vQE.ntrack<1 ||*/ abs(vQE.vz)>0.27 || vQE.eHCAL<c_.eHCAL_L || abs((vQE.ePS+vQE.eSH)/(vQE.trP)-1)>0.2 || vQE.ePS<0.2 ||
            (c_.W2_L>vQE.W2 || vQE.W2>c_.W2_H) || (c_.dy_L>vQE.dy || vQE.dy>c_.dy_H)) continue;

        if(vQE.fnucl == 0) {
            if(std::strcmp(kin_, "GEN3_He3") == 0){
                hQE_neutron.Fill(vQE.dx-0.1,vQE.weight);
            }
            else if(std::strcmp(kin_, "GEN4_He3") == 0){
                hQE_neutron.Fill(vQE.dx-0.07,vQE.weight);
            }
            else if(std::strcmp(kin_, "GEN4b_He3") == 0){
                hQE_neutron.Fill(vQE.dx-0.025,vQE.weight);
            }
            else{
                hQE_neutron.Fill(vQE.dx,vQE.weight);
            }

        }

        if(vQE.fnucl == 1) {
            if(std::strcmp(kin_, "GEN3_He3") == 0){
                hQE_proton.Fill(vQE.dx,vQE.weight);
            }
            else if(std::strcmp(kin_, "GEN4_He3") == 0){
                hQE_proton.Fill(vQE.dx-0.02,vQE.weight);
            }
            else{
                hQE_proton.Fill(vQE.dx,vQE.weight);
            }
        }//+0.05 is for GEN3 its hard coded for now

        // progress bar
        if (i % step == 0 || i == nentries_QE - 1) {
            double frac = double(i + 1) / nentries_QE;
            int barw = 42, pos = static_cast<int>(barw * frac);
            std::cout << '\r' << '[';
            for (int j = 0; j < barw; ++j)
                std::cout << (j < pos ? '=' : (j == pos ? '>' : ' '));
            std::cout << "] " << static_cast<int>(frac * 100) << " %" << std::flush;
        }
        
    }

    // inelastic sim
    Long64_t nentries_inel=ch_inel.GetEntries();
    std::cout << "\n"<<"[InelasticCorrection] looping over Inelastic sim " << nentries_inel << " events\n";
    for(Long64_t i=0;i<ch_inel.GetEntries();++i){ 
        ch_inel.GetEntry(i);
    
        if(/*vInel.ntrack<1 ||*/ abs(vInel.vz)>0.27 || vInel.eHCAL<c_.eHCAL_L || abs((vInel.ePS+vInel.eSH)/(vInel.trP)-1)>0.2 || vInel.ePS<0.2) continue;

        double dxi = dx_shifted_Inel(vInel);
        if (((pow((vInel.dy-0.0)/0.2,2)+pow((dxi-0.0)/0.2,2))<=1)) {
            hInelastic_W2.Fill(vInel.W2, vInel.weight);                    // see #3
        }

        if(/*vInel.ntrack<1 ||*/ abs(vInel.vz)>0.27 || vInel.eHCAL<c_.eHCAL_L || abs((vInel.ePS+vInel.eSH)/(vInel.trP)-1)>0.2 || vInel.ePS<0.2 ||
            (c_.W2_L>vInel.W2 || vInel.W2>c_.W2_H) || (c_.dy_L>vInel.dy || vInel.dy>c_.dy_H)) continue;    

        if(std::strcmp(kin_, "GEN3_He3") == 0){
            
            hInelastic.Fill(vInel.dx+0.4,vInel.weight);//+0.4 is for GEN3 its hard coded for now

            if (vInel.fnucl == 0) hInelastic_neutron.Fill(vInel.dx+0.4,vInel.weight);

            if (vInel.fnucl == 1) hInelastic_proton.Fill(vInel.dx+0.4,vInel.weight);
        }
        else if(std::strcmp(kin_, "GEN4_He3") == 0){
            
            hInelastic.Fill(vInel.dx+0.6,vInel.weight);//+0.4 is for GEN3 its hard coded for now

            if (vInel.fnucl == 0) hInelastic_neutron.Fill(vInel.dx+0.6,vInel.weight);

            if (vInel.fnucl == 1) hInelastic_proton.Fill(vInel.dx+0.6,vInel.weight);
        }
        else if(std::strcmp(kin_, "GEN4b_He3") == 0){
            
            //std::cout<<"debug here"<<'\n';

            hInelastic.Fill(vInel.dx+0.25,vInel.weight);//+0.4 is for GEN3 its hard coded for now

            if (vInel.fnucl == 0) hInelastic_neutron.Fill(vInel.dx+0.25,vInel.weight);

            if (vInel.fnucl == 1) hInelastic_proton.Fill(vInel.dx+0.25,vInel.weight);
        }
        else{
            //std::cout<<"debug here in else"<<'\n';

            hInelastic.Fill(vInel.dx,vInel.weight);

            if (vInel.fnucl == 0) hInelastic_neutron.Fill(vInel.dx,vInel.weight);

            if (vInel.fnucl == 1) hInelastic_proton.Fill(vInel.dx,vInel.weight);

        }

        // progress bar
        if (i % step == 0 || i == nentries_inel - 1) {
            double frac = double(i + 1) / nentries_inel;
            int barw = 42, pos = static_cast<int>(barw * frac);
            std::cout << '\r' << '[';
            for (int j = 0; j < barw; ++j)
                std::cout << (j < pos ? '=' : (j == pos ? '>' : ' '));
            std::cout << "] " << static_cast<int>(frac * 100) << " %" << std::flush;
        }

    }    

    double par0=1,par1=1,par2=1; 

    TH1D * h_combined =  performFit(&hData,&hInelastic,&hQE_proton,&hQE_neutron,par0,par1,par2);

    double parP0=1,parP1=1,parP2=1;

    TH1D * h_combined_pos =  performFit(&hData_pos,&hInelastic,&hQE_proton,&hQE_neutron,parP0,parP1,parP2);

    double parN0=1,parN1=1,parN2=1;

    TH1D * h_combined_neg =  performFit(&hData_neg,&hInelastic,&hQE_proton,&hQE_neutron,parN0,parN1,parN2);

    h_combined->Scale(hData.Integral());
    h_combined_pos->Scale(hData_pos.Integral());
    h_combined_neg->Scale(hData_neg.Integral());

    TH1D* hQE_proton_pos = (TH1D*) hQE_proton.Clone("hQE_proton_pos");
    TH1D* hQE_proton_neg = (TH1D*) hQE_proton.Clone("hQE_proton_neg");

    TH1D* hQE_neutron_pos = (TH1D*) hQE_neutron.Clone("hQE_neutron_pos");
    TH1D* hQE_neutron_neg = (TH1D*) hQE_neutron.Clone("hQE_neutron_neg");

    TH1D* hInelastic_pos = (TH1D*) hInelastic.Clone("hInelastic_pos");
    TH1D* hInelastic_neg = (TH1D*) hInelastic.Clone("hInelastic_neg");

    hQE_proton.Scale(par0*hData.Integral()*1/hQE_proton.Integral());
    hQE_neutron.Scale(par0*par1*hData.Integral()*1/hQE_neutron.Integral());
    hInelastic.Scale(par0*par2*hData.Integral()*1/hInelastic.Integral());

    hQE_proton_pos->Scale(parP0*hData_pos.Integral()*1/hQE_proton_pos->Integral());
    hQE_neutron_pos->Scale(parP0*parP1*hData_pos.Integral()*1/hQE_neutron_pos->Integral());
    hInelastic_pos->Scale(parP0*parP2*hData_pos.Integral()*1/hInelastic_pos->Integral());

    hQE_proton_neg->Scale(parN0*hData_neg.Integral()*1/hQE_proton_neg->Integral());
    hQE_neutron_neg->Scale(parN0*parN1*hData_neg.Integral()*1/hQE_neutron_neg->Integral());
    hInelastic_neg->Scale(parN0*parN2*hData_neg.Integral()*1/hInelastic_neg->Integral());

    //frac_  = wI;                 // N_inel / (N_inel+QE)
    //dfrac_ = std::sqrt(wI*wQ/(wI+wQ)/(wI+wQ)/(hD.Integral()));

    double QE_events = hData.Integral(hData.FindBin(c_.dx_L),hData.FindBin(c_.dx_H));
    double inelastic_events = hInelastic.Integral(hInelastic.FindBin(c_.dx_L),hInelastic.FindBin(c_.dx_H));
    double proton_events = hQE_proton.Integral(hQE_proton.FindBin(c_.dx_L),hQE_proton.FindBin(c_.dx_H));
    double neutron_events = hQE_neutron.Integral(hQE_neutron.FindBin(c_.dx_L),hQE_neutron.FindBin(c_.dx_H));

    double Rnop = neutron_events/proton_events;

    double inelastic_events_pos = hInelastic_pos->Integral();//(hInelastic_pos->FindBin(c_.dx_L),hInelastic_pos->FindBin(c_.dx_H));
    double inelastic_events_neg = hInelastic_neg->Integral();//(hInelastic_neg->FindBin(c_.dx_L),hInelastic_neg->FindBin(c_.dx_H));

    const double R     = (1 - facc - fN2 - fpi) / QE_events;
    const double F     = inelastic_events * R;              // inelastic_frac

    const double dN_in   = std::sqrt(inelastic_events);     // Poisson
    const double dN_QE   = std::sqrt(QE_events);            // Poisson
    const double dFacc   = errfacc;                       // your TXT value
    const double dFN2    = errfN2;
    const double dFpi    = errfpi;

    const double dF2 =
        std::pow(R * dN_in, 2) +
        std::pow(F / QE_events * dN_QE, 2) +
        std::pow(inelastic_events / QE_events * dFacc, 2) +
        std::pow(inelastic_events / QE_events * dFN2 , 2) +
        std::pow(inelastic_events / QE_events * dFpi , 2);

    const double dFin = std::sqrt(dF2);

    double background_frac = inelastic_events/QE_events; //before removing other fractions
    double inelastic_frac = inelastic_events * (1 - facc - fN2 - fpi)/QE_events;    
    double errinelastic_frac =  dFin;//(inelastic_events/QE_events)*sqrt((1/inelastic_events)+(1/QE_events)); 
    double proton_frac = proton_events/QE_events;
    double errproton_frac = (proton_events/QE_events)*sqrt((1/proton_events)+(1/QE_events));


    double A_in = (inelastic_events_pos - inelastic_events_neg)/(inelastic_events_pos + inelastic_events_neg);
    double err_A_in = 2.0 * sqrt(inelastic_events_pos * inelastic_events_neg * (inelastic_events_pos + inelastic_events_neg)) / ((inelastic_events_pos + inelastic_events_neg)*(inelastic_events_pos + inelastic_events_neg));

    // store
    std::ofstream txt(Form("corrections/%s/InelasticCorrection_%s.txt",kin_,kin_));
    //TNamed n("inel_fraction", (std::to_string(frac_)+","+std::to_string(dfrac_)).c_str());
    //n.Write(); fout.Close();
    txt<<"par0 = "<< par0 <<"\n";
    txt<<"par1 = "<< par1 <<"\n";
    txt<<"par2 = "<< par2 <<"\n";
    txt<<"background_fraction = "<<background_frac<<"\n";
    txt<<"err_inelastic_fraction = "<<errinelastic_frac<<"\n";
    txt<<"proton_fraction = "<<proton_frac<<"\n";
    txt<<"err_proton_fraction = "<<errproton_frac<<"\n";
    
    txt<<"f_in = "<<inelastic_frac<<"\n";
    txt<<"err_f_in = "<<errinelastic_frac<<"\n";
    txt<<"f_p = "<<proton_frac<<"\n";
    txt<<"err_f_p = "<<errproton_frac<<"\n";
    //below values should be calculated separately, for now they are set to zero
    txt<<"inelastic_events_pos = "<<inelastic_events_pos<<"\n";
    txt<<"inelastic_events_neg = "<<inelastic_events_neg<<"\n";
    txt<<"A_in = "<<A_in<<"\n";
    txt<<"err_A_in = "<<err_A_in<<"\n";
    txt<<"A_p = "<<0.0<<"\n";
    txt<<"err_A_p = "<<0.0<<"\n";


    txt.close();

    std::cout << "[InelasticCorrection] par0 = "<< par0 << std::endl;
    //     << " saved to "<< outFile_ <<"\n";


    ////////////////////W2 fitting /////////////////////////////////////

    std::cout<<"Rnop : "<< Rnop <<std::endl;

    double parW0=1,parW1=1; 

    TH1D * h_combined_W2 =  performFitW2(&hData_W2,&hInelastic_W2,&hQE_proton_W2,&hQE_neutron_W2,parW0,parW1,Rnop);

    h_combined_W2->Scale(hData_W2.Integral());

    hQE_proton_W2.Scale(parW0*hData_W2.Integral()*1/hQE_proton_W2.Integral());
    hQE_neutron_W2.Scale(parW0*Rnop*hData_W2.Integral()*1/hQE_neutron_W2.Integral());
    hInelastic_W2.Scale(parW0*parW1*hData_W2.Integral()*1/hInelastic_W2.Integral());

    std::cout<<"parW0 : "<<parW0<<std::endl;
    std::cout<<"parW1 : "<<parW1<<std::endl;
    //std::cout<<"parW2 : "<<parW2<<std::endl;





    ///////////////////plotting and printing////////////////////////////

    TCanvas *C = new TCanvas("c","c",2400,1500);
    TCanvas *C1 = new TCanvas("c1","c1",2400,1500);
    //TCanvas *C2 = new TCanvas("c2","c2",2400,1500);

    C->Divide(2,2);
    C->cd(1);
    
    h_combined->SetLineColor(3);
    hQE_proton.SetLineColor(6);
    hQE_neutron.SetLineColor(9);
    hInelastic.SetLineColor(7);
    hData.SetLineColor(kBlack);

    h_combined->SetLineWidth(4);
    hQE_proton.SetLineWidth(4);
    hQE_neutron.SetLineWidth(4);
    hInelastic.SetLineWidth(4);
    hData.SetLineWidth(4);

    h_combined->SetFillColorAlpha(19,0.1);
    h_combined->SetFillStyle(3009);
    hQE_proton.SetFillColorAlpha(6,0.5);
    hQE_proton.SetFillStyle(3004);
    hQE_neutron.SetFillColorAlpha(9,0.5);
    hQE_neutron.SetFillStyle(3005);
    hInelastic.SetFillColorAlpha(7,0.5);
    hInelastic.SetFillStyle(3009);

    hData.SetMarkerStyle(kFullCircle);

    hData.Draw("p");
    h_combined->Draw("hist same");
    hQE_proton.Draw("hist same");
    hQE_neutron.Draw("hist same");
    hInelastic.Draw("hist same");

    C->cd(2);
    h_combined_pos->SetLineColor(3);
    hQE_proton_pos->SetLineColor(6);
    hQE_neutron_pos->SetLineColor(9);
    hInelastic_pos->SetLineColor(7);
    hData_pos.SetLineColor(kBlack);

    h_combined_pos->SetLineWidth(4);
    hQE_proton_pos->SetLineWidth(4);
    hQE_neutron_pos->SetLineWidth(4);
    hInelastic_pos->SetLineWidth(4);
    hData_pos.SetLineWidth(4);

    h_combined_pos->SetFillColorAlpha(19,0.1);
    h_combined_pos->SetFillStyle(3009);
    hQE_proton_pos->SetFillColorAlpha(6,0.5);
    hQE_proton_pos->SetFillStyle(3004);
    hQE_neutron_pos->SetFillColorAlpha(9,0.5);
    hQE_neutron_pos->SetFillStyle(3005);
    hInelastic_pos->SetFillColorAlpha(7,0.5);
    hInelastic_pos->SetFillStyle(3009);

    hData_pos.SetMarkerStyle(kFullCircle);

    hData_pos.Draw("p");
    h_combined_pos->Draw("hist same");
    hQE_proton_pos->Draw("hist same");
    hQE_neutron_pos->Draw("hist same");
    hInelastic_pos->Draw("hist same");


    C->cd(3);
    h_combined_neg->SetLineColor(3);
    hQE_proton_neg->SetLineColor(6);
    hQE_neutron_neg->SetLineColor(9);
    hInelastic_neg->SetLineColor(7);
    hData_neg.SetLineColor(kBlack);

    h_combined_neg->SetLineWidth(4);
    hQE_proton_neg->SetLineWidth(4);
    hQE_neutron_neg->SetLineWidth(4);
    hInelastic_neg->SetLineWidth(4);
    hData_neg.SetLineWidth(4);

    h_combined_neg->SetFillColorAlpha(19,0.1);
    h_combined_neg->SetFillStyle(3009);
    hQE_proton_neg->SetFillColorAlpha(6,0.5);
    hQE_proton_neg->SetFillStyle(3004);
    hQE_neutron_neg->SetFillColorAlpha(9,0.5);
    hQE_neutron_neg->SetFillStyle(3005);
    hInelastic_neg->SetFillColorAlpha(7,0.5);
    hInelastic_neg->SetFillStyle(3009);

    hData_neg.SetMarkerStyle(kFullCircle);

    hData_neg.Draw("p");
    h_combined_neg->Draw("hist same");
    hQE_proton_neg->Draw("hist same");
    hQE_neutron_neg->Draw("hist same");
    hInelastic_neg->Draw("hist same");

    C->cd(4);
    auto leg = new TLegend(0.2, 0.2, 0.8, 0.8); // x1,y1,x2,y2 (NDC)
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    // pointer vs object: use h_combined (ptr) and &hX (object)
    leg->AddEntry(&hData,       "Data",                       "p");
    leg->AddEntry(h_combined,   "Fit par0*(QE_{p} + par1*QE_{n} + par2*Inel.)", "lf");
    leg->AddEntry(&hQE_proton,  "QE proton (sim)",            "f");
    leg->AddEntry(&hQE_neutron, "QE neutron (sim)",           "f");
    leg->AddEntry(&hInelastic,  "Inelastic (sim)",            "f");

    leg->Draw();

    C1->Divide(2,2);
    C1->cd(1);

    h_combined_W2->SetLineColor(3);
    hQE_proton_W2.SetLineColor(6);
    hQE_neutron_W2.SetLineColor(9);
    hInelastic_W2.SetLineColor(7);
    hData_W2.SetLineColor(kBlack);

    h_combined_W2->SetLineWidth(4);
    hQE_proton_W2.SetLineWidth(4);
    hQE_neutron_W2.SetLineWidth(4);
    hInelastic_W2.SetLineWidth(4);
    hData_W2.SetLineWidth(4);

    h_combined_W2->SetFillColorAlpha(19,0.1);
    h_combined_W2->SetFillStyle(3009);
    hQE_proton_W2.SetFillColorAlpha(6,0.5);
    hQE_proton_W2.SetFillStyle(3004);
    hQE_neutron_W2.SetFillColorAlpha(9,0.5);
    hQE_neutron_W2.SetFillStyle(3005);
    hInelastic_W2.SetFillColorAlpha(7,0.5);
    hInelastic_W2.SetFillStyle(3009);

    hData_W2.SetMarkerStyle(kFullCircle);

    hData_W2.Draw("p");
    h_combined_W2->Draw("hist same");
    hQE_proton_W2.Draw("hist same");
    hQE_neutron_W2.Draw("hist same");
    hInelastic_W2.Draw("hist same");    

    C->Print(Form("images/%s/InelasticCorrection_%s.png",kin_,kin_));
    C1->Print(Form("images/%s/InelasticCorrection_W2_%s.png",kin_,kin_));
}

