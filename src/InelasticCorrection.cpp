#include "InelasticCorrection.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TNamed.h>
#include <iostream>
#include <cmath>

using std::cout;

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

    TF1* f = new TF1("fInel",[&](double*x,double*par){
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

// --- main driver --------------------------------------------------------------
void InelasticCorrection::process(TChain& ch, TChain& ch_QE, TChain& ch_inel,
                                  BranchVars& v, BranchVarsSim& vQE, BranchVarsSim& vInel)
{
    //if(c_.inel_W2_L==0 && c_.inel_W2_H==0){
    //    std::cerr << "[InelCorrection] inel_W2_L/H not set\n"; return; }

    // histograms for fit (tight selection)
    TH1D hData ("hData",  "dx data; W2 (GeV^2)", 100, -4.0, 3.0);
    TH1D hQE_proton("hQE_proton", "dx QE sim protons",          100, -4.0, 3.0);
    TH1D hQE_neutron("hQE_neutron", "dx QE sim neutrons",          100, -4.0, 3.0);
    TH1D hInelastic ("hInelastic",  "dx inelastic sim",        100, -4.0, 3.0);
    TH1D hInelastic_proton ("hInelastic_proton",  "dx inelastic sim protons",        100, -4.0, 3.0);
    TH1D hInelastic_neutron ("hInelastic_neutron",  "dx inelastic sim neutrons",        100, -4.0, 3.0);

    // --- loop data
    Long64_t n=ch.GetEntries();
    std::cout << "[InelasticCorrection] looping over " << n << " events\n";
    const Long64_t step     = 100;
    for(Long64_t i=0;i<n;++i){ 
        
        ch.GetEntry(i);
        
        if(rq_ && (!rq_->helicityOK(v.runnum)||!rq_->mollerOK(v.runnum))) continue;
            
        if(v.ntrack<1 || abs(v.vz)>0.27 || v.eHCAL<c_.eHCAL_L || abs((v.ePS+v.eSH)/(v.trP)-1)>0.2 || v.ePS<0.2 ||
            (c_.coin_L>v.coin_time || v.coin_time>c_.coin_H) || (c_.W2_L>v.W2 || v.W2>c_.W2_H) || (c_.dy_L>v.dy || v.dy>c_.dy_H) || 
            abs(v.helicity)!=1) continue; // no dx cut since we are fitting it
        
        hData.Fill(v.dx);

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
    std::cout << "[InelasticCorrection] looping over QE sim " << nentries_QE << " events\n";
    for(Long64_t i=0;i<ch_QE.GetEntries();++i){ 

        ch_QE.GetEntry(i);

        if(/*vQE.ntrack<1 ||*/ abs(vQE.vz)>0.27 || vQE.eHCAL<c_.eHCAL_L || abs((vQE.ePS+vQE.eSH)/(vQE.trP)-1)>0.2 || vQE.ePS<0.2 ||
            (c_.W2_L>vQE.W2 || vQE.W2>c_.W2_H) || (c_.dy_L>vQE.dy || vQE.dy>c_.dy_H)) continue;

        if(vQE.fnucl == 0) hQE_neutron.Fill(vQE.dx,vQE.weight);

        if(vQE.fnucl == 1) hQE_proton.Fill(vQE.dx+0.05,vQE.weight);//+0.05 is for GEN3 its hard coded for now

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
    std::cout << "[InelasticCorrection] looping over Inelastic sim " << nentries_QE << " events\n";
    for(Long64_t i=0;i<ch_inel.GetEntries();++i){ 
        ch_inel.GetEntry(i);
    
        if(/*vInel.ntrack<1 ||*/ abs(vInel.vz)>0.27 || vInel.eHCAL<c_.eHCAL_L || abs((vInel.ePS+vInel.eSH)/(vInel.trP)-1)>0.2 || vInel.ePS<0.2 ||
            (c_.W2_L>vInel.W2 || vInel.W2>c_.W2_H) || (c_.dy_L>vInel.dy || vInel.dy>c_.dy_H)) continue;    

        hInelastic.Fill(vInel.dx+0.4,vInel.weight);//+0.4 is for GEN3 its hard coded for now

        if (vInel.fnucl == 0) hInelastic_neutron.Fill(vInel.dx+0.4,vInel.weight);

        if (vInel.fnucl == 1) hInelastic_proton.Fill(vInel.dx+0.4,vInel.weight);

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

    h_combined->Scale(hData.Integral());

    hQE_proton.Scale(par0*hData.Integral()*1/hQE_proton.Integral());
    hQE_neutron.Scale(par0*par1*hData.Integral()*1/hQE_neutron.Integral());
    hInelastic.Scale(par0*par2*hData.Integral()*1/hInelastic.Integral());

    //frac_  = wI;                 // N_inel / (N_inel+QE)
    //dfrac_ = std::sqrt(wI*wQ/(wI+wQ)/(wI+wQ)/(hD.Integral()));

    double QE_events = hData.Integral(hData.FindBin(c_.dx_L),hData.FindBin(c_.dx_H));
    double inelastic_events = hInelastic.Integral(hInelastic.FindBin(c_.dx_L),hInelastic.FindBin(c_.dx_H));
    double proton_events = hQE_proton.Integral(hQE_proton.FindBin(c_.dx_L),hQE_proton.FindBin(c_.dx_H));
    double inelastic_frac = inelastic_events/QE_events;// this is not correct should remove other fractions before doing this
    double errinelastic_frac = (inelastic_events/QE_events)*sqrt((1/inelastic_events)+(1/QE_events)); // this is not correct should remove other fractions before doing this
    double proton_frac = proton_events/QE_events;
    double errproton_frac = (proton_events/QE_events)*sqrt((1/proton_events)+(1/QE_events));

    // store
    std::ofstream txt(Form("corrections/InelasticCorrection_%s.txt",kin_));
    //TNamed n("inel_fraction", (std::to_string(frac_)+","+std::to_string(dfrac_)).c_str());
    //n.Write(); fout.Close();
    txt<<"par0 = "<< par0 <<"\n";
    txt<<"par1 = "<< par1 <<"\n";
    txt<<"par2 = "<< par2 <<"\n";
    txt<<"inelastic_fraction = "<<inelastic_frac<<"\n";
    txt<<"err_inelastic_fraction = "<<errinelastic_frac<<"\n";
    txt<<"proton_fraction = "<<proton_frac<<"\n";
    txt<<"err_proton_fraction = "<<errproton_frac<<"\n";
    
    txt<<"f_in = "<<inelastic_frac<<"\n";
    txt<<"err_f_in = "<<errinelastic_frac<<"\n";
    txt<<"f_p = "<<proton_frac<<"\n";
    txt<<"err_f_p = "<<errproton_frac<<"\n";
    //below values should be calculated separately, for now they are set to zero
    txt<<"A_in = "<<0.0<<"\n";
    txt<<"err_A_in = "<<0.0<<"\n";
    txt<<"A_p = "<<0.0<<"\n";
    txt<<"err_A_p = "<<0.0<<"\n";


    txt.close();

    std::cout << "[InelasticCorrection] par0 = "<< par0 << std::endl;
    //     << " saved to "<< outFile_ <<"\n";


    ///////////////////plotting and printing////////////////////////////

    TCanvas *C = new TCanvas("c","c",2400,1500);

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

    C->Print(Form("images/InelasticCorrection_%s.png",kin_));

}

