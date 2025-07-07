#include "PionCorrection.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TNamed.h>
#include <iostream>
#include <cmath>

PionCorrection::PionCorrection(const AnalysisCuts& cuts,
                               const RunQuality* rq,
                               const char* kin,
                               const char* rootFile)
    : c_(cuts), rq_(rq),kin_(kin), outFile_(rootFile) {}


TH1D* PionCorrection::performFit(TH1D* h_data, TH1D* h_pion, TH1D* h_QE,
                 double &pion_weight, double &qe_weight)
{
    // 1) Clone each histogram so we can re-scale them for shape fitting
    //    without affecting the originals.
    TH1D* h_data_norm = (TH1D*)h_data->Clone("h_data_norm");
    TH1D* h_pion_norm = (TH1D*)h_pion->Clone("h_pion_norm");
    TH1D* h_QE_norm   = (TH1D*)h_QE->Clone("h_QE_norm");

    // 2) Normalize each histogram to unity (so we are doing a shape-only fit).
    h_data_norm->Scale(1.0 / h_data_norm->Integral());
    h_pion_norm->Scale(1.0 / h_pion_norm->Integral());
    h_QE_norm->Scale(1.0 / h_QE_norm->Integral());

    // 3) Define a TF1 that returns pion_weight * (pion shape) + qe_weight * (QE shape).
    TF1* fitFunction = new TF1("fitFunction",
        [h_pion_norm, h_QE_norm](double *x, double *par) {
            double pion_w   = par[0];
            double qe_w     = par[1];
            int bin         = h_pion_norm->FindBin(x[0]);
            double pion_val = h_pion_norm->GetBinContent(bin);
            double qe_val   = h_QE_norm->GetBinContent(bin);
            return pion_w * pion_val + qe_w * qe_val;
        },
        0.01, 5, 2);

    fitFunction->SetParNames("Pion Weight", "QE Weight");
    fitFunction->SetParameters(0.5, 0.5); // initial guesses
    fitFunction->SetParLimits(0, 0, 1);
    fitFunction->SetParLimits(1, 0, 1);

    // 4) Fit the normalized data histogram using the shape model.
    h_data_norm->Fit(fitFunction, "RQ"); // R=Range-limit, Q=Quiet

    // 5) Retrieve the fitted weights.
    pion_weight = fitFunction->GetParameter(0);
    qe_weight   = fitFunction->GetParameter(1);

    // 6) Construct a shape-only combined fit histogram:
    //    combined_fit(bin) = pion_weight*pion_norm(bin) + qe_weight*QE_norm(bin).
    TH1D* h_combined_fit = (TH1D*)h_pion_norm->Clone("h_combined_fit");
    h_combined_fit->Reset();
    for (int bin = 1; bin <= h_combined_fit->GetNbinsX(); ++bin) {
        double val = pion_weight * h_pion_norm->GetBinContent(bin)
                   + qe_weight   * h_QE_norm->GetBinContent(bin);
        h_combined_fit->SetBinContent(bin, val);
    }

    // 7) At this point, h_combined_fit is also normalized to area = 1
    //    (just like h_data_norm). We return it for you to scale as needed.
    return h_combined_fit;
}


void PionCorrection::process(TChain& ch, TChain& ch_QE_sim, TChain& ch_pim_sim, BranchVars& v, 
    BranchVarsSim& vQE, BranchVarsSim& vPim)
{
    if (c_.pion_L == 0 && c_.pion_H == 0) {
        std::cerr << "[PionCorrection] Warning: pion_L/H not defined — skipping." << std::endl;
        return;
    }

    TH1D *h_PSe_pion = new TH1D("h_PSe_pion",  "Preshower Energy (pion sim) ; Energy(GeV)",  200, 0.01, 2);
    TH1D *h_PSe_QE   = new TH1D("h_PSe_QE",    "Preshower Energy (QE sim) ; Energy(GeV)",      200, 0.01, 2);
    TH1D *h_PSe_pion_loose_cuts = new TH1D("h_PSe_pion_loose_cuts",  "Preshower Energy (pion sim) ; Energy(GeV)",  200, 0.01, 2);
    TH1D *h_PSe_QE_loose_cuts   = new TH1D("h_PSe_QE_loose_cuts",    "Preshower Energy (QE sim) ; Energy(GeV)",      200, 0.01, 2);
    TH1D *h_PSe_data = new TH1D("h_PSe_data",  "Preshower Energy ; Energy (GeV)",    200, 0.01, 2);
    TH1D *h_PSe_data_pos = new TH1D("h_PSe_data_pos","Preshower Energy (Helicity +1) ; Energy (GeV)", 200,0.01,2);
    TH1D *h_PSe_data_neg = new TH1D("h_PSe_data_neg","Preshower Energy (Helicity -1) ; Energy (GeV)", 200,0.01,2);

    TH1D *h_PSe_data_grinch = new TH1D("h_PSe_data_grinch",  "Preshower Energy (grinch cuts) ; Energy (GeV)",    200, 0.01, 2);

    double Ngrinch_pos = 0;
    double Ngrinch_neg = 0;

    //////////////// data tree /////////////////////////

    Long64_t nentries = ch.GetEntries();
    const Long64_t step     = 100;
    std::cout << "[PionCorrection] looping over " << nentries << " events\n";

    for (Long64_t i=0;i<nentries;++i){
        ch.GetEntry(i);
        ////////////////////////run quality checks//////////////////////
        if (rq_ && (!rq_->helicityOK(v.runnum) || !rq_->mollerOK(v.runnum))) continue;

        ////////////////////////loose cuts for asymmetry calculation/////////////
        if (v.ntrack<1 || abs(v.vz)>0.27 ||/* v.eHCAL<c_.eHCAL_L || abs((v.ePS+v.eSH)/(v.trP)-1)>0.2 ||
            */(c_.coin_L>v.coin_time || v.coin_time>c_.coin_H) || (c_.W2_L>v.W2 || v.W2>c_.W2_H) || /*(c_.dx_L>v.dx || v.dx>c_.dx_H) || (c_.dy_L>v.dy || v.dy>c_.dy_H) ||*/ 
            abs(v.helicity)!=1) continue; //no ePS since we are looking at pions

        if(v.helicity==1) h_PSe_data_pos->Fill(v.ePS);
            
        if(v.helicity==-1) h_PSe_data_neg->Fill(v.ePS);

        //grinch for verification

        if(v.grinch_track!=0 && v.grinch_clus_size<2){ // this could be an or 
            h_PSe_data_grinch->Fill(v.ePS);

            if(v.helicity==1) Ngrinch_pos++;
            if(v.helicity==-1) Ngrinch_neg++;
        
        }


        ///////////////////////tight cuts for fraction calculation///////////////
        if (v.ntrack<1 || abs(v.vz)>0.27 || v.eHCAL<c_.eHCAL_L || abs((v.ePS+v.eSH)/(v.trP)-1)>0.2 ||
            (c_.coin_L>v.coin_time || v.coin_time>c_.coin_H) || (c_.W2_L>v.W2 || v.W2>c_.W2_H) || (c_.dx_L>v.dx || v.dx>c_.dx_H) || (c_.dy_L>v.dy || v.dy>c_.dy_H) || 
            abs(v.helicity)!=1) continue; //no ePS since we are looking at pions

        h_PSe_data->Fill(v.ePS);


        if (v.ePS<c_.pion_L){ 
            if (v.helicity > 0) ++Np_;
            else if (v.helicity < 0) ++Nm_;
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
    
    /////////////////// QE Sim ///////////////////////////////////////

    Long64_t nentries_QE = ch_QE_sim.GetEntries();
    std::cout << "[PionCorrection] looping over QE sim " << nentries_QE << " events\n";

    for (Long64_t i=0;i<nentries_QE;++i){
        ch_QE_sim.GetEntry(i);

        ///////////////loose cuts for asymmetry calculation/////////////////

        if (abs(vQE.vz)>0.27 || vQE.eHCAL<0.025 /*|| abs((vQE.ePS+vQE.eSH)/(vQE.trP)-1)>0.2*/ ||
           (c_.W2_L>vQE.W2 || vQE.W2>c_.W2_H) /*|| (c_.dx_L>vQE.dx || vQE.dx>c_.dx_H) || (c_.dy_L>vQE.dy || vQE.dy>c_.dy_H)*/) continue;

        h_PSe_QE_loose_cuts->Fill(vQE.ePS,vQE.weight);

        ///////////////tight cuts for fraction calculation/////////////////
        if (abs(vQE.vz)>0.27 || vQE.eHCAL<c_.eHCAL_L || abs((vQE.ePS+vQE.eSH)/(vQE.trP)-1)>0.2 ||
           (c_.W2_L>vQE.W2 || vQE.W2>c_.W2_H) || (c_.dx_L>vQE.dx || vQE.dx>c_.dx_H) || (c_.dy_L>vQE.dy || vQE.dy>c_.dy_H)) continue; //no ePS since we are looking at pions

        h_PSe_QE->Fill(vQE.ePS,vQE.weight);

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

    //////////////////////////////pim sim/////////////////////////////////////

    Long64_t nentries_pim = ch_pim_sim.GetEntries();
    std::cout << "[PionCorrection] looping over pim sim " << nentries_pim << " events\n";

    for (Long64_t i=0;i<nentries_pim;++i){
        ch_pim_sim.GetEntry(i);

        //////////////loose cuts for asymmetry calculation//////////////////
        if (abs(vPim.vz)>0.27 || /*vPim.eHCAL<c_.eHCAL_L || abs((vPim.ePS+vPim.eSH)/(vPim.trP)-1)>0.2 ||*/
           (c_.W2_L>vPim.W2 || vPim.W2>c_.W2_H) /*|| (c_.dx_L>vPim.dx || vPim.dx>c_.dx_H) || (c_.dy_L>vPim.dy || vPim.dy>c_.dy_H)*/) continue; 

        h_PSe_pion_loose_cuts->Fill(vPim.ePS,vPim.weight);

        ///////////////tight cuts for fraction calculation//////////////////
        if (abs(vPim.vz)>0.27 || /*vPim.eHCAL<c_.eHCAL_L ||*/ abs((vPim.ePS+vPim.eSH)/(vPim.trP)-1)>0.2 ||
           (c_.W2_L>vPim.W2 || vPim.W2>c_.W2_H) || (c_.dx_L>vPim.dx || vPim.dx>c_.dx_H) || (c_.dy_L>vPim.dy || vPim.dy>c_.dy_H)) continue; //no ePS since we are looking at pions

        h_PSe_pion->Fill(vPim.ePS,vPim.weight);

        // progress bar
        if (i % step == 0 || i == nentries_pim - 1) {
            double frac = double(i + 1) / nentries_pim;
            int barw = 42, pos = static_cast<int>(barw * frac);
            std::cout << '\r' << '[';
            for (int j = 0; j < barw; ++j)
                std::cout << (j < pos ? '=' : (j == pos ? '>' : ' '));
            std::cout << "] " << static_cast<int>(frac * 100) << " %" << std::flush;
        }
    }

    double pion_weight_all=0, qe_weight_all=0;
    TH1D* h_combined_all = performFit(h_PSe_data, h_PSe_pion, h_PSe_QE,
                                      pion_weight_all, qe_weight_all);

    double pion_weight_pos=0, qe_weight_pos=0;
    TH1D* h_combined_pos = performFit(h_PSe_data_pos, h_PSe_pion_loose_cuts, h_PSe_QE_loose_cuts,
                                      pion_weight_pos, qe_weight_pos);

    double pion_weight_neg=0, qe_weight_neg=0;
    TH1D* h_combined_neg = performFit(h_PSe_data_neg, h_PSe_pion_loose_cuts, h_PSe_QE_loose_cuts,
                                      pion_weight_neg, qe_weight_neg);

    // For "all data" 
    double dataIntegral_all = h_PSe_data->Integral();
    double dataIntegral_pos = h_PSe_data_pos->Integral();
    double dataIntegral_neg = h_PSe_data_neg->Integral();

    // Create scaled copies of the MC shape so we can superimpose them on data:
    TH1D* h_pion_scaled_all = (TH1D*)h_PSe_pion->Clone("h_pion_scaled_all");
    TH1D* h_QE_scaled_all   = (TH1D*)h_PSe_QE->Clone("h_QE_scaled_all");
    TH1D* h_combined_scaled_all = (TH1D*)h_combined_all->Clone("h_combined_scaled_all");

    TH1D* h_pion_scaled_pos = (TH1D*)h_PSe_pion_loose_cuts->Clone("h_pion_scaled_pos");
    TH1D* h_QE_scaled_pos   = (TH1D*)h_PSe_QE_loose_cuts->Clone("h_QE_scaled_pos");
    TH1D* h_combined_scaled_pos = (TH1D*)h_combined_pos->Clone("h_combined_scaled_pos");
    
    TH1D* h_pion_scaled_neg = (TH1D*)h_PSe_pion_loose_cuts->Clone("h_pion_scaled_neg");
    TH1D* h_QE_scaled_neg   = (TH1D*)h_PSe_QE_loose_cuts->Clone("h_QE_scaled_neg");
    TH1D* h_combined_scaled_neg = (TH1D*)h_combined_neg->Clone("h_combined_scaled_neg");

    // Each one is area=1 originally (from performFit), so multiply by dataIntegral_all:
    h_pion_scaled_all     ->Scale(pion_weight_all * dataIntegral_all/h_pion_scaled_all->Integral());
    h_QE_scaled_all       ->Scale(qe_weight_all   * dataIntegral_all/h_QE_scaled_all->Integral());
    h_combined_scaled_all ->Scale(dataIntegral_all); // The combined sum is 1, times total data.

    h_pion_scaled_pos     ->Scale(pion_weight_pos * dataIntegral_pos/h_pion_scaled_pos->Integral());
    h_QE_scaled_pos       ->Scale(qe_weight_pos   * dataIntegral_pos/h_QE_scaled_pos->Integral());
    h_combined_scaled_pos ->Scale(dataIntegral_pos); // The combined sum is 1, times total data.

    h_pion_scaled_neg     ->Scale(pion_weight_neg * dataIntegral_neg/h_pion_scaled_neg->Integral());
    h_QE_scaled_neg       ->Scale(qe_weight_neg   * dataIntegral_neg/h_QE_scaled_neg->Integral());
    h_combined_scaled_neg ->Scale(dataIntegral_neg); // The combined sum is 1, times total data.

    long long sum = Np_ + Nm_;
    asym_ = (sum>0)? double(Np_ - Nm_) / sum : 0.0;
    err_  = (sum>0)? 2.0*std::sqrt(double(Np_)*Nm_)/(sum*sum) : 0.0;

    double Cpi_pos = pion_weight_pos;
    double Cpi_neg = pion_weight_neg;

    double pos_integral = h_pion_scaled_pos->Integral();
    double neg_integral = h_pion_scaled_neg->Integral();

    double asymCpi = (Cpi_pos - Cpi_neg)/(Cpi_pos + Cpi_neg);
    double err_asymCpi = 2*sqrt(Cpi_pos*Cpi_neg/pow(Cpi_pos+Cpi_neg,3));

    double asymIntegral = (pos_integral - neg_integral)/(pos_integral + neg_integral);
    double err_asymIntegral = 2*sqrt(pos_integral*neg_integral/pow(pos_integral+neg_integral,3));

    double asymGrinch = (Ngrinch_pos - Ngrinch_neg)/(Ngrinch_pos + Ngrinch_neg);
    double err_asymGrinch = 2*sqrt(Ngrinch_pos*Ngrinch_neg/pow((Ngrinch_pos+Ngrinch_neg),3));

    double f_pi = 0.0;
    double err_f_pi = 0.0;

    double pi_events_all = h_pion_scaled_all->Integral(h_pion_scaled_all->FindBin(c_.pion_L),h_pion_scaled_all->FindBin(c_.pion_H));
    double QE_events_all = h_PSe_data->Integral(h_PSe_data->FindBin(c_.pion_L),h_PSe_data->FindBin(c_.pion_H));

    if( pi_events_all<1) {
        f_pi = 0.0;
        err_f_pi = 0.0;
    }
    else{
        f_pi = pi_events_all/QE_events_all;

        err_f_pi = (pi_events_all/QE_events_all)*sqrt((1/pi_events_all)+(1/QE_events_all));
    }

    // ROOT output: store as TNamed (name:value,error)
    std::ofstream f(Form("corrections/PionCorrection_%s.txt",kin_)) /*outFile_.c_str(), "RECREATE")*/;
    //std::string payload = std::to_string(asym_) + "," + std::to_string(err_);
    //TNamed n("pion_asym", payload.c_str());
    //n.Write();
    f<< "Cpi_pos = "<< Cpi_pos <<"\n";
    f<< "Cpi_neg = "<< Cpi_neg <<"\n";
    f<< "pos_integral = "<< pos_integral <<"\n";
    f<< "neg_integral = "<< neg_integral <<"\n";
    f<< "asymCpi = "<<asymCpi<<"\n";
    f<< "err_asymCpi = "<<err_asymCpi<<"\n";
    f<< "asymIntegral = "<<asymIntegral<<"\n";
    f<< "err_asymIntegral = "<<err_asymIntegral<<"\n";
    f<< "asymGrinch = "<<asymGrinch<<"\n";
    f<< "err_asymGrinch = "<<err_asymGrinch<<"\n";
    f<< "A_pi = "<<asymGrinch<<"\n";
    f<< "err_A_pi = "<<err_asymGrinch<<"\n";
    f<< "f_pi = "<<f_pi<<"\n";
    f<< "err_f_pi = "<<err_f_pi<<"\n";

    f.close();

    TFile f1(Form("rootfiles/pion_correction_%s.root",kin_) /*outFile_.c_str()*/, "RECREATE");
    h_PSe_data->Write();
    h_pion_scaled_all->Write();
    h_QE_scaled_all->Write();
    h_combined_scaled_all->Write();

    h_PSe_data_pos->Write();
    h_pion_scaled_pos->Write();
    h_QE_scaled_pos->Write();
    h_combined_scaled_pos->Write();

    h_PSe_data_neg->Write();
    h_pion_scaled_neg->Write();
    h_QE_scaled_neg->Write();
    h_combined_scaled_neg->Write();

    h_PSe_data_grinch->Write();

    f1.Close();

    std::cout << "[PionCorrection] raw pion asym = " << asym_
              << " ± " << err_ << " written to " << outFile_ << "\n";


    /////////////////// Canvas and Printing ////////////////////////
    TCanvas *C = new TCanvas("c","c",2400,1500);

    C->Divide(2,2);
    C->cd(1);
    h_combined_scaled_all->SetLineColor(kGreen);
    h_PSe_data->SetLineColor(kBlack);
    h_QE_scaled_all->SetLineColor(kBlue);
    h_pion_scaled_all->SetLineColor(kRed);

    h_combined_scaled_all->SetLineWidth(4);
    h_PSe_data->SetLineWidth(4);
    h_QE_scaled_all->SetLineWidth(4);
    h_pion_scaled_all->SetLineWidth(4);

    h_PSe_data->SetMarkerStyle(kFullCircle);

    h_combined_scaled_all->Draw("hist");
    h_PSe_data->Draw("same p");
    h_pion_scaled_all->Draw("same hist");
    h_QE_scaled_all->Draw("same hist");

    C->cd(2);
    h_combined_scaled_pos->SetLineColor(kGreen);
    h_PSe_data_pos->SetLineColor(kBlack);
    h_QE_scaled_pos->SetLineColor(kBlue);
    h_pion_scaled_pos->SetLineColor(kRed);

    h_combined_scaled_pos->SetLineWidth(4);
    h_PSe_data_pos->SetLineWidth(4);
    h_QE_scaled_pos->SetLineWidth(4);
    h_pion_scaled_pos->SetLineWidth(4);

    h_PSe_data_pos->SetMarkerStyle(kFullCircle);

    h_combined_scaled_pos->Draw("hist");
    h_PSe_data_pos->Draw("same p");
    h_pion_scaled_pos->Draw("same hist");
    h_QE_scaled_pos->Draw("same hist");

    C->cd(3);
    h_combined_scaled_neg->SetLineColor(kGreen);
    h_PSe_data_neg->SetLineColor(kBlack);
    h_QE_scaled_neg->SetLineColor(kBlue);
    h_pion_scaled_neg->SetLineColor(kRed);

    h_combined_scaled_neg->SetLineWidth(4);
    h_PSe_data_neg->SetLineWidth(4);
    h_QE_scaled_neg->SetLineWidth(4);
    h_pion_scaled_neg->SetLineWidth(4);

    h_PSe_data_neg->SetMarkerStyle(kFullCircle);

    h_combined_scaled_neg->Draw("hist");
    h_PSe_data_neg->Draw("same p");
    h_pion_scaled_neg->Draw("same hist");
    h_QE_scaled_neg->Draw("same hist");

    C->cd(4);
    h_PSe_data_grinch->Draw();

    C->Print(Form("images/PionPlots_%s.png",kin_));

}
