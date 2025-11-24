#include "PionCorrection.h"
#include "Utility.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TNamed.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLatex.h>
#include <string>

using std::string; using std::vector; using std::unordered_map;

PionCorrection::PionCorrection(const AnalysisCuts& cuts,
                               const RunQuality* rq,
                               const char* kin,
                               const char* rootFile)
    : c_(cuts), rq_(rq),kin_(kin), outFile_(rootFile) {}


TH1D* PionCorrection::performFit(TH1D* h_data, TH1D* h_pion, TH1D* h_QE,
                 double &pion_weight, double &qe_weight,
                 double &pion_weight_err, double &qe_weight_err)
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
    // Request the fit and quietly retrieve parameter values/errors.
    // Using option "RQ" (Range, Quiet) — TF1 stores parameter errors accessible
    // via GetParError after the fit.
    h_data_norm->Fit(fitFunction, "RQ"); // R=Range-limit, Q=Quiet

    // 5) Retrieve the fitted weights and their 1-sigma errors.
    pion_weight = fitFunction->GetParameter(0);
    qe_weight   = fitFunction->GetParameter(1);
    pion_weight_err = fitFunction->GetParError(0);
    qe_weight_err   = fitFunction->GetParError(1);

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
    Utility utility;

    auto   corrFile = [&](const char* stem){
        return string("corrections/")+kin_+"/" + stem + "Correction_" + kin_ + ".txt"; };

    const string accFile = corrFile("Accidental");
    const string nitFile = corrFile("Nitrogen");

    // -------- read TXT correction files directly -------------------------------
    const auto accMap  = utility.readSimpleKVGlobal(accFile);
    const auto nitMap  = utility.readSimpleKVGlobal(nitFile);

    auto V=[&](const auto& M,const char* k){ auto it=M.find(k); return it!=M.end()? it->second : 0.0; };

    // central values
    double Aacc = V(accMap ,"A_acc");           // not provided yet

    double facc = V(accMap ,"f_acc");
    double fN2  = V(nitMap ,"f_N2");

    // errors
    double errAacc = V(accMap ,"err_A_acc");

    double errfacc = V(accMap ,"err_f_acc");
    double errfN2  = V(nitMap ,"err_f_N2");


    if (c_.pion_L == 0 && c_.pion_H == 0) {
        std::cerr << "[PionCorrection] Warning: pion_L/H not defined — skipping." << std::endl;
        return;
    }

    std::cout<< "pion_L" << c_.pion_L<<"\n";
    std::cout<< "pion_H" << c_.pion_H<<"\n";

    TH1D *h_PSe_pion = new TH1D("h_PSe_pion",  "Preshower Energy (GeV) ; Energy(GeV)",  200, 0.01, 2.5);
    TH1D *h_PSe_QE   = new TH1D("h_PSe_QE",    "Preshower Energy (QE sim) ; Energy(GeV)",      200, 0.01, 2.5);
    TH1D *h_PSe_pion_loose_cuts = new TH1D("h_PSe_pion_loose_cuts",  "Preshower Energy (GeV) ; Energy(GeV)",  200, 0.01, 2.5);
    TH1D *h_PSe_QE_loose_cuts   = new TH1D("h_PSe_QE_loose_cuts",    "Preshower Energy (QE sim) ; Energy(GeV)",      200, 0.01, 2.5);
    TH1D *h_PSe_data_loose_cuts = new TH1D("h_PSe_data_loose_cuts",  "Preshower Energy (relaxed cuts) ; Energy (GeV)",    200, 0.01, 2.5);
    TH1D *h_PSe_data = new TH1D("h_PSe_data",  "Preshower Energy (analysis cuts); Energy (GeV)",    200, 0.01, 2.5);
    TH1D *h_PSe_data_pos = new TH1D("h_PSe_data_pos","Preshower Energy (Helicity +1) ; Energy (GeV)", 200,0.01,2.5);
    TH1D *h_PSe_data_neg = new TH1D("h_PSe_data_neg","Preshower Energy (Helicity -1) ; Energy (GeV)", 200,0.01,2.5);

    TH1D *h_PSe_data_grinch = new TH1D("h_PSe_data_grinch",  "Preshower Energy (grinch cuts) ; Energy (GeV)",    200, 0.01, 2.5);
    TH1D *h_PSe_data_antigrinch = new TH1D("h_PSe_data_antigrinch",  "Preshower Energy (anti-grinch cuts) ; Energy (GeV)",    200, 0.01, 2.5);

    TH1D *h_PSe_data_grinch_loose_cuts = new TH1D("h_PSe_data_grinch_loose_cuts",  "Preshower Energy (grinch cuts) ; Energy (GeV)",    200, 0.01, 2.5);
    TH1D *h_PSe_data_antigrinch_loose_cuts = new TH1D("h_PSe_data_antigrinch_loose_cuts",  "Preshower Energy (anti-grinch cuts) ; Energy (GeV)",    200, 0.01, 2.5);

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

        h_PSe_data_loose_cuts->Fill(v.ePS);

        if(v.helicity==1) h_PSe_data_pos->Fill(v.ePS);
            
        if(v.helicity==-1) h_PSe_data_neg->Fill(v.ePS);


        //grinch for verification
        if(v.grinch_track!=0 || v.grinch_clus_size<2){ // this could be an or 
            h_PSe_data_grinch_loose_cuts->Fill(v.ePS);

            if(v.helicity==1) Ngrinch_pos++;
            if(v.helicity==-1) Ngrinch_neg++;
        
        }else{
            h_PSe_data_antigrinch_loose_cuts->Fill(v.ePS);
        }

        ///////////////////////tight cuts for fraction calculation///////////////
        if (v.ntrack<1 || abs(v.vz)>0.27 || v.eHCAL<c_.eHCAL_L || abs((v.ePS+v.eSH)/(v.trP)-1)>0.2 ||
            (c_.coin_L>v.coin_time || v.coin_time>c_.coin_H) || (c_.W2_L>v.W2 || v.W2>c_.W2_H) || ((pow((v.dy-0.0)/0.4,2)+pow((v.dx-0.0)/0.4,2))>1) || 
            abs(v.helicity)!=1) continue; //no ePS since we are looking at pions

        h_PSe_data->Fill(v.ePS);

        //grinch for verification

        if(v.grinch_track!=0 || v.grinch_clus_size<2){ // this could be an or 
            h_PSe_data_grinch->Fill(v.ePS);
        
        }else{
            h_PSe_data_antigrinch->Fill(v.ePS);
        }


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
           (c_.W2_L>vQE.W2 || vQE.W2>c_.W2_H) || ((pow((vQE.dy-0.0)/0.4,2)+pow((vQE.dx-0.0)/0.4,2))>1)) continue;

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
           (c_.W2_L>vPim.W2 || vPim.W2>c_.W2_H) /*|| ((pow((v.dy-0.0)/0.4,2)+pow((v.dx-0.0)/0.4,2))>1)*/) continue;

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
    double pion_weight_all_err=0, qe_weight_all_err=0;
    TH1D* h_combined_all = performFit(h_PSe_data, h_PSe_pion, h_PSe_QE,
                                      pion_weight_all, qe_weight_all,
                                      pion_weight_all_err, qe_weight_all_err);

    double pion_weight_pos=0, qe_weight_pos=0;
    double pion_weight_pos_err=0, qe_weight_pos_err=0;
    TH1D* h_combined_pos = performFit(h_PSe_data_pos, h_PSe_pion_loose_cuts, h_PSe_QE_loose_cuts,
                                      pion_weight_pos, qe_weight_pos,
                                      pion_weight_pos_err, qe_weight_pos_err);

    double pion_weight_neg=0, qe_weight_neg=0;
    double pion_weight_neg_err=0, qe_weight_neg_err=0;
    TH1D* h_combined_neg = performFit(h_PSe_data_neg, h_PSe_pion_loose_cuts, h_PSe_QE_loose_cuts,
                                      pion_weight_neg, qe_weight_neg,
                                      pion_weight_neg_err, qe_weight_neg_err);

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
    if (sum > 0) {
        const double sum_d = double(sum);
        // Correct propagation for A = (N+ - N-)/N where N=N+ + N-
        // sigma_A = 2*sqrt(N+ * N-) / N^(3/2)
        err_ = 2.0 * std::sqrt(double(Np_) * double(Nm_)) / std::pow(sum_d, 1.5);
    } else {
        err_ = 0.0;
    }

    double Cpi_pos = pion_weight_pos;
    double Cpi_neg = pion_weight_neg;

    double pos_integral = h_pion_scaled_pos->Integral();
    double neg_integral = h_pion_scaled_neg->Integral();

    double asymCpi = 0.0;
    double err_asymCpi = 0.0;
    double denomC = (Cpi_pos + Cpi_neg);
    if (denomC > 0) {
        asymCpi = (Cpi_pos - Cpi_neg) / denomC;
        // propagate parameter errors from the fit. We neglect the covariance term
        // between the two fit parameters here (could be added if available).
        double dCpos = 2.0 * Cpi_neg / (denomC * denomC);
        double dCneg = -2.0 * Cpi_pos / (denomC * denomC);
        err_asymCpi = std::sqrt(dCpos*dCpos * pion_weight_pos_err*pion_weight_pos_err
                              + dCneg*dCneg * pion_weight_neg_err*pion_weight_neg_err);
    }

    double asymIntegral = (pos_integral - neg_integral)/(pos_integral + neg_integral);
    double err_asymIntegral = 2*sqrt(pos_integral*neg_integral/pow(pos_integral+neg_integral,3));

    double asymGrinch = (Ngrinch_pos - Ngrinch_neg)/(Ngrinch_pos + Ngrinch_neg);
    double err_asymGrinch = 2*sqrt(Ngrinch_pos*Ngrinch_neg/pow((Ngrinch_pos+Ngrinch_neg),3));

    double f_pi = 0.0;
    double err_f_pi = 0.0;

    double pi_events_all = h_pion_scaled_all->Integral(h_pion_scaled_all->FindBin(c_.pion_L),h_pion_scaled_all->FindBin(c_.pion_H));
    double QE_events_all = h_PSe_data->Integral(h_PSe_data->FindBin(c_.pion_L),h_PSe_data->FindBin(c_.pion_H));

    // full propagation for f_pi = (A * C) / B, with A = pi_events_all,
    // B = QE_events_all, C = (1 - facc - fN2).
    if (pi_events_all <= 0.0 || QE_events_all <= 0.0) {
        f_pi = 0.0;
        err_f_pi = 0.0;
    } else {
        const double A = pi_events_all;
        const double B = QE_events_all;
        const double C = 1.0 - facc - fN2;

        f_pi = (A * C) / B;

        // Uncertainties: sigma_A = sqrt(A) (Poisson), sigma_B = sqrt(B),
        // sigma_C = sqrt(err_facc^2 + err_fN2^2) (assume independent)
        const double sigmaA = std::sqrt(std::max(0.0, A));
        const double sigmaB = std::sqrt(std::max(0.0, B));
        const double sigmaC = std::sqrt(errfacc*errfacc + errfN2*errfN2);

        // Propagate using partial derivatives:
        // df/dA = C / B
        // df/dC = A / B
        // df/dB = -A*C / B^2
        const double termA = (C / B) * sigmaA;
        const double termC = (A / B) * sigmaC;
        const double termB = (A * C / (B*B)) * sigmaB;

        err_f_pi = std::sqrt(termA*termA + termC*termC + termB*termB);
    }

    // ROOT output: store as TNamed (name:value,error)
    std::ofstream f(Form("corrections/%s/PionCorrection_%s.txt",kin_,kin_)) /*outFile_.c_str(), "RECREATE")*/;
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
    f<< "pi_events_all = "<<pi_events_all<<"\n";
    f<< "QE_events_all = "<<QE_events_all<<"\n";
    f<< "A_pi = "<<asymGrinch<<"\n";
    f<< "err_A_pi = "<<err_asymGrinch<<"\n";
    f<< "f_pi = "<<f_pi<<"\n";
    f<< "err_f_pi = "<<err_f_pi<<"\n";

    f.close();

    TFile f1(Form("rootfiles/%s/pion_correction_%s.root",kin_,kin_) /*outFile_.c_str()*/, "RECREATE");
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

    // Diagnostic output: show integrals so an empty canvas can be investigated
    std::cout << "[PionCorrection] integrals: data=" << dataIntegral_all
              << " combined=" << h_combined_scaled_all->Integral()
              << " pion_scaled=" << h_pion_scaled_all->Integral()
              << " qe_scaled=" << h_QE_scaled_all->Integral() << "\n";
    if (h_combined_scaled_all->Integral() <= 0.0) {
        std::cerr << "[PionCorrection] Warning: combined fit histogram is empty. "
                  << "This usually means no events passed the tight selection or the fit failed." << std::endl;
    }

    gStyle->SetOptStat(0);
    /////////////////// Canvas and Printing ////////////////////////
    TCanvas *C = new TCanvas("c","c",2400,1500);
    TCanvas *C1 = new TCanvas("c1","c1",2400,1500);
    TCanvas *C2 = new TCanvas("c2","c2",2400,1500);

    C->Divide(1,1);
    C->cd(1);
    // If combined histogram is empty, draw a helpful message instead of blank canvas.
    if (h_combined_scaled_all->Integral() > 0.0) {
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
    } else {
        TH1D *h_dummy = new TH1D("h_dummy_no_entries", "No entries after selection", 10, 0, 1);
        h_dummy->SetStats(0);
        h_dummy->GetXaxis()->SetTitle("preshower energy");
        h_dummy->GetYaxis()->SetTitle("counts");
        h_dummy->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.03);
        t.DrawLatex(0.12, 0.8, Form("dataIntegral_all = %.0f", dataIntegral_all));
        t.DrawLatex(0.12, 0.75, Form("pion_weight = %.4g  qe_weight = %.4g", pion_weight_all, qe_weight_all));
        t.DrawLatex(0.12, 0.70, "Check tight selection windows (eHCAL, dx/dy, W2, coin_time) and fit status.");
    }

    TLegend* leg = new TLegend(0.55, 0.62, 0.88, 0.88); // x1,y1,x2,y2 (NDC)
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    leg->AddEntry(h_PSe_data,           "Data (PS e)",   "p"); // points
    leg->AddEntry(h_QE_scaled_all,      "QE ",   "l"); // blue line
    leg->AddEntry(h_pion_scaled_all,    "#pi^{-} background","l"); // red line
    leg->AddEntry(h_combined_scaled_all,"Fit (C_{#pi^{-}}S_{#pi^{-}} + C_{e}S_{e})","l"); // green line

    leg->Draw();
    gPad->Update();

    gPad->Update(); // make sure pad limits are known
    double yminPS = gPad->GetUymin();
    double ymaxPS = gPad->GetUymax();

    auto cutLinePS = new TLine(0.2, yminPS, 0.2, ymaxPS);
    cutLinePS->SetLineStyle(1); //1-solid 2-dashed
    cutLinePS->SetLineWidth(4);
    cutLinePS->SetLineColor(kMagenta);
    cutLinePS->Draw();
    gPad->RedrawAxis();

    C1->Divide(2,1);

    C1->cd(1);
    h_combined_scaled_pos->SetLineColor(kGreen);
    h_PSe_data_pos->SetLineColor(kBlack);
    h_QE_scaled_pos->SetLineColor(kBlue);
    h_pion_scaled_pos->SetLineColor(kRed);

    h_combined_scaled_pos->SetLineWidth(4);
    h_PSe_data_pos->SetLineWidth(4);
    h_QE_scaled_pos->SetLineWidth(4);
    h_pion_scaled_pos->SetLineWidth(4);

    h_PSe_data_pos->SetMarkerStyle(kFullCircle);

    h_PSe_data_pos->Draw("p");
    h_combined_scaled_pos->Draw("hist same");
    h_pion_scaled_pos->Draw("same hist");
    h_QE_scaled_pos->Draw("same hist");

    TLegend* leg_p = new TLegend(0.55, 0.62, 0.88, 0.88); // x1,y1,x2,y2 (NDC)
    leg_p->SetBorderSize(0);
    leg_p->SetFillStyle(0);
    leg_p->SetTextSize(0.035);

    leg_p->AddEntry(h_PSe_data_pos,           "Data (PS e)",   "p"); // points
    leg_p->AddEntry(h_QE_scaled_pos,      "QE ",   "l"); // blue line
    leg_p->AddEntry(h_pion_scaled_pos,    "#pi^{-} background","l"); // red line
    leg_p->AddEntry(h_combined_scaled_pos,"Fit (C_{#pi^{-}}^{+} S_{#pi^{-}} + C_{e}^{+} S_{e})","l"); // green line

    leg_p->Draw();
    gPad->Update();

    gPad->Update(); // make sure pad limits are known
    double yminPSp = gPad->GetUymin();
    double ymaxPSp = gPad->GetUymax();

    auto cutLinePSp = new TLine(0.2, yminPSp, 0.2, ymaxPSp);
    cutLinePSp->SetLineStyle(1); //1-solid 2-dashed
    cutLinePSp->SetLineWidth(4);
    cutLinePSp->SetLineColor(kMagenta);
    cutLinePSp->Draw();
    gPad->RedrawAxis();

    C1->cd(2);
    h_combined_scaled_neg->SetLineColor(kGreen);
    h_PSe_data_neg->SetLineColor(kBlack);
    h_QE_scaled_neg->SetLineColor(kBlue);
    h_pion_scaled_neg->SetLineColor(kRed);

    h_combined_scaled_neg->SetLineWidth(4);
    h_PSe_data_neg->SetLineWidth(4);
    h_QE_scaled_neg->SetLineWidth(4);
    h_pion_scaled_neg->SetLineWidth(4);

    h_PSe_data_neg->SetMarkerStyle(kFullCircle);

    h_PSe_data_neg->Draw("p");
    h_combined_scaled_neg->Draw("hist same");
    h_pion_scaled_neg->Draw("same hist");
    h_QE_scaled_neg->Draw("same hist");

    TLegend* leg_m = new TLegend(0.55, 0.62, 0.88, 0.88); // x1,y1,x2,y2 (NDC)
    leg_m->SetBorderSize(0);
    leg_m->SetFillStyle(0);
    leg_m->SetTextSize(0.035);

    leg_m->AddEntry(h_PSe_data_neg,           "Data (PS e)",   "p"); // points
    leg_m->AddEntry(h_QE_scaled_neg,      "QE ",   "l"); // blue line
    leg_m->AddEntry(h_pion_scaled_neg,    "#pi^{-} background","l"); // red line
    leg_m->AddEntry(h_combined_scaled_neg,"Fit (C_{#pi^{-}}^{-} S_{#pi^{-}} + C_{e}^{-} S_{e})","l"); // green line

    leg_m->Draw();
    gPad->Update();

    gPad->Update(); // make sure pad limits are known
    double yminPSm = gPad->GetUymin();
    double ymaxPSm = gPad->GetUymax();

    auto cutLinePSm = new TLine(0.2, yminPSm, 0.2, ymaxPSm);
    cutLinePSm->SetLineStyle(1); //1-solid 2-dashed
    cutLinePSm->SetLineWidth(4);
    cutLinePSm->SetLineColor(kMagenta);
    cutLinePSm->Draw();
    gPad->RedrawAxis();

    C2->Divide(2,1);

    C2->cd(1);

    h_PSe_data->SetLineColor(kBlue);
    h_PSe_data_grinch->SetLineColor(kBlack);
    h_PSe_data_antigrinch->SetLineColor(kOrange);
    
    h_PSe_data->SetLineWidth(4);
    h_PSe_data_grinch->SetLineWidth(4);
    h_PSe_data_antigrinch->SetLineWidth(4);

    h_PSe_data->Draw();
    h_PSe_data_grinch->Draw("same");
    h_PSe_data_antigrinch->Draw("same");

    gPad->Update(); // make sure pad limits are known
    double ymin = gPad->GetUymin();
    double ymax = gPad->GetUymax();

    auto legg = new TLegend(0.60, 0.70, 0.88, 0.88); // x1,y1,x2,y2 (NDC)
    legg->SetBorderSize(0);
    legg->SetFillStyle(0);
    legg->SetTextSize(0.035);
    legg->AddEntry(h_PSe_data,          "All data",        "l");
    legg->AddEntry(h_PSe_data_grinch,   "GRINCH", "l");
    legg->AddEntry(h_PSe_data_antigrinch,"anti-GRINCH",     "l");
    legg->Draw();

    auto cutLine = new TLine(0.2, ymin, 0.2, ymax);
    cutLine->SetLineStyle(1); //1-solid 2-dashed
    cutLine->SetLineWidth(4);
    cutLine->SetLineColor(kMagenta);
    cutLine->Draw();
    gPad->RedrawAxis();

    C2->cd(2);
    h_PSe_data_loose_cuts->SetLineColor(kBlue);
    h_PSe_data_grinch_loose_cuts->SetLineColor(kBlack);
    h_PSe_data_antigrinch_loose_cuts->SetLineColor(kOrange);
    
    h_PSe_data_loose_cuts->SetLineWidth(4);
    h_PSe_data_grinch_loose_cuts->SetLineWidth(4);
    h_PSe_data_antigrinch_loose_cuts->SetLineWidth(4);

    h_PSe_data_loose_cuts->Draw();
    h_PSe_data_grinch_loose_cuts->Draw("same");
    h_PSe_data_antigrinch_loose_cuts->Draw("same");

    gPad->Update(); // make sure pad limits are known
    double yminl = gPad->GetUymin();
    double ymaxl = gPad->GetUymax();

    auto leggl = new TLegend(0.60, 0.70, 0.88, 0.88); // x1,y1,x2,y2 (NDC)
    leggl->SetBorderSize(0);
    leggl->SetFillStyle(0);
    leggl->SetTextSize(0.035);
    leggl->AddEntry(h_PSe_data_loose_cuts,          "All data",        "l");
    leggl->AddEntry(h_PSe_data_grinch_loose_cuts,   "GRINCH", "l");
    leggl->AddEntry(h_PSe_data_antigrinch_loose_cuts,"anti-GRINCH",     "l");
    leggl->Draw();

    auto cutLinel = new TLine(0.2, yminl, 0.2, ymaxl);
    cutLinel->SetLineStyle(1); //1-solid 2-dashed
    cutLinel->SetLineWidth(4);
    cutLinel->SetLineColor(kMagenta);
    cutLinel->Draw();
    gPad->RedrawAxis();

    C->Print(Form("images/%s/PionPlots_%s.png",kin_,kin_));
    C1->Print(Form("images/%s/PionPlots_1_%s.png",kin_,kin_));
    C2->Print(Form("images/%s/PionPlots_GRINCH_%s.png",kin_,kin_));
}
