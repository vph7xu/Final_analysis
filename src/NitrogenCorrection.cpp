#include "NitrogenCorrection.h"

#include <TFile.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>

NitrogenCorrection::NitrogenCorrection(const AnalysisCuts& cuts,
		                       const char* kin,
		                       const char* rootFile
				      )	                      
	:c_(cuts),kin_(kin),outfile_(rootFile){}


void NitrogenCorrection::process(TChain& ch_He3_sim, TChain& ch_N2_sim, BranchVarsSim& vHe3, BranchVarsSim& vN2){

	TH1D *h_dx_He3 = new TH1D("h_dx_He3"," delta-x distribution (He3); dx(m)", 100,-4,3);
	TH1D *h_dx_He3_protons = new TH1D("h_dx_He3_protons"," delta-x distribution (He3); dx(m)", 100,-4,3);
	TH1D *h_dx_He3_neutrons = new TH1D("h_dx_He3_neutrons"," delta-x distribution (He3); dx(m)", 100,-4,3);
	TH1D *h_dx_N2  = new TH1D("h_dx_N2","delta-x distribution (N2); dx(m)", 100,-4,3);
	TH1D *h_dx_N2_protons  = new TH1D("h_dx_N2_protons","delta-x distribution (N2); dx(m)", 100,-4,3);
	TH1D *h_dx_N2_neutrons  = new TH1D("h_dx_N2_neutrons","delta-x distribution (N2); dx(m)", 100,-4,3);

	double N_He3_sim = 0.0;
	double N_N2_sim = 0.0;
    // counters for elliptical selection integrals (sum of weights) and
    // their sum-of-weights-squared for uncertainty estimation
    double N_He3_sim_integral_ellipse = 0.0;
    double N_He3_sim_integral_ellipse_var = 0.0;
    double N_N2_sim_integral_ellipse = 0.0;
    double N_N2_sim_integral_ellipse_var = 0.0;

    /////////////////// He3 Sim ///////////////////////////////////////

    Long64_t nentries_He3 = ch_He3_sim.GetEntries();
    const Long64_t step     = 100;
    std::cout << "[NitrogenCorrection] looping over He3 sim " << nentries_He3 << " events\n";

    for (Long64_t i=0;i<nentries_He3;++i){
        ch_He3_sim.GetEntry(i);

        if(/*vHe3.ntrack<1 ||*/ abs(vHe3.vz)>0.27 || vHe3.eHCAL<c_.eHCAL_L || abs((vHe3.ePS+vHe3.eSH)/(vHe3.trP)-1)>0.2 || vHe3.ePS<0.2 ||
            (c_.W2_L>vHe3.W2 || vHe3.W2>c_.W2_H) || (c_.dy_L>vHe3.dy || vHe3.dy>c_.dy_H)) continue;

        h_dx_He3->Fill(vHe3.dx,vHe3.weight);

        if(vHe3.fnucl == 0) {
        	h_dx_He3_neutrons->Fill(vHe3.dx,vHe3.weight);
        	N_He3_sim++;
            // elliptical selection in (dy,dx) plane centered at 0 with radii 0.4
            if ((pow((vHe3.dy-c_.dy_c)/c_.dy_r,2) + pow((vHe3.dx-c_.dx_c)/c_.dx_r,2)) < 1.0) {
                N_He3_sim_integral_ellipse += vHe3.weight;
                N_He3_sim_integral_ellipse_var += vHe3.weight * vHe3.weight;
            }
        }

        if(vHe3.fnucl == 1) h_dx_He3_protons->Fill(vHe3.dx+0.05,vHe3.weight);//+0.05 is for GEN3 its hard coded for now

        // progress bar
        if (i % step == 0 || i == nentries_He3 - 1) {
            double frac = double(i + 1) / nentries_He3;
            int barw = 42, pos = static_cast<int>(barw * frac);
            std::cout << '\r' << '[';
            for (int j = 0; j < barw; ++j)
                std::cout << (j < pos ? '=' : (j == pos ? '>' : ' '));
            std::cout << "] " << static_cast<int>(frac * 100) << " %" << std::flush;
        }

    }

    /////////////////// N2 Sim ///////////////////////////////////////

    Long64_t nentries_N2 = ch_N2_sim.GetEntries();
    //const Long64_t step     = 100;
    std::cout << "[NitrogenCorrection] looping over N2 sim " << nentries_N2 << " events\n";

    for (Long64_t i=0;i<nentries_N2;++i){
        ch_N2_sim.GetEntry(i);

        if(/*vHe3.ntrack<1 ||*/ abs(vN2.vz)>0.27 || vN2.eHCAL<c_.eHCAL_L || abs((vN2.ePS+vN2.eSH)/(vN2.trP)-1)>0.2 || vN2.ePS<0.2 ||
            (c_.W2_L>vN2.W2 || vN2.W2>c_.W2_H) || (c_.dy_L>vN2.dy || vN2.dy>c_.dy_H)) continue;

        h_dx_N2->Fill(vN2.dx,vN2.weight);

        if(vN2.fnucl == 0) {
        	h_dx_N2_neutrons->Fill(vN2.dx,vN2.weight);
        	N_N2_sim++;
            // elliptical selection in (dy,dx) plane centered at 0 with radii 0.4
            if ((pow((vN2.dy-c_.dy_c)/c_.dy_r,2) + pow((vN2.dx-c_.dx_c)/c_.dx_r,2)) < 1.0) {
                N_N2_sim_integral_ellipse += vN2.weight;
                N_N2_sim_integral_ellipse_var += vN2.weight * vN2.weight;
            }
        }

        if(vN2.fnucl == 1) h_dx_N2_protons->Fill(vN2.dx+0.05,vN2.weight);//+0.05 is for GEN3 its hard coded for now

        // progress bar
        if (i % step == 0 || i == nentries_N2 - 1) {
            double frac = double(i + 1) / nentries_N2;
            int barw = 42, pos = static_cast<int>(barw * frac);
            std::cout << '\r' << '[';
            for (int j = 0; j < barw; ++j)
                std::cout << (j < pos ? '=' : (j == pos ? '>' : ' '));
            std::cout << "] " << static_cast<int>(frac * 100) << " %" << std::flush;
        }
    }

    // Use the elliptical-selection integrals calculated above. As a fallback
    // (for comparison) keep the 1D dx-integrals commented out.
    // double N_N2_sim_integral = h_dx_N2_neutrons->Integral(h_dx_N2_neutrons->FindBin(c_.dx_L),h_dx_N2_neutrons->FindBin(c_.dx_H));
    // double N_He3_sim_integral = h_dx_He3_neutrons->Integral(h_dx_He3_neutrons->FindBin(c_.dx_L),h_dx_He3_neutrons->FindBin(c_.dx_H));

    double N_N2_sim_integral = N_N2_sim_integral_ellipse; 
    double N_He3_sim_integral = N_He3_sim_integral_ellipse;

    double S = (N_N2_sim_integral * nentries_He3 / nentries_N2) / N_He3_sim_integral / 14;

    double N_N2_fill = 0.015;//fill values defined at the fill
    double N_He3_fill = 0.985;//fill values defined at the fill ask Hunter for detailed values

    double frac_N2 = 14*N_N2_fill/(14*N_N2_fill + N_He3_fill)*S;

    ////////////////////////////////////error with the fraction //////////////////////////////////////////
        // --- inputs (already computed) ---------------------------------
    double A = N_N2_sim_integral;
    double B = N_He3_sim_integral;
    double C = nentries_He3;
    double D = nentries_N2;
    double E = N_N2_fill;
    double F = N_He3_fill;
    const double k = 14.0;

    // --- statistical uncertainties (examples) ----------------------
    // For the elliptical selection we don't have a single ROOT bin; use the
    // sum-of-weights-squared to estimate the variance (sigma^2 = sum w_i^2).
    double sA = (N_N2_sim_integral_ellipse_var > 0.0) ? std::sqrt(N_N2_sim_integral_ellipse_var) : 0.0;
    double sB = (N_He3_sim_integral_ellipse_var > 0.0) ? std::sqrt(N_He3_sim_integral_ellipse_var) : 0.0;

    // Poisson for raw entry counts
    double sC = sqrt(C);
    double sD = sqrt(D);

    // fill-fraction uncertainties (replace with survey errors if known)
    double sE = N_N2_fill*0.02;  //check these numbers // e.g. 0.015 * 0.02  (2 % relative)  
    double sF = N_He3_fill*0.02;  //check these numbers // e.g. 0.985 * 0.02

    // --- derived quantities ---------------------------------------
    double S_1 = (A*C)/(D*B*k);
    double g = (k*E)/(k*E + F);
    double frac_N2_1 = g * S_1;

    // partials
    double dSdA =  C / (D*B*k);
    double dSdB = -A*C / (D*B*B*k);
    double dSdC =  A / (D*B*k);
    double dSdD = -A*C / (D*D*B*k);

    double dgdE =  k*F / pow(k*E + F, 2);
    double dgdF = -k*E / pow(k*E + F, 2);

    // --- error propagation ----------------------------------------
    double var_S   = pow(dSdA*sA,2) + pow(dSdB*sB,2)
                   + pow(dSdC*sC,2) + pow(dSdD*sD,2);
    double var_g   = pow(dgdE*sE,2) + pow(dgdF*sF,2);

    double sigma_frac = sqrt( pow(g,2)*var_S + pow(S,2)*var_g );


    std::ofstream txt(Form("corrections/%s/NitrogenCorrection_%s.txt",kin_,kin_));
    txt<<"N_N2_sim = "<<N_N2_sim<<"\n";
    txt<<"N_He3_sim = "<<N_He3_sim<<"\n";
    txt<<"N_N2_sim_integral = "<<N_N2_sim_integral<<"\n";
    txt<<"N_He3_sim_integral = "<<N_He3_sim_integral<<"\n";
    txt<<"nentries_N2 = "<<nentries_N2<<"\n";
    txt<<"nentries_He3 = "<<nentries_He3<<"\n";
    txt<<"S = "<<S<<"\n";
    txt<<"Nitrogen fraction = "<<frac_N2<<"\n";
    txt<<"f_N2 = "<<frac_N2_1<<"\n";
    txt<<"err_f_N2 = "<<sigma_frac <<"\n";

    txt.close();

    ///////////////////plotting and printing////////////////////////////

    TCanvas *C1 = new TCanvas("c","c",2400,1500);
    h_dx_He3->SetLineWidth(4);
    h_dx_N2->SetLineWidth(4);

    C1->Divide(2,2);
    C1->cd(1);
    h_dx_He3->Draw();
    C1->cd(2);
    h_dx_N2->Draw();

    C1->Print(Form("images/%s/NitrogenCorrection_%s.png",kin_,kin_));

}