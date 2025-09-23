#include "InelasticCorrection.h"
#include "Utility.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TNamed.h>
#include <TLegend.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TLine.h>

using std::string; using std::vector; using std::unordered_map;

InelasticCorrection::InelasticCorrection(const AnalysisCuts& cuts,
                                         const RunQuality* rq,
                                         const char* kin,
                                         const char* rootFile)
    : c_(cuts), rq_(rq), kin_(kin), outFile_(rootFile) {}

// --- template shape fit -------------------------------------------------------
/*TH1D* InelasticCorrection::performFit(TH1D* hD, TH1D* hInel, TH1D* hQE_proton, TH1D* hQE_neutron,
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
*/

TH1D* InelasticCorrection::performFit(TH1D* hD, TH1D* hInel, TH1D* hQE_proton, TH1D* hQE_neutron,
                                      double& par0, double& par1, double& par2
                                      /* optional: add outputs for shifts if you like
                                       * , double& dx_p_out, double& dx_n_out, double& dx_inel_out */
                                      )
{
    // --- Clone & safe-normalize (avoid div by 0) ---
    auto clone_norm = [](TH1D* src, const char* newname){
        TH1D* h = (TH1D*)src->Clone(newname);
        double I = h->Integral();
        if (I > 0) h->Scale(1.0/I);
        return h;
    };

    TH1D *d    = clone_norm(hD,         Form("hD_%p",         hD));
    TH1D *inel = clone_norm(hInel,      Form("hInel_%p",      hInel));
    TH1D *qe_p = clone_norm(hQE_proton, Form("hQE_proton_%p", hQE_proton));
    TH1D *qe_n = clone_norm(hQE_neutron,Form("hQE_neutron_%p",hQE_neutron));

    const double xmin = d->GetXaxis()->GetXmin();
    const double xmax = d->GetXaxis()->GetXmax();

    // helper: robust interpolate with edge clamping
    auto interp = [](TH1D* h, double xx){
        const double hxmin = h->GetXaxis()->GetXmin();
        const double hxmax = h->GetXaxis()->GetXmax();
        if (xx <= hxmin) return h->GetBinContent(1);
        if (xx >= hxmax) return h->GetBinContent(h->GetNbinsX());
        return h->Interpolate(xx);
    };

    // --- Fit model with 6 params: [0]=A, [1]=rN/P, [2]=rInel/P, [3]=dx_p, [4]=dx_n, [5]=dx_inel
    TF1* f = new TF1("fdx", [&](double* x, double* p){
        const double xx = x[0];
        const double P  = interp(qe_p, xx - p[3]);   // shift QE-p by dx_p
        const double N  = interp(qe_n, xx - p[4]);   // shift QE-n by dx_n
        const double I  = interp(inel, xx - p[5]);   // shift Inelastic by dx_inel
        // same weight structure as your original: overall A times (P + rN*N + rI*I)
        return p[0] * ( P + p[1]*N + p[2]*I );
    }, xmin, xmax, 6);

    f->SetParNames("A", "rN_over_P", "rInel_over_P", "dx_p", "dx_n", "dx_inel");

    // sensible initials
    f->SetParameters(1.0, 0.3, 0.1, 0.0, 0.0, 0.0);

    // limits: tune as you like; keep shifts small to avoid edge-area loss
    f->SetParLimits(0, 0.0, 10.0); // amplitude
    f->SetParLimits(1, 0.0, 5.0);  // N/P ratio
    f->SetParLimits(2, 0.0, 5.0);  // Inel/P ratio
    const double maxShift = 0.02;  // units = your dx units; e.g., 0.02 m or 2 cm (adjust!)
    f->SetParLimits(3, -maxShift, +maxShift);
    f->SetParLimits(4, -maxShift, +maxShift);
    f->SetParLimits(5, -maxShift, +maxShift);

    // Fit quietly over the histogram range
    d->Fit(f, "RQ0"); // R: use range, Q: quiet, 0: don't draw

    const int kDxP = 3, kDxN = 4, kDxI = 5;

    double dx_p    = f->GetParameter(kDxP);
    double dx_n    = f->GetParameter(kDxN);
    double dx_inel = f->GetParameter(kDxI);

    double edx_p    = f->GetParError(kDxP);
    double edx_n    = f->GetParError(kDxN);
    double edx_inel = f->GetParError(kDxI);

    std::cout << std::fixed << std::setprecision(4)
              << "[InelasticCorrection] Shifts (dx, axis units):  "
              << "p = "     << dx_p    << " ± " << edx_p
              << ", n = "   << dx_n    << " ± " << edx_n
              << ", inel = "<< dx_inel << " ± " << edx_inel
              << std::endl;

    // output the three "physics" parameters as before
    par0 = f->GetParameter(0);           // overall A
    par1 = f->GetParameter(1);           // rN_over_P
    par2 = f->GetParameter(2);           // rInel_over_P
    // if you added refs for shifts, set them here:
    // dx_p_out    = f->GetParameter(3);
    // dx_n_out    = f->GetParameter(4);
    // dx_inel_out = f->GetParameter(5);

    // --- Build combined best-fit histogram (with shifts) ---
    TH1D* comb = (TH1D*)inel->Clone(Form("hComb_%p", inel));
    comb->Reset();

    for (int i = 1; i <= comb->GetNbinsX(); ++i){
        const double x = comb->GetXaxis()->GetBinCenter(i);
        const double P = interp(qe_p, x - f->GetParameter(3));
        const double N = interp(qe_n, x - f->GetParameter(4));
        const double I = interp(inel, x - f->GetParameter(5));
        comb->SetBinContent(i, par0 * ( P + par1*N + par2*I ));
    }

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
        return par[0]*(qe_p->GetBinContent(bin1) + Rnop*qe_n->GetBinContent(bin2)) + par[1]*inel->GetBinContent(bin3);
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
        comb->SetBinContent(i, par0*(qe_p->GetBinContent(i) + Rnop*qe_n->GetBinContent(i)) + par1*inel->GetBinContent(i));
    return comb;

}

// α = N_inel / (N_QE_p + N_QE_n), Δ = W2 shift of inelastic (GeV^2)
TH1D* InelasticCorrection::performFitW2_1(
    TH1D* hD, TH1D* hInel, TH1D* hQEp, TH1D* hQEn,
    double& alpha, double& delta, double Rnop)
{
    auto unit = [&](TH1D* h, const char* tag)->bool{
        double s=h->Integral(); if(s<=0){ std::cerr<<"[W2_1] Empty "<<tag<<"\n"; return false; }
        h->Scale(1.0/s); return true;
    };

    TH1D *d  =(TH1D*)hD->Clone("hDW2");        if(!unit(d,"data"))  return (TH1D*)hInel->Clone("hComb_empty");
    TH1D *pi =(TH1D*)hQEp->Clone("qepW2");     if(!unit(pi,"QEp"))  return (TH1D*)hInel->Clone("hComb_empty");
    TH1D *ni =(TH1D*)hQEn->Clone("qenW2");     if(!unit(ni,"QEn"))  return (TH1D*)hInel->Clone("hComb_empty");
    TH1D *ii =(TH1D*)hInel->Clone("inelW2");   if(!unit(ii,"Inel")) return (TH1D*)hInel->Clone("hComb_empty");

    const double xmin = d->GetXaxis()->GetXmin();
    const double xmax = d->GetXaxis()->GetXmax();

    auto inel_shift = [&](double x, double del)->double {
        const double xq = x; //- del;
        //const double xq = x - del;
        if (xq<=xmin || xq>=xmax) return 0.0;
        const int bi = ii->FindBin(xq);
        return ii->GetBinContent(bi);               // smooth shift (better than FindBin)
    };

    TF1* f = new TF1(Form("fW2mix_%p", d),
        [&](double* x, double* par){
            const double a  = par[0];             // alpha
            //const double dL = par[1];             // delta (GeV^2)
            
            //const double qep = pi->Interpolate(x[0]);      // was FindBin/GetBinContent
            //const double qen = ni->Interpolate(x[0]);      // was FindBin/GetBinContent
            //const double ine = ii->Interpolate(x[0] - dL); // your inel_shift does this too

            const int bp = pi->FindBin(x[0]);
            const int bn = ni->FindBin(x[0]);
            const double qep = pi->GetBinContent(bp);
            const double qen = ni->GetBinContent(bn);
            const double ine = ii->GetBinContent(ii->FindBin(x[0]));//inel_shift(x[0], dL);
            //const double ine = inel_shift(x[0], dL);
            // simple normalized mixture (good if Δ is small or range wide)
            

            return (qep + Rnop*qen + a*ine)/(1.0 + Rnop + a);
        }, xmin, xmax, 1);

    f->SetParameter(0, 0.5);     // alpha start
    f->SetParLimits(0, 0.0, 1000.0);
    //f->SetParameter(1, 0.0);     // delta start (GeV^2)
    //f->SetParLimits(1,-0.4, 0.4); // tune to your expected shift range

    d->Fit(f, "RQ");

    alpha = f->GetParameter(0);
    //delta = f->GetParameter(1);

    // Build combined unit-area PDF using fitted Δ
    TH1D* comb=(TH1D*)ii->Clone("hCombW2"); comb->Reset();
    for (int i=1;i<=comb->GetNbinsX();++i){
        const double x = comb->GetXaxis()->GetBinCenter(i);
        const double qep = pi->GetBinContent(pi->FindBin(x));
        const double qen = ni->GetBinContent(ni->FindBin(x));
        const double ine = ii->GetBinContent(ii->FindBin(x));//inel_shift(x, delta);
        //const double ine = inel_shift(x, delta);
        // double qep = pi->Interpolate(x);           // was GetBinContent(FindBin(x))
        // double qen = ni->Interpolate(x);
        // double ine = ii->Interpolate(x - delta);   // or inel_shift(x, delta)
        comb->SetBinContent(i, (qep + Rnop*qen + alpha*ine)/(1.0 + Rnop + alpha));
    }
    return comb;
}



// par0 = QE scale (fraction or yield proxy), par1 = Inel scale, par2 = Δ shift (GeV^2)
TH1D* InelasticCorrection::performFitW2_2(
    TH1D* hD, TH1D* hInel, TH1D* hQE, double& par0, double& par1, double& delta)
{
    auto unit = [&](TH1D* h, const char* tag)->bool{
        double s=h->Integral(); if(s<=0){ std::cerr<<"[W2_2] Empty "<<tag<<"\n"; return false; }
        h->Scale(1.0/s); return true;
    };
    TH1D *d  =(TH1D*)hD->Clone("hDW2_2");       if(!unit(d,"data")) return (TH1D*)hInel->Clone("hComb_empty");
    TH1D *qe =(TH1D*)hQE->Clone("qeW2");        if(!unit(qe,"QE"))  return (TH1D*)hInel->Clone("hComb_empty");
    TH1D *ii =(TH1D*)hInel->Clone("inelW2_2");  if(!unit(ii,"Inel"))return (TH1D*)hInel->Clone("hComb_empty");

    const double xmin = d->GetXaxis()->GetXmin();
    const double xmax = d->GetXaxis()->GetXmax();

    auto inel_shift = [&](double x, double del)->double {
        const double xq = x; //- del;
        //const double xq = x- del;
        if (xq<=xmin || xq>=xmax) return 0.0;
        const int bi = ii->FindBin(xq);
        return ii->GetBinContent(bi); 
    };

    TF1* f = new TF1(Form("fW2qeIn_%p", d),
        [&](double* x, double* par){
            const double aQ = par[0];  // QE weight
            const double aI = par[1];  // Inel weight
            //const double dL = par[2];  // Δ

            // const double q  = qe->Interpolate(x[0]);
            //const double ine= ii->Interpolate(x[0] - dL);
            const double q  = qe->GetBinContent(qe->FindBin(x[0]));
            const double ine= ii->GetBinContent(qe->FindBin(x[0]));//inel_shift(x[0], dL);
            //const double ine= inel_shift(x[0], dL);
            return aQ*q + aI*ine;      // not forced to unit area; fit handles norm
        }, xmin, xmax, 2);

    f->SetParameters(0.5, 0.5, 0.0);
    f->SetParLimits(0, 0.0, 100.0);
    f->SetParLimits(1, 0.0, 100.0);
    //f->SetParLimits(2, -0.4, 0.4);

    d->Fit(f,"RQ");

    par0  = f->GetParameter(0);
    par1  = f->GetParameter(1);
    //delta = f->GetParameter(2);

    // Combined model (same form as fit)
    TH1D* comb=(TH1D*)ii->Clone("hCombW2_2"); comb->Reset();
    for(int i=1;i<=comb->GetNbinsX();++i){
        const double x = comb->GetXaxis()->GetBinCenter(i);
        // const double q  = qe->Interpolate(x);
        // const double ine = ii->Interpolate(x - delta);
        const double q  = qe->GetBinContent(qe->FindBin(x));
        const double ine= ii->GetBinContent(ii->FindBin(x));//inel_shift(x, delta);
        //const double ine= inel_shift(x, delta);
        comb->SetBinContent(i, par0*q + par1*ine);
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

        // Build a unit-area shifted copy of a 1D template
    auto make_shifted_unit = [&](TH1D* src, const char* name, double delta)->TH1D* {
        TH1D* h = (TH1D*)src->Clone(name);
        h->Reset();
        TAxis* ax = h->GetXaxis();
        const double xmin = ax->GetXmin(), xmax = ax->GetXmax();

        for (int i=1; i<=h->GetNbinsX(); ++i) {
            double x    = ax->GetBinCenter(i);
            double xsrc = x - delta;
            double y    = (xsrc<=xmin || xsrc>=xmax) ? 0.0 : src->GetBinContent(src->FindBin(xsrc)); //src->Interpolate(xsrc);
            h->SetBinContent(i, y);
        }
        double s = h->Integral();
        if (s>0) h->Scale(1.0/s); // unit area
        return h;
    };


    auto dx_shifted_QE = [&](const BranchVarsSim& s){
        double d = s.dx;
        if (std::strcmp(kin_, "GEN3_He3") == 0 && s.fnucl==0) d -= 0.0;//0.1;      // n
        else if (std::strcmp(kin_, "GEN4_He3") == 0 && s.fnucl==1) d -= 0.0;//0.02;     // p
        else if (std::strcmp(kin_, "GEN4_He3")  == 0 && s.fnucl==0) d -= 0.0;//0.07;
        else if (std::strcmp(kin_, "GEN4b_He3")== 0 && s.fnucl==0) d -= 0.0;//0.025;    // n
        return d;
    };

    auto dx_shifted_Inel = [&](const BranchVarsSim& s){
        double d = s.dx;
        if (std::strcmp(kin_, "GEN3_He3") == 0)  d += 0.0;//0.4;
        else if (std::strcmp(kin_, "GEN4_He3")==0)  d += 0.0;//0.6;
        else if (std::strcmp(kin_, "GEN4b_He3")==0) d += 0.0;//0.25;
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


    auto makeAsymGraph = [&](const TH1D& hpos, const TH1D& hneg,
                             const char* name, int mcolor, int mstyle)->TGraphErrors* {
      auto* g = new TGraphErrors(); g->SetName(name);
      g->SetTitle("Asymmetry vs W^{2};W^{2} (GeV^{2});Asymmetry (%)");
      int ip=0;
      for (int i=1;i<=hpos.GetNbinsX();++i){
        const double Np = hpos.GetBinContent(i);
        const double Nn = hneg.GetBinContent(i);
        const double N  = Np + Nn;
        if (N<=0) continue;
        const double A   = (Np - Nn)/N;
        const double sA  = std::sqrt(std::max(0.0, (1.0 - A*A)/N)); // σ_A = sqrt((1-A^2)/N)
        const double x   = hpos.GetXaxis()->GetBinCenter(i);
        g->SetPoint(ip, x, 100.0*A);          // percent
        g->SetPointError(ip, 0.0, 100.0*sA);  // only vertical error
        ++ip;
      }
      g->SetMarkerStyle(mstyle);
      g->SetMarkerSize(1.2);
      g->SetMarkerColor(mcolor);
      g->SetLineColor(mcolor);
      return g;
    };



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

    double W2_hist_upper_limit = 0.0;

    if(std::strcmp(kin_, "GEN3_He3") == 0){
        W2_hist_upper_limit = 2.0;
    }
    else if(std::strcmp(kin_, "GEN4b_He3") == 0 || std::strcmp(kin_, "GEN4_He3") == 0){
        W2_hist_upper_limit = 2.0;
    }
    else{
        W2_hist_upper_limit = 2.0;
    }

    const int NBW2 = 100; // match your other W² binning

    TH1D hData_W2("hData_W2",  "W2 data; W^{2} (GeV^{2})", NBW2, -1, W2_hist_upper_limit);
    TH1D hData_W2_Neutrons("hData_W2_Neutrons",  "W2 data; W^{2} (GeV^{2})", NBW2, -1, W2_hist_upper_limit);
    TH1D hQE_W2("hQE_W2", "W^{2} QE sim",          NBW2, -1, W2_hist_upper_limit);
    TH1D hQE_W2_Neutrons("hQE_W2_Neutrons", "W^{2} QE sim",          NBW2, -1, W2_hist_upper_limit);
    TH1D hQE_proton_W2("hQE_proton_W2", "W^{2} QE sim protons",          NBW2, -1, W2_hist_upper_limit);
    TH1D hQE_proton_W2_Neutrons("hQE_proton_W2_Neutrons", "W^{2} QE sim protons",          NBW2, -1, W2_hist_upper_limit);
    TH1D hQE_neutron_W2("hQE_neutron_W2", "W^{2} QE sim neutrons",          NBW2, -1, W2_hist_upper_limit);
    TH1D hQE_neutron_W2_Neutrons("hQE_neutron_W2_Neutrons", "W^{2} QE sim neutrons",          NBW2, -1, W2_hist_upper_limit);
    TH1D hInelastic_W2("hInelastic_W2",  "W^{2} inelastic sim",        NBW2, -1, W2_hist_upper_limit);
    TH1D hInelastic_W2_Neutrons("hInelastic_W2_Neutrons",  "W^{2} inelastic sim",        NBW2, -1, W2_hist_upper_limit);
    TH1D hInelastic_W2_2("hInelastic_W2_2",  "W^{2} inelastic sim",        NBW2, -1, W2_hist_upper_limit);
    TH1D hInelastic_W2_2_Neutrons("hInelastic_W2_2_Neutrons",  "W^{2} inelastic sim",        NBW2, -1, W2_hist_upper_limit);

    TH2D hDxdy("hDxdy", "dxdy distribution ; dx(m)",100,-4,3,100,-4,3);
    TH2D hDxdy_cut("hDxdy_cut", "dxdy distribution ; dx(m)",100,-4,3,100,-4,3);

    TH2D hDxdy_inelastic("hDxdy_inelastic", "dxdy distribution ; dy(m); dx(m)",100,-4,3,100,-4,3);
    TH1D hDx_inelastic("hDx_inelastic", "dx distribution ; dx(m)",100,-4,3);
    TH1D hDy_inelastic("hDy_inelastic", "dy distribution ; dy(m)",100,-4,3);
    TH1D hDx_elastic("hDx_elastic", "dx distribution ; dx(m)",100,-4,3);
    TH1D hDy_elastic("hDy_elastic", "dy distribution ; dy(m)",100,-4,3);
    TH1D hDx_both("hDx_both", "dx distribution ; dx(m)",100,-4,3);
    TH1D hDy_both("hDy_both", "dy distribution ; dy(m)",100,-4,3);

    const int NBW2_low = 100;

    double W2hist_low = -2;
    double W2hist_high = 8;

    double dxhist_low = -4;
    double dxhist_high = 3;

    double dyhist_low = -2;
    double dyhist_high = 2;

    if (std::strcmp(kin_, "GEN2_He3") == 0){W2hist_low = -0.5; W2hist_high = 2;}
    else if (std::strcmp(kin_, "GEN3_He3") == 0){W2hist_low = -1.5; W2hist_high = 4;} 
    else if (std::strcmp(kin_, "GEN4_He3") == 0){W2hist_low = -2; W2hist_high = 6;}
    else if (std::strcmp(kin_, "GEN4b_He3") == 0){W2hist_low = -2; W2hist_high = 6;}

    if (std::strcmp(kin_, "GEN2_He3") == 0){dxhist_low = -4; dxhist_high = 3;}
    else if (std::strcmp(kin_, "GEN3_He3") == 0){dxhist_low = -4; dxhist_high = 3;} 
    else if (std::strcmp(kin_, "GEN4_He3") == 0){dxhist_low = -3; dxhist_high = 2;}
    else if (std::strcmp(kin_, "GEN4b_He3") == 0){dxhist_low = -3; dxhist_high = 2;}

    if (std::strcmp(kin_, "GEN2_He3") == 0){dyhist_low = -2; dyhist_high = 2;}
    else if (std::strcmp(kin_, "GEN3_He3") == 0){dyhist_low = -2; dyhist_high = 2;} 
    else if (std::strcmp(kin_, "GEN4_He3") == 0){dyhist_low = -1.5; dyhist_high = 1.5;}
    else if (std::strcmp(kin_, "GEN4b_He3") == 0){dyhist_low = -1.5; dyhist_high = 1.5;}

    // W² distributions per cut (no helicity split) for the right panel
    TH1D hW2_dxonly     ("hW2_dxonly",     "W^{2};W^{2} (GeV^{2})", NBW2, W2hist_low, W2hist_high);
    TH1D hW2_dxdy       ("hW2_dxdy",       "W^{2};W^{2} (GeV^{2})", NBW2, W2hist_low, W2hist_high);
    TH1D hW2_dxAntiDy   ("hW2_dxAntiDy",   "W^{2};W^{2} (GeV^{2})", NBW2, W2hist_low, W2hist_high);
    TH1D hW2_dyAntiDx   ("hW2_dyAntiDx",   "W^{2};W^{2} (GeV^{2})", NBW2, W2hist_low, W2hist_high);

    // W² split by helicity for asymmetry vs W² (left panel)
    TH1D hW2_dxonly_pos   ("hW2_dxonly_pos",   "W^{2};W^{2} (GeV^{2})", NBW2_low, W2hist_low, W2hist_high);
    TH1D hW2_dxonly_neg   ("hW2_dxonly_neg",   "W^{2};W^{2} (GeV^{2})", NBW2_low, W2hist_low, W2hist_high);

    TH1D hW2_dxdy_pos     ("hW2_dxdy_pos",     "W^{2};W^{2} (GeV^{2})", NBW2_low, W2hist_low, W2hist_high);
    TH1D hW2_dxdy_neg     ("hW2_dxdy_neg",     "W^{2};W^{2} (GeV^{2})", NBW2_low, W2hist_low, W2hist_high);

    TH1D hW2_dxAntiDy_pos ("hW2_dxAntiDy_pos", "W^{2};W^{2} (GeV^{2})", NBW2_low, W2hist_low, W2hist_high);
    TH1D hW2_dxAntiDy_neg ("hW2_dxAntiDy_neg", "W^{2};W^{2} (GeV^{2})", NBW2_low, W2hist_low, W2hist_high);

    TH1D hW2_dyAntiDx_pos ("hW2_dyAntiDx_pos", "W^{2};W^{2} (GeV^{2})", NBW2_low, W2hist_low, W2hist_high);
    TH1D hW2_dyAntiDx_neg ("hW2_dyAntiDx_neg", "W^{2};W^{2} (GeV^{2})", NBW2_low, W2hist_low, W2hist_high);

    // (optional) get ROOT to store per-bin variances
    hW2_dxonly.Sumw2(); hW2_dxdy.Sumw2(); hW2_dxAntiDy.Sumw2(); hW2_dyAntiDx.Sumw2();
    hW2_dxonly_pos.Sumw2(); hW2_dxonly_neg.Sumw2();
    hW2_dxdy_pos.Sumw2();   hW2_dxdy_neg.Sumw2();
    hW2_dxAntiDy_pos.Sumw2(); hW2_dxAntiDy_neg.Sumw2();
    hW2_dyAntiDx_pos.Sumw2(); hW2_dyAntiDx_neg.Sumw2();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    // dx distributions per cut (no helicity split) for the right panel
    TH1D hDx_W2only     ("hDx_W2only",     "dx;W^{2} (m)", NBW2, dxhist_low, dxhist_high);
    TH1D hDx_W2dy       ("hDx_W2dy",       "dx;W^{2} (m)", NBW2, dxhist_low, dxhist_high);
    TH1D hDx_W2AntiDy   ("hDx_W2AntiDy",   "dx;W^{2} (m)", NBW2, dxhist_low, dxhist_high);
    TH1D hDx_dyAntiW2   ("hDx_dyAntiW2",   "dx;W^{2} (m)", NBW2, dxhist_low, dxhist_high);

    // dx split by helicity for asymmetry vs dx (left panel)
    TH1D hDx_W2only_pos   ("hDx_W2only_pos",   "dx;W^{2} (m)", NBW2_low, dxhist_low, dxhist_high);
    TH1D hDx_W2only_neg   ("hDx_W2only_neg",   "dx;W^{2} (m)", NBW2_low, dxhist_low, dxhist_high);

    TH1D hDx_W2dy_pos     ("hDx_W2dy_pos",     "dx;W^{2} (m)", NBW2_low, dxhist_low, dxhist_high);
    TH1D hDx_W2dy_neg     ("hDx_W2dy_neg",     "dx;W^{2} (m)", NBW2_low, dxhist_low, dxhist_high);

    TH1D hDx_W2AntiDy_pos ("hDx_W2AntiDy_pos", "dx;W^{2} (m)", NBW2_low, dxhist_low, dxhist_high);
    TH1D hDx_W2AntiDy_neg ("hDx_W2AntiDy_neg", "dx;W^{2} (m)", NBW2_low, dxhist_low, dxhist_high);

    TH1D hDx_dyAntiW2_pos ("hDx_dyAntiW2_pos", "dx;W^{2} (m)", NBW2_low, dxhist_low, dxhist_high);
    TH1D hDx_dyAntiW2_neg ("hDx_dyAntiW2_neg", "dx;W^{2} (m)", NBW2_low, dxhist_low, dxhist_high);

    // (optional) get ROOT to store per-bin variances
    hDx_W2only.Sumw2(); hDx_W2dy.Sumw2(); hDx_W2AntiDy.Sumw2(); hDx_dyAntiW2.Sumw2();
    hDx_W2only_pos.Sumw2(); hDx_W2only_neg.Sumw2();
    hDx_W2dy_pos.Sumw2();   hDx_W2dy_neg.Sumw2();
    hDx_W2AntiDy_pos.Sumw2(); hDx_W2AntiDy_neg.Sumw2();
    hDx_dyAntiW2_pos.Sumw2(); hDx_dyAntiW2_neg.Sumw2();

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    TH1D hDy_W2dx     ("hDy_W2dx",     "dy;W^{2} (m)", NBW2, dyhist_low, dyhist_high);


    TH1D hDy_W2dx_pos     ("hDy_W2dx_pos",     "dy;W^{2} (m)", NBW2, dyhist_low, dyhist_high);
    TH1D hDy_W2dx_neg     ("hDy_W2dx_neg",     "dy;W^{2} (m)", NBW2, dyhist_low, dyhist_high);

    hDy_W2dx.Sumw2();
    hDy_W2dx_pos.Sumw2();   hDy_W2dx_neg.Sumw2();

    ////////////////////////// --- loop data   ///////////////////////////////////////////////////////
    //                                                                                              //
    //                                                                                              //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    Long64_t n=ch.GetEntries();
    std::cout << "\n"<<"[InelasticCorrection] looping over data " << n << " events\n";
    const Long64_t step     = 100;
    for(Long64_t i=0;i<n;++i){ 
        
        int helCorr = -1*v.helicity*v.IHWP*c_.Pkin_L;

        ch.GetEntry(i);
        
        if(rq_ && (!rq_->helicityOK(v.runnum)||!rq_->mollerOK(v.runnum))) continue;
   
        if( v.runnum<c_.runnum_L || v.runnum>c_.runnum_H ||
            v.ntrack<1 || abs(v.vz)>0.27 || v.eHCAL<c_.eHCAL_L || abs((v.ePS+v.eSH)/(v.trP)-1)>0.2 || v.ePS<0.2 ||
            (c_.coin_L>v.coin_time || v.coin_time>c_.coin_H) || 
            abs(v.helicity)!=1) continue;

        
        if(c_.W2_L<v.W2 && v.W2<c_.W2_H){
            hDxdy.Fill(v.dy,v.dx);
        }

        if(2.5<v.W2 && v.W2<5){
            hDxdy_inelastic.Fill(v.dy,v.dx);
            hDx_inelastic.Fill(v.dx);
            hDy_inelastic.Fill(v.dy);
        }

        if(-1<v.W2 && v.W2<2.5){
            hDx_elastic.Fill(v.dx);
            hDy_elastic.Fill(v.dy);
        }

        if(-1<v.W2 && v.W2<5){
            hDx_both.Fill(v.dx);
            hDy_both.Fill(v.dy);
        }

        if((c_.dy_L<v.dy && v.dy<c_.dy_H) && (c_.dx_L<v.dx && v.dx<c_.dx_H) /*((pow((v.dy-0.0)/0.4,2)+pow((v.dx-0.0)/0.3,2))<=1)*/) {
            hData_W2_Neutrons.Fill(v.W2);
        }

        if(((pow((v.dy-0.0)/0.4,2)+pow((v.dx-0.0)/0.3,2))<=1) || ((pow((v.dy-0.0)/0.4,2)+pow((v.dx-(c_.dx_P_L+c_.dx_P_H)/2)/0.3,2))<=1)){
            hData_W2.Fill(v.W2);
            if(c_.W2_L<v.W2 && v.W2<c_.W2_H){
                hDxdy_cut.Fill(v.dy,v.dx);
            }
        }

        // define the 3 cut flavors you want to compare
        const bool pass_dx      = (c_.dx_L < v.dx && v.dx < c_.dx_H);
        const bool pass_dy      = (c_.dy_L < v.dy && v.dy < c_.dy_H);
        const bool pass_anti_dy = (-1.5>v.dy || 1.5<v.dy); // outside the dy window
        const bool pass_anti_dx = (c_.dx_P_L-0.3>v.dx || c_.dx_H+0.3<v.dx); // outside the dy window

        // --- W² shapes (right panel) ---
        if (pass_dx)                 hW2_dxonly.Fill(v.W2);
        if (pass_dx && pass_dy)      hW2_dxdy.Fill(v.W2);
        if (pass_dx && pass_anti_dy) hW2_dxAntiDy.Fill(v.W2);
        if (pass_dy && pass_anti_dx) hW2_dyAntiDx.Fill(v.W2);

        // --- helicity-split W² (for A(W²), left panel) ---
        if (pass_dx){
          if (helCorr==+1) hW2_dxonly_pos.Fill(v.W2);
          else if (helCorr==-1) hW2_dxonly_neg.Fill(v.W2);
        }
        if (pass_dx && pass_dy){
          if (helCorr==+1) hW2_dxdy_pos.Fill(v.W2);
          else if (helCorr==-1) hW2_dxdy_neg.Fill(v.W2);
        }
        if (pass_dx && pass_anti_dy){
          if (helCorr==+1) hW2_dxAntiDy_pos.Fill(v.W2);
          else if (helCorr==-1) hW2_dxAntiDy_neg.Fill(v.W2);
        }
        if (pass_dy && pass_anti_dx){
          if (helCorr==+1) hW2_dyAntiDx_pos.Fill(v.W2);
          else if (helCorr==-1) hW2_dyAntiDx_neg.Fill(v.W2);
        }

        // define the 3 cut flavors you want to compare
        const bool pass_W2      = (c_.W2_L < v.W2 && v.W2 < c_.W2_H);
        //const bool pass_dy      = (c_.dy_L < v.dy && v.dy < c_.dy_H);
        //const bool pass_anti_dy = (-1.5>v.dy || 1.5<v.dy); // outside the dy window
        const bool pass_anti_W2 = (2.5<v.W2); // outside the dy window

        // --- dx shapes (right panel) ---
        if (pass_W2)                 hDx_W2only.Fill(v.dx);
        if (pass_W2 && pass_dy)      hDx_W2dy.Fill(v.dx);
        if (pass_W2 && pass_anti_dy) hDx_W2AntiDy.Fill(v.dx);
        if (pass_dy && pass_anti_W2) hDx_dyAntiW2.Fill(v.dx);

        // --- helicity-split dx (for A(dx), left panel) ---
        if (pass_W2){
          if (helCorr==+1) hDx_W2only_pos.Fill(v.dx);
          else if (helCorr==-1) hDx_W2only_neg.Fill(v.dx);
        }
        if (pass_W2 && pass_dy){
          if (helCorr==+1) hDx_W2dy_pos.Fill(v.dx);
          else if (helCorr==-1) hDx_W2dy_neg.Fill(v.dx);
        }
        if (pass_W2 && pass_anti_dy){
          if (helCorr==+1) hDx_W2AntiDy_pos.Fill(v.dx);
          else if (helCorr==-1) hDx_W2AntiDy_neg.Fill(v.dx);
        }
        if (pass_dy && pass_anti_W2){
          if (helCorr==+1) hDx_dyAntiW2_pos.Fill(v.dx);
          else if (helCorr==-1) hDx_dyAntiW2_neg.Fill(v.dx);
        }


        //const bool pass_dx      = (c_.dx_L < v.dx && v.dx < c_.dx_H);

        if (pass_dx && pass_W2) hDy_W2dx.Fill(v.dy);

        if (pass_dx && pass_W2) {
            if (helCorr==+1) hDy_W2dx_pos.Fill(v.dy);
            if (helCorr==-1) hDy_W2dx_neg.Fill(v.dy);
        }  


        if( v.runnum<c_.runnum_L || v.runnum>c_.runnum_H ||
            v.ntrack<1 || abs(v.vz)>0.27 || v.eHCAL<c_.eHCAL_L || abs((v.ePS+v.eSH)/(v.trP)-1)>0.2 || v.ePS<0.2 ||
            (c_.coin_L>v.coin_time || v.coin_time>c_.coin_H) || (c_.W2_L>v.W2 || v.W2>c_.W2_H) || (c_.dy_L>v.dy || v.dy>c_.dy_H) || 
            abs(v.helicity)!=1) continue; // no dx cut since we are fitting it
        
        hData.Fill(v.dx);
        
        if (helCorr == 1){
            hData_pos.Fill(v.dx);
        }
        else if(helCorr == -1){
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
    
    ////////////////////////// --- loop QE sim   /////////////////////////////////////////////////////
    //                                                                                              //
    //                                                                                              //
    //////////////////////////////////////////////////////////////////////////////////////////////////


    Long64_t nentries_QE=ch_QE.GetEntries();
    std::cout << "\n"<< "[InelasticCorrection] looping over QE sim " << nentries_QE << " events\n";
    for(Long64_t i=0;i<ch_QE.GetEntries();++i){ 

        ch_QE.GetEntry(i);

        if(/*vQE.ntrack<1 ||*/ abs(vQE.vz)>0.27 || vQE.eHCAL<c_.eHCAL_L || abs((vQE.ePS+vQE.eSH)/(vQE.trP)-1)>0.2 || vQE.ePS<0.2) continue;

        double dxq = vQE.dx;//dx_shifted_QE(vQE);

        if ((c_.dy_L<vQE.dy && vQE.dy<c_.dy_H) && (c_.dx_L<dxq && dxq<c_.dx_H) /*((pow((vQE.dy-0.0)/0.4,2)+pow((dxq-0.0)/0.3,2))<=1)*/) {
            hQE_W2_Neutrons.Fill(vQE.W2, vQE.weight);
            if (vQE.fnucl == 0) hQE_neutron_W2_Neutrons.Fill(vQE.W2, vQE.weight);  // see #3
            if (vQE.fnucl == 1) hQE_proton_W2_Neutrons.Fill(vQE.W2, vQE.weight);
        }

        if (((pow((vQE.dy-0.0)/0.4,2)+pow((dxq-0.0)/0.3,2))<=1) || ((pow((vQE.dy-0.0)/0.4,2)+pow((dxq-(c_.dx_P_L+c_.dx_P_H)/2)/0.3,2))<=1)) {
            hQE_W2.Fill(vQE.W2, vQE.weight);
            if (vQE.fnucl == 0) hQE_neutron_W2.Fill(vQE.W2, vQE.weight);  // see #3
            if (vQE.fnucl == 1) hQE_proton_W2.Fill(vQE.W2, vQE.weight);
        }


        if(/*vQE.ntrack<1 ||*/ abs(vQE.vz)>0.27 || vQE.eHCAL<c_.eHCAL_L || abs((vQE.ePS+vQE.eSH)/(vQE.trP)-1)>0.2 || vQE.ePS<0.2 ||
            (c_.W2_L>vQE.W2 || vQE.W2>c_.W2_H) || (c_.dy_L>vQE.dy || vQE.dy>c_.dy_H)) continue;

        if(vQE.fnucl == 0) {
            if(std::strcmp(kin_, "GEN3_He3") == 0){
                hQE_neutron.Fill(vQE.dx/*-0.1*/,vQE.weight);
            }
            else if(std::strcmp(kin_, "GEN4_He3") == 0){
                hQE_neutron.Fill(vQE.dx/*-0.07*/,vQE.weight);
            }
            else if(std::strcmp(kin_, "GEN4b_He3") == 0){
                hQE_neutron.Fill(vQE.dx/*-0.025*/,vQE.weight);
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
                hQE_proton.Fill(vQE.dx/*-0.02*/,vQE.weight);
            }
            else if(std::strcmp(kin_, "GEN4b_He3") == 0){
                hQE_proton.Fill(vQE.dx/*+0.05*/,vQE.weight);
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

    ////////////////////////// --- loop inelastic sim ////////////////////////////////////////////////
    //                                                                                              //
    //                                                                                              //
    //////////////////////////////////////////////////////////////////////////////////////////////////


    Long64_t nentries_inel=ch_inel.GetEntries();
    std::cout << "\n"<<"[InelasticCorrection] looping over Inelastic sim " << nentries_inel << " events\n";
    for(Long64_t i=0;i<ch_inel.GetEntries();++i){ 
        ch_inel.GetEntry(i);
    
        if(/*vInel.ntrack<1 ||*/ abs(vInel.vz)>0.27 || vInel.eHCAL<c_.eHCAL_L || abs((vInel.ePS+vInel.eSH)/(vInel.trP)-1)>0.2 || vInel.ePS<0.2) continue;

        double dxi = vInel.dx; //dx_shifted_Inel(vInel);

        if ( (c_.dy_L<vInel.dy && vInel.dy<c_.dy_H) && (c_.dx_L<dxi && dxi<c_.dx_H) /*((pow((vInel.dy-0.0)/0.4,2)+pow((dxi-0.0)/0.3,2))<=1)*/) {
            hInelastic_W2_Neutrons.Fill(vInel.W2, vInel.weight);
            hInelastic_W2_2_Neutrons.Fill(vInel.W2, vInel.weight);                    // see #3
        }

        if (((pow((vInel.dy-0.0)/0.4,2)+pow((dxi-0.0)/0.3,2))<=1) || ((pow((vInel.dy-0)/0.4,2)+pow((dxi-(c_.dx_P_L+c_.dx_P_H)/2)/0.3,2))<=1) ) {
            hInelastic_W2.Fill(vInel.W2, vInel.weight);
            hInelastic_W2_2.Fill(vInel.W2, vInel.weight);                    // see #3
        }

        if(/*vInel.ntrack<1 ||*/ abs(vInel.vz)>0.27 || vInel.eHCAL<c_.eHCAL_L || abs((vInel.ePS+vInel.eSH)/(vInel.trP)-1)>0.2 || vInel.ePS<0.2 ||
            (c_.W2_L>vInel.W2 || vInel.W2>c_.W2_H) || (c_.dy_L>vInel.dy || vInel.dy>c_.dy_H)) continue;    

        if(std::strcmp(kin_, "GEN3_He3") == 0){
            
            hInelastic.Fill(vInel.dx/*+0.4*/,vInel.weight);//+0.4 is for GEN3 its hard coded for now

            if (vInel.fnucl == 0) hInelastic_neutron.Fill(vInel.dx/*+0.4*/,vInel.weight);

            if (vInel.fnucl == 1) hInelastic_proton.Fill(vInel.dx/*+0.4*/,vInel.weight);
        }
        else if(std::strcmp(kin_, "GEN4_He3") == 0){
            
            hInelastic.Fill(vInel.dx/*+0.6*/,vInel.weight);//+0.4 is for GEN3 its hard coded for now

            if (vInel.fnucl == 0) hInelastic_neutron.Fill(vInel.dx/*+0.6*/,vInel.weight);

            if (vInel.fnucl == 1) hInelastic_proton.Fill(vInel.dx/*+0.6*/,vInel.weight);
        }
        else if(std::strcmp(kin_, "GEN4b_He3") == 0){
            
            //std::cout<<"debug here"<<'\n';

            hInelastic.Fill(vInel.dx+0.25,vInel.weight);//+0.4 is for GEN3 its hard coded for now

            if (vInel.fnucl == 0) hInelastic_neutron.Fill(vInel.dx/*+0.25*/,vInel.weight);

            if (vInel.fnucl == 1) hInelastic_proton.Fill(vInel.dx/*+0.25*/,vInel.weight);
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

    double Rn = neutron_events/(neutron_events+proton_events);
    double Rp = proton_events/(neutron_events+proton_events);

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

    //TH1D * h_combined_W2 =  performFitW2(&hData_W2,&hInelastic_W2,&hQE_proton_W2,&hQE_neutron_W2,parW0,parW1,Rnop);

    //h_combined_W2->Scale(hData_W2.Integral());

    //hQE_proton_W2.Scale(parW0*hData_W2.Integral()*1/hQE_proton_W2.Integral());
    //hQE_neutron_W2.Scale(parW0*Rnop*hData_W2.Integral()*1/hQE_neutron_W2.Integral());
    //hInelastic_W2.Scale(parW1*hData_W2.Integral()*1/hInelastic_W2.Integral());

    std::cout<<"parW0 : "<<parW0<<std::endl;
    std::cout<<"parW1 : "<<parW1<<std::endl;
    //std::cout<<"parW2 : "<<parW2<<std::endl;


    double alpha = 0.0;
    double delta_1 = 0.1;
    double delta_2 = 0.1;

    TH1D * h_combined_W2 =  performFitW2_1(&hData_W2,&hInelastic_W2,&hQE_proton_W2,&hQE_neutron_W2,alpha, delta_1 ,Rnop);

    double par0_2 = 0.5; double par1_2 = 0.5;

    TH1D * h_combined_W2_2 = performFitW2_2(&hData_W2,&hInelastic_W2_2,&hQE_W2,par0_2, par1_2, delta_2);

    // --- build shifted inelastic templates ---
    TH1D* hInelastic_W2_shift1 = make_shifted_unit(&hInelastic_W2,   "hInelastic_W2_shift1", delta_1);
    TH1D* hInelastic_W2_shift2 = make_shifted_unit(&hInelastic_W2_2, "hInelastic_W2_shift2", delta_2);

    // --- scale components to data counts ---
    const double N = hData_W2.Integral();
    const double denom = (1.0 + Rnop + alpha);

    // branch 1: (QEp + Rnop*QEn + alpha*Inel_shift)/denom
    h_combined_W2->Scale(N);
    hQE_proton_W2.Scale( N * (1.0/hQE_proton_W2.Integral()) / denom );
    hQE_neutron_W2.Scale( N * Rnop * (1.0/hQE_neutron_W2.Integral()) / denom );
    hInelastic_W2_shift1->Scale( N * alpha * (1.0/hInelastic_W2_shift1->Integral()) / denom );

    // branch 2: par0_2*QE + par1_2*Inel_shift
    h_combined_W2_2->Scale(N);
    hQE_W2.Scale(          N * (1.0/hQE_W2.Integral())           * par0_2 );
    hInelastic_W2_shift2->Scale( N * (1.0/hInelastic_W2_shift2->Integral())* par1_2 );

    std::cout<<"alpha : "<<alpha<<std::endl;
    std::cout<<"denom : "<<denom<<std::endl;
    std::cout<<"N : "<<N<<std::endl;

    std::cout<<"par0_2 : "<<par0_2<<std::endl;
    std::cout<<"par1_2 : "<<par1_2<<std::endl;

    std::cout<<"delta_1: "<<delta_1<<std::endl;
    std::cout<<"delta_2: "<<delta_2<<std::endl;

    /////////////////////////////////////W2 Fitting Neutrons only/////////////////////////////////////////////////////////

    double alpha_Neutrons = 0.0;
    double delta_1_Neutrons = 0.0;//0.1;
    double delta_2_Neutrons = 0.0;//0.1;

    TH1D * h_combined_W2_Neutrons =  performFitW2_1(&hData_W2_Neutrons,&hInelastic_W2_Neutrons,
        &hQE_proton_W2_Neutrons,&hQE_neutron_W2_Neutrons,alpha_Neutrons,delta_1_Neutrons,Rnop);

    double par0_2_Neutrons = 0.5; double par1_2_Neutrons = 0.5;

    TH1D * h_combined_W2_2_Neutrons = performFitW2_2(&hData_W2_Neutrons,&hInelastic_W2_2_Neutrons,&hQE_W2_Neutrons,
        par0_2_Neutrons,par1_2_Neutrons,delta_2_Neutrons);


    // --- build shifted inelastic templates ---
    TH1D* hInelastic_W2_shift1_Neutrons = make_shifted_unit(&hInelastic_W2_Neutrons,   "hInelastic_W2_shift1_Neutrons", 0);//delta_1_Neutrons);
    TH1D* hInelastic_W2_shift2_Neutrons = make_shifted_unit(&hInelastic_W2_2_Neutrons, "hInelastic_W2_shift2_Neutrons", 0);//delta_2_Neutrons);

    // --- scale components to data counts ---
    const double N_Neutrons = hData_W2_Neutrons.Integral();
    const double denom_Neutrons = (1.0 + Rnop + alpha_Neutrons);

    // branch 1: (QEp + Rnop*QEn + alpha*Inel_shift)/denom
    h_combined_W2_Neutrons->Scale(N_Neutrons);
    hQE_proton_W2_Neutrons.Scale( N_Neutrons * (1.0/hQE_proton_W2_Neutrons.Integral()) / denom_Neutrons );
    hQE_neutron_W2_Neutrons.Scale( N_Neutrons * Rnop * (1.0/hQE_neutron_W2_Neutrons.Integral()) / denom_Neutrons );
    hInelastic_W2_shift1_Neutrons->Scale( N_Neutrons * alpha_Neutrons /** (1.0/hInelastic_W2_shift1_Neutrons->Integral())*/ / denom_Neutrons );

    // branch 2: par0_2*QE + par1_2*Inel_shift
    h_combined_W2_2_Neutrons->Scale(N_Neutrons);
    hQE_W2_Neutrons.Scale(          N_Neutrons * (1.0/hQE_W2_Neutrons.Integral())           * par0_2_Neutrons );
    hInelastic_W2_shift2_Neutrons->Scale( N_Neutrons /*(1.0/hInelastic_W2_shift2_Neutrons->Integral())*/* par1_2_Neutrons );

    std::cout<<"alpha_Neutrons : "<<alpha_Neutrons<<std::endl;
    std::cout<<"denom_Neutrons : "<<denom_Neutrons<<std::endl;
    std::cout<<"N_Neutrons : "<<N_Neutrons<<std::endl;

    std::cout<<"par0_2_Neutrons : "<<par0_2_Neutrons<<std::endl;
    std::cout<<"par1_2_Neutrons : "<<par1_2_Neutrons<<std::endl;

    std::cout<<"delta_1_Neutrons : "<<delta_1_Neutrons<<std::endl;
    std::cout<<"delta_2_Neutrons : "<<delta_2_Neutrons<<std::endl;


    TH1D *h_QE_W2_split_N  =(TH1D*)hQE_W2_Neutrons.Clone("h_QE_W2_split_N"); 
    TH1D *h_QE_W2_split_P  =(TH1D*)hQE_W2_Neutrons.Clone("h_QE_W2_split_P");

    h_QE_W2_split_N->Scale(Rn);
    h_QE_W2_split_P->Scale(Rp); 


    ///////////////////plotting and printing////////////////////////////

    TCanvas *C = new TCanvas("c","c",2400,1500);
    TCanvas *C1 = new TCanvas("c1","c1",2400,1500);
    TCanvas *C2 = new TCanvas("c2","c2",2400,1500);
    TCanvas *C3 = new TCanvas("c3","c3",2400,1500);
    TCanvas *C4 = new TCanvas("c4","c4",2400,1500);
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
    //gPad->SetLogy();

    h_combined_W2->SetLineColor(3);
    hQE_proton_W2.SetLineColor(6);
    hQE_neutron_W2.SetLineColor(9);
    hInelastic_W2_shift1->SetLineColor(7);
    hData_W2.SetLineColor(kBlack);

    h_combined_W2->SetLineWidth(4);
    hQE_proton_W2.SetLineWidth(4);
    hQE_neutron_W2.SetLineWidth(4);
    hInelastic_W2_shift1->SetLineWidth(4);
    hData_W2.SetLineWidth(4);

    h_combined_W2->SetFillColorAlpha(19,0.1);
    h_combined_W2->SetFillStyle(3009);
    hQE_proton_W2.SetFillColorAlpha(6,0.5);
    hQE_proton_W2.SetFillStyle(3004);
    hQE_neutron_W2.SetFillColorAlpha(9,0.5);
    hQE_neutron_W2.SetFillStyle(3005);
    hInelastic_W2_shift1->SetFillColorAlpha(7,0.5);
    hInelastic_W2_shift1->SetFillStyle(3009);

    hData_W2.SetMarkerStyle(kFullSquare);

    hData_W2.Draw("PE1");
    hQE_proton_W2.Draw("hist same");
    hQE_neutron_W2.Draw("hist same");
    hInelastic_W2_shift1->Draw("hist same");   
    h_combined_W2->Draw("hist same"); 

    C1->cd(2);

    //gPad->SetLogy();

    h_combined_W2_2->SetLineColor(3);
    hQE_W2.SetLineColor(207);
    hInelastic_W2_shift2->SetLineColor(7);
    hData_W2.SetLineColor(kBlack);

    h_combined_W2_2->SetLineWidth(4);
    hQE_W2.SetLineWidth(4);
    hInelastic_W2_shift2->SetLineWidth(4);
    hData_W2.SetLineWidth(4);

    h_combined_W2_2->SetFillColorAlpha(19,0.1);
    h_combined_W2_2->SetFillStyle(3009);
    hQE_W2.SetFillColorAlpha(207,0.5);
    hQE_W2.SetFillStyle(3004);
    hInelastic_W2_shift2->SetFillColorAlpha(7,0.5);
    hInelastic_W2_shift2->SetFillStyle(3009);

    hData_W2.SetMarkerStyle(kFullSquare);

    hData_W2.Draw("PE1");
    hQE_W2.Draw("hist same");
    hInelastic_W2_shift2->Draw("hist same");   
    h_combined_W2_2->Draw("same");   

    C1->cd(3);
    hDxdy.Draw("COLZ");

    C1->cd(4);
    hDxdy_cut.Draw("COLZ");


    C2->Divide(2,2);
    C2->cd(1);
    h_combined_W2_Neutrons->SetLineColor(3);
    hQE_proton_W2_Neutrons.SetLineColor(6);
    hQE_neutron_W2_Neutrons.SetLineColor(9);
    hInelastic_W2_shift1_Neutrons->SetLineColor(7);
    hData_W2_Neutrons.SetLineColor(kBlack);

    h_combined_W2_Neutrons->SetLineWidth(4);
    hQE_proton_W2_Neutrons.SetLineWidth(4);
    hQE_neutron_W2_Neutrons.SetLineWidth(4);
    hInelastic_W2_shift1_Neutrons->SetLineWidth(4);
    hData_W2_Neutrons.SetLineWidth(4);

    h_combined_W2_Neutrons->SetFillColorAlpha(19,0.1);
    h_combined_W2_Neutrons->SetFillStyle(3009);
    hQE_proton_W2_Neutrons.SetFillColorAlpha(6,0.5);
    hQE_proton_W2_Neutrons.SetFillStyle(3004);
    hQE_neutron_W2_Neutrons.SetFillColorAlpha(9,0.5);
    hQE_neutron_W2_Neutrons.SetFillStyle(3005);
    hInelastic_W2_shift1_Neutrons->SetFillColorAlpha(7,0.5);
    hInelastic_W2_shift1_Neutrons->SetFillStyle(3009);

    hData_W2_Neutrons.SetMarkerStyle(kFullSquare);

    hData_W2_Neutrons.Draw("PE1");

    hQE_proton_W2_Neutrons.Draw("hist same");
    hQE_neutron_W2_Neutrons.Draw("hist same");
    hInelastic_W2_shift1_Neutrons->Draw("hist same");    
    h_combined_W2_Neutrons->Draw("hist same");

    C2->cd(2);

    //gPad->SetLogy();

    h_combined_W2_2_Neutrons->SetLineColor(3);
    hQE_W2_Neutrons.SetLineColor(207);
    hInelastic_W2_shift2_Neutrons->SetLineColor(7);
    hData_W2_Neutrons.SetLineColor(kBlack);
    h_QE_W2_split_N->SetLineColor(9);
    h_QE_W2_split_P->SetLineColor(6);

    h_combined_W2_2_Neutrons->SetLineWidth(4);
    hQE_W2_Neutrons.SetLineWidth(4);
    hInelastic_W2_shift2_Neutrons->SetLineWidth(4);
    hData_W2_Neutrons.SetLineWidth(4);
    h_QE_W2_split_N->SetLineWidth(4);
    h_QE_W2_split_P->SetLineWidth(4);

    h_combined_W2_2_Neutrons->SetFillColorAlpha(19,0.1);
    h_combined_W2_2_Neutrons->SetFillStyle(3009);
    hQE_W2_Neutrons.SetFillColorAlpha(207,0.5);
    hQE_W2_Neutrons.SetFillStyle(3004);
    hInelastic_W2_shift2_Neutrons->SetFillColorAlpha(7,0.5);
    hInelastic_W2_shift2_Neutrons->SetFillStyle(3009);
    h_QE_W2_split_N->SetFillColorAlpha(9,0.5);
    h_QE_W2_split_N->SetFillStyle(3005);
    h_QE_W2_split_P->SetFillColorAlpha(6,0.5);
    h_QE_W2_split_P->SetFillStyle(3005);

    hData_W2_Neutrons.SetMarkerStyle(kFullSquare);

    hData_W2_Neutrons.Draw("PE1");
    hQE_W2_Neutrons.Draw("hist same");
    hInelastic_W2_shift2_Neutrons->Draw("hist same");   
    h_combined_W2_2_Neutrons->Draw("hist same");
    h_QE_W2_split_N->Draw("hist same");
    h_QE_W2_split_P->Draw("hist same");

    C2->cd(3);
    hInelastic_W2_Neutrons.Draw("E");

    C2->cd(4);
    hInelastic_W2_2_Neutrons.Draw("E");

    C3->Divide(2,2);

    // --- Pad (1): Neutrons, (QEp + Rnp*QEn + α*Inel^Δ)/(1+Rnp+α)
    C3->cd(1);
    // ... your Draw() calls above ...
    auto legN1 = new TLegend(0.28, 0.28, 0.88, 0.88);  // x1,y1,x2,y2 in NDC
    legN1->SetBorderSize(0);
    legN1->SetFillStyle(0);
    legN1->SetTextSize(0.032);
    legN1->AddEntry(&hData_W2_Neutrons,             "Data", "p");
    legN1->AddEntry(h_combined_W2_Neutrons,         "Fit: (QEp + R_{np} QEn + #alpha Inel{#Delta})/(1+R_{np}+ #alpha)", "lf");
    legN1->AddEntry(&hQE_proton_W2_Neutrons,        "QE p (sim)", "f");
    legN1->AddEntry(&hQE_neutron_W2_Neutrons,       "QE n (sim)", "f");
    legN1->AddEntry(hInelastic_W2_shift1_Neutrons,  "Inelastic (shifted)", "f");
    legN1->Draw();

    // --- Pad (2): Neutrons, a_QE*QE + a_Inel*Inel^Δ
    C3->cd(2);
    // ... your Draw() calls above ...
    auto legN2 = new TLegend(0.28, 0.28, 0.88, 0.88);
    legN2->SetBorderSize(0);
    legN2->SetFillStyle(0);
    legN2->SetTextSize(0.032);
    legN2->AddEntry(&hData_W2_Neutrons,               "Data", "p");
    legN2->AddEntry(h_combined_W2_2_Neutrons,         "Fit: a_{QE} QE + a_{inel} Inel{#Delta}", "lf");
    legN2->AddEntry(&hQE_W2_Neutrons,                 "QE (sim)", "f");
    legN2->AddEntry(hInelastic_W2_shift2_Neutrons,    "Inelastic (sim shifted)", "f");
    legN2->AddEntry(h_QE_W2_split_N,    "QE Neutrons (sim)", "f");
    legN2->AddEntry(h_QE_W2_split_P,    "QE Protons (sim)", "f");

    legN2->Draw();


    C4->Divide(2,2);
    hDx_both.SetLineColor(kBlack);
    hDx_elastic.SetLineColor(kBlue);
    hDx_inelastic.SetLineColor(kRed);

    hDy_both.SetLineColor(kBlack);
    hDy_elastic.SetLineColor(kBlue);
    hDy_inelastic.SetLineColor(kRed);

    C4->cd(1);
    hDxdy_inelastic.Draw("COLZ");
    C4->cd(2);
    hDx_both.Draw();
    hDx_inelastic.Draw("same");
    hDx_elastic.Draw("same");
    C4->cd(3);
    hDy_both.Draw();
    hDy_inelastic.Draw("same");
    hDy_elastic.Draw("same");




    ///////////////////////// Build asymmetry graphs  /////////////////////////////////////////
    //                                                                                       //
    //                                                                                       //
    ///////////////////////////////////////////////////////////////////////////////////////////
     

    auto gA_dxonly   = makeAsymGraph(hW2_dxonly_pos,   hW2_dxonly_neg,   "gA_dxonly",   kBlack, kFullCircle);
    auto gA_dxdy     = makeAsymGraph(hW2_dxdy_pos,     hW2_dxdy_neg,     "gA_dxdy",     kRed, kFullCircle);
    auto gA_dxAntiDy = makeAsymGraph(hW2_dxAntiDy_pos, hW2_dxAntiDy_neg, "gA_dxAntiDy", kAzure+2, kOpenSquare);
    auto gA_dyAntiDx = makeAsymGraph(hW2_dyAntiDx_pos, hW2_dyAntiDx_neg, "gA_dyAntiDx", kGreen+2, kOpenSquare);
    

    // Optionally normalize the right-panel shapes to show just shape differences
    auto norm = [](TH1D& h){ double s=h.Integral(); if(s>0) h.Scale(1.0/s); };
    //norm(hW2_dxonly); norm(hW2_dxdy); norm(hW2_dxAntiDy);

    gA_dxonly->SetMinimum(-20);
    gA_dxonly->SetMaximum(20);
    gA_dxdy->SetMinimum(-20);
    gA_dxdy->SetMaximum(20);
    gA_dxAntiDy->SetMinimum(-20);
    gA_dxAntiDy->SetMaximum(20);
    gA_dyAntiDx->SetMinimum(-20);
    gA_dyAntiDx->SetMaximum(20);

    // Style for right panel
    hW2_dxonly.SetLineColor(kBlack);
    hW2_dxAntiDy.SetLineColor(kAzure+2);
    hW2_dyAntiDx.SetLineColor(kGreen+2);
    hW2_dxdy.SetLineColor(kRed);
    hW2_dxonly.SetLineWidth(3);
    hW2_dxAntiDy.SetLineWidth(3);
    hW2_dyAntiDx.SetLineWidth(3);
    hW2_dxdy.SetLineWidth(3);


    // Canvas like your example
    TCanvas* Casym = new TCanvas("Casym","Asym vs W2 and W2 shapes", 2400, 1500);
    Casym->Divide(2,2);

    // Left: A(W²)
    Casym->cd(1);

    Casym->SetGridy(true);    // turn ON horizontal grid
    Casym->SetGridx(false);   // keep vertical grid OFF

    gStyle->SetGridStyle(2);   // 2 = dashed (1=solid, 3=dotted, …)
    gStyle->SetGridWidth(1);   // thin lines; raise to 2-3 for thicker

    double yMin=-15, yMax=15; // adjust as you like
    TH2F* frame = new TH2F("Aframe","Asymmetry vs W^{2};W^{2} (GeV^{2});Asymmetry (%)",
                           20, -2, 8, 20, yMin, yMax);
    frame->Draw();


    gA_dxonly->SetMarkerSize(1);
    gA_dxAntiDy->SetMarkerSize(1);

    //gA_dxonly->Draw("AP");
    //gA_dxAntiDy->Draw("P SAME"); // include if you want three curves
    gA_dxdy->Draw("P SAME");  
    //gA_dyAntiDx->Draw("P SAME");

    // add/remove depending on what you want to compare
    
    auto legA = new TLegend(0.55,0.1,0.88,0.30);
    legA->SetBorderSize(0); legA->SetFillStyle(0);
    legA->AddEntry(gA_dxonly,   "|dx| cut only", "p");
    legA->AddEntry(gA_dxAntiDy, "|dx| & anti-dy", "p");
    legA->AddEntry(gA_dxdy,     "|dx| & |dy| cut", "p");
    legA->AddEntry(gA_dyAntiDx, "|dy| & anti-dx", "p");
    //legA->Draw();

    // Right: W² distributions
    Casym->cd(2);
    hW2_dxonly.Draw("HIST");
    hW2_dxAntiDy.Draw("HIST SAME");
    hW2_dxdy.Draw("HIST SAME");
    hW2_dyAntiDx.Draw("HIST SAME");
    auto legW = new TLegend(0.1,0.65,0.3,0.88);
    legW->SetBorderSize(0); legW->SetFillStyle(0);
    legW->AddEntry(&hW2_dxonly,   "dx cut but no dy cut", "l");
    legW->AddEntry(&hW2_dxAntiDy, "dx cut and anti dy cut", "l");
    legW->AddEntry(&hW2_dyAntiDx, "dy cut and anti dx cut", "l");
    legW->AddEntry(&hW2_dxdy,     "dx and dy cut", "l");
    //legW->Draw();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    auto gA_W2only   = makeAsymGraph(hDx_W2only_pos,   hDx_W2only_neg,   "gA_W2only",   kBlack, kFullCircle);
    auto gA_W2dy     = makeAsymGraph(hDx_W2dy_pos,     hDx_W2dy_neg,     "gA_W2dy",     kRed, kFullCircle);
    auto gA_W2AntiDy = makeAsymGraph(hDx_W2AntiDy_pos, hDx_W2AntiDy_neg, "gA_W2AntiDy", kAzure+2, kOpenSquare);
    auto gA_dyAntiW2 = makeAsymGraph(hDx_dyAntiW2_pos, hDx_dyAntiW2_neg, "gA_dyAntiW2", kGreen+2, kOpenSquare);
    

    // Optionally normalize the right-panel shapes to show just shape differences
    auto norm1 = [](TH1D& h){ double s=h.Integral(); if(s>0) h.Scale(1.0/s); };
    //norm(hW2_dxonly); norm(hW2_dxdy); norm(hW2_dxAntiDy);

    gA_W2only->SetMinimum(-20);
    gA_W2only->SetMaximum(20);
    gA_W2dy->SetMinimum(-20);
    gA_W2dy->SetMaximum(20);
    gA_W2AntiDy->SetMinimum(-20);
    gA_W2AntiDy->SetMaximum(20);
    gA_dyAntiW2->SetMinimum(-20);
    gA_dyAntiW2->SetMaximum(20);

    // Style for right panel
    hDx_W2only.SetLineColor(kBlack);
    hDx_W2AntiDy.SetLineColor(kAzure+2);
    hDx_dyAntiW2.SetLineColor(kGreen+2);
    hDx_W2dy.SetLineColor(kRed);
    hDx_W2only.SetLineWidth(3);
    hDx_W2AntiDy.SetLineWidth(3);
    hDx_dyAntiW2.SetLineWidth(3);
    hDx_W2dy.SetLineWidth(3);


    // Canvas like your example
    TCanvas* Casym1 = new TCanvas("Casym1","Asym vs Dx and Dx shapes", 2400, 1500);
    Casym1->Divide(2,2);

    // Left: A(W²)
    Casym1->cd(1);

    Casym1->SetGridy(true);    // turn ON horizontal grid
    Casym1->SetGridx(false);   // keep vertical grid OFF

    gStyle->SetGridStyle(2);   // 2 = dashed (1=solid, 3=dotted, …)
    gStyle->SetGridWidth(1);   // thin lines; raise to 2-3 for thicker

    //double yMin=-15, yMax=15; // adjust as you like
    TH2F* frame1 = new TH2F("Aframe","Asymmetry vs dx;dx (m);Asymmetry (%)",
                           20, -4, 3, 20, yMin, yMax);
    frame1->Draw();


    gA_W2only->SetMarkerSize(1);
    gA_W2AntiDy->SetMarkerSize(1);

    //gA_W2only->Draw("AP");
    //gA_W2AntiDy->Draw("P SAME"); // include if you want three curves
    gA_W2dy->Draw("P SAME");  
    //gA_dyAntiW2->Draw("P SAME");

    // add/remove depending on what you want to compare
    
    auto legA1 = new TLegend(0.55,0.1,0.88,0.30);
    legA1->SetBorderSize(0); legA->SetFillStyle(0);
    legA1->AddEntry(gA_W2only,   "W2 cut only", "p");
    legA1->AddEntry(gA_W2AntiDy, "W2 & anti-dy", "p");
    legA1->AddEntry(gA_W2dy,     "W2 & |dy| cut", "p");
    legA1->AddEntry(gA_dyAntiW2, "|dy| & anti-W2", "p");
    //legA1->Draw();

    // Right: W² distributions
    Casym1->cd(2);
    //hDx_W2only.Draw("HIST");
    //hDx_W2AntiDy.Draw("HIST SAME");
    hDx_W2dy.Draw("HIST SAME");
    //hDx_dyAntiW2.Draw("HIST SAME");
    auto legW1 = new TLegend(0.1,0.65,0.3,0.88);
    legW1->SetBorderSize(0); legW->SetFillStyle(0);
    legW1->AddEntry(&hDx_W2only,   "W2 cut but no dy cut", "l");
    legW1->AddEntry(&hDx_W2AntiDy, "W2 cut and anti dy cut", "l");
    legW1->AddEntry(&hDx_dyAntiW2, "dy cut and anti W2 cut", "l");
    legW1->AddEntry(&hDx_W2dy,     "W2 and dy cut", "l");
    //legW1->Draw();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // --- Stacked, x-aligned dx panels (TOP: dx shape, BOTTOM: A(dx)) ---
    gStyle->SetOptStat(0);

    const double xMin_dx = dxhist_low;
    const double xMax_dx = dxhist_high;
    const double yMin_A  = -15.0;   // adjust if you like
    const double yMax_A  =  15.0;

    TCanvas* Cdx = new TCanvas("Cdx","dx: shape + asym",1800,1200);

    // pads: same left/right margins so axes line up exactly
    double Left = 0.12, Right = 0.03, Ttop = 0.06, Btop = 0.02, Tbot = 0.02, Bbot = 0.18;
    TPad* padTop = new TPad("padTop","",0.0,0.36,1.0,1.0);
    TPad* padBot = new TPad("padBot","",0.0,0.00,1.0,0.36);
    padTop->SetLeftMargin(Left); padTop->SetRightMargin(Right); padTop->SetTopMargin(Ttop); padTop->SetBottomMargin(Btop);
    padBot->SetLeftMargin(Left); padBot->SetRightMargin(Right); padBot->SetTopMargin(Tbot); padBot->SetBottomMargin(Bbot);
    padTop->SetTickx(); padTop->SetTicky();
    padBot->SetTickx(); padBot->SetTicky();
    padTop->Draw(); padBot->Draw();

    // ---------------- TOP: dx histogram(s) ----------------
    padTop->cd();

    // pick the shape you want on top; here I use your "W2 & dy cut" version
    hDx_W2dy.GetXaxis()->SetRangeUser(xMin_dx, xMax_dx);

    // hide top x labels/titles so the pads can butt together cleanly
    hDx_W2dy.GetXaxis()->SetLabelSize(0);
    hDx_W2dy.GetXaxis()->SetTitleSize(0);
    hDx_W2dy.GetYaxis()->SetTitle("Counts");
    hDx_W2dy.GetYaxis()->SetTitleSize(0.05);
    hDx_W2dy.GetYaxis()->SetLabelSize(0.045);

    hDx_W2dy.SetLineWidth(2);
    hDx_W2dy.SetLineColor(kBlack);
    hDx_W2dy.SetFillStyle(1001);
    hDx_W2dy.SetFillColorAlpha(kViolet+2,0.3);
    hDx_W2dy.Draw("hist");

    // (optional) overlay other dx shapes, sideband shading, etc.
    // Example sideband shading:
    // TBox sb1(xMin_dx,0,-1.5, hDx_W2dy.GetMaximum()); sb1.SetFillColorAlpha(kRed,0.1); sb1.SetLineColor(0); sb1.Draw("same");
    // TBox sb2( 1.5,0, xMax_dx,hDx_W2dy.GetMaximum()); sb2.SetFillColorAlpha(kRed,0.1); sb2.SetLineColor(0); sb2.Draw("same");

    // ---------------- BOTTOM: A(dx) ----------------
    padBot->cd();

    // frame sets the x-range to MATCH the top pad
    TH2F* frame_dx = new TH2F("frame_dx",";dx (m);Asymmetry (%)",
                              100, xMin_dx, xMax_dx, 100, yMin_A, yMax_A);
    frame_dx->GetXaxis()->SetTitleSize(0.06);
    frame_dx->GetXaxis()->SetLabelSize(0.05);
    frame_dx->GetYaxis()->SetTitleSize(0.06);
    frame_dx->GetYaxis()->SetLabelSize(0.05);
    frame_dx->Draw();

    // your asymmetry graph (built earlier) goes on the frame
    gA_W2dy->SetMarkerSize(1.0);
    gA_W2dy->SetLineColor(kBlack);
    gA_W2dy->SetMarkerColor(kBlack);
    gA_W2dy->Draw("P SAME");

    // zero line helps the eye
    TLine* L0 = new TLine(xMin_dx,0,xMax_dx,0);
    L0->SetLineStyle(2);
    L0->Draw("same");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gA_W2dx     = makeAsymGraph(hDy_W2dx_pos,     hDy_W2dx_neg,     "gA_W2dx",     kRed, kFullCircle);


        // --- Stacked, x-aligned dy panels (TOP: dy shape, BOTTOM: A(dy)) ---
    gStyle->SetOptStat(0);

    const double xMin_dy = dyhist_low;
    const double xMax_dy =  dyhist_high;
    //const double yMin_A  = -15.0;   // adjust if you like
    //const double yMax_A  =  15.0;

    TCanvas* Cdy = new TCanvas("Cdy","dy: shape + asym",1800,1200);

    // pads: same left/right margins so axes line up exactly
    //double Left = 0.12, Right = 0.03, Ttop = 0.06, Btop = 0.02, Tbot = 0.02, Bbot = 0.18;
    TPad* padTop_dy = new TPad("padTop_dy","",0.0,0.36,1.0,1.0);
    TPad* padBot_dy = new TPad("padBot_dy","",0.0,0.00,1.0,0.36);
    padTop_dy->SetLeftMargin(Left); padTop_dy->SetRightMargin(Right); padTop_dy->SetTopMargin(Ttop); padTop_dy->SetBottomMargin(Btop);
    padBot_dy->SetLeftMargin(Left); padBot_dy->SetRightMargin(Right); padBot_dy->SetTopMargin(Tbot); padBot_dy->SetBottomMargin(Bbot);
    padTop_dy->SetTickx(); padTop_dy->SetTicky();
    padBot_dy->SetTickx(); padBot_dy->SetTicky();
    padTop_dy->Draw(); padBot_dy->Draw();

    // ---------------- TOP: dx histogram(s) ----------------
    padTop_dy->cd();

    // pick the shape you want on top; here I use your "W2 & dy cut" version
    hDy_W2dx.GetXaxis()->SetRangeUser(xMin_dy, xMax_dy);

    // hide top x labels/titles so the pads can butt together cleanly
    hDy_W2dx.SetLineColor(kBlack);
    hDy_W2dx.GetXaxis()->SetLabelSize(0);
    hDy_W2dx.GetXaxis()->SetTitleSize(0);
    hDy_W2dx.GetYaxis()->SetTitle("Counts");
    hDy_W2dx.GetYaxis()->SetTitleSize(0.05);
    hDy_W2dx.GetYaxis()->SetLabelSize(0.045);

    hDy_W2dx.SetLineWidth(2);
    hDy_W2dx.SetLineColor(kBlack);
    hDy_W2dx.SetFillStyle(1001);
    hDy_W2dx.SetFillColorAlpha(kViolet+2,0.3);
    hDy_W2dx.Draw("hist");

    // (optional) overlay other dx shapes, sideband shading, etc.
    // Example sideband shading:
    // TBox sb1(xMin_dx,0,-1.5, hDx_W2dy.GetMaximum()); sb1.SetFillColorAlpha(kRed,0.1); sb1.SetLineColor(0); sb1.Draw("same");
    // TBox sb2( 1.5,0, xMax_dx,hDx_W2dy.GetMaximum()); sb2.SetFillColorAlpha(kRed,0.1); sb2.SetLineColor(0); sb2.Draw("same");

    // ---------------- BOTTOM: A(dx) ----------------
    padBot_dy->cd();

    // frame sets the x-range to MATCH the top pad
    TH2F* frame_dy = new TH2F("frame_dy",";dy (m);Asymmetry (%)",
                              100, xMin_dy, xMax_dy, 100, yMin_A, yMax_A);
    frame_dy->GetXaxis()->SetTitleSize(0.06);
    frame_dy->GetXaxis()->SetLabelSize(0.05);
    frame_dy->GetYaxis()->SetTitleSize(0.06);
    frame_dy->GetYaxis()->SetLabelSize(0.05);
    frame_dy->Draw();

    // your asymmetry graph (built earlier) goes on the frame
    gA_W2dx->SetMarkerSize(1.0);
    gA_W2dx->SetLineColor(kBlack);
    gA_W2dx->SetMarkerColor(kBlack);
    gA_W2dx->Draw("P SAME");


    // zero line helps the eye
    TLine* Ldy0 = new TLine(xMin_dy,0,xMax_dy,0);
    Ldy0->SetLineStyle(2);
    Ldy0->Draw("same");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    gStyle->SetOptStat(0);

    const double xMin_W2 =  W2hist_low;
    const double xMax_W2 =  W2hist_high;
    //const double yMin_A  = -15.0;   // adjust if you like
    //const double yMax_A  =  15.0;

    TCanvas* CW2 = new TCanvas("CW2","W^{2}: shape + asym",1800,1200);

    // pads: same left/right margins so axes line up exactly
    //double Left = 0.12, Right = 0.03, Ttop = 0.06, Btop = 0.02, Tbot = 0.02, Bbot = 0.18;
    TPad* padTop_W2 = new TPad("padTop_W2","",0.0,0.36,1.0,1.0);
    TPad* padBot_W2 = new TPad("padBot_W2","",0.0,0.00,1.0,0.36);
    padTop_W2->SetLeftMargin(Left); padTop_W2->SetRightMargin(Right); padTop_W2->SetTopMargin(Ttop); padTop_W2->SetBottomMargin(Btop);
    padBot_W2->SetLeftMargin(Left); padBot_W2->SetRightMargin(Right); padBot_W2->SetTopMargin(Tbot); padBot_W2->SetBottomMargin(Bbot);
    padTop_W2->SetTickx(); padTop_W2->SetTicky();
    padBot_W2->SetTickx(); padBot_W2->SetTicky();
    padTop_W2->Draw(); padBot_W2->Draw();

    // ---------------- TOP: dx histogram(s) ----------------
    padTop_W2->cd();

    // pick the shape you want on top; here I use your "W2 & dy cut" version
    hW2_dxdy.GetXaxis()->SetRangeUser(xMin_W2, xMax_W2);

    // hide top x labels/titles so the pads can butt together cleanly
    hW2_dxdy.GetXaxis()->SetLabelSize(0);
    hW2_dxdy.GetXaxis()->SetTitleSize(0);
    hW2_dxdy.GetYaxis()->SetTitle("Counts");
    hW2_dxdy.GetYaxis()->SetTitleSize(0.05);
    hW2_dxdy.GetYaxis()->SetLabelSize(0.045);

    hW2_dxdy.SetLineWidth(2);
    hW2_dxdy.SetLineColor(kBlack);
    hW2_dxdy.SetFillStyle(1001);
    hW2_dxdy.SetFillColorAlpha(kViolet+2,0.3);
    hW2_dxdy.Draw("hist");

    // (optional) overlay other dx shapes, sideband shading, etc.
    // Example sideband shading:
    // TBox sb1(xMin_dx,0,-1.5, hDx_W2dy.GetMaximum()); sb1.SetFillColorAlpha(kRed,0.1); sb1.SetLineColor(0); sb1.Draw("same");
    // TBox sb2( 1.5,0, xMax_dx,hDx_W2dy.GetMaximum()); sb2.SetFillColorAlpha(kRed,0.1); sb2.SetLineColor(0); sb2.Draw("same");

    // ---------------- BOTTOM: A(dx) ----------------
    padBot_W2->cd();

    // frame sets the x-range to MATCH the top pad
    TH2F* frame_W2 = new TH2F("frame_W2",";W^{2} (GeV^{2});Asymmetry (%)",
                              100, xMin_W2, xMax_W2, 100, yMin_A, yMax_A);
    frame_W2->GetXaxis()->SetTitleSize(0.06);
    frame_W2->GetXaxis()->SetLabelSize(0.05);
    frame_W2->GetYaxis()->SetTitleSize(0.06);
    frame_W2->GetYaxis()->SetLabelSize(0.05);
    frame_W2->Draw();

    // your asymmetry graph (built earlier) goes on the frame
    gA_dxdy->SetMarkerSize(1.0);
    gA_dxdy->SetLineColor(kBlack);
    gA_dxdy->SetMarkerColor(kBlack);
    gA_dxdy->Draw("P SAME");

    // zero line helps the eye
    TLine* LW20 = new TLine(xMin_W2,0,xMax_W2,0);
    LW20->SetLineStyle(2);
    LW20->Draw("same");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Save
    Cdx->Print(Form("images/%s/Asymmetry_vs_Dx_stacked_%s.png", kin_, kin_));
    Cdy->Print(Form("images/%s/Asymmetry_vs_Dy_stacked_%s.png", kin_, kin_));
    CW2->Print(Form("images/%s/Asymmetry_vs_W2_stacked_%s.png", kin_, kin_));

    // Save

    C->Print(Form("images/%s/InelasticCorrection_%s.png",kin_,kin_));
    C1->Print(Form("images/%s/InelasticCorrection_W2_%s.png",kin_,kin_));
    C2->Print(Form("images/%s/InelasticCorrection_W2_Neutrons_%s.png",kin_,kin_));
    C3->Print(Form("images/%s/InelasticCorrection_W2_Neutrons_legend%s.png",kin_,kin_));
    C4->Print(Form("images/%s/InelasticDxDy%s.png",kin_,kin_));
    Casym->Print(Form("images/%s/Asymmetry_vs_W2_%s.png", kin_, kin_));
    Casym1->Print(Form("images/%s/Asymmetry_vs_Dx_%s.png", kin_, kin_));

}

