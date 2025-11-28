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
#include <TH3F.h>
#include <TH3D.h>
#include <TF2.h>
#include <TF3.h>
#include <TStyle.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TGaxis.h>

using std::string; using std::vector; using std::unordered_map;

#include <limits>

// ---------- SSR helpers (W2 neutrons) ----------------------------------------
namespace {
  // Unweighted SSR over [binLo, binHi] inclusive.
  double SumSquaredResiduals(const TH1* data, const TH1* model,
                             int binLo = 1, int binHi = -1)
  {
    if (!data || !model) return std::numeric_limits<double>::quiet_NaN();
    if (data->GetNbinsX() != model->GetNbinsX()) return std::numeric_limits<double>::quiet_NaN();

    const int nb = data->GetNbinsX();
    if (binHi < 0) binHi = nb;

    // (Optional) axis compatibility check
    for (int i = 1; i <= nb + 1; ++i) {
      if (std::abs(data->GetBinLowEdge(i) - model->GetBinLowEdge(i)) > 1e-9 ||
          std::abs(data->GetBinWidth(i)   - model->GetBinWidth(i))   > 1e-12) {
        return std::numeric_limits<double>::quiet_NaN();
      }
    }

    binLo = std::max(1,  binLo);
    binHi = std::min(nb, binHi);

    double ssr = 0.0;
    for (int i = binLo; i <= binHi; ++i) {
      const double r = data->GetBinContent(i) - model->GetBinContent(i);
      ssr += r * r;
    }
    return ssr;
  }

  // Convenience: SSR over axis coords [W2_low, W2_high].
  double SumSquaredResidualsInRange(const TH1* data, const TH1* model,
                                    double W2_low, double W2_high)
  {
    if (!data) return std::numeric_limits<double>::quiet_NaN();
    const TAxis* ax = data->GetXaxis();
    int iLo = std::max(1,            ax->FindBin(W2_low));
    int iHi = std::min(ax->GetNbins(), ax->FindBin(W2_high));
    return SumSquaredResiduals(data, model, iLo, iHi);
  }

  double Chi2InRange(const TH1* data, const TH1* model, double W2_low, double W2_high) {
      if (!data || !model) return std::numeric_limits<double>::quiet_NaN();
      const TAxis* ax = data->GetXaxis();
      int iLo = std::max(1, ax->FindBin(W2_low));
      int iHi = std::min(ax->GetNbins(), ax->FindBin(W2_high));

      double chi2 = 0.0;
      for (int i = iLo; i <= iHi; ++i) {
        const double d = data->GetBinContent(i);
        const double m = model->GetBinContent(i);
        // Use model or data for variance; pick one consistently
        const double var = std::max(m, 1.0); // or std::max(d,1.0)
        const double r = d - m;
        chi2 += (r * r) / var;
      }
      return chi2;
    }

} // anonymous namespace


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
                                      /* optional: add outputs for shifts if you like*/
                                       , double& dx_p_out, double& dx_n_out, double& dx_inel_out 
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
    const double maxShift = 0.05;  // units = your dx units; e.g., 0.02 m or 2 cm (adjust!)
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
    dx_p_out    = f->GetParameter(3);
    dx_n_out    = f->GetParameter(4);
    dx_inel_out = f->GetParameter(5);

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


// --- 2D template shape fit on (dx,dy) ---------------------------------------
// Model:  A * [  P(x-Δx_p, y-Δy_p) + rNP * N(x-Δx_n, y-Δy_n) + rI * I(x-Δx_i, y-Δy_i)  ]
// Histos should be unit-normalized PDFs (method does that defensively).
TH2D* InelasticCorrection::performFitDxDy(
    TH2D* hD2, TH2D* hInel2, TH2D* hQEp2, TH2D* hQEn2,
    double& A, double& rNP, double& rI,
    double& dx_p, double& dy_p,
    double& dx_n, double& dy_n,
    double& dx_i, double& dy_i
){
    auto clone_norm2D = [](TH2D* src, const char* name){
        TH2D* h = (TH2D*)src->Clone(name);
        double S = h->Integral();
        if (S > 0) h->Scale(1.0/S);
        return h;
    };

    TH2D* d    = clone_norm2D(hD2,    Form("hD2_%p",    hD2));
    TH2D* inel = clone_norm2D(hInel2, Form("hInel2_%p", hInel2));
    TH2D* qep  = clone_norm2D(hQEp2,  Form("hQEp2_%p",  hQEp2));
    TH2D* qen  = clone_norm2D(hQEn2,  Form("hQEn2_%p",  hQEn2));

    const double xmin = d->GetXaxis()->GetXmin();
    const double xmax = d->GetXaxis()->GetXmax();
    const double ymin = d->GetYaxis()->GetXmin();
    const double ymax = d->GetYaxis()->GetXmax();

    auto interp2D = [&](TH2D* h, double x, double y){
        if (x <= xmin) x = xmin + 1e-9;
        if (x >= xmax) x = xmax - 1e-9;
        if (y <= ymin) y = ymin + 1e-9;
        if (y >= ymax) y = ymax - 1e-9;
        return h->Interpolate(x,y); // bilinear
    };

    // Params:
    // [0]=A, [1]=rNP, [2]=rI, [3..4]=dx,dy for P; [5..6]=dx,dy for N; [7..8]=dx,dy for Inel
    TF2* f2 = new TF2(Form("fDXDY_%p", d),
        [&](double* xx, double* p){
            const double x = xx[0], y = xx[1];
            const double P = interp2D(qep,  x - p[3], y - p[4]);
            const double N = interp2D(qen,  x - p[5], y - p[6]);
            const double I = interp2D(inel, x - p[7], y - p[8]);
            return p[0] * ( P + p[1]*N + p[2]*I );
        },
        xmin, xmax, ymin, ymax, 9);

    f2->SetParNames("A","rNP","rI","dx_p","dy_p","dx_n","dy_n","dx_i","dy_i");
    f2->SetParameters(1.0, 0.3, 0.1, 0,0, 0,0, 0,0);

    f2->SetParLimits(0, 0.0,  50.0);   // A
    f2->SetParLimits(1, 0.0,   5.0);   // rNP
    f2->SetParLimits(2, 0.0,   5.0);   // rI
    const double maxShiftX = 0.1;     // tune for your dx units
    const double maxShiftY = 0.2;     // tune for your dy units
    f2->SetParLimits(3, -maxShiftX, +maxShiftX); // dx_p
    f2->SetParLimits(4, -maxShiftY, +maxShiftY); // dy_p
    f2->SetParLimits(5, -maxShiftX, +maxShiftX); // dx_n
    f2->SetParLimits(6, -maxShiftY, +maxShiftY); // dy_n
    f2->SetParLimits(7, -maxShiftX, +maxShiftX); // dx_i
    f2->SetParLimits(8, -maxShiftY, +maxShiftY); // dy_i

    // Fit quietly (least-squares over 2D bins)
    d->Fit(f2, "RQ0");

    // Extract
    A    = f2->GetParameter(0);
    rNP  = f2->GetParameter(1);
    rI   = f2->GetParameter(2);
    dx_p = f2->GetParameter(3); dy_p = f2->GetParameter(4);
    dx_n = f2->GetParameter(5); dy_n = f2->GetParameter(6);
    dx_i = f2->GetParameter(7); dy_i = f2->GetParameter(8);

    std::cout << std::fixed << std::setprecision(4)
              << "[DxDyFit] A="<<A
              << " rNP="<<rNP<<" rI="<<rI
              << "  (dx,dy)_p=("<<dx_p<<","<<dy_p<<")"
              << " (dx,dy)_n=("<<dx_n<<","<<dy_n<<")"
              << " (dx,dy)_inel=("<<dx_i<<","<<dy_i<<")\n";

    // Build best-fit combined TH2
    TH2D* comb = (TH2D*)inel->Clone(Form("hComb2D_%p", d));
    comb->Reset();

    for (int ix=1; ix<=comb->GetNbinsX(); ++ix){
        for (int iy=1; iy<=comb->GetNbinsY(); ++iy){
            const double x = comb->GetXaxis()->GetBinCenter(ix);
            const double y = comb->GetYaxis()->GetBinCenter(iy);
            const double P = interp2D(qep,  x - dx_p, y - dy_p);
            const double N = interp2D(qen,  x - dx_n, y - dy_n);
            const double I = interp2D(inel, x - dx_i, y - dy_i);
            comb->SetBinContent(ix, iy, A*(P + rNP*N + rI*I));
        }
    }

    return comb;
}


// --- 3D template shape fit on (dx,dy,W2) ------------------------------------
// Model:  A * [ P(dx-Δx_p,dy-Δy_p,W2-ΔW2_p)
//             + rNP * N(dx-Δx_n,dy-Δy_n,W2-ΔW2_n)
//             + rI  * I(dx-Δx_i,dy-Δy_i,W2-ΔW2_i) ]
// Histos are defensively unit-normalized here.
TH3D* InelasticCorrection::performFitDxDyW2(
    TH3D* hD3, TH3D* hInel3, TH3D* hQEp3, TH3D* hQEn3,
    double& A, double& rNP, double& rI,
    double& dx_p, double& dy_p, double& dW2_p,
    double& dx_n, double& dy_n, double& dW2_n,
    double& dx_i, double& dy_i, double& dW2_i
){
    auto clone_norm3D = [](TH3D* src, const char* nm){
        TH3D* h=(TH3D*)src->Clone(nm);
        double S=h->Integral();
        if(S>0) h->Scale(1.0/S);
        return h;
    };
    TH3D* d    = clone_norm3D(hD3,    Form("hD3_%p",    hD3));
    TH3D* inel = clone_norm3D(hInel3, Form("hInel3_%p", hInel3));
    TH3D* qep  = clone_norm3D(hQEp3,  Form("hQEp3_%p",  hQEp3));
    TH3D* qen  = clone_norm3D(hQEn3,  Form("hQEn3_%p",  hQEn3));

    const double xmin=d->GetXaxis()->GetXmin(), xmax=d->GetXaxis()->GetXmax();
    const double ymin=d->GetYaxis()->GetXmin(), ymax=d->GetYaxis()->GetXmax();
    const double zmin=d->GetZaxis()->GetXmin(), zmax=d->GetZaxis()->GetXmax();

    auto interp3D = [&](TH3D* h, double x, double y, double z){
        // clamp into range to avoid Interpolate() edge issues
        if (x<=xmin) x=xmin+1e-9; if (x>=xmax) x=xmax-1e-9;
        if (y<=ymin) y=ymin+1e-9; if (y>=ymax) y=ymax-1e-9;
        if (z<=zmin) z=zmin+1e-9; if (z>=zmax) z=zmax-1e-9;
        return h->Interpolate(x,y,z);
    };

    // Params:
    // [0]=A, [1]=rNP, [2]=rI,
    // [3..5]=dx,dy,dW2 for P; [6..8]=dx,dy,dW2 for N; [9..11]=dx,dy,dW2 for Inel
    TF3* f3 = new TF3(Form("fDXDYW2_%p", d),
        [&](double* xx, double* p){
            const double x=xx[0], y=xx[1], z=xx[2];
            const double P = interp3D(qep,  x - p[3], y - p[4], z - p[5]);
            const double N = interp3D(qen,  x - p[6], y - p[7], z - p[8]);
            const double I = interp3D(inel, x - p[9], y - p[10], z - p[11]);
            return p[0]*( P + p[1]*N + p[2]*I );
        },
        xmin,xmax, ymin,ymax, zmin,zmax, 12);

    // ---- Set parameter names (12 params) ----
    const char* pnames[12] = {
      "A","rNP","rI",
      "dx_p","dy_p","dW2_p",
      "dx_n","dy_n","dW2_n",
      "dx_i","dy_i","dW2_i"
    };
    for (int i = 0; i < 12; ++i) f3->SetParName(i, pnames[i]);

    // ---- Initial guesses (use array form for 12 params) ----
    double pinits[12] = {
      1.0, 0.3, 0.1,   // A, rNP, rI
      0.0, 0.0, 0.0,   // dx_p, dy_p, dW2_p
      0.0, 0.0, 0.0,   // dx_n, dy_n, dW2_n
      0.0, 0.0, 0.0    // dx_i, dy_i, dW2_i
    };
    f3->SetParameters(pinits);

    // limits (tune to your experiment units)
    const double maxShiftX = 0.2; // dx (m)
    const double maxShiftY = 0.06; // dy (m)
    const double maxShiftW = 0.05; // W2 (GeV^2)
    f3->SetParLimits(0, 0.0,  100.0);   // A
    f3->SetParLimits(1, 0.0,  100.0);   // rNP
    f3->SetParLimits(2, 0.0,  100.0);   // rI

    f3->SetParLimits(3, -maxShiftX, +maxShiftX);
    f3->SetParLimits(4, -maxShiftY, +maxShiftY);
    f3->SetParLimits(5, -maxShiftW, +maxShiftW);

    f3->SetParLimits(6, -maxShiftX, +maxShiftX);
    f3->SetParLimits(7, -maxShiftY, +maxShiftY);
    f3->SetParLimits(8, -maxShiftW, +maxShiftW);

    f3->SetParLimits(9,  -maxShiftX, +maxShiftX);
    f3->SetParLimits(10, -maxShiftY, +maxShiftY);
    f3->SetParLimits(11, -maxShiftW, +maxShiftW);

    // Quiet least-squares fit over 3D bins
    d->Fit(f3,"RQ0");

    // output
    A    = f3->GetParameter(0);
    rNP  = f3->GetParameter(1);
    rI   = f3->GetParameter(2);

    dx_p = f3->GetParameter(3);  dy_p = f3->GetParameter(4);  dW2_p = f3->GetParameter(5);
    dx_n = f3->GetParameter(6);  dy_n = f3->GetParameter(7);  dW2_n = f3->GetParameter(8);
    dx_i = f3->GetParameter(9);  dy_i = f3->GetParameter(10); dW2_i = f3->GetParameter(11);

    std::cout<<std::fixed<<std::setprecision(4)
             <<"[DxDyW2Fit] A="<<A<<" rNP="<<rNP<<" rI="<<rI
             <<"  (dx,dy,dW2)_p=("<<dx_p<<","<<dy_p<<","<<dW2_p<<")"
             <<" (dx,dy,dW2)_n=("<<dx_n<<","<<dy_n<<","<<dW2_n<<")"
             <<" (dx,dy,dW2)_i=("<<dx_i<<","<<dy_i<<","<<dW2_i<<")\n";

    // Build combined TH3 with shifts
    TH3D* comb=(TH3D*)inel->Clone(Form("hComb3D_%p", d));
    comb->Reset();
    for(int ix=1; ix<=comb->GetNbinsX(); ++ix){
      const double x=comb->GetXaxis()->GetBinCenter(ix);
      for(int iy=1; iy<=comb->GetNbinsY(); ++iy){
        const double y=comb->GetYaxis()->GetBinCenter(iy);
        for(int iz=1; iz<=comb->GetNbinsZ(); ++iz){
          const double z=comb->GetZaxis()->GetBinCenter(iz);
          const double P = interp3D(qep,  x - dx_p, y - dy_p, z - dW2_p);
          const double N = interp3D(qen,  x - dx_n, y - dy_n, z - dW2_n);
          const double I = interp3D(inel, x - dx_i, y - dy_i, z - dW2_i);
          comb->SetBinContent(ix,iy,iz, A*(P + rNP*N + rI*I));
        }
      }
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
            
            const double qep = pi->Interpolate(x[0]);      // was FindBin/GetBinContent
            const double qen = ni->Interpolate(x[0]);      // was FindBin/GetBinContent
            const double ine = ii->Interpolate(x[0]); //- dL); // your inel_shift does this too

            //const int bp = pi->FindBin(x[0]);
            //const int bn = ni->FindBin(x[0]);
            //const double qep = pi->GetBinContent(bp);
            //const double qen = ni->GetBinContent(bn);
            //const double ine = ii->GetBinContent(ii->FindBin(x[0]));//inel_shift(x[0], dL);
            //const double ine = inel_shift(x[0], dL);
            // simple normalized mixture (good if Δ is small or range wide)
            

            return (qep + Rnop*qen + a*ine)/(1.0 + Rnop + a);
        }, xmin, xmax, 1);

    f->SetParameter(0, 0.5);     // alpha start
    f->SetParLimits(0, 0.0, 40.0); // tune to your expected alpha range
    //f->SetParameter(1, 0.0);     // delta start (GeV^2)
    //f->SetParLimits(1,-0.4, 0.4); // tune to your expected shift range

    d->Fit(f, "RQ");

    alpha = f->GetParameter(0);
    //delta = f->GetParameter(1);

    // Build combined unit-area PDF using fitted Δ
    TH1D* comb=(TH1D*)ii->Clone("hCombW2"); comb->Reset();
    for (int i=1;i<=comb->GetNbinsX();++i){
        const double x = comb->GetXaxis()->GetBinCenter(i);
        //const double qep = pi->GetBinContent(pi->FindBin(x));
        //const double qen = ni->GetBinContent(ni->FindBin(x));
        //const double ine = ii->GetBinContent(ii->FindBin(x));//inel_shift(x, delta);
        //const double ine = inel_shift(x, delta);
        double qep = pi->Interpolate(x);           // was GetBinContent(FindBin(x))
        double qen = ni->Interpolate(x);
        double ine = ii->Interpolate(x); //- delta);   // or inel_shift(x, delta)
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
    f->SetParLimits(0, 0.0, 1.0);
    f->SetParLimits(1, 0.0, 1.0);
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
    double Api  = V(pionMap,"A_pi");

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

    double W2_hist_upper_limit = -0.05;//c_.W2_H;//0.0;
    double W2_hist_lower_limit = 1.3;//c_.W2_L;//0.0;

    // if(std::strcmp(kin_, "GEN3_He3") == 0){
    //     W2_hist_upper_limit = 1.45;
    //     W2_hist_lower_limit = -1; 
    // }
    // else if( std::strcmp(kin_, "GEN4_He3") == 0){
    //     W2_hist_upper_limit = 1.45;
    //     W2_hist_lower_limit = -1; 
    // }
    // else if(std::strcmp(kin_, "GEN4b_He3") == 0 ){
    //     W2_hist_upper_limit = 1.45;
    //     W2_hist_lower_limit = -1; 
    // }
    // else{
    //     W2_hist_upper_limit = 1.45;
    //     W2_hist_lower_limit = -1; 
    // }

    const int NBW2 = 100; // match your other W² binning

    TH1D hData_W2("hData_W2",  "W2 data; W^{2} (GeV^{2})", NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hData_W2_Neutrons("hData_W2_Neutrons",  "W2 data; W^{2} (GeV^{2})", NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hQE_W2("hQE_W2", "W^{2} QE sim",          NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hQE_W2_Neutrons("hQE_W2_Neutrons", "W^{2} QE sim",          NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hQE_proton_W2("hQE_proton_W2", "W^{2} QE sim protons",          NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hQE_proton_W2_Neutrons("hQE_proton_W2_Neutrons", "W^{2} QE sim protons",          NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hQE_neutron_W2("hQE_neutron_W2", "W^{2} QE sim neutrons",          NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hQE_neutron_W2_Neutrons("hQE_neutron_W2_Neutrons", "W^{2} QE sim neutrons",          NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hInelastic_W2("hInelastic_W2",  "W^{2} inelastic sim",        NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hInelastic_W2_Neutrons("hInelastic_W2_Neutrons",  "W^{2} inelastic sim",        NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hInelastic_W2_Neutrons_eHCALcut_1("hInelastic_W2_Neutrons_eHCALcut_1",  "W^{2} inelastic sim",        NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hInelastic_W2_Neutrons_eHCALcut_2("hInelastic_W2_Neutrons_eHCALcut_2",  "W^{2} inelastic sim",        NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hInelastic_W2_Neutrons_eHCALcut_3("hInelastic_W2_Neutrons_eHCALcut_3",  "W^{2} inelastic sim",        NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hInelastic_W2_Neutrons_eHCALcut_4("hInelastic_W2_Neutrons_eHCALcut_4",  "W^{2} inelastic sim",        NBW2, W2_hist_lower_limit, W2_hist_upper_limit);


    TH1D hInelastic_W2_2("hInelastic_W2_2",  "W^{2} inelastic sim",        NBW2, W2_hist_lower_limit, W2_hist_upper_limit);
    TH1D hInelastic_W2_2_Neutrons("hInelastic_W2_2_Neutrons",  "W^{2} inelastic sim",        NBW2, W2_hist_lower_limit, W2_hist_upper_limit);

    TH2D hDxdy("hDxdy", "dxdy distribution ; dx(m)",100, c_.dy_L-0.3, c_.dy_H+0.3,100,-4,3);
    TH2D hDxdy_cut("hDxdy_cut", "dxdy distribution ; dx(m)",100,-4,3,100,-4,3);

    TH2D hDxdy_inelastic("hDxdy_inelastic", "dxdy distribution ; dy(m); dx(m)",100,-4,3,100,-4,3);

    // 2D templates for QE p/n and inelastic (same binning as hDxdy)
    TH2D hQE_dxdy_p ("hQE_dxdy_p",  "QE p: dy vs dx;dy (m);dx (m)", 100, c_.dy_L-0.3, c_.dy_H+0.3, 100, -4, 3);
    TH2D hQE_dxdy_n ("hQE_dxdy_n",  "QE n: dy vs dx;dy (m);dx (m)", 100,  c_.dy_L-0.3, c_.dy_H+0.3, 100, -4, 3);
    TH2D hInel_dxdy ("hInel_dxdy",  "Inelastic: dy vs dx;dy (m);dx (m)", 100,  c_.dy_L-0.3, c_.dy_H+0.3, 100, -4, 3);

    // Bin ranges: reuse your dx, dy, W2 windows
    const int NBX = 100, NBY = 100, NBZ = 100; // tune rebinning as needed
    const double dyhist_low3 = c_.dy_L-0.3, dyhist_high3 = c_.dy_H+0.3;
    const double dxhist_low3 = -0.7 , dxhist_high3 = 0.7;
    const double W2hist_low3 = c_.W2_L - 0.05 , W2hist_high3 = c_.W2_H + 0.05;  

    TH3D hD3   ("hD3",   "Data;dy (m);dx (m);W^{2} (GeV^{2})", NBX, dyhist_low3, dyhist_high3,
                                                      NBY, dxhist_low3, dxhist_high3,
                                                      NBZ, W2hist_low3, W2hist_high3);
    TH3D hP3   ("hP3",   "QE p (sim);dy;dx;W^{2}", NBX, dyhist_low3, dyhist_high3,
                                               NBY, dxhist_low3, dxhist_high3,
                                               NBZ, W2hist_low3, W2hist_high3);
    TH3D hN3   ("hN3",   "QE n (sim);dy;dx;W^{2}", NBX, dyhist_low3, dyhist_high3,
                                               NBY, dxhist_low3, dxhist_high3,
                                               NBZ, W2hist_low3, W2hist_high3);
    TH3D hI3   ("hI3",   "Inelastic (sim);dy;dx;W^{2}", NBX, dyhist_low3, dyhist_high3,
                                                   NBY, dxhist_low3, dxhist_high3,
                                                   NBZ, W2hist_low3, W2hist_high3);

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
    
    TH1D hW2_for_asym   ("hW2_for_asym",       "W^{2};W^{2} (GeV^{2})", NBW2, -2, 6);
    TH1D hW2_for_asym_pos   ("hW2_for_asym_pos",       "W^{2};W^{2} (GeV^{2})", NBW2, -2, 6);
    TH1D hW2_for_asym_neg   ("hW2_for_asym_neg",       "W^{2};W^{2} (GeV^{2})", NBW2, -2, 6);

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


    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    TH1D h_eratio ("h_eratio","eratio distribution ; eratio ; number of secondary clusters",100,0,2);
    TH1D h_eratio_test ("h_eratio_test","h_eratio_test",100,0,2);
    TH1D h_eratio_wide_bins ("h_eratio_wide_bins","eratio distribution ; eratio ; number of clusters",10,0.0,1.0);
    TH1D h_tdiff ("h_tdiff","tdiff distribution; tdiff (ns);number of clusters",100,-100,100);
    TH1D hdist ("hdist","distance from the primary cluster to the secondaries (QE + tdiff cut); dist (m)", 100, 0 , 5);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                      LOOP 1                                                    //
    //                                                                                                //
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////// --- loop data   ///////////////////////////////////////////////////////
    //                                                                                              //
    //                                                                                              //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    Long64_t n=ch.GetEntries();

    double A = 0.287186374;
    double B = 0.798543608;
    double m = 0.938;

    std::cout << "\n"<<"[InelasticCorrection] looping over data " << n << " events\n";
    const Long64_t step     = 100;
    for(Long64_t i=0;i<n;++i){ 

        ch.GetEntry(i);
        
        int helCorr = -1*v.helicity*v.IHWP*c_.Pkin_L;

        double ee = A*(1+B*v.thtgt)/v.thetabend;
        double Q2_new = 2*v.ebeam*ee*(1-cos(v.etheta));
        double nu_new = (v.ebeam - ee);
        double W2_new = m*m+2*m*nu_new-Q2_new;
        

        if(rq_ && (!rq_->helicityOK(v.runnum)||!rq_->mollerOK(v.runnum))) continue;
   
        if( v.runnum<c_.runnum_L || v.runnum>c_.runnum_H ||
            v.ntrack<1 || abs(v.vz)>0.27 || v.eHCAL<c_.eHCAL_L || abs(((v.ePS+v.eSH)/v.trP)-1)>0.2 || v.ePS<0.2 ||
            (c_.coin_L>v.coin_time || v.coin_time>c_.coin_H) || abs(v.helicity)!=1 /*|| (v.ntrack_sbs>0) || abs(v.vz_sbs)<0.27*/) continue; //change here to remove sbs track cut

        if(((v.ePS+v.eSH)/v.trP)<0.8 || ((v.ePS+v.eSH)/v.trP)>1.2) continue; //if abs for E/P malfunctions

        // filling W2 for asymmetry calculation
        if(v.W2>-1 && v.W2<6){    
            hW2_for_asym.Fill(v.W2);
            if(helCorr==+1) hW2_for_asym_pos.Fill(v.W2);
            else if(helCorr==-1) hW2_for_asym_neg.Fill(v.W2);
        }
        if(c_.W2_L<v.W2 && v.W2<c_.W2_H && (c_.dy_L-0.1<v.dy && v.dy<c_.dy_H+0.1)){
            if(std::strcmp(kin_, "GEN3_He3") == 0){
                hDxdy.Fill(v.dy,v.dx);
            }
            else{
                hDxdy.Fill(v.dy,v.dx);
            }

        }

        if(c_.W2_L-0.05<v.W2 && v.W2<c_.W2_H+0.05 && (c_.dy_L-0.1<v.dy && v.dy<c_.dy_H+0.1) && (c_.dx_L-0.1<v.dx && v.dx<c_.dx_H+0.1)){
            if(std::strcmp(kin_, "GEN3_He3") == 0){
                hD3.Fill(v.dy,v.dx,v.W2);
            }
            else{
                hD3.Fill(v.dy,v.dx,v.W2);
            }
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

        if(std::strcmp(kin_, "GEN4_He3") == 0){
                W2_new = W2_new;
        }
        
        else{
                W2_new = v.W2;
        }

        if((c_.dy_L<v.dy && v.dy<c_.dy_H) && (c_.dx_L<v.dx && v.dx<c_.dx_H) && (W2_hist_lower_limit<W2_new && W2_new<W2_hist_upper_limit) 
        &&((pow((v.dy-c_.dy_c)/c_.dy_r,2)+pow((v.dx-c_.dx_c)/c_.dx_r,2))<=1)) {
            if(std::strcmp(kin_, "GEN4_He3") == 0){
                hData_W2_Neutrons.Fill(W2_new);
            }
            else{
                hData_W2_Neutrons.Fill(v.W2);
            }
        }

        if(((pow((v.dy-c_.dy_c)/c_.dy_r,2)+pow((v.dx-c_.dx_c)/c_.dx_r,2))<=1) || ((pow((v.dy-c_.dy_P_c)/c_.dy_P_r,2)+pow((v.dx-c_.dx_P_c)/c_.dx_P_r,2))<=1)
            &&(W2_hist_lower_limit<W2_new && W2_new<W2_hist_upper_limit)){
            
            if(std::strcmp(kin_, "GEN4_He3") == 0){
                hData_W2.Fill(W2_new);
            }
            else{
                hData_W2.Fill(v.W2);
            }

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
        
        if(v.hcal_nclus<2) continue;

        std::vector<int> order(v.hcal_nclus);
        std::iota(order.begin(), order.end(), 0);          // 0,1,2,…

        std::sort(order.begin(), order.end(),
            [&](int a, int b){ return v.hcal_clus_e[a] > v.hcal_clus_e[b]; });

        const int i0 = order[0];
        const int i1 = order[1];   

        double  e0   = v.hcal_clus_e[i0];                     
        double  t0   = v.hcal_clus_atime[i0];
        double  x0   = v.hcal_clus_x[i0];
        double  y0   = v.hcal_clus_y[i0];
        
        int nstrong = 0;
        int nstrong2 = 0;
        int nstrong4 = 0;
        int nstrong6 = 0;
        int nstrong8 = 0;

        hData.Fill(v.dx);

        for (int k = 1; k < v.hcal_nclus; ++k) {                   // start at 1 ⇒ “secondary”
            int    ik   = order[k];                         // kth-highest-E cluster
            double er   = v.hcal_clus_e[ik]   / e0;           // E_secondary / E_primary
            double dt   = t0 - v.hcal_clus_atime[ik];         // Δt = t_primary – t_sec
            double er0  = v.hcal_clus_e[ik]   / v.eHCAL;        // optional cross-check
            double dist = std::hypot(v.hcal_clus_x[ik] - x0, v.hcal_clus_y[ik] - y0);

            //double dxi = dx_all_clus[ik];
            //double dyi = dy_all_clus[ik];

            h_eratio.Fill(er);                        // all secondary ratios
            h_eratio_wide_bins.Fill(er);
            h_tdiff.Fill(dt);                        // all secondary Δt
            h_eratio_test.Fill(er0);                       // sec / total-cluster-E
            
            hdist.Fill(dist);
            //h_secondary_dx_dist->Fill(dxi,dist);
            //h_secondary_dy_dist->Fill(dyi,dist);

            bool gooddt = abs(dt+3)<2;


        }
        
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
    std::cout << "\n"<< "[InelasticCorrection] looping over QE sim for dx fit " << nentries_QE << " events\n";
    for(Long64_t i=0;i<ch_QE.GetEntries();++i){ 

        ch_QE.GetEntry(i);

        if(/*vQE.ntrack<1 ||*/ abs(vQE.vz)>0.27 || vQE.eHCAL<c_.eHCAL_L || abs((vQE.ePS+vQE.eSH)/(vQE.trP)-1)>0.2 || vQE.ePS<0.2) continue;

        double dxq = vQE.dx;//-0.02;//dx_shifted_QE(vQE);

        if(/*vQE.ntrack<1 ||*/ abs(vQE.vz)>0.27 || vQE.eHCAL<c_.eHCAL_L || abs((vQE.ePS+vQE.eSH)/(vQE.trP)-1)>0.2 || vQE.ePS<0.2 ||
            (c_.W2_L-0.05>vQE.W2 || vQE.W2>c_.W2_H+0.05) || (c_.dy_L-0.1>vQE.dy || vQE.dy>c_.dy_H+0.1))continue;

        if (c_.dx_L-0.1<dxq && dxq<c_.dx_H+0.1){
            if (vQE.fnucl==1) hP3.Fill(vQE.dy, dxq, vQE.W2, vQE.weight);
            if (vQE.fnucl==0) hN3.Fill(vQE.dy, dxq, vQE.W2, vQE.weight);
        }

        if(/*vQE.ntrack<1 ||*/ abs(vQE.vz)>0.27 || vQE.eHCAL<c_.eHCAL_L || abs((vQE.ePS+vQE.eSH)/(vQE.trP)-1)>0.2 || vQE.ePS<0.2 ||
            (c_.W2_L>vQE.W2 || vQE.W2>c_.W2_H) || (c_.dy_L-0.1>vQE.dy || vQE.dy>c_.dy_H+0.1))continue;

        // after the "continue" filters and inside the kept region for the dx fit:
        if (vQE.fnucl == 1) hQE_dxdy_p.Fill(vQE.dy, dxq, vQE.weight);
        if (vQE.fnucl == 0) hQE_dxdy_n.Fill(vQE.dy, dxq, vQE.weight);

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
    std::cout << "\n"<<"[InelasticCorrection] looping over Inelastic sim for dx fit" << nentries_inel << " events\n";
    for(Long64_t i=0;i<ch_inel.GetEntries();++i){ 
        ch_inel.GetEntry(i);
    
        if(/*vInel.ntrack<1 ||*/ abs(vInel.vz)>0.27 || vInel.eHCAL<c_.eHCAL_L || abs((vInel.ePS+vInel.eSH)/(vInel.trP)-1)>0.2 || vInel.ePS<0.2) continue;

        double dxi = vInel.dx; /*-0.02*/ //dx_shifted_Inel(vInel);

        if(/*vInel.ntrack<1 ||*/ abs(vInel.vz)>0.27 || vInel.eHCAL<c_.eHCAL_L || abs((vInel.ePS+vInel.eSH)/(vInel.trP)-1)>0.2 || vInel.ePS<0.2 ||
            (c_.W2_L-0.05>vInel.W2 || vInel.W2>c_.W2_H+0.05) || (c_.dy_L-0.1>vInel.dy || vInel.dy>c_.dy_H+0.1)) continue; 

        if(c_.dx_L-0.1<dxi && dxi<c_.dx_H+0.1){
            hI3.Fill(vInel.dy,dxi,vInel.W2,vInel.weight);
        }   

        if(/*vInel.ntrack<1 ||*/ abs(vInel.vz)>0.27 || vInel.eHCAL<c_.eHCAL_L || abs((vInel.ePS+vInel.eSH)/(vInel.trP)-1)>0.2 || vInel.ePS<0.2 ||
            (c_.W2_L>vInel.W2 || vInel.W2>c_.W2_H) || (c_.dy_L-0.1>vInel.dy || vInel.dy>c_.dy_H+0.1)) continue; 

        // after the "continue" filters and inside the kept region for the dx fit:
        hInel_dxdy.Fill(vInel.dy, dxi, vInel.weight);


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

            hInelastic.Fill(vInel.dx/*+0.25*/,vInel.weight);//+0.4 is for GEN3 its hard coded for now

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

    double par0=1,par1=1,par2=1,dx_p_out=0,dx_n_out=0,dx_inel_out=0; 

    TH1D * h_combined =  performFit(&hData,&hInelastic,&hQE_proton,&hQE_neutron,par0,par1,par2,dx_p_out,dx_n_out,dx_inel_out);

    double parP0=1,parP1=1,parP2=1,dx_p_outP=0,dx_n_outP=0,dx_inel_outP=0;

    TH1D * h_combined_pos =  performFit(&hData_pos,&hInelastic,&hQE_proton,&hQE_neutron,parP0,parP1,parP2,dx_p_outP,dx_n_outP,dx_inel_outP);

    double parN0=1,parN1=1,parN2=1,dx_p_outN=0,dx_n_outN=0,dx_inel_outN=0;

    TH1D * h_combined_neg =  performFit(&hData_neg,&hInelastic,&hQE_proton,&hQE_neutron,parN0,parN1,parN2,dx_p_outN,dx_n_outN,dx_inel_outN);

    std::cout<<"dx fitting method"<<std::endl;
    std::cout<<"par0 : "<<par0<<std::endl;
    std::cout<<"par1 : "<<par1<<std::endl;
    std::cout<<"par2 : "<<par2<<std::endl;
    std::cout<<"dx_p_out : "<<dx_p_out<<std::endl;
    std::cout<<"dx_n_out : "<<dx_n_out<<std::endl;
    std::cout<<"dx_inel_out : "<<dx_inel_out<<std::endl;

    h_combined->Scale(hData.Integral());
    h_combined_pos->Scale(hData_pos.Integral());
    h_combined_neg->Scale(hData_neg.Integral());

    TH1D* hQE_proton_dx_shifted = make_shifted_unit(&hQE_proton,   "hQE_proton_dx_shifted", dx_p_out);
    TH1D* hQE_neutron_dx_shifted = make_shifted_unit(&hQE_neutron,   "hQE_neutron_dx_shifted", dx_n_out);
    TH1D* hInelastic_dx_shifted = make_shifted_unit(&hInelastic,   "hinelastic_dx_shifted", dx_inel_out);

    TH1D* hQE_proton_pos = (TH1D*) hQE_proton.Clone("hQE_proton_pos");
    TH1D* hQE_proton_neg = (TH1D*) hQE_proton.Clone("hQE_proton_neg");

    TH1D* hQE_neutron_pos = (TH1D*) hQE_neutron.Clone("hQE_neutron_pos");
    TH1D* hQE_neutron_neg = (TH1D*) hQE_neutron.Clone("hQE_neutron_neg");

    TH1D* hInelastic_pos = (TH1D*) hInelastic.Clone("hInelastic_pos");
    TH1D* hInelastic_neg = (TH1D*) hInelastic.Clone("hInelastic_neg");

    hQE_proton_dx_shifted->Scale(par0*hData.Integral()*1/hQE_proton_dx_shifted->Integral());
    hQE_neutron_dx_shifted->Scale(par0*par1*hData.Integral()*1/hQE_neutron_dx_shifted->Integral());
    hInelastic_dx_shifted->Scale(par0*par2*hData.Integral()*1/hInelastic_dx_shifted->Integral());

    hQE_proton_pos->Scale(parP0*hData_pos.Integral()*1/hQE_proton_pos->Integral());
    hQE_neutron_pos->Scale(parP0*parP1*hData_pos.Integral()*1/hQE_neutron_pos->Integral());
    hInelastic_pos->Scale(parP0*parP2*hData_pos.Integral()*1/hInelastic_pos->Integral());

    hQE_proton_neg->Scale(parN0*hData_neg.Integral()*1/hQE_proton_neg->Integral());
    hQE_neutron_neg->Scale(parN0*parN1*hData_neg.Integral()*1/hQE_neutron_neg->Integral());
    hInelastic_neg->Scale(parN0*parN2*hData_neg.Integral()*1/hInelastic_neg->Integral());

    //frac_  = wI;                 // N_inel / (N_inel+QE)
    //dfrac_ = std::sqrt(wI*wQ/(wI+wQ)/(wI+wQ)/(hD.Integral()));

    double QE_events = hData.Integral(hData.FindBin(c_.dx_L),hData.FindBin(c_.dx_H));
    double inelastic_events = hInelastic_dx_shifted->Integral(hInelastic_dx_shifted->FindBin(c_.dx_L),hInelastic_dx_shifted->FindBin(c_.dx_H));
    double proton_events = hQE_proton_dx_shifted->Integral(hQE_proton_dx_shifted->FindBin(c_.dx_L),hQE_proton_dx_shifted->FindBin(c_.dx_H));
    double neutron_events = hQE_neutron_dx_shifted->Integral(hQE_neutron_dx_shifted->FindBin(c_.dx_L),hQE_neutron_dx_shifted->FindBin(c_.dx_H));

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
    double err_background_frac = (inelastic_events/QE_events)*sqrt((1/inelastic_events)+(1/QE_events));
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
    txt<<"data events = "<<QE_events<<"\n";
    txt<<"background events = "<<inelastic_events<<"\n";
    txt<<"background_fraction = "<<background_frac<<"\n";
    txt<<"err_background_fraction = "<<err_background_frac<<"\n";
    //txt<<"err_inelastic_fraction = "<<errinelastic_frac<<"\n";
    txt<<"proton_fraction = "<<proton_frac<<"\n";
    txt<<"err_proton_fraction = "<<errproton_frac<<"\n";
    
    txt<<"f_in_dx = "<<inelastic_frac<<"\n";
    txt<<"err_f_in_dx = "<<errinelastic_frac<<"\n";
    txt<<"f_p = "<<proton_frac<<"\n";
    txt<<"err_f_p = "<<errproton_frac<<"\n";
    //below values should be calculated separately, for now they are set to zero
    txt<<"inelastic_events_pos = "<<inelastic_events_pos<<"\n";
    txt<<"inelastic_events_neg = "<<inelastic_events_neg<<"\n";
    txt<<"A_in_dx_method = "<<A_in<<"\n";
    txt<<"err_A_in_dx_method = "<<err_A_in<<"\n";
    txt<<"A_p = "<<0.0<<"\n";
    txt<<"err_A_p = "<<0.0<<"\n";


    //txt.close();

    std::cout << "[InelasticCorrection] par0 = "<< par0 << std::endl;
    //     << " saved to "<< outFile_ <<"\n";



    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////        2D fit //////////////////////////////////////////////////////////////////////////////////////

    gStyle->SetOptStat(0);

    // ---------------- 2D (dx,dy) template fit ----------------
    double A2=1.0, rNP2=0.3, rI2=0.1;
    double dxp=0, dyp=0, dxn=0, dyn=0, dxi=0, dyi=0;

    TH2D* hComb2D = performFitDxDy(
        &hDxdy, &hInel_dxdy, &hQE_dxdy_p, &hQE_dxdy_n,
        A2, rNP2, rI2, dxp, dyp, dxn, dyn, dxi, dyi);

    const double Ndata2D  = hDxdy.Integral();
    const double Nmodel2D = hComb2D->Integral();
    if (Nmodel2D > 0) hComb2D->Scale(Ndata2D / Nmodel2D);

    //TH2D* hComb2D_scaled = (TH2D*)hDxdy.Clone("hComb2D_scaled");

    // Make residual and pull maps for QA
    TH2D* hResid2D = (TH2D*)hDxdy.Clone("hResid2D");
    hResid2D->Reset();
    hResid2D->Add(&hDxdy, hComb2D, 1.0, -1.0);

    TH2D* hPull2D  = (TH2D*)hDxdy.Clone("hPull2D");
    for (int ix=1; ix<=hPull2D->GetNbinsX(); ++ix){
        for (int iy=1; iy<=hPull2D->GetNbinsY(); ++iy){
            const double D = hDxdy.GetBinContent(ix,iy);
            const double M = hComb2D->GetBinContent(ix,iy);
            double pull = 0.0;
            if (M > 1e-12) pull = (D - M)/std::sqrt(M);
            hPull2D->SetBinContent(ix,iy, pull);
        }
    }

    // Quick visual
    TCanvas* Cdxdy = new TCanvas("Cdxdy","2D dx-dy fit", 2400, 1500);
    Cdxdy->Divide(3,1);
    Cdxdy->cd(1); gPad->SetRightMargin(0.13); hDxdy.Draw("COLZ");     hDxdy.SetTitle("Data: dy vs dx");
    Cdxdy->cd(2); gPad->SetRightMargin(0.13); hComb2D->Draw("COLZ");  hComb2D->SetTitle("Fit Model: A[P + rNP*N + rI*I]");
    Cdxdy->cd(3); gPad->SetRightMargin(0.13); hResid2D->Draw("COLZ"); hResid2D->SetTitle("Residual: Data - Model");

    // Optional: pull map
    TCanvas* CdxdyPull = new TCanvas("CdxdyPull","2D pulls", 2400, 1500);
    gPad->SetRightMargin(0.13);
    hPull2D->SetTitle("Pulls: (Data-Model)/#sqrt{Model}");
    hPull2D->Draw("COLZ");

    // Save
    Cdxdy->Print(Form("images/%s/InelasticCorrection_dxdy_%s.png", kin_, kin_));
    CdxdyPull->Print(Form("images/%s/InelasticCorrection_dxdy_pulls_%s.png", kin_, kin_));

    // Report
    std::cout << "[DxDyFit] A="<<A2<<"  rNP="<<rNP2<<"  rI="<<rI2<<"\n"
              << "         (dx,dy)_p=("<<dxp<<","<<dyp<<")  (dx,dy)_n=("<<dxn<<","<<dyn<<")  (dx,dy)_inel=("<<dxi<<","<<dyi<<")\n";

    // =================== Projections & Overlays ===================

    // Make unit PDFs with the fitted shifts
    auto make_shifted_unit2D = [&](TH2D* src, const char* name, double dx, double dy){
        TH2D* h = (TH2D*)src->Clone(name);
        h->Reset();
        for (int ix=1; ix<=h->GetNbinsX(); ++ix){
            const double x = h->GetXaxis()->GetBinCenter(ix);
            for (int iy=1; iy<=h->GetNbinsY(); ++iy){
                const double y = h->GetYaxis()->GetBinCenter(iy);
                h->SetBinContent(ix, iy, src->Interpolate(x - dx, y - dy));
            }
        }
        double S = h->Integral();
        if (S > 0) h->Scale(1.0/S);
        return h;
    };

    TH2D* hP2 = make_shifted_unit2D(&hQE_dxdy_p, "hP2",  dxp, dyp);
    TH2D* hN2 = make_shifted_unit2D(&hQE_dxdy_n, "hN2",  dxn, dyn);
    TH2D* hI2 = make_shifted_unit2D(&hInel_dxdy, "hI2",  dxi, dyi);

    // Now scale so that P + rNP*N + rI*I has the same total as data
    const double Ndata2 = hDxdy.Integral();
    const double denom2  = 1.0 + rNP2 + rI2;

    // Each component becomes the **absolute** contribution
    hP2->Scale(Ndata2 * (1.0/denom2));      // QE p
    hN2->Scale(Ndata2 * (rNP2/denom2));     // QE n
    hI2->Scale(Ndata2 * (rI2/denom2));      // Inelastic

    // Build a perfectly-normalized combined model from components
    TH2D* hComb2D_scaled = (TH2D*)hP2->Clone("hComb2D_scaled");
    hComb2D_scaled->Add(hN2);
    hComb2D_scaled->Add(hI2);


    // Data projections
    TH1D* data_dx = hDxdy.ProjectionY("data_dx");   // Y = dx
    TH1D* data_dy = hDxdy.ProjectionX("data_dy");   // X = dy

    // Model projections (combined)
    TH1D* model_dx = hComb2D_scaled->ProjectionY("model_dx");
    TH1D* model_dy = hComb2D_scaled->ProjectionX("model_dy");

    // Component projections (already absolute counts if you used Option B)
    TH1D* p_dx = hP2->ProjectionY("p_dx");
    TH1D* n_dx = hN2->ProjectionY("n_dx");
    TH1D* i_dx = hI2->ProjectionY("i_dx");

    TH1D* p_dy = hP2->ProjectionX("p_dy");
    TH1D* n_dy = hN2->ProjectionX("n_dy");
    TH1D* i_dy = hI2->ProjectionX("i_dy");

    // --- Styling: dx ---
    data_dx->SetMarkerStyle(kFullCircle);
    data_dx->SetLineColor(kBlack);

    model_dx->SetLineColor(kGreen+2);
    model_dx->SetLineWidth(3);

    p_dx->SetLineColor(6);        // magenta-ish (same as your 1D colors)
    p_dx->SetLineWidth(3);
    p_dx->SetFillStyle(3004);
    p_dx->SetFillColorAlpha(6,0.35);

    n_dx->SetLineColor(9);        // blue-ish
    n_dx->SetLineWidth(3);
    n_dx->SetFillStyle(3005);
    n_dx->SetFillColorAlpha(9,0.35);

    i_dx->SetLineColor(7);        // olive-ish
    i_dx->SetLineWidth(3);
    i_dx->SetFillStyle(3009);
    i_dx->SetFillColorAlpha(7,0.35);

    // --- Styling: dy ---
    data_dy->SetMarkerStyle(kFullSquare);
    data_dy->SetLineColor(kBlack);

    model_dy->SetLineColor(kGreen+2);
    model_dy->SetLineWidth(3);

    p_dy->SetLineColor(6);
    p_dy->SetLineWidth(3);
    p_dy->SetFillStyle(3004);
    p_dy->SetFillColorAlpha(6,0.35);

    n_dy->SetLineColor(9);
    n_dy->SetLineWidth(3);
    n_dy->SetFillStyle(3005);
    n_dy->SetFillColorAlpha(9,0.35);

    i_dy->SetLineColor(7);
    i_dy->SetLineWidth(3);
    i_dy->SetFillStyle(3009);
    i_dy->SetFillColorAlpha(7,0.35);

    // --- Draw overlays: dx (Y projection) ---
    TCanvas* Cproj_dx = new TCanvas("Cproj_dx","dx projection: data vs model",2400,1500);
    Cproj_dx->cd();
    data_dx->SetTitle("dx projection;dx (m);Counts");
    data_dx->Draw("E1");                // data with errors
    p_dx->Draw("HIST SAME");
    n_dx->Draw("HIST SAME");
    i_dx->Draw("HIST SAME");
    model_dx->Draw("HIST SAME");

    auto leg_dx = new TLegend(0.62,0.62,0.88,0.88);
    leg_dx->SetBorderSize(0); leg_dx->SetFillStyle(0);
    leg_dx->AddEntry(data_dx,   "Data", "lep");
    leg_dx->AddEntry(model_dx,  "Model (P+N+Inel)", "l");
    leg_dx->AddEntry(p_dx,      "QE p (model comp.)", "f");
    leg_dx->AddEntry(n_dx,      "QE n (model comp.)", "f");
    leg_dx->AddEntry(i_dx,      "Inelastic (model comp.)", "f");
    leg_dx->Draw();

    // --- Draw overlays: dy (X projection) ---
    TCanvas* Cproj_dy = new TCanvas("Cproj_dy","dy projection: data vs model",2400,1500);
    Cproj_dy->cd();
    data_dy->SetTitle("dy projection;dy (m);Counts");
    data_dy->Draw("E1");
    p_dy->Draw("HIST SAME");
    n_dy->Draw("HIST SAME");
    i_dy->Draw("HIST SAME");
    model_dy->Draw("HIST SAME");

    auto leg_dy = new TLegend(0.62,0.62,0.88,0.88);
    leg_dy->SetBorderSize(0); leg_dy->SetFillStyle(0);
    leg_dy->AddEntry(data_dy,   "Data", "lep");
    leg_dy->AddEntry(model_dy,  "Model (P+N+Inel)", "l");
    leg_dy->AddEntry(p_dy,      "QE p (model comp.)", "f");
    leg_dy->AddEntry(n_dy,      "QE n (model comp.)", "f");
    leg_dy->AddEntry(i_dy,      "Inelastic (model comp.)", "f");
    leg_dy->Draw();

    // (Optional) save
    Cproj_dx->Print(Form("images/%s/Projection_dx_%s.png", kin_, kin_));
    Cproj_dy->Print(Form("images/%s/Projection_dy_%s.png", kin_, kin_));


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////fraction calculation /////////////////////////////////////////////////////////////////////////////////////

        // ---------- Setup ROI + helpers ----------
    auto ensure_sumw2 = [](TH2D* h){ if(!h->GetSumw2N()) h->Sumw2(); };

    ensure_sumw2(&hDxdy);
    ensure_sumw2(hP2);
    ensure_sumw2(hN2);
    ensure_sumw2(hI2);

    // Elliptical neutron ROI centered at (0,0):
    struct ROI { double a_dx=0.30; double b_dy=0.40; }; // semi-axes in dx and dy
    ROI roi; // adjust if needed

    // Integrate inside (dx/a)^2 + (dy/b)^2 <= 1.
    // NOTE: X axis = dy, Y axis = dx in your hists.
    auto integrate_ellipse = [&](TH2D* h, double a_dx, double b_dy){
        double sum = 0.0, sumw2 = 0.0;
        auto* axX = h->GetXaxis(); // dy
        auto* axY = h->GetYaxis(); // dx
        const int nx = h->GetNbinsX(), ny = h->GetNbinsY();

        for(int ix=1; ix<=nx; ++ix){
            const double dy = axX->GetBinCenter(ix);
            for(int iy=1; iy<=ny; ++iy){
                const double dx = axY->GetBinCenter(iy);
                const double val = h->GetBinContent(ix,iy);
                const double w2  = h->GetBinError(ix,iy); // this is sigma; square it
                if ( (dx*dx)/(a_dx*a_dx) + (dy*dy)/(b_dy*b_dy) <= 1.0 ){
                    sum   += val;
                    sumw2 += (w2*w2);
                }
            }
        }
        double err = (sumw2>0.0) ? std::sqrt(sumw2) : std::sqrt(std::max(0.0,sum)); // safe fallback
        return std::pair<double,double>(sum, err);
    };

    // ---------- Do the integrals ----------
    auto [Ndata,  eData]  = integrate_ellipse(&hDxdy, roi.a_dx, roi.b_dy);
    auto [Nprot,  eProt]  = integrate_ellipse(hP2,    roi.a_dx, roi.b_dy);
    auto [Nneut,  eNeut]  = integrate_ellipse(hN2,    roi.a_dx, roi.b_dy);
    auto [Ninel,  eInel]  = integrate_ellipse(hI2,    roi.a_dx, roi.b_dy);

    // ---------- Ratios x / data with errors ----------
    auto ratio_err = [](double x, double ex, double d, double ed){
        if (d<=0 || x<0) return std::pair<double,double>(0.0, 0.0);
        double r = x/d;
        double dr = r * std::sqrt( (ex>0? (ex*ex)/(x*x) : 0.0) + (ed>0? (ed*ed)/(d*d) : 0.0) );
        return std::pair<double,double>(r, dr);
    };

    auto [Rprot, eRprot] = ratio_err(Nprot, eProt, Ndata, eData);
    auto [Rneut, eRneut] = ratio_err(Nneut, eNeut, Ndata, eData);
    auto [Rinel, eRinel] = ratio_err(Ninel, eInel, Ndata, eData);

    // ---------- Print + (optional) save ----------
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\n[Neutron-spot ROI] ellipse: (dx/"<<roi.a_dx<<")^2 + (dy/"<<roi.b_dy<<")^2 <= 1\n";
    std::cout << "Data     : " << Ndata << " ± " << eData << "\n";
    std::cout << "Proton   : " << Nprot << " ± " << eProt << "   => Proton/Data = " << Rprot << " ± " << eRprot << "\n";
    std::cout << "Neutron  : " << Nneut << " ± " << eNeut << "   => Neutron/Data = " << Rneut << " ± " << eRneut << "\n";
    std::cout << "Inelastic: " << Ninel << " ± " << eInel << "   => Inel/Data    = " << Rinel << " ± " << eRinel << "\n";

    // Write to your correction file too (append)
    {
        std::ofstream out(Form("corrections/%s/InelasticCorrection_2D_fit_%s.txt", kin_, kin_), std::ios::app);
        out << std::setprecision(6);
        out << "\n# --- Neutron-spot ROI results (dx,dy centered at 0,0) ---\n";
        out << "roi_a_dx = " << roi.a_dx << "\n";
        out << "roi_b_dy = " << roi.b_dy << "\n";
        out << "ROI_Data_counts = " << Ndata  << "  err = " << eData  << "\n";
        out << "ROI_Proton      = " << Nprot  << "  err = " << eProt  << "  ratio = " << Rprot << "  err_ratio = " << eRprot << "\n";
        out << "ROI_Neutron     = " << Nneut  << "  err = " << eNeut  << "  ratio = " << Rneut << "  err_ratio = " << eRneut << "\n";
        out << "ROI_Inelastic   = " << Ninel  << "  err = " << eInel  << "  ratio = " << Rinel << "  err_ratio = " << eRinel << "\n";
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////  3D  Fit ///////////////////////////////////////////////////////////////////////////////


    // ---------------- 3D fit ----------------
    double A3=1.0, rNP3=0.3, rI3=0.1;
    double dWp=0, dWn=0, dWi=0;
    dxp=0, dyp=0 ,dxn=0, dyn=0, dxi=0, dyi=0;


    TH3D* hComb3D = performFitDxDyW2(&hD3, &hI3, &hP3, &hN3,
        A3, rNP3, rI3,
        dxp,dyp,dWp, dxn,dyn,dWn, dxi,dyi,dWi);

    // --- Option B-style component builds in absolute counts (recommended) ---
    auto make_shifted_unit3D = [&](TH3D* src, const char* nm, double ddx, double ddy, double dW){
        TH3D* h=(TH3D*)src->Clone(nm); h->Reset();
        for(int ix=1; ix<=h->GetNbinsX(); ++ix){
          const double x=h->GetXaxis()->GetBinCenter(ix);
          for(int iy=1; iy<=h->GetNbinsY(); ++iy){
            const double y=h->GetYaxis()->GetBinCenter(iy);
            for(int iz=1; iz<=h->GetNbinsZ(); ++iz){
              const double z=h->GetZaxis()->GetBinCenter(iz);
              h->SetBinContent(ix,iy,iz, src->Interpolate(x - ddx, y - ddy, z - dW));
            }
          }
        }
        double S=h->Integral(); if(S>0) h->Scale(1.0/S);
        return h;
    };

    TH3D* hP3u = make_shifted_unit3D(&hP3, "hP3u", dxp,dyp,dWp);
    TH3D* hN3u = make_shifted_unit3D(&hN3, "hN3u", dxn,dyn,dWn);
    TH3D* hI3u = make_shifted_unit3D(&hI3, "hI3u", dxi,dyi,dWi);

    // scale components -> absolute counts that sum to data
    const double Ndata3 = hD3.Integral();
    const double denom3 = 1.0 + rNP3 + rI3;
    hP3u->Scale(Ndata3 * (1.0/denom3));
    hN3u->Scale(Ndata3 * (rNP3/denom3));
    hI3u->Scale(Ndata3 * (rI3 /denom3));

    // build perfectly normalized combined model
    TH3D* hComb3D_scaled=(TH3D*)hP3u->Clone("hComb3D_scaled");
    hComb3D_scaled->Add(hN3u);
    hComb3D_scaled->Add(hI3u);

    // residual cube (optional)
    TH3D* hResid3D=(TH3D*)hD3.Clone("hResid3D");
    hResid3D->Add(hComb3D_scaled,-1.0);

    // // dy projection (X axis)
    // TH1D* data_dy_3 = hD3.ProjectionX("data_dy_3");
    // TH1D* modl_dy_3 = hComb3D_scaled->ProjectionX("modl_dy_3");
    // TH1D* P_dy_3 = hP3u->ProjectionX("P_dy_3");
    // TH1D* N_dy_3 = hN3u->ProjectionX("N_dy_3");
    // TH1D* I_dy_3 = hI3u->ProjectionX("I_dy_3");

    // // dx projection (Y axis)
    // TH1D* data_dx_3 = hD3.ProjectionY("data_dx_3");
    // TH1D* modl_dx_3 = hComb3D_scaled->ProjectionY("modl_dx_3");
    // TH1D* P_dx_3 = hP3u->ProjectionY("P_dx_3");
    // TH1D* N_dx_3 = hN3u->ProjectionY("N_dx_3");
    // TH1D* I_dx_3 = hI3u->ProjectionY("I_dx_3");

    // // W2 projection (Z axis)
    // TH1D* data_w2_3 = hD3.ProjectionZ("data_w2_3");
    // TH1D* modl_w2_3 = hComb3D_scaled->ProjectionZ("modl_w2_3");
    // TH1D* P_w2_3 = hP3u->ProjectionZ("P_w2_3");
    // TH1D* N_w2_3 = hN3u->ProjectionZ("N_w2_3");
    // TH1D* I_w2_3 = hI3u->ProjectionZ("I_w2_3");

    auto make_pulls_1D = [&](TH1D* hData, TH1D* hModel, const char* name){
        TH1D* p = (TH1D*)hData->Clone(name);
        p->Reset("ICES"); // keep axis & sumw2
        const int nb = hData->GetNbinsX();
        for (int i=1;i<=nb;++i){
            const double D = hData->GetBinContent(i);
            const double M = hModel->GetBinContent(i);
            const double s = std::sqrt(std::max(M,1.0));
            p->SetBinContent(i, (s>0)? (D - M)/s : 0.0);
            p->SetBinError(i, 0.0);
        }
        return p;
    };

    auto style_fill = [&](TH1* h, Color_t lc, Color_t fc, Int_t lw, Int_t fs, double fa){
        h->SetLineColor(lc); h->SetLineWidth(lw);
        h->SetFillColorAlpha(fc, fa); h->SetFillStyle(fs);
    };

    // Project a TH3D into a TH2D (X vs Y) over a Z-bin window [zb1, zb2].
    // Works on older ROOT via TH3::Project3D("xy").
    auto projXY_in_Z = [&](TH3D* h3, int zb1, int zb2, const char* name)->TH2D* {
        // Save current Z range
        TAxis* zax = h3->GetZaxis();
        const int oldFirst = zax->GetFirst();
        const int oldLast  = zax->GetLast();

        // Limit to requested Z slice and project
        zax->SetRange(zb1, zb2);
        // "xy" means output has X = original X, Y = original Y.
        // Your convention is X=dy, Y=dx, so "xy" is exactly dy vs dx.
        TH2D* h2 = (TH2D*)h3->Project3D("xy");
        h2->SetName(name);

        // Restore previous Z range
        zax->SetRange(oldFirst, oldLast);
        return h2;
    };


    auto style_data = [&](TH1* h, Int_t ms){ h->SetMarkerStyle(ms); h->SetLineColor(kBlack); };

    auto draw_proj_canvas = [&](TH1D* data, TH1D* model, TH1D* hp, TH1D* hn, TH1D* hi,
                                const char* title, const char* xlab,
                                const char* cname, const char* fname_png){

        // Styles
        style_data(data, kFullCircle);
        model->SetLineColor(kGreen+2); model->SetLineWidth(3);

        style_fill(hp, kMagenta+2, kMagenta+1, 3, 3004, 0.35);
        style_fill(hn, kAzure+2,   kAzure+1,   3, 3005, 0.35);
        style_fill(hi, kCyan+2,   kCyan+1,   3, 3009, 0.35);

        // Layout: stacked pads
        TCanvas* C = new TCanvas(cname, title, 2400, 1500);
        double L=0.12, R=0.03, Tt=0.06, Bt=0.02, Tb=0.02, Bb=0.28;
        TPad* top = new TPad(Form("%s_top",cname),"",0, Bb, 1, 1);
        TPad* bot = new TPad(Form("%s_bot",cname),"",0, 0,  1, Bb);
        for(TPad* p : {top,bot}){ p->SetLeftMargin(L); p->SetRightMargin(R); p->SetTickx(); p->SetTicky(); }
        top->SetTopMargin(Tt); top->SetBottomMargin(Bt);
        bot->SetTopMargin(Tb); bot->SetBottomMargin(0.22);
        top->Draw(); bot->Draw();

        // ----- TOP: overlay -----
        top->cd();
        data->SetTitle(Form("%s;%s;Counts", title, xlab));
        data->Draw("E1");
        hp->Draw("HIST SAME");
        hn->Draw("HIST SAME");
        hi->Draw("HIST SAME");
        model->Draw("HIST SAME");

        auto leg = new TLegend(0.62,0.62,0.88,0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->AddEntry(data, "Data", "lep");
        leg->AddEntry(model,"Model (P+N+Inel)", "l");
        leg->AddEntry(hp,   "QE p (model comp.)", "f");
        leg->AddEntry(hn,   "QE n (model comp.)", "f");
        leg->AddEntry(hi,   "Inelastic (model comp.)", "f");
        leg->Draw();

        // ----- BOTTOM: pulls -----
        bot->cd();
        TH1D* pulls = make_pulls_1D(data, model, Form("%s_pulls", cname));
        pulls->SetTitle(Form("; %s; (Data-Model)/#sqrt{Model}", xlab));
        pulls->SetLineColor(kBlack);
        pulls->GetYaxis()->SetNdivisions(505);
        pulls->GetYaxis()->SetTitleSize(0.08);
        pulls->GetYaxis()->SetLabelSize(0.07);
        pulls->GetXaxis()->SetTitleSize(0.08);
        pulls->GetXaxis()->SetLabelSize(0.07);
        pulls->GetYaxis()->SetRangeUser(-5,5);
        pulls->Draw("HIST");

        TLine* z0 = new TLine(pulls->GetXaxis()->GetXmin(),0, pulls->GetXaxis()->GetXmax(),0);
        z0->SetLineStyle(2); z0->Draw("same");

        C->Print(fname_png);
    };

    // ===== Build projections =====
    TH1D* data_dy_3 = hD3.ProjectionX("data_dy_3");
    TH1D* modl_dy_3 = hComb3D_scaled->ProjectionX("modl_dy_3");
    TH1D* P_dy_3    = hP3u->ProjectionX("P_dy_3");
    TH1D* N_dy_3    = hN3u->ProjectionX("N_dy_3");
    TH1D* I_dy_3    = hI3u->ProjectionX("I_dy_3");

    TH1D* data_dx_3 = hD3.ProjectionY("data_dx_3");
    TH1D* modl_dx_3 = hComb3D_scaled->ProjectionY("modl_dx_3");
    TH1D* P_dx_3    = hP3u->ProjectionY("P_dx_3");
    TH1D* N_dx_3    = hN3u->ProjectionY("N_dx_3");
    TH1D* I_dx_3    = hI3u->ProjectionY("I_dx_3");

    TH1D* data_w2_3 = hD3.ProjectionZ("data_w2_3");
    TH1D* modl_w2_3 = hComb3D_scaled->ProjectionZ("modl_w2_3");
    TH1D* P_w2_3    = hP3u->ProjectionZ("P_w2_3");
    TH1D* N_w2_3    = hN3u->ProjectionZ("N_w2_3");
    TH1D* I_w2_3    = hI3u->ProjectionZ("I_w2_3");

    // ===== Draw canvases =====
    draw_proj_canvas(data_dy_3, modl_dy_3, P_dy_3, N_dy_3, I_dy_3,
                     "3D model: dy projection", "dy (m)",
                     "C3D_dy", Form("images/%s/3D_dy_projection_%s.png", kin_, kin_));

    draw_proj_canvas(data_dx_3, modl_dx_3, P_dx_3, N_dx_3, I_dx_3,
                     "3D model: dx projection", "dx (m)",
                     "C3D_dx", Form("images/%s/3D_dx_projection_%s.png", kin_, kin_));

    draw_proj_canvas(data_w2_3, modl_w2_3, P_w2_3, N_w2_3, I_w2_3,
                     "3D model: W^{2} projection", "W^{2} (GeV^{2})",
                     "C3D_W2", Form("images/%s/3D_W2_projection_%s.png", kin_, kin_));


    auto draw_dxdy_slice = [&](double W2_lo, double W2_hi,
                               const char* tag_out){

        // locate Z bins
        int zb1 = hD3.GetZaxis()->FindBin(W2_lo);
        int zb2 = hD3.GetZaxis()->FindBin(W2_hi);

        // Slice: ProjectionXY over the W² window
        TH2D* D_xy  = projXY_in_Z(&hD3,                 zb1, zb2, Form("Dxy_%s",tag_out));
        TH2D* M_xy  = projXY_in_Z(hComb3D_scaled,       zb1, zb2, Form("Mxy_%s",tag_out));
        TH2D* R_xy  = (TH2D*)D_xy->Clone(Form("Rxy_%s",tag_out));
        R_xy->Add(M_xy, -1.0);

        // Style
        D_xy->SetTitle(Form("Data: dx vs dy  (%.2f < W^{2} < %.2f);dy (m);dx (m)", W2_lo, W2_hi));
        M_xy->SetTitle(Form("Model: dx vs dy  (%.2f < W^{2} < %.2f);dy (m);dx (m)", W2_lo, W2_hi));
        R_xy->SetTitle(Form("Residual: Data - Model  (%.2f < W^{2} < %.2f);dy (m);dx (m)", W2_lo, W2_hi));

        // Canvas
        TCanvas* C = new TCanvas(Form("Cslice_%s",tag_out),
                                 Form("dx–dy slices @ %.2f<W^{2}<%.2f",W2_lo,W2_hi), 2400, 1500);
        C->Divide(3,1);
        C->cd(1); D_xy->Draw("COLZ");
        C->cd(2); M_xy->Draw("COLZ");
        C->cd(3); R_xy->Draw("COLZ");

        C->Print(Form("images/%s/3D_slice_dxdy_%s_%s.png", kin_, tag_out, kin_));
    };

    // Examples: elastic-ish window and a higher-W² band (tune to your analysis)
    draw_dxdy_slice(/*W2_lo=*/c_.W2_L, /*W2_hi=*/c_.W2_H, "signalW2");
    draw_dxdy_slice(/*W2_lo=*/2.5,     /*W2_hi=*/5.0,     "inelW2");

    auto draw_component_maps = [&](double W2_lo, double W2_hi, const char* tag_out){
        int zb1 = hD3.GetZaxis()->FindBin(W2_lo);
        int zb2 = hD3.GetZaxis()->FindBin(W2_hi);

        TH2D* P_xy = projXY_in_Z(hP3u, zb1, zb2, Form("Pxy_%s",tag_out));
        TH2D* N_xy = projXY_in_Z(hN3u, zb1, zb2, Form("Nxy_%s",tag_out));
        TH2D* I_xy = projXY_in_Z(hI3u, zb1, zb2, Form("Ixy_%s",tag_out));

        P_xy->SetTitle(Form("QE p (model comp.)  %.2f<W^{2}<%.2f;dy (m);dx (m)", W2_lo, W2_hi));
        N_xy->SetTitle(Form("QE n (model comp.)  %.2f<W^{2}<%.2f;dy (m);dx (m)", W2_lo, W2_hi));
        I_xy->SetTitle(Form("Inelastic (model comp.)  %.2f<W^{2}<%.2f;dy (m);dx (m)", W2_lo, W2_hi));

        TCanvas* C = new TCanvas(Form("Ccomp_%s",tag_out),
                                 Form("Components in dx–dy (%.2f<W^{2}<%.2f)",W2_lo,W2_hi),
                                 2400, 1500);
        C->Divide(3,1);
        C->cd(1); P_xy->Draw("COLZ");
        C->cd(2); N_xy->Draw("COLZ");
        C->cd(3); I_xy->Draw("COLZ");
        C->Print(Form("images/%s/3D_slice_components_%s_%s.png", kin_, tag_out, kin_));
    };

    draw_component_maps(c_.W2_L, c_.W2_H, "signalW2");

    // Over all W²
    // Use full Z range: pass the entire Z span
    int zAll1 = hD3.GetZaxis()->GetFirst();
    int zAll2 = hD3.GetZaxis()->GetLast();

    TH2D* D_xy_all = projXY_in_Z(&hD3,           zAll1, zAll2, "D_xy_all");
    TH2D* M_xy_all = projXY_in_Z(hComb3D_scaled, zAll1, zAll2, "M_xy_all");

    TH2D* R_xy_all = (TH2D*)D_xy_all->Clone("R_xy_all"); R_xy_all->Add(M_xy_all, -1.0);

    TCanvas* Cres = new TCanvas("Cres_3D","Global residuals (dx–dy)", 2400, 1500);
    Cres->Divide(3,1);
    Cres->cd(1); D_xy_all->SetTitle("Data;dy (m);dx (m)"); D_xy_all->Draw("COLZ");
    Cres->cd(2); M_xy_all->SetTitle("Model;dy (m);dx (m)"); M_xy_all->Draw("COLZ");
    Cres->cd(3); R_xy_all->SetTitle("Residual: Data - Model;dy (m);dx (m)"); R_xy_all->Draw("COLZ");
    Cres->Print(Form("images/%s/3D_dxdy_global_residuals_%s.png", kin_, kin_));



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                     LOOP 2                                                                     //
    //                                                                                                                                                //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////// --- loop QE sim   /////////////////////////////////////////////////////
    //                                                                                              //
    //                                                                                              //
    //////////////////////////////////////////////////////////////////////////////////////////////////


    std::cout << "\n"<< "[InelasticCorrection] looping over QE sim for W2 fit " << nentries_QE << " events\n";
    for(Long64_t i=0;i<ch_QE.GetEntries();++i){ 

        ch_QE.GetEntry(i);

        if(/*vQE.ntrack<1 ||*/ abs(vQE.vz)>0.27 || vQE.eHCAL<c_.eHCAL_L || abs((vQE.ePS+vQE.eSH)/(vQE.trP)-1)>0.2 || vQE.ePS<0.2) continue;

        double dxq = vQE.dx;//-0.02;//dx_shifted_QE(vQE);
        double dyq = vQE.dy;

        if(vQE.fnucl==1){
            dxq = vQE.dx - dxp;//dx_p_out;
            dyq = vQE.dy - dyp;
        }
        else if(vQE.fnucl==0){
            dxq = vQE.dx - dxn;//dx_n_out;
            dyq = vQE.dy - dyn;
        }

        if ((c_.dy_L<dyq && dyq<c_.dy_H) && (c_.dx_L<dxq && dxq<c_.dx_H) && (W2_hist_lower_limit<vQE.W2 && vQE.W2<W2_hist_upper_limit)
         && ((pow((dyq-c_.dy_c)/c_.dy_r,2)+pow((dxq-c_.dx_c)/c_.dx_r,2))<=1)) {
            hQE_W2_Neutrons.Fill(vQE.W2, vQE.weight);
            if (vQE.fnucl == 0) hQE_neutron_W2_Neutrons.Fill(vQE.W2, vQE.weight);  // see #3
            if (vQE.fnucl == 1) hQE_proton_W2_Neutrons.Fill(vQE.W2, vQE.weight);
        }

        if (((pow((dyq-c_.dy_c)/c_.dy_r,2)+pow((dxq-c_.dx_c)/c_.dx_r,2))<=1) || ((pow((dyq-c_.dy_P_c)/c_.dy_P_r,2)+pow((dxq-c_.dx_P_r)/c_.dx_P_r,2))<=1)
            && (W2_hist_lower_limit<vQE.W2 && vQE.W2<W2_hist_upper_limit)) {
            hQE_W2.Fill(vQE.W2, vQE.weight);
            if (vQE.fnucl == 0) hQE_neutron_W2.Fill(vQE.W2, vQE.weight);  // see #3
            if (vQE.fnucl == 1) hQE_proton_W2.Fill(vQE.W2, vQE.weight);
        }

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


    std::cout << "\n"<<"[InelasticCorrection] looping over Inelastic sim for W2 fit" << nentries_inel << " events\n";
    for(Long64_t i=0;i<ch_inel.GetEntries();++i){ 
        ch_inel.GetEntry(i);
    
        if(/*vInel.ntrack<1 ||*/ abs(vInel.vz)>0.27 || vInel.eHCAL<c_.eHCAL_L || abs((vInel.ePS+vInel.eSH)/(vInel.trP)-1)>0.2 || vInel.ePS<0.2) continue;

        double dxii = vInel.dx - dxi;//dx_inel_out;
        double dyii = vInel.dy - dyi;

        if ( (c_.dy_L<dyii && dyii<c_.dy_H) && (c_.dx_L<dxii && dxii<c_.dx_H) && (W2_hist_lower_limit<vInel.W2 && vInel.W2<W2_hist_upper_limit) 
        && ((pow((dyii-c_.dy_c)/c_.dy_r,2)+pow((dxii-c_.dx_c)/c_.dx_r,2))<=1)) {
        
            hInelastic_W2_Neutrons.Fill(vInel.W2, vInel.weight);
            hInelastic_W2_2_Neutrons.Fill(vInel.W2, vInel.weight);                    // see #3

            if (vInel.eHCAL>0.250) {
                hInelastic_W2_Neutrons_eHCALcut_1.Fill(vInel.W2,vInel.weight);
            }

            if (vInel.eHCAL>0.275) {
                hInelastic_W2_Neutrons_eHCALcut_2.Fill(vInel.W2,vInel.weight);
            }

            if (vInel.eHCAL>0.300) {
                hInelastic_W2_Neutrons_eHCALcut_3.Fill(vInel.W2,vInel.weight);
            }
            if (vInel.eHCAL>0.325) {
                hInelastic_W2_Neutrons_eHCALcut_4.Fill(vInel.W2,vInel.weight);
            }

        }

        if (((pow((dyii-c_.dy_c)/c_.dy_r,2)+pow((dxii-c_.dx_c)/c_.dx_r,2))<=1) || ((pow((vInel.dy-c_.dy_P_c)/c_.dy_P_r,2)+pow((dxii-c_.dx_P_c)/c_.dx_P_r,2))<=1) 
            && (W2_hist_lower_limit<vInel.W2 && vInel.W2<W2_hist_upper_limit) ) {
            hInelastic_W2.Fill(vInel.W2, vInel.weight);
            hInelastic_W2_2.Fill(vInel.W2, vInel.weight);                    // see #3
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


    ////////////////////W2 fitting /////////////////////////////////////

    std::cout<<"Rnop 1D : "<< Rnop <<std::endl;

    double parW0=1,parW1=1; 

    //TH1D * h_combined_W2 =  performFitW2(&hData_W2,&hInelastic_W2,&hQE_proton_W2,&hQE_neutron_W2,parW0,parW1,Rnop);

    //h_combined_W2->Scale(hData_W2.Integral());

    //hQE_proton_W2.Scale(parW0*hData_W2.Integral()*1/hQE_proton_W2.Integral());
    //hQE_neutron_W2.Scale(parW0*Rnop*hData_W2.Integral()*1/hQE_neutron_W2.Integral());
    //hInelastic_W2.Scale(parW1*hData_W2.Integral()*1/hInelastic_W2.Integral());

    //std::cout<<"parW0 : "<<parW0<<std::endl;
    //std::cout<<"parW1 : "<<parW1<<std::endl;
    //std::cout<<"parW2 : "<<parW2<<std::endl;


    double alpha = 0.0;
    double delta_1 = 0.1;
    double delta_2 = 0.1;

    Rnop = Rneut/Rprot;

    std::cout<<"Rnop 2D : "<< Rnop <<std::endl;

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


    std::cout<<"single paramater method (protons+neutrons)"<<std::endl;
    std::cout<<"alpha : "<<alpha<<std::endl;
    std::cout<<"denom : "<<denom<<std::endl;
    std::cout<<"N : "<<N<<std::endl;


    std::cout<<"double paramater method (protons+neutrons)"<<std::endl;
    std::cout<<"par0_2 : "<<par0_2<<std::endl;
    std::cout<<"par1_2 : "<<par1_2<<std::endl;


    std::cout<<"delta shifts for W^{2} (protons+neutrons) (not used anymore)"<<std::endl;
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

    double ssr_neutrons_1 =
      SumSquaredResidualsInRange(&hData_W2_Neutrons,
                                  h_combined_W2_Neutrons,
                                  c_.W2_L, c_.W2_H);

    double ssr_neutrons_2 =
      SumSquaredResidualsInRange(&hData_W2_Neutrons,
                                  h_combined_W2_2_Neutrons,
                                  c_.W2_L, c_.W2_H);

    std::cout << std::fixed << std::setprecision(6)
              << "[SSR] neutrons (single-param)  W2["<<c_.W2_L<<","<<c_.W2_H<<"] = " << ssr_neutrons_1 << "\n"
              << "[SSR] neutrons (two-param)     W2["<<c_.W2_L<<","<<c_.W2_H<<"] = " << ssr_neutrons_2 << "\n";



    // branch 2: par0_2*QE + par1_2*Inel_shift
    h_combined_W2_2_Neutrons->Scale(N_Neutrons);
    hQE_W2_Neutrons.Scale(          N_Neutrons * (1.0/hQE_W2_Neutrons.Integral())           * par0_2_Neutrons );
    hInelastic_W2_shift2_Neutrons->Scale( N_Neutrons /*(1.0/hInelastic_W2_shift2_Neutrons->Integral())*/* par1_2_Neutrons );

    double W2_data_events = hData_W2_Neutrons.Integral(hData_W2_Neutrons.FindBin(c_.W2_L),hData_W2_Neutrons.FindBin(c_.W2_H));
    double W2_fit_events = h_combined_W2_Neutrons->Integral(h_combined_W2_Neutrons->FindBin(c_.W2_L),h_combined_W2_Neutrons->FindBin(c_.W2_H));
    double W2_neutron_events = hQE_neutron_W2_Neutrons.Integral(hQE_neutron_W2_Neutrons.FindBin(c_.W2_L),hQE_neutron_W2_Neutrons.FindBin(c_.W2_H));
    double W2_bkg_events = hInelastic_W2_shift1_Neutrons->Integral(hInelastic_W2_shift1_Neutrons->FindBin(c_.W2_L),hInelastic_W2_shift1_Neutrons->FindBin(c_.W2_H));
    double W2_proton_events = hQE_proton_W2_Neutrons.Integral(hQE_proton_W2_Neutrons.FindBin(c_.W2_L),hQE_proton_W2_Neutrons.FindBin(c_.W2_H));

    double W2_fit_bkg_fraction = W2_bkg_events/W2_data_events;
    double W2_fit_err_background_frac = (W2_bkg_events/W2_data_events)*sqrt((1/W2_bkg_events)+(1/W2_data_events));

    const double R_W2     = (1 - facc - fN2 - fpi) / W2_data_events;
    const double F_W2     = W2_bkg_events * R_W2;              // inelastic_frac

    const double dN_in_W2   = std::sqrt(W2_bkg_events);     // Poisson
    const double dN_QE_W2   = std::sqrt(W2_data_events);            // Poisson
    const double dFacc_W2   = errfacc;                       // your TXT value
    const double dFN2_W2    = errfN2;
    const double dFpi_W2    = errfpi;

    const double dF2_W2 =
        std::pow(R_W2 * dN_in_W2, 2) +
        std::pow(F_W2 / W2_data_events * dN_QE_W2, 2) +
        std::pow(W2_bkg_events / W2_data_events * dFacc_W2, 2) +
        std::pow(W2_bkg_events / W2_data_events * dFN2_W2 , 2) +
        std::pow(W2_bkg_events / W2_data_events * dFpi_W2 , 2);

    const double dFin_W2 = std::sqrt(dF2_W2);


    double W2_inelastic_frac = W2_bkg_events * (1 - facc - fN2 - fpi)/W2_data_events;    
    double W2_errinelastic_frac =  dFin_W2;
    
    txt<<"W2_data_events = "<<W2_data_events<<"\n";
    txt<<"W2_fit_events = "<<W2_fit_events<<"\n";
    txt<<"W2_neutron_events = "<<W2_neutron_events<<"\n";
    txt<<"W2_bkg_events = "<<W2_bkg_events<<"\n";
    txt<<"W2_proton_events = "<<W2_proton_events<<"\n";
    txt<<"background_fraction_W2_fit = "<<W2_fit_bkg_fraction<<"\n";
    txt<<"err_background_fraction_W2_fit = "<<W2_fit_err_background_frac<<"\n";
    txt<<"inelastic_fraction_W2_fit = "<<W2_inelastic_frac<<"\n";
    txt<<"err_inelastic_fraction_W2_fit = "<<W2_errinelastic_frac<<"\n";
    txt<<"f_in = "<<W2_inelastic_frac<<"\n";
    txt<<"err_f_in = "<<W2_errinelastic_frac<<"\n";
    txt << "SSR_W2_single_param = " << ssr_neutrons_1 << "\n";
    txt << "SSR_W2_two_param    = " << ssr_neutrons_2 << "\n";
    //txt.close();


    std::cout<<"single paramater method (neutrons)"<<std::endl;
    std::cout<<"alpha_Neutrons : "<<alpha_Neutrons<<std::endl;
    std::cout<<"denom_Neutrons : "<<denom_Neutrons<<std::endl;
    std::cout<<"N_Neutrons : "<<N_Neutrons<<std::endl;

    std::cout<<"double paramater method (neutrons)"<<std::endl;
    std::cout<<"par0_2_Neutrons : "<<par0_2_Neutrons<<std::endl;
    std::cout<<"par1_2_Neutrons : "<<par1_2_Neutrons<<std::endl;

    std::cout<<"delta shifts for W^{2} (neutrons) (not used anymore)"<<std::endl;
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
    TCanvas *C5 = new TCanvas("c5","c5",2400,1500);
    //TCanvas *C2 = new TCanvas("c2","c2",2400,1500);

    C->Divide(2,2);
    C->cd(1);
    
    h_combined->SetLineColor(3);
    hQE_proton_dx_shifted->SetLineColor(6);
    hQE_neutron_dx_shifted->SetLineColor(9);
    hInelastic_dx_shifted->SetLineColor(7);
    hData.SetLineColor(kBlack);

    h_combined->SetLineWidth(4);
    hQE_proton_dx_shifted->SetLineWidth(4);
    hQE_neutron_dx_shifted->SetLineWidth(4);
    hInelastic_dx_shifted->SetLineWidth(4);
    hData.SetLineWidth(4);

    h_combined->SetFillColorAlpha(19,0.1);
    h_combined->SetFillStyle(3009);
    hQE_proton_dx_shifted->SetFillColorAlpha(6,0.5);
    hQE_proton_dx_shifted->SetFillStyle(3004);
    hQE_neutron_dx_shifted->SetFillColorAlpha(9,0.5);
    hQE_neutron_dx_shifted->SetFillStyle(3005);
    hInelastic_dx_shifted->SetFillColorAlpha(7,0.5);
    hInelastic_dx_shifted->SetFillStyle(3009);

    hData.SetMarkerStyle(kFullCircle);

    hData.Draw("p");
    h_combined->Draw("hist same");
    hQE_proton_dx_shifted->Draw("hist same");
    hQE_neutron_dx_shifted->Draw("hist same");
    hInelastic_dx_shifted->Draw("hist same");

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
    hData_W2_Neutrons.SetTitle("single parameter method");
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
    hData_W2_Neutrons.SetTitle("double parameter method");
    hQE_W2_Neutrons.Draw("hist same");
    hInelastic_W2_shift2_Neutrons->Draw("hist same");   
    h_combined_W2_2_Neutrons->Draw("hist same");
    h_QE_W2_split_N->Draw("hist same");
    h_QE_W2_split_P->Draw("hist same");

    C2->cd(3);

    hInelastic_W2_Neutrons.SetLineColor(kBlack);
    hInelastic_W2_Neutrons_eHCALcut_1.SetLineColor(kBlue);
    hInelastic_W2_Neutrons_eHCALcut_2.SetLineColor(kRed);
    hInelastic_W2_Neutrons_eHCALcut_3.SetLineColor(kOrange);
    hInelastic_W2_Neutrons_eHCALcut_4.SetLineColor(kGreen);

    hInelastic_W2_Neutrons.Draw("hist");
    hInelastic_W2_Neutrons_eHCALcut_1.Draw("hist same");
    hInelastic_W2_Neutrons_eHCALcut_2.Draw("hist same");
    hInelastic_W2_Neutrons_eHCALcut_3.Draw("hist same");
    hInelastic_W2_Neutrons_eHCALcut_4.Draw("hist same");

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
    legN1->AddEntry(h_combined_W2_Neutrons,         "Fit: (QEp + R_{np} QEn + #alpha Inel)/(1+R_{np}+ #alpha)", "lf");
    legN1->AddEntry(&hQE_proton_W2_Neutrons,        "QE p (sim)", "f");
    legN1->AddEntry(&hQE_neutron_W2_Neutrons,       "QE n (sim)", "f");
    legN1->AddEntry(hInelastic_W2_shift1_Neutrons,  "Inelastic (sim)", "f");
    legN1->Draw();

    // --- Pad (2): Neutrons, a_QE*QE + a_Inel*Inel^Δ
    C3->cd(2);
    // ... your Draw() calls above ...
    auto legN2 = new TLegend(0.28, 0.28, 0.88, 0.88);
    legN2->SetBorderSize(0);
    legN2->SetFillStyle(0);
    legN2->SetTextSize(0.032);
    legN2->AddEntry(&hData_W2_Neutrons,               "Data", "p");
    legN2->AddEntry(h_combined_W2_2_Neutrons,         "Fit: a_{QE} QE + a_{inel} Inel", "lf");
    legN2->AddEntry(&hQE_W2_Neutrons,                 "QE (sim)", "f");
    legN2->AddEntry(hInelastic_W2_shift2_Neutrons,    "Inelastic (sim)", "f");
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

    C5->Divide(1,1);
    C5->cd(1);
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
    hData_W2_Neutrons.SetTitle("W^{2} distribution fit for Neutrons spot on the #Delta-x #Delta-y");
    hQE_proton_W2_Neutrons.Draw("hist same");
    hQE_neutron_W2_Neutrons.Draw("hist same");
    hInelastic_W2_shift1_Neutrons->Draw("hist same");    
    h_combined_W2_Neutrons->Draw("hist same");

    auto legN1_c5 = new TLegend(0.28, 0.58, 0.58, 0.88);  // x1,y1,x2,y2 in NDC
    legN1_c5->SetBorderSize(0);
    legN1_c5->SetFillStyle(0);
    legN1_c5->SetTextSize(0.032);
    legN1_c5->AddEntry(&hData_W2_Neutrons,             "Data", "p");
    legN1_c5->AddEntry(h_combined_W2_Neutrons,         "Fit: (QEp + R_{np} QEn + #alpha Inel)/(1+R_{np}+ #alpha)", "lf");
    legN1_c5->AddEntry(&hQE_proton_W2_Neutrons,        "QE p (sim)", "f");
    legN1_c5->AddEntry(&hQE_neutron_W2_Neutrons,       "QE n (sim)", "f");
    legN1_c5->AddEntry(hInelastic_W2_shift1_Neutrons,  "Inelastic (sim)", "f");
    legN1_c5->Draw();


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
    //gStyle->SetOptStat(0);

    const double xMin_dx = dxhist_low;
    const double xMax_dx = dxhist_high;
    const double yMin_A  = -15.0;   // adjust if you like
    const double yMax_A  =  15.0;

    TCanvas* Cdx = new TCanvas("Cdx","dx: shape + asym",2400,1500);

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
    //gStyle->SetOptStat(0);

    const double xMin_dy = dyhist_low;
    const double xMax_dy =  dyhist_high;
    //const double yMin_A  = -15.0;   // adjust if you like
    //const double yMax_A  =  15.0;

    TCanvas* Cdy = new TCanvas("Cdy","dy: shape + asym",2400,1500);

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


    //gStyle->SetOptStat(0);

    const double xMin_W2 =  W2hist_low;
    const double xMax_W2 =  W2hist_high;
    //const double yMin_A  = -15.0;   // adjust if you like
    //const double yMax_A  =  15.0;

    TCanvas* CW2 = new TCanvas("CW2","W^{2}: shape + asym",2400,1500);

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
    // W² asymmetry extrapolation: fit inelastic region (W2>2) and extrapolate to elastic (W2<1.6)
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto gA_for_asym     = makeAsymGraph(hW2_for_asym_pos,     hW2_for_asym_neg,     "gA_for_asym",     kRed, kFullCircle);

    // Extract asymmetry vs W2 from the existing dxdy graph
    // Separate into inelastic (W2 > 2) and elastic (W2 < 1.6) regions
    TGraphErrors* gA_inelastic = new TGraphErrors();
    TGraphErrors* gA_elastic   = new TGraphErrors();
    TGraphErrors* gA_all       = (TGraphErrors*)gA_for_asym->Clone("gA_all");

    int nInelastic = 0, nElastic = 0;
    for (int i = 0; i < gA_for_asym->GetN(); ++i) {
        double x, y;
        gA_for_asym->GetPoint(i, x, y);
        double ey = gA_for_asym->GetErrorY(i);

        if (x > 2.0) {
            // Inelastic region (W2 > 2)
            gA_inelastic->SetPoint(nInelastic, x, y);
            gA_inelastic->SetPointError(nInelastic, 0.0, ey);
            ++nInelastic;
        } else if (x < c_.W2_H) {
            // Elastic region (W2 < 1.6)
            gA_elastic->SetPoint(nElastic, x, y);
            gA_elastic->SetPointError(nElastic, 0.0, ey);
            ++nElastic;
        }
    }

    gA_inelastic->SetName("gA_inelastic");
    gA_elastic->SetName("gA_elastic");
    gA_inelastic->SetTitle("Asymmetry (inelastic region, W2 > 2);W^{2} (GeV^{2});Asymmetry (%)");
    gA_elastic->SetTitle("Asymmetry (elastic region, W2 < W2_H);W^{2} (GeV^{2});Asymmetry (%)");

    std::cout << "[InelasticCorrection] Inelastic region (W2 > 2): " << nInelastic << " points\n";
    std::cout << "[InelasticCorrection] Elastic region (W2 < W2_H): " << nElastic << " points\n";

    // Fit a linear model to the inelastic region: A(W2) = a0 + a1*W2
    TF1* fLinear = new TF1("fLinear", "[0] + [1]*x", 2.0, 6.0);
    fLinear->SetParameters(0.0, 0.01);
    if (nInelastic > 1) {
        gA_inelastic->Fit(fLinear, "RQ");
    }

    // TF1* fPL = new TF1("fPL", "[0] + [1]*pow(x,[2])", 2.0, 6.0);
    // fPL->SetParameters(0.0, 1.0, 1.0);
    // if (nInelastic > 1) {
    //     gA_inelastic->Fit(fPL, "RQ");
    // }

    double a0 = fLinear->GetParameter(0);
    double a1 = fLinear->GetParameter(1);
    double ea0 = fLinear->GetParError(0);
    double ea1 = fLinear->GetParError(1);

    std::cout << std::fixed << std::setprecision(6)
              << "[InelasticCorrection] Linear fit: A(W2) = " << a0 << " + " << a1 << "*W2\n"
              << "  Errors: a0 = " << ea0 << ", a1 = " << ea1 << "\n";

    // Create extrapolation graph in elastic region using the fit parameters
    TGraphErrors* gA_extrapolated = new TGraphErrors();
    gA_extrapolated->SetName("gA_extrapolated");
    gA_extrapolated->SetTitle("Asymmetry (extrapolated into elastic region);W^{2} (GeV^{2});Asymmetry (%)");

    // Extrapolate at the same W2 points where we have elastic data
    for (int i = 0; i < gA_elastic->GetN(); ++i) {
        double W2, A_meas;
        gA_elastic->GetPoint(i, W2, A_meas);
        double A_extrap = a0 + a1 * W2;
        // Error propagation: σ_extrap = sqrt((σa0)^2 + (W2*σa1)^2)
        double sigma_extrap = std::sqrt(ea0*ea0 + (W2*ea1)*(W2*ea1));

        gA_extrapolated->SetPoint(i, W2, A_extrap);
        gA_extrapolated->SetPointError(i, 0.0, sigma_extrap);
    }

    // Style the graphs
    gA_inelastic->SetMarkerStyle(kFullCircle);
    gA_inelastic->SetMarkerColor(kRed);
    gA_inelastic->SetLineColor(kRed);
    gA_inelastic->SetMarkerSize(1.2);

    gA_elastic->SetMarkerStyle(kFullSquare);
    gA_elastic->SetMarkerColor(kBlue);
    gA_elastic->SetLineColor(kBlue);
    gA_elastic->SetMarkerSize(1.2);

    gA_extrapolated->SetLineStyle(2);  // dashed
    gA_extrapolated->SetLineColor(kGreen+2);
    gA_extrapolated->SetLineWidth(3);

    fLinear->SetLineStyle(1);
    fLinear->SetLineColor(kMagenta);
    fLinear->SetLineWidth(3);

    // fPL->SetLineStyle(1);
    // fPL->SetLineColor(kGreen+4);
    // fPL->SetLineWidth(3);

    // ======================================================================
    // Canvas split: top = W^2 distribution (counts), bottom = asymmetry vs W^2
    // Both pads share the same X range (xMin..xMax)
    // ======================================================================
    TCanvas* CW2_extrap = new TCanvas("CW2_extrap", "W2 Asymmetry Extrapolation with Distribution", 2400, 1500);
    CW2_extrap->cd();

    // X-range already defined above for both pads

    // Determine histogram Y range (counts)
    double yminW2 = 0.0;
    double ymaxW2 = hW2_for_asym.GetMaximum();
    double xMin = -1.0;
    double xMax = 6.0;
    if (ymaxW2 <= 0) ymaxW2 = 1.0;

    // --- Top pad: W^2 distribution ---------------------------------------
    TPad* padTop_asym = new TPad("padTop_asym", "", 0.0, 0.55, 1.0, 1.0);
    padTop_asym->SetBottomMargin(0.02);
    padTop_asym->SetTopMargin(0.06);
    padTop_asym->SetLeftMargin(0.12);
    padTop_asym->SetRightMargin(0.12);
    padTop_asym->SetTicks(1,1);
    padTop_asym->Draw();
    padTop_asym->cd();

    TH2F* frame_top = new TH2F("frame_top", ";W^{2} (GeV^{2});Counts",
                               100, xMin, xMax, 100, yminW2, 1.05*ymaxW2);
    frame_top->GetXaxis()->SetLabelSize(0); // hide x labels on top pad
    frame_top->GetYaxis()->SetTitleSize(0.045);
    frame_top->GetYaxis()->SetLabelSize(0.04);
    frame_top->Draw();

    TH1D* hW2_draw = (TH1D*)hW2_for_asym.Clone("hW2_draw");
    hW2_draw->GetXaxis()->SetRangeUser(xMin, xMax);
    hW2_draw->SetLineColor(kOrange+2);
    hW2_draw->SetLineWidth(2);
    hW2_draw->SetFillStyle(3004);
    hW2_draw->SetFillColorAlpha(kOrange+2, 0.15);
    hW2_draw->Draw("HIST SAME");

    // --- Bottom pad: asymmetry -------------------------------------------
    CW2_extrap->cd();
    TPad* padBot_asym = new TPad("padBot_asym", "", 0.0, 0.0, 1.0, 0.55);
    padBot_asym->SetTopMargin(0.02);
    padBot_asym->SetBottomMargin(0.12);
    padBot_asym->SetLeftMargin(0.12);
    padBot_asym->SetRightMargin(0.12);
    padBot_asym->SetTicks(1,1);
    padBot_asym->Draw();
    padBot_asym->cd();

    // bottom pad will be drawn by the existing asymmetry frame and plots below
    // if(gA_inelastic->GetN()+gA_elastic->GetN()+gA_extrapolated->GetN() > 0){
    //     double ymin =  1e9, ymax = -1e9;
    //     if(gA_inelastic->GetN()){ ymin = std::min(ymin, getYmin(gA_inelastic)); ymax = std::max(ymax, getYmax(gA_inelastic)); }
    //     if(gA_elastic->GetN()){   ymin = std::min(ymin, getYmin(gA_elastic));   ymax = std::max(ymax, getYmax(gA_elastic)); }
    //     if(gA_extrapolated->GetN()){ ymin = std::min(ymin, getYmin(gA_extrapolated)); ymax = std::max(ymax, getYmax(gA_extrapolated)); }
    //     // add 10% headroom
    //     double pad = 0.10*(ymax - ymin + 1e-9);
    //     ylo = ymin - pad; yhi = ymax + pad;
    // }

    // X-range unified for both pads
    //const double xMin = c_.W2_L;
    //const double xMax =  6.0;

    double ylo = -15.0;
    double yhi =  15.0;

    TH2F* frame_extrap = new TH2F("frame_extrap",
        "Asymmetry Extrapolation with W^{2} Distribution;W^{2} (GeV^{2});Asymmetry (%)",
        100, xMin, xMax, 100, ylo, yhi);
    frame_extrap->GetXaxis()->SetLabelSize(0.035);
    frame_extrap->GetXaxis()->SetTitleSize(0.040);
    // Restore asymmetry axis on the left (default)
    frame_extrap->GetYaxis()->SetLabelSize(0.035);
    frame_extrap->GetYaxis()->SetTitleSize(0.040);
    frame_extrap->GetYaxis()->SetTitleOffset(0.8);
    frame_extrap->Draw();

    gPad->SetGridy(true);
    gStyle->SetGridStyle(2);
    gStyle->SetGridWidth(1);

    // Draw asymmetry elements on the left-axis pad
    gA_inelastic->Draw("P SAME");
    gA_elastic->Draw("P SAME");
    //gA_extrapolated->Draw("P SAME");

    // Extend and draw the fit line across the full X-range
    fLinear->SetRange(xMin, xMax);
    fLinear->Draw("SAME");

    //fPL->SetRange(xMin, xMax);
    //fPL->Draw("SAME");

    // Shaded W^2 bands and vertical edges
    TBox* box_elastic = new TBox(xMin, ylo, c_.W2_H, yhi);
    box_elastic->SetFillColorAlpha(kBlue, 0.08);
    box_elastic->SetLineStyle(0);
    box_elastic->Draw("SAME");

    TBox* box_inelastic = new TBox(2.0, ylo, xMax, yhi);
    box_inelastic->SetFillColorAlpha(kRed, 0.08);
    box_inelastic->SetLineStyle(0);
    box_inelastic->Draw("SAME");

    TLine* L_edge1 = new TLine(c_.W2_H, ylo, c_.W2_H, yhi);
    L_edge1->SetLineStyle(3);
    L_edge1->SetLineColor(kGray+2);
    L_edge1->Draw("SAME");

    TLine* L_edge2 = new TLine(2.0, ylo, 2.0, yhi);
    L_edge2->SetLineStyle(3);
    L_edge2->SetLineColor(kGray+2);
    L_edge2->Draw("SAME");

    // switch back to bottom pad to draw legend/annotations
    padBot->cd();

    auto legExtrap = new TLegend(0.55, 0.58, 0.88, 0.88);
    legExtrap->SetBorderSize(0);
    legExtrap->SetFillStyle(0);
    legExtrap->SetTextSize(0.032);
    legExtrap->AddEntry(gA_inelastic, "Measured (W^{2} > 2)", "p");
    legExtrap->AddEntry(gA_elastic,   "Measured (W^{2} < 1.6)", "p");
    legExtrap->AddEntry(gA_extrapolated, "Extrapolated from fit", "p");
    legExtrap->AddEntry(fLinear, Form("Linear fit: A = %.4f + %.4f#timesW^{2}", a0, a1), "l");
    legExtrap->AddEntry(hW2_draw, "W^{2} distribution", "f");
    legExtrap->Draw();

    TPaveText* pave = new TPaveText(0.12, 0.12, 0.55, 0.42, "NDC");
    pave->SetBorderSize(1);
    pave->SetFillColor(kWhite);
    pave->SetTextAlign(12);
    pave->SetTextSize(0.032);
    pave->AddText("Linear fit parameters:");
    pave->AddText(Form("a_{0} = %.6f #pm %.6f", a0, ea0));
    pave->AddText(Form("a_{1} = %.6f #pm %.6f", a1, ea1));
    pave->AddText(" ");
    pave->AddText(Form("Inelastic region: %d points (W^{2} > 2)", nInelastic));
    pave->AddText(Form("Elastic region: %d points (W^{2} < 1.6)", nElastic));
    //pave->Draw();

    // Compute average asymmetry and error in elastic region based on extrapolation
    if (gA_extrapolated->GetN() > 0) {
        double sum_A = 0.0, sum_sigma2 = 0.0;
        int n_elastic = gA_extrapolated->GetN();
        
        for (int i = 0; i < n_elastic; ++i) {
            double W2, A_extrap;
            gA_extrapolated->GetPoint(i, W2, A_extrap);
            double sigma_extrap = gA_extrapolated->GetErrorY(i);
            
            // Accumulate weighted by inverse variance (1/sigma^2)
            if (sigma_extrap > 0) {
                double weight = 1.0 / (sigma_extrap * sigma_extrap);
                sum_A += A_extrap * weight;
                sum_sigma2 += weight;
            }
        }
        
        double avg_A = sum_A / sum_sigma2;
        double avg_sigma = 1.0 / std::sqrt(sum_sigma2);
        
        std::cout << std::fixed << std::setprecision(6)
                  << "\n[InelasticCorrection] Average Asymmetry in QE Region (Extrapolated):\n"
                  << "  Average A = " << avg_A << " % ± " << avg_sigma << " %\n"
                  << "  Average A (fraction) = " << avg_A/100.0 << " ± " << avg_sigma/100.0 << "\n"
                  << "  Computed from " << n_elastic << " elastic points (extrapolated using inelastics) using inverse-variance weighting\n";

        // Write averaged value to the main corrections txt file (same file opened earlier)
        // Check that the output file stream 'txt' is open before writing.
        if (txt.good()) {
            txt << "avg_extrap_A_percent = " << std::fixed << std::setprecision(6) << avg_A << "\n";
            txt << "err_avg_extrap_A_percent = " << std::fixed << std::setprecision(6) << avg_sigma << "\n";
            txt << "A_in = " << std::fixed << std::setprecision(8) << (avg_A/100.0) << "\n";
            txt << "err_A_in = " << std::fixed << std::setprecision(8) << (avg_sigma/100.0) << "\n";
            txt << "num_elastic_points_used = " << n_elastic << "\n";
            txt.flush();
            txt.close();
        } else {
            std::cerr << "[InelasticCorrection] Warning: output file not open, cannot write averaged extrapolated asymmetry.\n";
        }
    }

    // Save
    CW2_extrap->Print(Form("images/%s/W2_Asymmetry_Extrapolation_%s.png", kin_, kin_));

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas* Chcal = new TCanvas("Chcal","secondary clusters",2400,1500);

    Chcal->Divide(2,2);
    Chcal->cd(1);
    h_eratio.Draw();
    Chcal->cd(2);
    h_eratio_test.Draw();
    Chcal->cd(3);
    h_eratio_wide_bins.Draw();
    Chcal->cd(4);
    h_tdiff.Draw();

    TCanvas* Chcal1 = new TCanvas("Chcal1","secondary clusters",2400,1500);
    Chcal1->Divide(2,2);
    Chcal1->cd(1);
    hdist.Draw();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Save
    Cdx->Print(Form("images/%s/Asymmetry_vs_Dx_stacked_%s.png", kin_, kin_));
    Cdy->Print(Form("images/%s/Asymmetry_vs_Dy_stacked_%s.png", kin_, kin_));
    CW2->Print(Form("images/%s/Asymmetry_vs_W2_stacked_%s.png", kin_, kin_));

    // Save

    C->Print(Form("images/%s/InelasticCorrection_%s.png",kin_,kin_));
    C1->Print(Form("images/%s/InelasticCorrection_W2_%s.png",kin_,kin_));
    C2->Print(Form("images/%s/InelasticCorrection_W2_Neutrons_comparison_%s.png",kin_,kin_));
    C3->Print(Form("images/%s/InelasticCorrection_W2_Neutrons_legend%s.png",kin_,kin_));
    C4->Print(Form("images/%s/InelasticDxDy%s.png",kin_,kin_));
    C5->Print(Form("images/%s/InelasticCorrection_W2_Neutrons_%s.png",kin_,kin_));
    Casym->Print(Form("images/%s/Asymmetry_vs_W2_%s.png", kin_, kin_));
    Casym1->Print(Form("images/%s/Asymmetry_vs_Dx_%s.png", kin_, kin_));
    Chcal->Print(Form("images/%s/Secondary_cluster_%s.png",kin_,kin_));
    Chcal1->Print(Form("images/%s/Secondary_cluster_1_%s.png",kin_,kin_));

}
