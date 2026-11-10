// Standalone GEn-II GEN3 2D dx-dy background model comparison.
// From Final_analysis/: root -l -b -q 'compare_inelastic_dxdy_models_ABC.C+()'
// Requires ROOT 6.x + Minuit2. No external JSON library.
//
// F(dx,dy) = A * [ P(dx,dy) + B*N(dx,dy) + C*I(dx,dy) ].
// Unit-normalized P,N,I over fit grid. A=N_protons, B=N_neutrons/N_protons,
// C=N_inelastic/N_protons.  I = simulation or 2D Bernstein polynomial
// of degree n in BOTH dx and dy with nonnegative coefficients.
// Objective: Poisson deviance, all 2D cells including zero-count cells.
// Fit diagnostics: Pearson chi2/ndf (approximate for sparse grids), dev/ndf,
// AIC; neutron elliptical ROI uses (dx/dx_r)^2+(dy/dy_r)^2 <= 1.
// Unlike InelasticCorrection.cpp's TH2 histograms, here ROOT X=dx and Y=dy.
// NEW: data stays on the original fit grid, while simulation is filled on
// larger (dx,dy) support. Translations integrate bin contents via a bilinear
// 2D cumulative histogram and renormalize within the FIT grid. This removes
// artificial shifts/discontinuities near the fit plot boundaries.

#include <TCanvas.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TPad.h>
#include <TString.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dxdyabc {
static const char* kData =
  "/media/vimukthi/Sata/Research_JLab/outfiles/pass3/intimehighest/"
  "QE_data_GEN3_sbs100p_nucleon_np_model2_sbstrackingon_intimehighest.root";
static const char* kQE =
  "/media/vimukthi/Sata/Research_JLab/outfiles/SIMC/noNRad_correct_raster/"
  "QE_sim_simc_He3eN_GEN3_sbs100p_nucleon_np_model2_elastic_simc.root";
static const char* kInel =
  "/media/vimukthi/Sata/Research_JLab/outfiles/inelastic_simulation_2025/"
  "QE_sim_GEN3_sbs100p_nucleon_np_model2_inelastic.root";

struct Cuts {
    double eHCAL_L, EoverP_L, EoverP_H, W2_L, W2_H;
    double dx_L, dx_H, dy_L, dy_H, dx_r, dy_r;
    double coin_L, coin_H, runnum_L, runnum_H;
};
std::string contents(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot read " + path);
    std::ostringstream s; s<<f.rdbuf(); return s.str();
}
double number(const std::string& text, const std::string& key) {
    std::regex exp("\\\""+key+"\\\"\\s*:\\s*([-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?)");
    std::smatch match;
    if (!std::regex_search(text,match,exp)) throw std::runtime_error("Missing numeric cut '"+key+"'");
    return std::stod(match[1].str());
}
Cuts readCuts(const std::string& path) {
    const auto j=contents(path);
    Cuts c;
    c.eHCAL_L=number(j,"eHCAL_L"); c.EoverP_L=number(j,"EoverP_L");c.EoverP_H=number(j,"EoverP_H");
    c.W2_L=number(j,"W2_L");c.W2_H=number(j,"W2_H");
    c.dx_L=number(j,"dx_L");c.dx_H=number(j,"dx_H");
    c.dy_L=number(j,"dy_L");c.dy_H=number(j,"dy_H");
    c.dx_r=number(j,"dx_r");c.dy_r=number(j,"dy_r");
    c.coin_L=number(j,"coin_L");c.coin_H=number(j,"coin_H");
    c.runnum_L=number(j,"runnum_L");c.runnum_H=number(j,"runnum_H");
    if (!(c.dx_r>0 && c.dy_r>0 && c.dy_L<c.dy_H)) throw std::runtime_error("Invalid dx_r/dy_r or dy interval");
    return c;
}
std::set<int> readGoodRuns(const std::string& path) {
    const auto text=contents(path);
    std::regex line("\\\"([0-9]+)\\\"\\s*:\\s*\\{([^{}]*)\\}");
    std::regex hel("\\\"helicity\\\"\\s*:\\s*(true|false|1|0)",std::regex_constants::icase);
    std::regex mol("\\\"moller\\\"\\s*:\\s*(true|false|1|0)",std::regex_constants::icase);
    std::set<int> runs;
    for (auto it=std::sregex_iterator(text.begin(),text.end(),line);it!=std::sregex_iterator();++it) {
        const std::string body=(*it)[2].str();std::smatch m1,m2;
        if (!std::regex_search(body,m1,hel)||!std::regex_search(body,m2,mol))continue;
        auto okay=[](std::string s){return s=="1"||s=="true"||s=="True"||s=="TRUE";};
        if (okay(m1[1].str())&&okay(m2[1].str())) runs.insert(std::stoi((*it)[1].str()));
    }
    if (runs.empty())throw std::runtime_error("No good runs in "+path+"; use applyRunQuality=false if appropriate");
    return runs;
}

struct HistogramInput {
    int nx=0,ny=0;
    double xmin=0,xmax=0,ymin=0,ymax=0;
    std::vector<double> data,p,n,I;
    int index(int ix,int iy)const {return ix*ny+iy;}
    double wx()const{return (xmax-xmin)/nx;}
    double wy()const{return (ymax-ymin)/ny;}
    double cx(int i)const{return xmin+(i+0.5)*wx();}
    double cy(int j)const{return ymin+(j+0.5)*wy();}
    int getIndex(double dx,double dy) const {
        int ix=(int)std::floor((dx-xmin)/wx()),iy=(int)std::floor((dy-ymin)/wy());
        return ix>=0&&ix<nx&&iy>=0&&iy<ny ? index(ix,iy) : -1;
    }
};
double sum(const std::vector<double>& a) {return std::accumulate(a.begin(),a.end(),0.0);}
void normalize(std::vector<double>& a) {
    double s=sum(a);if (!(s>0))throw std::runtime_error("Empty template");
    for (auto& v:a) {
        if (!std::isfinite(v)||v<-1e-9*s) throw std::runtime_error("Nonfinite or negative weighted MC bin");
        v=std::max(0.,v)/s;
    }
}
std::unique_ptr<TFile> openFile(const std::string& path) {
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(),"READ"));
    if (!f||f->IsZombie())throw std::runtime_error("Cannot open "+path);
    return f;
}
TTree* tree(TFile* f) {
    for (const char* name:{"Tout","T"}) {
        if (auto* t=dynamic_cast<TTree*>(f->Get(name)))return t;
    }
    throw std::runtime_error(std::string("Tout or T not found in ")+f->GetName());
}
using Leaves=std::map<std::string,TLeaf*>;
Leaves leavesFor(TTree* t, const std::vector<std::string>& names) {
    t->SetBranchStatus("*",0);
    Leaves m;
    for (const auto& name:names) {
        auto* leaf=t->GetLeaf(name.c_str());
        if (!leaf)throw std::runtime_error(std::string(t->GetName())+" missing leaf "+name);
        t->SetBranchStatus(leaf->GetBranch()->GetName(),1);
        m[name]=leaf;
    }
    return m;
}
double value(const Leaves& l,const char* name){return l.at(name)->GetValue(0);}
void fill(HistogramInput& in,const std::string& filename,const Cuts& c,
          const std::set<int>& runs,bool quality,int kind,Long64_t maxEntries,
          double preshiftP,double preshiftN,double preshiftI) {
    auto file=openFile(filename);auto* t=tree(file.get());
    std::vector<std::string> names={"dx","dy","vz","W2","eHCAL","ePS","eSH","trP"};
    if(kind==0){for(const char* v:{"runnum","ntrack","coin_time","helicity"})names.emplace_back(v);}
    else names.emplace_back("weight");
    if(kind==1)names.emplace_back("fnucl");
    Leaves l=leavesFor(t,names);
    Long64_t ne=t->GetEntries();if(maxEntries>=0)ne=std::min(ne,maxEntries);
    std::cout<<"[input] "<<(kind==0?"data":kind==1?"QE":"inel")<<" entries="<<ne<<" file="<<filename<<std::endl;
    Long64_t accepted=0;
    for(Long64_t i=0;i<ne;++i) {
        t->GetEntry(i);
        double trp=value(l,"trP"),dx=value(l,"dx"),dy=value(l,"dy"),w2=value(l,"W2");
        if(!(trp>0)||!std::isfinite(dx)||!std::isfinite(dy)||!std::isfinite(w2))continue;
        const double ep=(value(l,"ePS")+value(l,"eSH"))/trp;
        if(!std::isfinite(ep)||std::abs(value(l,"vz"))>=0.27||
           value(l,"eHCAL")<=c.eHCAL_L||value(l,"ePS")<=0.2||
           ep<=c.EoverP_L||ep>=c.EoverP_H||w2<=c.W2_L||w2>=c.W2_H||
           dy<=in.ymin||dy>=in.ymax)continue;
        if(kind==0) {
            int run=(int)std::llround(value(l,"runnum"));
            if(run<c.runnum_L||run>c.runnum_H||value(l,"ntrack")<1||
               std::abs(value(l,"helicity"))!=1||
               value(l,"coin_time")<c.coin_L||value(l,"coin_time")>c.coin_H||
               (quality&&!runs.count(run)))continue;
        }
        double weight=kind==0?1.0:value(l,"weight");
        if(!std::isfinite(weight))continue;
        std::vector<double>* dst=&in.data;
        if(kind==1) {
            int nucl=(int)std::llround(value(l,"fnucl"));
            if(nucl==1){dst=&in.p;dx+=preshiftP;}
            else if(nucl==0){dst=&in.n;dx+=preshiftN;}
            else continue;
        } else if(kind==2){dst=&in.I;dx+=preshiftI;}
        int index=in.getIndex(dx,dy);
        if(index<0)continue;
        (*dst)[index]+=weight;
        ++accepted;
        if(i>0&&i%5000000==0)std::cout<<"   "<<i<<"/"<<ne<<std::endl;
    }
    std::cout<<"[input] accepted="<<accepted<<"/"<<ne<<std::endl;
}

// A 2D cumulative histogram supports exact shifted rectangular BIN integrals
// of the piecewise-constant MC template. This uses extended MC support, not
// the clipped fit histogram. Positive dx shift = P_fit(x,y)=P_MC(x-dx,y-dy).
struct ShiftedShape2D {
    int nx,ny;
    double xmin,ymin,wx,wy;
    std::vector<double> cdf; // (nx+1)*(ny+1), at CELL EDGES

    ShiftedShape2D(const std::vector<double>& raw,const HistogramInput& source)
      :nx(source.nx),ny(source.ny),xmin(source.xmin),ymin(source.ymin),
       wx(source.wx()),wy(source.wy()),cdf((nx+1)*(ny+1),0.) {
        if((int)raw.size()!=nx*ny)throw std::runtime_error("Invalid 2D MC source dimensions");
        std::vector<double> arr(raw);
        normalize(arr);
        for(int ix=1;ix<=nx;++ix)for(int iy=1;iy<=ny;++iy) {
            cdf[index(ix,iy)] = arr[(ix-1)*ny+iy-1] +
                cdf[index(ix-1,iy)]+cdf[index(ix,iy-1)]-cdf[index(ix-1,iy-1)];
        }
    }
    int index(int i,int j) const{return i*(ny+1)+j;}
    double integralTo(double x,double y)const {
        const double u=std::max(0.,std::min(double(nx),(x-xmin)/wx));
        const double v=std::max(0.,std::min(double(ny),(y-ymin)/wy));
        const int i=std::min(nx-1,int(std::floor(u)));
        const int j=std::min(ny-1,int(std::floor(v)));
        const double a=u-i,b=v-j;
        return (1-a)*((1-b)*cdf[index(i,j)]+b*cdf[index(i,j+1)])+
                   a *((1-b)*cdf[index(i+1,j)]+b*cdf[index(i+1,j+1)]);
    }
    std::vector<double> shifted(const HistogramInput& fit,double dx,double dy)const {
        std::vector<double> corners((fit.nx+1)*(fit.ny+1),0.);
        const int stride=fit.ny+1;
        for(int ix=0;ix<=fit.nx;++ix)for(int iy=0;iy<=fit.ny;++iy){
            const double x=fit.xmin+ix*fit.wx()-dx;
            const double y=fit.ymin+iy*fit.wy()-dy;
            corners[ix*stride+iy]=integralTo(x,y);
        }
        std::vector<double> result(fit.data.size(),0.);
        for(int ix=0;ix<fit.nx;++ix)for(int iy=0;iy<fit.ny;++iy) {
            // Integral over [xlo-dx,xhi-dx] x [ylo-dy,yhi-dy].
            const double v=corners[(ix+1)*stride+iy+1]-corners[ix*stride+iy+1]
                         -corners[(ix+1)*stride+iy]+corners[ix*stride+iy];
            result[fit.index(ix,iy)]=std::max(0.,v);
        }
        normalize(result); // Component yields A, AB, AC apply within FIT grid.
        return result;
    }
};

double choose(int n,int k) {
    double r=1;for(int i=1;i<=k;++i)r*=double(n-i+1)/i;return r;
}
std::vector<std::vector<double>> bernstein(int degree,int nb) {
    std::vector<std::vector<double>> b(degree+1,std::vector<double>(nb));
    for(int k=0;k<=degree;++k) {
        for(int i=0;i<nb;++i) {
            double t=(i+0.5)/nb;
            b[k][i]=choose(degree,k)*std::pow(t,k)*std::pow(1-t,degree-k);
        }
        normalize(b[k]);
    }
    return b;
}

struct FitProblem {
    const HistogramInput& h;
    ShiftedShape2D p,n,I;
    int degree,nShape,nShifts,nParams;
    bool floatShifts;
    std::vector<std::vector<double>> bx,by;
    explicit FitProblem(const HistogramInput& hh,const HistogramInput& source,int d,bool shift)
      :h(hh),p(source.p,source),n(source.n,source),I(source.I,source),degree(d),
       nShape(d<0?0:(d+1)*(d+1)-1),nShifts(!shift?0:(d<0?6:4)),
       nParams(3+nShape+nShifts),floatShifts(shift) {
        if(d>=0){bx=bernstein(d,h.nx);by=bernstein(d,h.ny);}
    }
    std::vector<double> polynomial(const double* pars)const {
        const int k=degree+1;
        std::vector<double> w(k*k,1.0);
        double maxz=0.;for(int i=0;i<nShape;++i)maxz=std::max(maxz,pars[3+i]);
        for(int i=0;i<nShape;++i)w[i]=std::exp(pars[3+i]-maxz);
        w.back()=std::exp(-maxz);
        const double sw=sum(w);
        for(double& a:w)a/=sw;
        std::vector<double> out(h.data.size(),0.);
        for(int i=0;i<h.nx;++i)for(int j=0;j<h.ny;++j) {
            double value=0.;
            for(int a=0;a<k;++a)for(int b=0;b<k;++b) {
                value+=w[a*k+b]*bx[a][i]*by[b][j];
            }
            out[h.index(i,j)]=value;
        }
        return out; // unit sum by construction
    }
    void components(const double* pars,std::vector<double>& vp,
                    std::vector<double>& vn,std::vector<double>& vi) const {
        const double A=pars[0],B=pars[1],C=pars[2];
        const int q=3+nShape;
        const double px=floatShifts?pars[q]:0.,py=floatShifts?pars[q+1]:0.;
        const double nx=floatShifts?pars[q+2]:0.,ny=floatShifts?pars[q+3]:0.;
        const double xi=(floatShifts&&degree<0)?pars[q+4]:0.;
        const double yi=(floatShifts&&degree<0)?pars[q+5]:0.;
        vp=p.shifted(h,px,py);
        vn=n.shifted(h,nx,ny);
        vi=degree<0?I.shifted(h,xi,yi):polynomial(pars);
        for(size_t k=0;k<vp.size();++k){vp[k]*=A;vn[k]*=A*B;vi[k]*=A*C;}
    }
    double deviance(const double* pars)const {
        std::vector<double> pp,nn,ii;components(pars,pp,nn,ii);
        double score=0.;
        for(size_t k=0;k<h.data.size();++k) {
            const double observed=h.data[k];
            const double predicted=std::max(1e-10,pp[k]+nn[k]+ii[k]);
            score+=2*(predicted-observed+(observed>0?observed*std::log(observed/predicted):0));
        }
        return score;
    }
};

struct Result {
    std::string name;
    bool success=false;
    int npars=0,ndf=0;
    double chi2=0,reduced=0,dev=0,devRed=0,pvalue=0,aic=0;
    double bgEll=0,dataEll=0,bgDx=0,dataDx=0;
    std::vector<double> pars,p,n,I,model,resid,pull;
};
Result fitOne(const std::string& name,const FitProblem& pr,const Cuts& cuts,
              double sxMax,double syMax) {
    const int npar=pr.nParams;
    if((int)pr.h.data.size()<=npar)throw std::runtime_error("Number of bins <= number of fit parameters");
    ROOT::Math::Functor fcn([&](const double* x){return pr.deviance(x);},npar);
    std::unique_ptr<ROOT::Math::Minimizer> min(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));
    if(!min)throw std::runtime_error("Minuit2 unavailable");
    min->SetFunction(fcn);
    min->SetMaxFunctionCalls(700000);
    min->SetMaxIterations(15000);
    min->SetTolerance(1e-5);
    min->SetPrintLevel(0);
    const double Ndata=sum(pr.h.data),A0=Ndata/(1+0.2+0.3);
    min->SetLimitedVariable(0,"A_protons",A0,std::max(0.5,A0*0.002),0.,3.*Ndata);
    min->SetLimitedVariable(1,"B_n_over_p",0.20,0.005,0.,20.);
    min->SetLimitedVariable(2,"C_inel_over_p",0.30,0.005,0.,20.);
    for(int k=0;k<pr.nShape;++k)min->SetLimitedVariable(3+k,"poly_log_weight_"+std::to_string(k),0.,0.05,-9.,9.);
    if(pr.floatShifts) {
        for(int j=0;j<pr.nShifts;++j) {
            const bool isX=j%2==0;
            const double limit=isX?sxMax:syMax;
            if(!(limit>0))throw std::runtime_error("Shift limit must be positive");
            min->SetLimitedVariable(3+pr.nShape+j,"shift_"+std::to_string(j),0.,
                                    std::min(0.002,limit*0.1),-limit,limit);
        }
    }
    const bool converged=min->Minimize();
    const double* params=min->X();
    if(!params||!std::isfinite(min->MinValue()))throw std::runtime_error("No finite fit result for "+name);
    Result r;r.name=name;r.success=converged;r.npars=npar;
    r.pars.assign(params,params+npar);
    pr.components(r.pars.data(),r.p,r.n,r.I);
    const size_t nb=pr.h.data.size();
    r.model.resize(nb);r.resid.resize(nb);r.pull.resize(nb);
    r.ndf=(int)nb-npar;
    for(int i=0;i<pr.h.nx;++i)for(int j=0;j<pr.h.ny;++j) {
        const int k=pr.h.index(i,j);
        const double observed=pr.h.data[k],mu=std::max(1e-10,r.p[k]+r.n[k]+r.I[k]);
        r.model[k]=mu;r.resid[k]=observed-mu;r.pull[k]=(observed-mu)/std::sqrt(mu);
        r.chi2+=(observed-mu)*(observed-mu)/mu;
        r.dev+=2*(mu-observed+(observed>0?observed*std::log(observed/mu):0));
        const double x=pr.h.cx(i),y=pr.h.cy(j);
        if((x/cuts.dx_r)*(x/cuts.dx_r)+(y/cuts.dy_r)*(y/cuts.dy_r)<=1.0) {
            r.bgEll+=r.I[k];r.dataEll+=observed;
        }
        if(x>=cuts.dx_L&&x<=cuts.dx_H){r.bgDx+=r.I[k];r.dataDx+=observed;}
    }
    r.reduced=r.chi2/r.ndf;r.devRed=r.dev/r.ndf;
    r.pvalue=TMath::Prob(r.chi2,r.ndf);r.aic=r.dev+2*npar;
    return r;
}

std::unique_ptr<TH2D> hist2(const HistogramInput& h,const std::vector<double>& bins,
                             const std::string& name) {
    auto hist=std::make_unique<TH2D>(name.c_str(),"",h.nx,h.xmin,h.xmax,h.ny,h.ymin,h.ymax);
    hist->SetDirectory(nullptr);
    for(int i=0;i<h.nx;++i)for(int j=0;j<h.ny;++j)
        hist->SetBinContent(i+1,j+1,bins[h.index(i,j)]);
    hist->GetXaxis()->SetTitle("dx (m)");hist->GetYaxis()->SetTitle("dy (m)");
    return hist;
}
std::unique_ptr<TH1D> projection(const HistogramInput& h,const std::vector<double>& a,
                                  const std::string& name,bool alongDx) {
    const int nb=alongDx?h.nx:h.ny;
    auto out=std::make_unique<TH1D>(name.c_str(),"",nb,alongDx?h.xmin:h.ymin,alongDx?h.xmax:h.ymax);
    out->SetDirectory(nullptr);
    for(int i=0;i<h.nx;++i)for(int j=0;j<h.ny;++j) {
        int ib=alongDx?i:j;
        out->AddBinContent(ib+1,a[h.index(i,j)]);
    }
    return out;
}
void save(TCanvas& c,const std::string& file) {
    c.SaveAs((file+".pdf").c_str());c.SaveAs((file+".png").c_str());
}
void plotOne(const HistogramInput& h,const Result& r,const std::string& outdir,TFile& root) {
    auto hd=hist2(h,h.data,"data2d_"+r.name);
    auto hp=hist2(h,r.p,"proton2d_"+r.name);
    auto hn=hist2(h,r.n,"neutron2d_"+r.name);
    auto hi=hist2(h,r.I,"background2d_"+r.name);
    auto hm=hist2(h,r.model,"fit2d_"+r.name);
    auto hr=hist2(h,r.resid,"residual2d_"+r.name);
    auto hpl=hist2(h,r.pull,"pull2d_"+r.name);
    root.cd();for(TH2D* object:{hd.get(),hp.get(),hn.get(),hi.get(),hm.get(),hr.get(),hpl.get()})object->Write();
    TCanvas c2(Form("c2_%s",r.name.c_str()),"2D comparison",1200,900);
    c2.Divide(2,2);
    for(int q=0;q<4;++q) {
        c2.cd(q+1)->SetRightMargin(0.16);
        TH2D* hh= q==0?hd.get():q==1?hm.get():q==2?hr.get():hpl.get();
        const char* title=q==0?"Data":q==1?"Total model":q==2?"Data - model":"Pearson pull";
        hh->SetTitle(Form("%s: %s;dx (m);dy (m)",r.name.c_str(),title));
        hh->SetStats(0);hh->Draw("COLZ");
    }
    save(c2,outdir+"/fit2d_"+r.name);

    TCanvas cproj(Form("cproj_%s",r.name.c_str()),"2D fit projections",1250,860);
    cproj.Divide(2,2);
    std::vector<std::unique_ptr<TH1D>> keep;
    std::vector<std::unique_ptr<TLegend>> legends;
    std::vector<std::unique_ptr<TLine>> guideLines;
    for(int coord=0;coord<2;++coord) {
        bool dx=coord==0;
        auto d=projection(h,h.data,r.name+"_data_"+std::to_string(coord),dx);
        auto m=projection(h,r.model,r.name+"_model_"+std::to_string(coord),dx);
        auto p=projection(h,r.p,r.name+"_proton_"+std::to_string(coord),dx);
        auto n=projection(h,r.n,r.name+"_neutron_"+std::to_string(coord),dx);
        auto I=projection(h,r.I,r.name+"_background_"+std::to_string(coord),dx);
        auto pulls=std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(d->Clone((r.name+"_pull_"+std::to_string(coord)).c_str())));
        pulls->SetDirectory(nullptr);
        for(int i=1;i<=d->GetNbinsX();++i) {
            const double y=d->GetBinContent(i),mu=std::max(1e-10,m->GetBinContent(i));
            d->SetBinError(i,std::sqrt(std::max(0.,y)));
            pulls->SetBinContent(i,(y-mu)/std::sqrt(mu));
        }
        cproj.cd(coord+1);gPad->SetBottomMargin(0.12);gPad->SetLeftMargin(.10);
        d->SetTitle(Form("%s projection; %s (m); Counts / bin",r.name.c_str(),dx?"dx":"dy"));
        d->SetMarkerStyle(20);d->SetMarkerSize(0.5);d->SetLineColor(kBlack);
        double maxD=d->GetMaximum(),maxM=m->GetMaximum();
        d->SetMaximum(1.26*std::max(maxD,maxM));d->Draw("E1");
        for(auto combo:{std::pair<TH1D*,int>{m.get(),kGreen+2},
                        {p.get(),kMagenta+1},{n.get(),kBlue+1},{I.get(),kCyan+2}}){
            combo.first->SetLineColor(combo.second);combo.first->SetLineWidth(2);
            combo.first->Draw("HIST SAME");
        }
        d->Draw("E1 SAME");
        auto leg=std::make_unique<TLegend>(0.62,0.58,0.93,0.88);
        leg->SetBorderSize(0);leg->SetFillStyle(0);leg->SetTextSize(.034);
        leg->AddEntry(d.get(),"Data","lep");leg->AddEntry(m.get(),"Fit","l");
        leg->AddEntry(p.get(),"QE proton","l");leg->AddEntry(n.get(),"QE neutron","l");
        leg->AddEntry(I.get(),"Inelastic background","l");leg->Draw();
        legends.emplace_back(std::move(leg));
        cproj.cd(coord+3);gPad->SetLeftMargin(.10);gPad->SetBottomMargin(.15);
        pulls->SetTitle(Form("Projection pull; %s (m); (Data-Fit)/#sqrt{Fit}",dx?"dx":"dy"));
        pulls->SetMarkerStyle(20);pulls->SetMarkerSize(.45);pulls->SetLineColor(kBlack);
        pulls->Draw("P");
        auto line=std::make_unique<TLine>(dx?h.xmin:h.ymin,0,dx?h.xmax:h.ymax,0);
        line->SetLineStyle(2);line->Draw();
        guideLines.emplace_back(std::move(line));
        root.cd();for(TH1D* ob:{d.get(),m.get(),p.get(),n.get(),I.get(),pulls.get()})ob->Write();
        keep.emplace_back(std::move(d));keep.emplace_back(std::move(m));
        keep.emplace_back(std::move(p));keep.emplace_back(std::move(n));
        keep.emplace_back(std::move(I));keep.emplace_back(std::move(pulls));
    }
    save(cproj,outdir+"/projections_"+r.name);
}
void compare(const std::vector<Result>& results,const std::string& outdir) {
    TCanvas c("cComparison","2D background model comparison",1080,650);
    c.SetBottomMargin(.15);
    TH1D h("chi2summary","GEN3 2D background models;Background model;Pearson #chi^{2} / ndf",
           (int)results.size(),0,(int)results.size());
    h.SetStats(0);h.SetFillColor(kAzure-9);
    double maxchi=1.;
    for(size_t i=0;i<results.size();++i){
        h.SetBinContent(i+1,results[i].reduced);
        h.GetXaxis()->SetBinLabel(i+1,results[i].name.c_str());
        maxchi=std::max(maxchi,results[i].reduced);
    }
    h.SetMinimum(0);h.SetMaximum(maxchi*1.24);h.Draw("HIST TEXT0");
    TLine one(0,1,results.size(),1);one.SetLineStyle(2);one.Draw();
    save(c,outdir+"/reduced_chi2_comparison");
}
} // namespace dxdyabc

// Empty path arguments use the GEN3 defaults declared above.
// Choose --equivalent here by changing arguments in the ROOT invocation.
// Example:
// root -l -b -q 'compare_inelastic_dxdy_models_ABC.C+("","","","cuts/GEN3.json","DB/run_quality.json","out_plots/test2d",60,40,true,"2,3,4,5",true,-1)'
void compare_inelastic_dxdy_models_ABC(
    const char* dataFile="",
    const char* qeFile="",
    const char* inelFile="",
    const char* cutsFile="cuts/GEN3.json",
    const char* runQualityFile="DB/run_quality.json",
    const char* outdir="out_plots/inelastic_dxdy_background_models_GEN3_ABC_cpp",
    int binsDx=60,
    int binsDy=40,
    bool floatShifts=false,
    const char* degreeList="2,3,4,5",
    bool applyRunQuality=true,
    Long64_t maxEntries=-1,
    double xmin=-2.2,
    double xmax=0.8,
    double dyPadding=0.05,
    double shiftLimitX=0.1,
    double shiftLimitY=0.1,
    double preShiftP=0.075,
    double preShiftN=0.0,
    double preShiftI=0.0,
    double templateGuardX=0.45,
    double templateGuardY=0.15) {
    using namespace dxdyabc;
    try {
        gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);
        if(binsDx<8||binsDy<8||xmin>=xmax)throw std::runtime_error("Invalid 2D histogram binning/range");
        if(dyPadding<0||templateGuardX<0||templateGuardY<0||shiftLimitX<0||shiftLimitY<0)
            throw std::runtime_error("Negative dy padding, template guard, or shift bound");
        Cuts cuts=readCuts(cutsFile);
        std::set<int> runs=applyRunQuality?readGoodRuns(runQualityFile):std::set<int>{};
        HistogramInput h;
        h.nx=binsDx;h.ny=binsDy;h.xmin=xmin;h.xmax=xmax;
        h.ymin=cuts.dy_L-dyPadding;h.ymax=cuts.dy_H+dyPadding;
        const int nb=binsDx*binsDy;
        h.data.assign(nb,0.);h.p.assign(nb,0.);h.n.assign(nb,0.);h.I.assign(nb,0.);
        // Keep the DATA selection/fitting/display grid completely unchanged.
        // MC gets an extended grid with the SAME bin widths and integer guards.
        const double wx=h.wx(),wy=h.wy();
        const int guardX=(int)std::ceil(std::max(templateGuardX,shiftLimitX+3*wx)/wx);
        const int guardY=(int)std::ceil(std::max(templateGuardY,shiftLimitY+3*wy)/wy);
        HistogramInput source;
        source.nx=binsDx+2*guardX;source.ny=binsDy+2*guardY;
        source.xmin=xmin-guardX*wx;source.xmax=xmax+guardX*wx;
        source.ymin=h.ymin-guardY*wy;source.ymax=h.ymax+guardY*wy;
        const int nsource=source.nx*source.ny;
        source.data.assign(nsource,0.);source.p.assign(nsource,0.);
        source.n.assign(nsource,0.);source.I.assign(nsource,0.);
        std::cout<<"[support] 2D fit dx=["<<h.xmin<<","<<h.xmax<<"] dy=["
                 <<h.ymin<<","<<h.ymax<<"] ; MC source dx=["<<source.xmin
                 <<","<<source.xmax<<"] dy=["<<source.ymin<<","<<source.ymax
                 <<"] (guard bins dx="<<guardX<<" dy="<<guardY<<")"<<std::endl;
        fill(h,(dataFile&&*dataFile)?dataFile:kData,cuts,runs,applyRunQuality,0,maxEntries,0,0,0);
        fill(source,(qeFile&&*qeFile)?qeFile:kQE,cuts,runs,applyRunQuality,1,maxEntries,preShiftP,preShiftN,0);
        fill(source,(inelFile&&*inelFile)?inelFile:kInel,cuts,runs,applyRunQuality,2,maxEntries,0,0,preShiftI);
        for(auto v:std::vector<std::pair<std::string,double>>{{"data",sum(h.data)},{"proton",sum(source.p)},
                  {"neutron",sum(source.n)},{"inelastic",sum(source.I)}}){
            if(!(v.second>0))throw std::runtime_error("Empty "+v.first+" histogram: check cuts or paths");
            std::cout<<"[hists] "<<v.first<<"="<<v.second<<std::endl;
        }
        gSystem->mkdir(outdir,true);
        std::set<int> degrees;
        std::stringstream iss(degreeList?degreeList:"");std::string token;
        while(std::getline(iss,token,',')){
            if(token.empty())continue;
            const int d=std::stoi(token);
            if(d<0||d>8)throw std::runtime_error("Supported polynomial degrees are 0..8");
            degrees.insert(d);
        }
        if(degrees.empty())throw std::runtime_error("Empty polynomial-degree list");
        std::vector<std::pair<std::string,int>> models={{"Simulation",-1}};
        for(int d:degrees)models.emplace_back("Pol"+std::to_string(d),d);
        std::vector<Result> results;
        for(const auto& m:models){
            FitProblem prob(h,source,m.second,floatShifts);
            std::cout<<"[fit] "<<m.first<<" with "<<prob.nParams<<" parameters and "<<nb<<" bins"<<std::endl;
            auto r=fitOne(m.first,prob,cuts,shiftLimitX,shiftLimitY);
            std::cout<<std::setprecision(5)<<"     converged="<<r.success
              <<" A="<<r.pars[0]<<" B="<<r.pars[1]<<" C="<<r.pars[2]
              <<" chi2/ndf="<<r.reduced<<" dev/ndf="<<r.devRed
              <<" f_inel(ellipse)="<<(r.dataEll>0?r.bgEll/r.dataEll:-1.)<<std::endl;
            if(!r.success)std::cerr<<"[warning] Minuit may not have fully converged: "<<m.first<<std::endl;
            results.emplace_back(std::move(r));
        }
        TFile output((std::string(outdir)+"/fit_histograms.root").c_str(),"RECREATE");
        if(output.IsZombie())throw std::runtime_error("Cannot create ROOT output file");
        for(const auto& r:results)plotOne(h,r,outdir,output);
        compare(results,outdir);output.Close();
        std::ofstream f(std::string(outdir)+"/fit_summary.csv");
        f<<"model,converged,A,B,C,n_parameters,ndf,pearson_chi2,reduced_chi2,deviance,deviance_per_ndf,p_value_approx,AIC,inelastic_full,inelastic_ellipse,data_ellipse,bg_frac_ellipse,inelastic_neutron_dx,data_neutron_dx,bg_frac_neutron_dx\n";
        f<<std::setprecision(12);
        for(const auto& r:results){
            f<<r.name<<","<<r.success<<","<<r.pars[0]<<","<<r.pars[1]<<","<<r.pars[2]<<","<<r.npars<<","<<r.ndf<<","<<r.chi2<<","<<r.reduced<<","<<r.dev<<","<<r.devRed<<","<<r.pvalue<<","<<r.aic<<","<<sum(r.I)<<","<<r.bgEll<<","<<r.dataEll<<","<<(r.dataEll>0?r.bgEll/r.dataEll:std::numeric_limits<double>::quiet_NaN())<<","<<r.bgDx<<","<<r.dataDx<<","<<(r.dataDx>0?r.bgDx/r.dataDx:std::numeric_limits<double>::quiet_NaN())<<"\n";
        }
        std::ofstream pars(std::string(outdir)+"/fit_parameters.csv");
        pars<<"model,parameter_index,value\n";pars<<std::setprecision(12);
        for(const auto& r:results)for(size_t i=0;i<r.pars.size();++i)
            pars<<r.name<<","<<i<<","<<r.pars[i]<<"\n";
        std::ofstream config(std::string(outdir)+"/configuration.txt");
        config<<"Model: A * (P + B*N + C*I)\n"
              <<"I polynomial: 2D tensor product Bernstein degree n in each axis.\n"
              <<"P,N,I have unit sum over 2D fit grid.\n"
              <<"Shift method: integrated 2D MC CDF on extended source, then normalize on fit grid.\n"
              <<"MC source dx=["<<source.xmin<<","<<source.xmax<<"] dy=["
              <<source.ymin<<","<<source.ymax<<"]\n"
              <<"templateGuardX="<<templateGuardX<<" templateGuardY="<<templateGuardY<<"\n"
              <<"Poisson deviance minimization, all cells, including empty.\n"
              <<"ROOT X=dx; Y=dy.\n"
              <<"Ellipse: (dx/"<<cuts.dx_r<<")^2 + (dy/"<<cuts.dy_r<<")^2 <= 1\n"
              <<"dx window: ["<<cuts.dx_L<<","<<cuts.dx_H<<"] (full dy fit range)\n"
              <<"preShiftP="<<preShiftP<<" preShiftN="<<preShiftN<<" preShiftI="<<preShiftI<<"\n"
              <<"floatShifts="<<floatShifts<<" dxShiftLimit="<<shiftLimitX<<" dyShiftLimit="<<shiftLimitY<<"\n"
              <<"data="<<((dataFile&&*dataFile)?dataFile:kData)<<"\n"
              <<"qe="<<((qeFile&&*qeFile)?qeFile:kQE)<<"\n"
              <<"inel="<<((inelFile&&*inelFile)?inelFile:kInel)<<"\n";
        std::cout<<"[done] "<<outdir<<std::endl;
    }catch(const std::exception& e){std::cerr<<"[error] "<<e.what()<<std::endl;throw;}
}
