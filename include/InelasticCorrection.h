// -----------------------------------------------------------------------------
// InelasticCorrection.h  —  correction factor for inelastic background
// -----------------------------------------------------------------------------
// * Idea: compute f_inel = N_inel / (N_QE + N_inel) using a model template
//   taken from simulation, then propagate its uncertainty to the physics
//   asymmetry.
// * Region-of-interest (inelastic peak) is defined by AnalysisCuts.W2.
//   For example, W2_H_inel and W2_L_inel thresholds can be added to GEN3.json.
// * The class works exactly like PionCorrection: give it data chain, a QE
//   simulation chain, and an inelastic simulation chain.
// -----------------------------------------------------------------------------
#ifndef INELASTIC_CORRECTION_H
#define INELASTIC_CORRECTION_H

#include "BranchVars.h"
#include "BranchVarsSim.h"
#include "AnalysisCuts.h"
#include "RunQuality.h"

#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <string>
#include <functional>

class InelasticCorrection {
public:
    InelasticCorrection(const AnalysisCuts& cuts,
                        const RunQuality*   rq   = nullptr,
                        const char*         kin  = "GEN3_He3",
                        const char*         rootFile = "inel_correction.root");

    // process chains (data, QE sim, inelastic sim)
    void process(TChain&        ch_data,
                 TChain&        ch_QE_sim,
                 TChain&        ch_inel_sim,
                 BranchVars&    vData,
                 BranchVarsSim& vQE,
                 BranchVarsSim& vInel);

    // accessors
    double fInel()  const { return frac_; }
    double dfInel() const { return dfrac_; }

private:
    const AnalysisCuts& c_;
    const RunQuality*   rq_;
    const char*         kin_;
    std::string         outFile_;

    // fractions
    double frac_  = 0.0;   // N_inel / (N_inel + N_QE)
    double dfrac_ = 0.0;

    // helpers
    TH1D* performFit(TH1D* h_data, TH1D* h_inel, TH1D* h_QE_proton, TH1D* h_QE_neutron,
                     double& par0, double& par1, double& par2, double& dx_p_out, double& dx_n_out, double& dx_inel_out);
                    
    TH2D* performFitDxDy(TH2D* hD2, TH2D* hInel2, TH2D* hQEp2, TH2D* hQEn2,
	    double& A, double& rNP, double& rI,
	    double& dx_p, double& dy_p,double& dx_n, double& dy_n,double& dx_i, double& dy_i);

    TH3D* performFitDxDyW2(TH3D* hD3, TH3D* hInel3, TH3D* hQEp3, TH3D* hQEn3,
    			    double& A, double& rNP, double& rI,
    			    double& dx_p, double& dy_p, double& dW2_p, double& dx_n, double& dy_n, double& dW2_n,
    			    double& dx_i, double& dy_i, double& dW2_i);
                     
    TH1D* performFitW2(TH1D* h_data, TH1D* h_inel, TH1D* h_QE_proton, TH1D* h_QE_neutron,
                     double& par0, double& par1, double Rnop);
                     
    TH1D* performFitW2_1(TH1D* hD, TH1D* hInel, TH1D* hQEp, TH1D* hQEn,
		      double& alpha, double& delta , double Rnop);
		      
    TH1D* performFitW2_2(TH1D* hD, TH1D* hInel, TH1D* hQE, 
    		      double& par0, double& par1, double& delta);

};

#endif // INELASTIC_CORRECTION_H

