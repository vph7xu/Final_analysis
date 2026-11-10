GEN3 HELICITY-SUMMED 2D dx-dy INELASTIC BACKGROUND MODEL COMPARISON
===================================================================

Files
-----
  compare_inelastic_dxdy_models_ABC.py    standalone Python (uproot, numpy, scipy, matplotlib)
  compare_inelastic_dxdy_models_ABC.C     standalone ROOT C++ macro (ROOT 6, Minuit2)

Run from the Final_analysis repository root so cuts/GEN3.json and DB/run_quality.json resolve.
The GEN3 experimental, QE SIMC, and inelastic ROOT paths from the previous raw_asymmetry
command are built into BOTH scripts. Use arguments if the paths change.

Default fits
------------
  M(dx,dy) = A * [ P(dx,dy) + B*N(dx,dy) + C*I(dx,dy) ]

  P: unit-normalized SIMC QE proton 2D histogram, fnucl=1.
  N: unit-normalized SIMC QE neutron 2D histogram, fnucl=0.
  I: unit-normalized simulated inelastic histogram OR unit-normalized positive
     bivariate Bernstein polynomial.

  A = expected QE proton count in the full 2D fit domain.
  B = expected QE neutron count divided by proton count.
  C = expected inelastic count divided by proton count.
  Predicted totals: P=A, N=A*B, I=A*C (within the fit domain).

  Simulation: 3 parameters (A,B,C), or 9 with all dx/dy shifts floated.
  Pol2: 2D Bernstein degree (2,2), 9 positive weights (8 independent)
        + A,B,C = 11 parameters (15 with P,N shifts).
  Pol3: 16 weights (15 independent) + A,B,C = 18 parameters.
  Pol4: 25 weights (24 independent) + A,B,C = 27 parameters.
  Pol5: 36 weights (35 independent) + A,B,C = 38 parameters.

  IMPORTANT: PolN is tensor-product degree N in *each* coordinate (bidegree
  (N,N)), NOT an unrestricted total-degree-N polynomial. All polynomial
  basis columns have unit bin sum over the full fitted dx-dy grid. Weights
  are >=0 (implemented with log ratios and one reference weight fixed).
  This is the two-dimensional counterpart of the previous 1D Bernstein
  polynomial fits and guarantees that I(dx,dy) is nonnegative.

The fit minimizes the Poisson deviance in all selected bins, including empty
bins. A Pearson chi-square divided by (N_dx*N_dy - N_free_parameters) is
reported for comparison; p-values are approximate if many bins are sparse.
AIC = deviance + 2*N_free_parameters also reported.

Python
------
  python3 -m pip install numpy scipy matplotlib uproot
  python3 compare_inelastic_dxdy_models_ABC.py

  Example with a fitted 2D offset for P,N and simulated I:
    python3 compare_inelastic_dxdy_models_ABC.py --float-shifts --bins-dx 50 --bins-dy 36

  To increase MC support further (without changing the fit interval):
    python3 compare_inelastic_dxdy_models_ABC.py --template-guard-dx 0.60 --template-guard-dy 0.20

  Example check fewer polynomials before including Pol4/Pol5:
    python3 compare_inelastic_dxdy_models_ABC.py --degrees 2 3 --bins-dx 50 --bins-dy 36

  To apply the manual shifts from the currently uploaded GEN3 2D filling:
    python3 compare_inelastic_dxdy_models_ABC.py --pre-shift-p-dx 0.05 --pre-shift-i-dx -0.10

  By default, manual shifts are ZERO; there are NO hidden hard-coded GEN3
  pre-shifts. This makes the different background-shape models comparable.

C++ ROOT
--------
  root -l -b -q 'compare_inelastic_dxdy_models_ABC.C+()'

  Example with floated offsets, 50 by 36 bins:
    root -l -b -q 'compare_inelastic_dxdy_models_ABC.C+("","","","cuts/GEN3.json","DB/run_quality.json","out_plots/inelastic_dxdy_background_models_GEN3_ABC_cpp",50,36,true)'

  The extra support guard parameters are optional final arguments to the
  ROOT entry function: templateGuardX=0.45 and templateGuardY=0.15 [m].
  The Python counterpart uses --template-guard-dx/--template-guard-dy.

  The C++ macro's top-level function is named exactly like its file, so the
  zero-argument ROOT +() syntax works. Change the default file paths at top
  of the macro or provide absolute paths as its first 3 arguments.

Selections and coordinates
--------------------------
  Read Tout (or T) trees from the original data and simulation files.
  Data: good run quality (unless explicitly disabled), abs(helicity)=1,
        ntrack>=1, coin_time and run range, abs(vz)<0.27, ePS>0.2,
        HCAL energy, E/P and W2 from cuts/GEN3.json.
  Simulation: the corresponding calorimeter, vertex, E/P and W2 cuts,
              weighted by the tree's weight branch.
  DATA retain the ORIGINAL fit dy window expanded by 0.05 m by default:
     dy_L - 0.05 < dy < dy_H + 0.05.
  DATA dx fit range stays -2.2 to +0.8 m by default.
  SIMULATION is now filled into a WIDER 2D support grid in both coordinates
  before shifting: the default minimum padding is 0.45 m on each dx side,
  and 0.15 m on each dy side. Padding is enlarged automatically as needed
  for allowed template shifts plus 3 extra original-width bins.
  NOTE: the additional MC dy range is an interpolation reservoir, not an
  additional fit region. This replaces the former preliminary MC dy cut.
  Each shifted 2D MC template is integrated into the UNCHANGED data fit bins
  using a bilinear 2D cumulative distribution (piecewise-constant MC bin
  density), then normalized to unit sum over the fit domain. No MC events
  outside the fit domain are counted in A, AB or AC unless translated into it.
  That prevents the shifted proton template being forcibly zero at the plot
  boundary simply because the source MC histogram ended there.
  Python arrays are indexed [dx_bin, dy_bin]. The C++ ROOT histograms
  have ROOT X=dx and ROOT Y=dy (the ORIGINAL analysis plots them swapped).
  All 2D templates are normalized over the SAME FIT dx-dy histogram limits.
  Source support is wider for the MC only; the data and fitting bins are NOT
  enlarged, and the neutron ellipse definition is unchanged.

Output
------
  Python default output: out_plots/inelastic_dxdy_background_models_GEN3_ABC_py/
  ROOT default output:   out_plots/inelastic_dxdy_background_models_GEN3_ABC_cpp/

  fit_summary.csv                  A, B, C, reduced chi2, deviance/ndf, AIC,
                                   fitted inelastic counts and ROI fractions.
  fit_parameters.csv               Full fitted parameter list.
  fit2d_<Model>.pdf/.png           Data, model, 2D residual, 2D pull heatmaps.
  projections_<Model>.pdf/.png     dx/dy fits with QE p/n and background,
                                   plus projection pull panels.
  reduced_chi2_comparison.pdf      Reduced chi2 for all models.
  Python: background_fraction_ellipse.pdf, fit_arrays.npz, configuration.json.
  C++: fit_histograms.root, configuration.txt.

ROI / fitted contamination
--------------------------
  Neutron ellipse is CENTERED AT dx=0, dy=0:
       (dx/dx_r)^2 + (dy/dy_r)^2 <= 1,
  using dx_r and dy_r from cuts/GEN3.json, and selecting bins whose centers
  fall inside it. The reported bg_frac_ellipse is the fitted inelastic count
  sum inside the ellipse divided by observed data counts inside the ellipse.

  A second fraction integrates only the configured neutron dx interval
  dx_L <= dx <= dx_H over the entire dy fit range.

  These are not equivalent selections. Additional cuts used later in the
  full asymmetry analysis must be included separately if applicable.

Cautions
--------
  - Polynomial Pol5 has many more shape parameters than the simulated
    template; chi2/ndf alone is not enough to prefer it. Compare AIC and
    the extracted inelastic fraction in the neutron ellipse.
  - The polynomial fit is a flexible empirical background; it need not
    reproduce physical dx-dy correlations from the inelastic simulation.
  - Poisson likelihood assumes unweighted observed DATA; simulation
    statistical fluctuations and generator-weight uncertainty are NOT
    profiled here. Consider template-statistics nuisance parameters
    if MC templates are sparse.
  - As with the prior 1D scripts, run_quality.json is assumed to contain
    per-run {"helicity": true, "moller": true} flags; if your schema differs,
    adjust the parser rather than silently skipping data quality checks.
  - The updated 2D cumulative-shift algorithm was tested against analytic
    Gaussian bin integrals and matched between Python and separately compiled
    C++ code at floating-point precision. Synthetic-data Python fit closure
    and projection plotting were also exercised. The FULL ROOT macro still
    could not be compiled/run because ROOT is not installed here.
  - The Python fitter internally rescales A by the observed event count for
    numerical optimization, but reports the usual physical count A; the
    model is exactly A*(P+B*N+C*I) with all shapes normalized in fit bins.
