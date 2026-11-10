Inelastic dx background comparisons — standalone ROOT C++ macro
================================================================

File: compare_inelastic_dx_models.C

Place the macro in your Final_analysis repository directory and run from
there (so cuts/GEN3.json and DB/run_quality.json resolve correctly).
ROOT 6.x with Minuit2 and a C++17 compiler are required.
No Python or nlohmann/json dependency is needed.

Basic command (the GEN3 data, SIMC QE, and inelastic simulation paths from
your raw_asymmetry invocation are embedded as defaults):

    root -l -b -q 'compare_inelastic_dx_models.C+()'

Models: simulation inelastic template, and positive Bernstein polynomials
of degrees 2,3,4,5.  All fits use the same data and SIMC QE p/n templates.

Use coarser 100-bin histograms, identical padded dy selection for data
and simulations, fixed template shifts, and polynomial degrees 2–6:

    root -l -b -q 'compare_inelastic_dx_models.C+("","","","cuts/GEN3.json","DB/run_quality.json","out_plots/inelastic_dx_background_models_GEN3_cpp",100,true,true,"2,3,4,5,6")'

To turn off run-quality filtering, change the 11th argument to false:

    root -l -b -q 'compare_inelastic_dx_models.C+("","","","cuts/GEN3.json","DB/run_quality.json","out_plots/inelastic_dx_background_models_GEN3_cpp",100,true,true,"2,3,4,5",false)'

To cap entries per file for debugging, supply the next (12th) argument:

    root -l -b -q 'compare_inelastic_dx_models.C+("","","","cuts/GEN3.json","DB/run_quality.json","out_plots/inelastic_dx_background_models_GEN3_cpp",100,true,true,"2,3,4,5",false,200000)'

Macro argument order:

  1 dataFile               "" uses default experimental ROOT file
  2 qeFile                 "" uses default SIMC QE p/n ROOT file
  3 inelFile               "" uses default inelastic ROOT file
  4 cutsFile               cuts/GEN3.json
  5 runQualityFile         DB/run_quality.json
  6 outdir                 out_plots/inelastic_dx_background_models_GEN3_cpp
  7 bins                   200
  8 matchedDy              false (data dy padded ±0.05m, sim tight dy)
                           true  (same padded dy selection for all)
  9 fixShifts              false (fit dx shifts within ±0.10m)
                           true  (all dx shifts fixed at zero)
 10 degreeList             "2,3,4,5"
 11 applyRunQuality        true
 12 maxEntriesPerFile      -1 (all entries)
 13 xmin                   -2.2 m
 14 xmax                    0.8 m
 15 shiftLimit              0.10 m
 16 fitxmin                 0.0 (disabled unless fitxmax>fitxmin)
 17 fitxmax                 0.0

Expected input trees: Tout or T. Required flat leaves:
  Data: dx dy W2 vz eHCAL ePS eSH trP runnum ntrack coin_time helicity
  SIMC QE: dx dy W2 vz eHCAL ePS eSH trP weight fnucl
  Inelastic: dx dy W2 vz eHCAL ePS eSH trP weight
Run-quality JSON is expected to contain run-number keys and nested
"helicity" / "moller" flags as true/false or 1/0, as in the Python script.
If the repository uses a different format, disable this flag or adapt the
readGoodRuns() parser.

Output directory (default):

    out_plots/inelastic_dx_background_models_GEN3_cpp/

Files:
    fit_simulation.pdf/.png
    fit_pol2.pdf/.png ... fit_pol5.pdf/.png
       individual fit: upper data/model/components, middle Data-Fit,
       lower (Data-Fit)/sqrt(model)
    all_models_overlay.pdf/.png
    inelastic_background_overlay.pdf/.png
    reduced_chi2_comparison.pdf/.png
    fit_summary.csv
       Pearson chi2, NDF, chi2/NDF, approximate p-value, Poisson
       deviance, AIC, fitted component yields, and background fraction
       inside the neutron dx_L..dx_H window
    fit_parameters.csv
       fitted coefficients and shifts
    fit_histograms.root
       experimental data, fitted total and components, residuals, pulls
    run_configuration.txt

All models minimize the Poisson deviance using unweighted DATA counts:
  2 sum_i [mu_i - n_i + n_i log(n_i/mu_i)]
The plotted chi2/ndf is Pearson sum_i (n_i-mu_i)^2/mu_i divided by
(number of fitted bins - number of fitted parameters).
High-degree polynomial fits use more parameters; reduced chi2 and AIC
should be considered together. These fits do not account for finite MC
template-statistics uncertainty or systematic uncertainty in SIMC shapes.

The "PolN" implementation uses nonnegative Bernstein coefficients instead
of unconstrained ROOT TF1 polN coefficients. This is a degree-N polynomial
but constrained to a nonnegative background. The background fraction quoted
in the CSV is for the 1D dx_L..dx_H window, not the full neutron 2D ellipse.

Testing status: Standalone C++ shape-shifting and Bernstein basis code was
compiled and tested with g++17. The ROOT macro was NOT executed against your
GEN3 ROOT files; ROOT and your data mount are unavailable in this container.
