Inelastic dx-background model comparison — GEn-II GEN3
=======================================================

Run from the Final_analysis repository root (where cuts/GEN3.json and DB/run_quality.json live):

    python3 -m pip install --user numpy scipy matplotlib uproot
    python3 compare_inelastic_dx_models.py

The script has your exact GEN3 data, SIMC QE, and inelastic simulation file paths
embedded as defaults, so no path arguments are required for the original files.
Use --data, --qe, --inel to override any of them. The pion, N2, and G4SBS
elastic files from raw_asymmetry are NOT needed for this particular dx fit.

Examples
--------

    # Original 1D selection convention from InelasticCorrection.cpp:
    python3 compare_inelastic_dx_models.py

    # Recommended control: identical padded dy range for data and simulation:
    python3 compare_inelastic_dx_models.py --selection-mode matched

    # Compare coarser binning and add higher polynomial orders:
    python3 compare_inelastic_dx_models.py --bins 100 --degrees 2 3 4 5 6 7

    # Remove fitted template shifts to isolate background models:
    python3 compare_inelastic_dx_models.py --fix-shifts

    # Restrict all fits to the same x region:
    python3 compare_inelastic_dx_models.py --fit-range -2.1 0.7

    # Quick input/branch check with a cap per ROOT tree:
    python3 compare_inelastic_dx_models.py --max-events 200000

Output directory (by default)
-----------------------------

    out_plots/inelastic_dx_background_models_GEN3/

Files:
  fit_simulation.pdf / .png         simulated inelastic component, 2 panels
  fit_pol2.pdf / .png               degree-2 Bernstein background
  fit_pol3.pdf / .png               degree-3 Bernstein background
  fit_pol4.pdf / .png               degree-4 Bernstein background
  fit_pol5.pdf / .png               degree-5 Bernstein background
  all_models_overlay.pdf / .png     totals and pull overlay
  reduced_chi2_comparison.pdf       side-by-side reduced chi-square
  inelastic_background_overlay.pdf compare the inferred backgrounds
  fit_summary.csv                   chi-square/ndf, p-value, AIC,
                                    background yield/fraction in neutron dx cut
  fit_parameters.json               best-fit parameters
  fit_histograms.npz                full arrays (data/templates/models/pulls)
  run_configuration.json            command options and exact JSON cut values

Statistics
----------
All fits minimize the Poisson deviance in unweighted DATA counts. The
quoted chi2 is Pearson chi2 = sum((data - model)^2 / model); degrees of
freedom = number of fitted bins - number of fitted parameters. "Pull"
is (data - model) / sqrt(model). This does not include finite weighted-MC
uncertainty. p-values are approximate, especially for sparse bins/boundary
parameters. Compare AIC as well as reduced chi2 because high-order
backgrounds have additional parameters.

"PolN" means a positive Bernstein polynomial of degree N, NOT an
unconstrained ROOT TF1 polN. It is still a degree-N polynomial, but all
coefficients are constrained nonnegative to avoid negative backgrounds.

By default, "current" selection mode reproduces the existing 1D convention:
the data dx histogram has dy padded by 0.05 m, while the 1D simulation
templates have the tight dy cut. "matched" mode uses the SAME padded dy
window for data and simulations. Other cuts mimic InelasticCorrection.cpp.
The simulated template shifts are performed smoothly using bin-integrated
histogram shifts; the fitting statistic is intentionally improved over the
original unweighted or histogram chi-square approach.

No experimental fit has been run here because the ROOT files are on your
local /media/vimukthi mount, not available in this environment. The fitting
and plotting functions were tested on synthetic count histograms.
