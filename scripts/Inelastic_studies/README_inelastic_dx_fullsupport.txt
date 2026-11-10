Full-support update for 1D GEN3 ROOT background-model comparison
==================================================================

ROOT macro: compare_inelastic_dx_models_ABC_fullsupport.C

Run from Final_analysis repository root:

    root -l -b -q 'compare_inelastic_dx_models_ABC_fullsupport.C+()'

Reason for this revision
------------------------
The previous macro filled the simulation templates only over the same dx
range used for fitting/plotting (default [-2.2, 0.8] m). When a proton template
is shifted by +delta, its leftmost fitted bins request the unshifted template
at x-delta, outside the histogram. The CDF returns zero outside the histogram,
producing a sharp, artificial step near the left plot boundary. The step can
change the minimized likelihood and the reduced chi-square.

This version separates the two ranges:

  * DATA FIT/PLOT RANGE: default [-2.2, 0.8] m, 200 bins.
  * MC TEMPLATE SOURCE RANGE: extends both sides by at least 0.45 m, and by
    at least (max absolute allowed shift + three fit bins).

The data and all three simulated components are histogrammed on this extended
support. Only the requested data fit interval is passed to the optimizer,
while the component template interpolator retains the entire extended MC
histogram. For each bin in the fit window the shifted component is obtained by
integrating the unshifted MC template over the shifted bin edges:

    P_shifted_bin_i = CDF_P(x_right_i - delta_p)
                      - CDF_P(x_left_i - delta_p)

Likewise for N and, in the simulated background model, I.
Each shifted component is normalized over the *fit* region, ensuring that

    M_i = A * [ P_i(dx-dp) + B*N_i(dx-dn) + C*I_i(dx-di) ]

and the component yields in that region are A, A*B, and A*C.
For polynomial models I_i is the same normalized nonnegative Bernstein
polynomial as before (not dx-shifted separately).

The wider source range does NOT add extra event counts to the likelihood:
fit data and fit bins remain the same by default. Only the MC template
support is wider, preventing a hard edge inside the fit window.

OPTIONAL: Increase the dx fit range as well, e.g. to [-2.5, 0.8] m, by
passing xmin=-2.5 and xmax=0.8 in the full macro argument list. Choose
nbins to preserve the desired bin width, and keep enough source guard on
each side.

Optional argument "templateGuard" is appended at the end of the macro
signature (default 0.45 m). For example, setting it to 0.8 m can test
sensitivity to the source-histogram boundaries. The effective guard may be
increased automatically to accommodate shiftLimit.

Outputs are unchanged (fit_*.png, fit_*.pdf, fit_summary.csv,
fit_parameters.csv, fit_histograms.root), but the new default output directory
is out_plots/inelastic_dx_background_models_GEN3_ABC_fullsupport/
so the earlier fits are not overwritten.

Important
---------
This removes a *histogram-boundary artifact*. It will not fix a real mismatch
between QE simulation and data, or Monte Carlo statistics fluctuations. If the
leftmost plots still show an abrupt cutoff, inspect the original source dx
histogram and verify events extend across the full source range.

ROOT is not installed in this execution environment, so the complete macro
cannot be compiled and run here; the isolated CDF/bin-shifting logic has been
tested with C++ test distributions.
