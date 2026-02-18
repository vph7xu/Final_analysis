#!/usr/bin/env python3
# overlay_eHCAL_and_ratio_multipage_with_gausfit_profile.py
#
# One loop over data.
# Produces TWO multipage PDFs:
#   (1) eHCAL: 1D overlay + 2D vs y/x with Gaussian-fit mean overlay
#   (2) R = eHCAL*2*0.938/Q2: same plots, separate PDF
#
# PLUS: per-bin fit PDFs (4x4 pads/page) that show ProjectionY + Gaussian fit
#       for each x-bin, for each 2D histogram.
#
# ROOT 6.26-friendly: manual loop + ASCII progress bar

import argparse
import sys
import time
import math
import ROOT

MP = 0.938  # GeV


# ------------------------- styling helpers -------------------------
def style_hist1d(h, color, width=3):
    h.SetLineColor(color)
    h.SetLineWidth(width)
    h.SetFillStyle(0)

def setup_pad(mleft=0.12, mright=0.04, mbot=0.12, mtop=0.08):
    ROOT.gPad.SetLeftMargin(mleft)
    ROOT.gPad.SetRightMargin(mright)
    ROOT.gPad.SetBottomMargin(mbot)
    ROOT.gPad.SetTopMargin(mtop)

def label_ndc(text, x=0.14, y=0.93, size=0.045, color=ROOT.kBlack):
    tx = ROOT.TLatex()
    tx.SetNDC(True)
    tx.SetTextFont(42)
    tx.SetTextSize(size)
    tx.SetTextColor(color)
    tx.DrawLatex(x, y, text)

def render_bar(frac, width=40):
    frac = max(0.0, min(1.0, frac))
    nfill = int(frac * width)
    return "[" + "#" * nfill + "-" * (width - nfill) + f"] {100.0*frac:6.2f}%"

def style_graph_mean(g, color=ROOT.kBlack, mstyle=20, msize=1.6, lwidth=2):
    g.SetMarkerStyle(mstyle)
    g.SetMarkerSize(msize)
    g.SetMarkerColor(color)
    g.SetLineColor(color)
    g.SetLineWidth(lwidth)

def style_graph_sigma(g, color=ROOT.kBlack, mstyle=24, msize=1.4, lwidth=2):
    g.SetMarkerStyle(mstyle)
    g.SetMarkerSize(msize)
    g.SetMarkerColor(color)
    g.SetLineColor(color)
    g.SetLineWidth(lwidth)


# ------------------------- gaussian fit per 1D -------------------------
def fit_gaussian_robust(h, fit_nsigma=2.0, min_sigma=1e-6):
    """
    Fit a Gaussian to TH1 'h' in a dynamic range around its mean.
    Returns:
      (ok, mu, mu_err, sig, sig_err, chi2, ndf, f_gaus)
    """
    entries = h.GetEntries()
    if entries <= 0:
        return (False, float("nan"), float("nan"), float("nan"), float("nan"), float("nan"), 0, None)

    mu0 = h.GetMean()
    sig0 = h.GetRMS()
    if not math.isfinite(mu0) or not math.isfinite(sig0) or sig0 < min_sigma:
        # fallback: try to estimate from bin content
        sig0 = max(min_sigma, 0.1 * (h.GetXaxis().GetXmax() - h.GetXaxis().GetXmin()))

    xmin = h.GetXaxis().GetXmin()
    xmax = h.GetXaxis().GetXmax()
    lo = max(xmin, mu0 - fit_nsigma * sig0)
    hi = min(xmax, mu0 + fit_nsigma * sig0)
    if not (hi > lo):
        lo, hi = xmin, xmax

    # Create a unique TF1 name to avoid collisions
    fname = f"fgaus_{h.GetName()}"
    fgaus = ROOT.TF1(fname, "gaus", lo, hi)

    # Reasonable initial params: amplitude ~ max bin, mean ~ mu0, sigma ~ sig0
    fgaus.SetParameters(h.GetMaximum(), mu0, sig0)

    # Quiet fit, no drawing; "R" uses range; "Q" quiet; "0" no draw
    # "S" returns fit result
    fitopt = "Q0RS"
    r = h.Fit(fgaus, fitopt)

    # In PyROOT, FitResult can be None depending on build; still read TF1 params
    mu = fgaus.GetParameter(1)
    sig = abs(fgaus.GetParameter(2))
    mu_err = fgaus.GetParError(1)
    sig_err = fgaus.GetParError(2)
    chi2 = fgaus.GetChisquare()
    ndf = fgaus.GetNDF()

    ok = (math.isfinite(mu) and math.isfinite(sig) and sig > 0.0)
    return (ok, mu, mu_err, sig, sig_err, chi2, ndf, fgaus)


# ------------------------- build gaussian-fit "profiles" from TH2 -------------------------
def gausfit_profile_from_th2(
    h2,
    fit_nsigma=2.0,
    min_entries_bin=50,
    max_bins=None,
    make_fit_pages_pdf=None,
    fit_pages_title="",
    pads_nx=4,
    pads_ny=4
):
    """
    For each X-bin of TH2 (X is the histogram x-axis), project Y distribution and fit Gaussian.
    Returns:
      (g_mean, g_sigma, fit_summaries)
    where
      g_mean  : TGraphErrors of mu(x) with mu_err
      g_sigma : TGraphErrors of sigma(x) with sigma_err
      fit_summaries: list of dict per bin (xbin, xcenter, entries, ok, mu, sig, etc.)
    Optionally writes a multipage PDF of the per-bin fits using 4x4 pads/page if make_fit_pages_pdf is set.
    """
    nbx = h2.GetNbinsX()
    if max_bins is not None:
        nbx = min(nbx, int(max_bins))

    g_mean = ROOT.TGraphErrors()
    g_sigma = ROOT.TGraphErrors()
    g_mean.SetName(f"gmean_{h2.GetName()}")
    g_sigma.SetName(f"gsig_{h2.GetName()}")

    summaries = []

    # Optional: multipage PDF for bin fits
    cfit = None
    per_page = pads_nx * pads_ny
    pad_index = 0
    page_open = False
    if make_fit_pages_pdf:
        cfit = ROOT.TCanvas(f"cfit_{h2.GetName()}", "bin fits", 1400, 1000)
        cfit.Print(make_fit_pages_pdf + "[")
        page_open = True

    # Iterate X bins (1..nbx, excluding under/overflow)
    for ibx in range(1, nbx + 1):
        # Project Y for this X bin
        proj_name = f"py_{h2.GetName()}_{ibx}"
        py = h2.ProjectionY(proj_name, ibx, ibx)
        py.SetDirectory(0)

        ent = int(py.GetEntries())
        xcenter = h2.GetXaxis().GetBinCenter(ibx)
        xlow = h2.GetXaxis().GetBinLowEdge(ibx)
        xhigh = h2.GetXaxis().GetBinUpEdge(ibx)

        info = {
            "xbin": ibx,
            "xcenter": xcenter,
            "xlow": xlow,
            "xhigh": xhigh,
            "entries": ent,
            "ok": False,
            "mu": float("nan"),
            "mu_err": float("nan"),
            "sig": float("nan"),
            "sig_err": float("nan"),
            "chi2": float("nan"),
            "ndf": 0
        }

        fgaus = None
        if ent >= min_entries_bin:
            ok, mu, mu_err, sig, sig_err, chi2, ndf, fgaus = fit_gaussian_robust(py, fit_nsigma=fit_nsigma)
            info.update({
                "ok": ok,
                "mu": mu,
                "mu_err": mu_err,
                "sig": sig,
                "sig_err": sig_err,
                "chi2": chi2,
                "ndf": ndf
            })
            if ok:
                ip = g_mean.GetN()
                g_mean.SetPoint(ip, xcenter, mu)
                g_mean.SetPointError(ip, 0.0, mu_err)

                ip2 = g_sigma.GetN()
                g_sigma.SetPoint(ip2, xcenter, sig)
                g_sigma.SetPointError(ip2, 0.0, sig_err)

        summaries.append(info)

        # Draw per-bin fit pages if requested
        if make_fit_pages_pdf and cfit:
            if (pad_index % per_page) == 0:
                cfit.Clear()
                cfit.Divide(pads_nx, pads_ny, 0.001, 0.001)
                pad_index = 0

            pad_index += 1
            cfit.cd(pad_index)
            ROOT.gPad.SetLeftMargin(0.14)
            ROOT.gPad.SetRightMargin(0.05)
            ROOT.gPad.SetBottomMargin(0.14)
            ROOT.gPad.SetTopMargin(0.10)

            py.SetMarkerStyle(20)
            py.SetMarkerSize(0.8)
            py.SetLineWidth(2)
            py.SetTitle("")

            py.GetXaxis().SetTitle("Y")
            py.GetYaxis().SetTitle("Counts")
            py.GetXaxis().SetTitleSize(0.07)
            py.GetYaxis().SetTitleSize(0.07)
            py.GetXaxis().SetLabelSize(0.06)
            py.GetYaxis().SetLabelSize(0.06)

            py.Draw("E1")

            if fgaus is not None and info["ok"]:
                fgaus.SetLineWidth(2)
                fgaus.Draw("SAME")

            # small label
            tx = ROOT.TLatex()
            tx.SetNDC(True)
            tx.SetTextFont(42)
            tx.SetTextSize(0.07)
            tx.DrawLatex(0.12, 0.92, f"x[{ibx}]={xlow:.3g}..{xhigh:.3g}")
            tx.SetTextSize(0.06)
            tx.DrawLatex(0.12, 0.84, f"N={ent}")
            if info["ok"]:
                tx.DrawLatex(0.12, 0.76, f"#mu={info['mu']:.3g}  #sigma={info['sig']:.3g}")
            else:
                tx.DrawLatex(0.12, 0.76, "fit: (skip/fail)")

            # If last pad of page or last bin, print page
            if (pad_index % per_page) == 0 or (ibx == nbx):
                # Add a page title on the canvas itself (in the first pad area)
                cfit.cd(1)
                ttop = ROOT.TLatex()
                ttop.SetNDC(True)
                ttop.SetTextFont(42)
                ttop.SetTextSize(0.035)
                ttop.DrawLatex(0.02, 0.98, fit_pages_title)
                cfit.Print(make_fit_pages_pdf)

        # cleanup
        # (py will be garbage-collected; SetDirectory(0) prevents ROOT ownership issues)

    if make_fit_pages_pdf and cfit and page_open:
        cfit.Print(make_fit_pages_pdf + "]")

    return g_mean, g_sigma, summaries


# ------------------------- main multipage PDFs -------------------------
def make_main_multipage_pdf(
    outpdf,
    var_label,
    xlab_1d, ylab_1d,
    xlab_y2d, xlab_x2d,
    ylab_2d,
    labels,
    counts,
    h1_list, h2y_list, h2x_list,
    gmean_y_list, gmean_x_list,
    gsigma_y_list=None, gsigma_x_list=None,
    prof_color=ROOT.kBlack
):
    """
    Main PDF: 1D overlay + 2D pages with Gaussian-fit mean overlay.
    Also adds summary pages for mean and sigma graphs (one page per 2D hist).
    """
    keep = []  # keep graphs alive

    c = ROOT.TCanvas("c_main", "c_main", 1000, 800)
    c.Print(outpdf + "[")

    # Page 1: 1D overlay
    c.Clear()
    setup_pad(mright=0.04)
    h1_list[0].GetXaxis().SetTitle(xlab_1d)
    h1_list[0].GetYaxis().SetTitle(ylab_1d)
    h1_list[0].Draw("hist")
    for hh in h1_list[1:]:
        hh.Draw("hist same")

    leg = ROOT.TLegend(0.55, 0.68, 0.94, 0.92)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for hh, lab, n in zip(h1_list, labels, counts):
        leg.AddEntry(hh, f"{lab} (N={n})", "l")
    leg.Draw()

    label_ndc(f"{var_label}", 0.14, 0.93, 0.045)
    c.Print(outpdf)

    # Pages 2-4: 2D vs y with gaussian-mean overlay
    for j in range(3):
        c.Clear()
        setup_pad(mright=0.14)
        h2y_list[j].GetXaxis().SetTitle(xlab_y2d)
        h2y_list[j].GetYaxis().SetTitle(ylab_2d)
        h2y_list[j].Draw("COLZ")

        g = gmean_y_list[j]
        style_graph_mean(g, color=prof_color, mstyle=20, msize=1.8, lwidth=2)
        g.Draw("SAME PE1")
        keep.append(g)

        label_ndc(f"{var_label} vs y | {labels[j]}  (Gaussian-fit #mu per x-bin)", 0.14, 0.93, 0.040)
        c.Print(outpdf)

    # Pages 5-7: 2D vs x with gaussian-mean overlay
    for j in range(3):
        c.Clear()
        setup_pad(mright=0.14)
        h2x_list[j].GetXaxis().SetTitle(xlab_x2d)
        h2x_list[j].GetYaxis().SetTitle(ylab_2d)
        h2x_list[j].Draw("COLZ")

        g = gmean_x_list[j]
        style_graph_mean(g, color=prof_color, mstyle=20, msize=1.8, lwidth=2)
        g.Draw("SAME PE1")
        keep.append(g)

        label_ndc(f"{var_label} vs x | {labels[j]}  (Gaussian-fit #mu per x-bin)", 0.14, 0.93, 0.040)
        c.Print(outpdf)

    # Optional: add summary pages (mean & sigma vs axis)
    # (One page per cut per axis type)
    if gsigma_y_list is not None and gsigma_x_list is not None:
        for j in range(3):
            c.Clear()
            c.Divide(1, 2, 0.001, 0.001)

            # top: mean vs y
            c.cd(1)
            setup_pad(mright=0.05)
            gmu = gmean_y_list[j]
            style_graph_mean(gmu, color=prof_color, mstyle=20, msize=1.6, lwidth=2)
            gmu.SetTitle("")
            gmu.GetXaxis().SetTitle(xlab_y2d)
            gmu.GetYaxis().SetTitle("#mu (Gaussian fit)")
            gmu.Draw("AP")
            label_ndc(f"{var_label} vs y | {labels[j]} : #mu(x)", 0.14, 0.92, 0.05)

            # bottom: sigma vs y
            c.cd(2)
            setup_pad(mright=0.05)
            gs = gsigma_y_list[j]
            style_graph_sigma(gs, color=prof_color, mstyle=24, msize=1.4, lwidth=2)
            gs.SetTitle("")
            gs.GetXaxis().SetTitle(xlab_y2d)
            gs.GetYaxis().SetTitle("#sigma (Gaussian fit)")
            gs.Draw("AP")
            label_ndc(f"{var_label} vs y | {labels[j]} : #sigma(x)", 0.14, 0.92, 0.05)

            c.Print(outpdf)

        for j in range(3):
            c.Clear()
            c.Divide(1, 2, 0.001, 0.001)

            # top: mean vs x
            c.cd(1)
            setup_pad(mright=0.05)
            gmu = gmean_x_list[j]
            style_graph_mean(gmu, color=prof_color, mstyle=20, msize=1.6, lwidth=2)
            gmu.SetTitle("")
            gmu.GetXaxis().SetTitle(xlab_x2d)
            gmu.GetYaxis().SetTitle("#mu (Gaussian fit)")
            gmu.Draw("AP")
            label_ndc(f"{var_label} vs x | {labels[j]} : #mu(x)", 0.14, 0.92, 0.05)

            # bottom: sigma vs x
            c.cd(2)
            setup_pad(mright=0.05)
            gs = gsigma_x_list[j]
            style_graph_sigma(gs, color=prof_color, mstyle=24, msize=1.4, lwidth=2)
            gs.SetTitle("")
            gs.GetXaxis().SetTitle(xlab_x2d)
            gs.GetYaxis().SetTitle("#sigma (Gaussian fit)")
            gs.Draw("AP")
            label_ndc(f"{var_label} vs x | {labels[j]} : #sigma(x)", 0.14, 0.92, 0.05)

            c.Print(outpdf)

    c.Print(outpdf + "]")


# ------------------------- main -------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="ROOT file")
    ap.add_argument("--tree", default="Tout", help="Tree name (default: Tout)")

    ap.add_argument("--evar", default="eHCAL", help="Energy variable/expression (default: eHCAL)")
    ap.add_argument("--xvar", default="xHCAL", help="HCAL X variable/expression (default: xHCAL)")
    ap.add_argument("--yvar", default="yHCAL", help="HCAL Y variable/expression (default: yHCAL)")
    ap.add_argument("--q2var", default="Q2",   help="Q2 branch/expression (default: Q2)")

    # eHCAL binning
    ap.add_argument("--nbinsE", type=int, default=100)
    ap.add_argument("--emin", type=float, default=0.0)
    ap.add_argument("--emax", type=float, default=1.0)

    # ratio binning
    ap.add_argument("--nbinsR", type=int, default=120)
    ap.add_argument("--rmin", type=float, default=0.0)
    ap.add_argument("--rmax", type=float, default=0.3)

    # y/x ranges for 2D
    ap.add_argument("--nbinsY", type=int, default=100)
    ap.add_argument("--ymin", type=float, default=-1.0)
    ap.add_argument("--ymax", type=float, default=1.0)

    ap.add_argument("--nbinsX", type=int, default=100)
    ap.add_argument("--xmin", type=float, default=-2.7)
    ap.add_argument("--xmax", type=float, default=1.2)

    ap.add_argument("--cut1", required=True)
    ap.add_argument("--cut2", required=True)
    ap.add_argument("--cut3", required=True)
    ap.add_argument("--label1", default="All")
    ap.add_argument("--label2", default="Coin")
    ap.add_argument("--label3", default="Singles")

    ap.add_argument("--norm1d", action="store_true", help="Normalize 1D hists to area=1")

    ap.add_argument("--out_e", default="eHCAL_multipage.pdf", help="Output PDF for eHCAL plots")
    ap.add_argument("--out_r", default="ratio_multipage.pdf", help="Output PDF for ratio plots")

    # Gaussian-fit profile controls
    ap.add_argument("--fit_nsigma", type=float, default=2.0, help="Fit range = mean +/- nsigma*RMS (per x-bin)")
    ap.add_argument("--min_entries_bin", type=int, default=80, help="Min entries in each x-bin projection to fit")
    ap.add_argument("--profile_color", default="white", help="black/white/red/blue")

    # Fit-pages PDFs controls
    ap.add_argument("--make_fit_pdfs", action="store_true", help="Write 4x4-per-page PDFs with per-bin fits")
    ap.add_argument("--pads_nx", type=int, default=4)
    ap.add_argument("--pads_ny", type=int, default=4)

    ap.add_argument("--print_every", type=int, default=20000)
    args = ap.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetNumberContours(80)

    prof_color = {
        "black": ROOT.kBlack,
        "white": ROOT.kWhite,
        "red":   ROOT.kRed+1,
        "blue":  ROOT.kBlue+1,
    }.get(args.profile_color.lower(), ROOT.kWhite)

    f = ROOT.TFile.Open(args.input)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open file: {args.input}")
    t = f.Get(args.tree)
    if not t:
        raise RuntimeError(f"Could not find tree '{args.tree}' in {args.input}")

    cuts   = [args.cut1, args.cut2, args.cut3]
    labels = [args.label1, args.label2, args.label3]
    colors = [ROOT.kBlack, ROOT.kRed+1, ROOT.kBlue+1]

    # Formulas (allow expressions)
    fE  = ROOT.TTreeFormula("fE",  args.evar,  t)
    fX  = ROOT.TTreeFormula("fX",  args.xvar,  t)
    fY  = ROOT.TTreeFormula("fY",  args.yvar,  t)
    fQ2 = ROOT.TTreeFormula("fQ2", args.q2var, t)
    fCuts = [ROOT.TTreeFormula(f"fCut{i+1}", cuts[i], t) for i in range(3)]

    # --- Histograms per cut for eHCAL ---
    h1_e, h2y_e, h2x_e = [], [], []
    # --- Histograms per cut for ratio R ---
    h1_r, h2y_r, h2x_r = [], [], []

    for i in range(3):
        # 1D e
        he = ROOT.TH1D(f"h_e_{i+1}", "", args.nbinsE, args.emin, args.emax)
        he.Sumw2()
        style_hist1d(he, colors[i])
        h1_e.append(he)

        # 2D e vs y/x (X=y/x, Y=e)
        hy_e = ROOT.TH2D(f"h2_e_vs_y_{i+1}", "", args.nbinsY, args.ymin, args.ymax,
                         args.nbinsE, args.emin, args.emax)
        hx_e = ROOT.TH2D(f"h2_e_vs_x_{i+1}", "", args.nbinsX, args.xmin, args.xmax,
                         args.nbinsE, args.emin, args.emax)
        hy_e.Sumw2(); hx_e.Sumw2()
        h2y_e.append(hy_e); h2x_e.append(hx_e)

        # 1D ratio
        hr = ROOT.TH1D(f"h_r_{i+1}", "", args.nbinsR, args.rmin, args.rmax)
        hr.Sumw2()
        style_hist1d(hr, colors[i])
        h1_r.append(hr)

        # 2D ratio vs y/x (X=y/x, Y=R)
        hy_r = ROOT.TH2D(f"h2_r_vs_y_{i+1}", "", args.nbinsY, args.ymin, args.ymax,
                         args.nbinsR, args.rmin, args.rmax)
        hx_r = ROOT.TH2D(f"h2_r_vs_x_{i+1}", "", args.nbinsX, args.xmin, args.xmax,
                         args.nbinsR, args.rmin, args.rmax)
        hy_r.Sumw2(); hx_r.Sumw2()
        h2y_r.append(hy_r); h2x_r.append(hx_r)

    # --- Event loop (single pass) ---
    nentries = t.GetEntries()
    acc = [0, 0, 0]
    t0 = time.time()

    if nentries > 0:
        t.GetEntry(0)
        fE.GetNdata(); fX.GetNdata(); fY.GetNdata(); fQ2.GetNdata()
        for fc in fCuts: fc.GetNdata()

    for i in range(nentries):
        t.GetEntry(i)

        e  = float(fE.EvalInstance())
        x  = float(fX.EvalInstance())
        y  = float(fY.EvalInstance())
        q2 = float(fQ2.EvalInstance())

        has_ratio = (math.isfinite(q2) and q2 != 0.0 and math.isfinite(e))
        r = (e * 2.0 * MP / q2) if has_ratio else float("nan")

        for ic, fc in enumerate(fCuts):
            if fc.EvalInstance():
                acc[ic] += 1

                # eHCAL
                h1_e[ic].Fill(e)
                h2y_e[ic].Fill(y, e)
                h2x_e[ic].Fill(x, e)

                # ratio (only if valid)
                if has_ratio and math.isfinite(r):
                    h1_r[ic].Fill(r)
                    h2y_r[ic].Fill(y, r)
                    h2x_r[ic].Fill(x, r)

        if (i % args.print_every) == 0 or (i == nentries - 1):
            frac = (i + 1) / float(nentries) if nentries > 0 else 1.0
            dt = time.time() - t0
            rate = (i + 1) / dt if dt > 0 else 0.0
            sys.stdout.write(f"\r{render_bar(frac)}  {i+1}/{nentries}  ({rate:,.0f} ev/s)")
            sys.stdout.flush()

    print("\nDone filling.\n")

    # Normalize 1D if requested (both sets)
    if args.norm1d:
        for hh in h1_e:
            if hh.Integral() > 0:
                hh.Scale(1.0 / hh.Integral())
        for hh in h1_r:
            if hh.Integral() > 0:
                hh.Scale(1.0 / hh.Integral())

    # ---- Build Gaussian-fit "profiles" (mean & sigma graphs) + optional per-bin fit PDFs ----
    # eHCAL vs y/x
    gmean_y_e, gsigma_y_e = [], []
    gmean_x_e, gsigma_x_e = [], []

    # ratio vs y/x
    gmean_y_r, gsigma_y_r = [], []
    gmean_x_r, gsigma_x_r = [], []

    for ic in range(3):
        # filenames for fit-pages
        fitpdf_e_y = f"fits_eHCAL_vs_y_cut{ic+1}.pdf" if args.make_fit_pdfs else None
        fitpdf_e_x = f"fits_eHCAL_vs_x_cut{ic+1}.pdf" if args.make_fit_pdfs else None
        fitpdf_r_y = f"fits_ratio_vs_y_cut{ic+1}.pdf" if args.make_fit_pdfs else None
        fitpdf_r_x = f"fits_ratio_vs_x_cut{ic+1}.pdf" if args.make_fit_pdfs else None

        gm, gs, _ = gausfit_profile_from_th2(
            h2y_e[ic],
            fit_nsigma=args.fit_nsigma,
            min_entries_bin=args.min_entries_bin,
            make_fit_pages_pdf=fitpdf_e_y,
            fit_pages_title=f"eHCAL vs {args.yvar} | {labels[ic]} | per-xbin ProjectionY + Gaussian fit",
            pads_nx=args.pads_nx,
            pads_ny=args.pads_ny
        )
        gmean_y_e.append(gm); gsigma_y_e.append(gs)

        gm, gs, _ = gausfit_profile_from_th2(
            h2x_e[ic],
            fit_nsigma=args.fit_nsigma,
            min_entries_bin=args.min_entries_bin,
            make_fit_pages_pdf=fitpdf_e_x,
            fit_pages_title=f"eHCAL vs {args.xvar} | {labels[ic]} | per-xbin ProjectionY + Gaussian fit",
            pads_nx=args.pads_nx,
            pads_ny=args.pads_ny
        )
        gmean_x_e.append(gm); gsigma_x_e.append(gs)

        gm, gs, _ = gausfit_profile_from_th2(
            h2y_r[ic],
            fit_nsigma=args.fit_nsigma,
            min_entries_bin=args.min_entries_bin,
            make_fit_pages_pdf=fitpdf_r_y,
            fit_pages_title=f"R = eHCAL*2*0.938/Q2 vs {args.yvar} | {labels[ic]} | per-xbin ProjectionY + Gaussian fit",
            pads_nx=args.pads_nx,
            pads_ny=args.pads_ny
        )
        gmean_y_r.append(gm); gsigma_y_r.append(gs)

        gm, gs, _ = gausfit_profile_from_th2(
            h2x_r[ic],
            fit_nsigma=args.fit_nsigma,
            min_entries_bin=args.min_entries_bin,
            make_fit_pages_pdf=fitpdf_r_x,
            fit_pages_title=f"R = eHCAL*2*0.938/Q2 vs {args.xvar} | {labels[ic]} | per-xbin ProjectionY + Gaussian fit",
            pads_nx=args.pads_nx,
            pads_ny=args.pads_ny
        )
        gmean_x_r.append(gm); gsigma_x_r.append(gs)

    # ---- Write main PDFs (mean overlay + summary pages for sigma too) ----
    make_main_multipage_pdf(
        outpdf=args.out_e,
        var_label="eHCAL",
        xlab_1d="eHCAL",
        ylab_1d=("Normalized counts" if args.norm1d else "Counts"),
        xlab_y2d=f"{args.yvar}  [{args.ymin:g},{args.ymax:g}]",
        xlab_x2d=f"{args.xvar}  [{args.xmin:g},{args.xmax:g}]",
        ylab_2d=f"eHCAL  [{args.emin:g},{args.emax:g}]",
        labels=labels,
        counts=acc,
        h1_list=h1_e, h2y_list=h2y_e, h2x_list=h2x_e,
        gmean_y_list=gmean_y_e, gmean_x_list=gmean_x_e,
        gsigma_y_list=gsigma_y_e, gsigma_x_list=gsigma_x_e,
        prof_color=prof_color
    )

    make_main_multipage_pdf(
        outpdf=args.out_r,
        var_label="R = eHCAL*2*0.938 / Q2",
        xlab_1d="R = eHCAL*2*0.938/Q2",
        ylab_1d=("Normalized counts" if args.norm1d else "Counts"),
        xlab_y2d=f"{args.yvar}  [{args.ymin:g},{args.ymax:g}]",
        xlab_x2d=f"{args.xvar}  [{args.xmin:g},{args.xmax:g}]",
        ylab_2d=f"R  [{args.rmin:g},{args.rmax:g}]",
        labels=labels,
        counts=acc,
        h1_list=h1_r, h2y_list=h2y_r, h2x_list=h2x_r,
        gmean_y_list=gmean_y_r, gmean_x_list=gmean_x_r,
        gsigma_y_list=gsigma_y_r, gsigma_x_list=gsigma_x_r,
        prof_color=prof_color
    )

    print("=== Accepted counts (cuts) ===")
    for lab, n in zip(labels, acc):
        print(f"{lab:>10s}: N = {n}")
    print(f"Saved main: {args.out_e}")
    print(f"Saved main: {args.out_r}")
    if args.make_fit_pdfs:
        print("Saved fit-pages PDFs: fits_*_cut*.pdf (4x4 pads/page)")

if __name__ == "__main__":
    main()

