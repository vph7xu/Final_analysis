#!/usr/bin/env python3
# overlay_eHCAL_and_ratio_multipage_with_profile.py
#
# One loop over data.
# Produces TWO multipage PDFs:
#   (1) eHCAL: 1D overlay + 2D vs y/x with ProfileX overlay
#   (2) R = eHCAL*2*0.938/Q2: same plots, separate PDF
#
# PLUS (optional):
#   (3) Trigger plots PDF:
#       - 1D overlay xHCAL for All + each trigger category
#       - 1D overlay yHCAL for All + each trigger category
#       - 2D yHCAL vs xHCAL per category
#       - (NEW) Ratio pages:
#           Coin/All vs x, Singles/All vs x
#           Coin/All vs y, Singles/All vs y
#
# ROOT 6.26-friendly: manual loop + ASCII progress bar

import argparse
import sys
import time
import math
import re
import ROOT

MP = 0.938  # GeV

def safe_name(s: str) -> str:
    """Make a ROOT-safe object name from an arbitrary string."""
    return re.sub(r"[^A-Za-z0-9_]+", "_", s)

def style_hist1d(h, color, width=3):
    h.SetLineColor(color)
    h.SetLineWidth(width)
    h.SetFillStyle(0)

def style_ratio_hist(h, color, mstyle=20, msize=1.2, lwidth=2):
    h.SetMarkerColor(color)
    h.SetLineColor(color)
    h.SetLineWidth(lwidth)
    h.SetMarkerStyle(mstyle)
    h.SetMarkerSize(msize)

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

def style_profile(p, color=ROOT.kBlack, mstyle=20, msize=1.5, lwidth=2):
    p.SetMarkerStyle(mstyle)
    p.SetMarkerSize(msize)
    p.SetMarkerColor(color)
    p.SetLineColor(color)
    p.SetLineWidth(lwidth)

def set_overlay_max(hlist, scale=1.25):
    m = 0.0
    for h in hlist:
        m = max(m, h.GetMaximum())
    if m > 0:
        for h in hlist:
            h.SetMaximum(scale * m)

def find_index_by_key(labels, key, prefer_exact=False):
    """
    Find index of a label by matching a substring key (case-insensitive).
    prefer_exact=True tries exact match first, then substring.
    Returns None if not found.
    """
    if not labels:
        return None
    key_l = key.strip().lower()
    if key_l == "":
        return None

    if prefer_exact:
        for i, lab in enumerate(labels):
            if lab.strip().lower() == key_l:
                return i

    for i, lab in enumerate(labels):
        if key_l in lab.strip().lower():
            return i
    return None

def parse_trigger_items(trigs_str, trigvar, trig_mode, trig_basecut):
    """
    Returns list of (default_label, expression) trigger categories.
    Supports tokens like "all,1,4" (case-insensitive 'all').
    If trigs_str is empty -> defaults to ["all"].
    """
    if trigs_str is None or trigs_str.strip() == "":
        toks = ["all"]
    else:
        toks = [s.strip() for s in trigs_str.split(",") if s.strip() != ""]

    want_all = any(s.lower() == "all" for s in toks)

    trig_vals = []
    for s in toks:
        if s.lower() == "all":
            continue
        trig_vals.append(int(s))

    items = []
    if want_all:
        items.append(("All", f"({trig_basecut})"))

    for tv in trig_vals:
        if trig_mode == "eq":
            expr = f"({trig_basecut}) && ({trigvar}=={tv})"
            lab  = f"{trigvar}=={tv}"
        else:
            expr = f"({trig_basecut}) && (({trigvar} & {tv}) != 0)"
            lab  = f"({trigvar}&{tv})!=0"
        items.append((lab, expr))

    return items

def apply_user_labels(default_labels, triglabels_str):
    """Override default labels with user-provided comma list, if provided."""
    if triglabels_str is None or triglabels_str.strip() == "":
        return default_labels
    user = [s.strip() for s in triglabels_str.split(",")]
    out = []
    for i, d in enumerate(default_labels):
        if i < len(user) and user[i] != "":
            out.append(user[i])
        else:
            out.append(d)
    return out

def make_multipage_pdf(
    outpdf,
    prefix,             # unique prefix for object names
    var_label,
    xlab_1d, ylab_1d,
    xlab_y2d, xlab_x2d,
    ylab_2d,
    labels,
    counts,
    h1_list, h2y_list, h2x_list,
    prof_color=ROOT.kBlack
):
    prof_keep = []

    c = ROOT.TCanvas(f"c_{safe_name(prefix)}", f"c_{safe_name(prefix)}", 1000, 800)
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

    # Pages 2-4: 2D vs y with profile
    for j in range(3):
        c.Clear()
        setup_pad(mright=0.14)
        h2y_list[j].GetXaxis().SetTitle(xlab_y2d)
        h2y_list[j].GetYaxis().SetTitle(ylab_2d)
        h2y_list[j].Draw("COLZ")

        pname = safe_name(f"pY_{prefix}_{j}")
        p = h2y_list[j].ProfileX(pname, 1, -1, "")
        style_profile(p, color=prof_color, mstyle=20, msize=1.6, lwidth=2)
        p.Draw("SAME PE1")
        prof_keep.append(p)

        label_ndc(f"{var_label} vs y | {labels[j]}  (ProfileX = mean per x-bin)", 0.14, 0.93, 0.040)
        c.Print(outpdf)

    # Pages 5-7: 2D vs x with profile
    for j in range(3):
        c.Clear()
        setup_pad(mright=0.14)
        h2x_list[j].GetXaxis().SetTitle(xlab_x2d)
        h2x_list[j].GetYaxis().SetTitle(ylab_2d)
        h2x_list[j].Draw("COLZ")

        pname = safe_name(f"pX_{prefix}_{j}")
        p = h2x_list[j].ProfileX(pname, 1, -1, "")
        style_profile(p, color=prof_color, mstyle=20, msize=1.6, lwidth=2)
        p.Draw("SAME PE1")
        prof_keep.append(p)

        label_ndc(f"{var_label} vs x | {labels[j]}  (ProfileX = mean per x-bin)", 0.14, 0.93, 0.040)
        c.Print(outpdf)

    c.Print(outpdf + "]")

def make_xy_trigger_pdf(
    outpdf,
    xlab, ylab,
    labels, counts, colors,
    hx_list, hy_list, hxy_list,
    norm1d=False,
    all_key="all",
    coin_key="coin",
    singles_key="single"
):
    """
    Multipage PDF:
      page 1: xHCAL overlay by category (All + per trigger)
      page 2: yHCAL overlay by category
      page 3..: y vs x 2D per category
      plus (NEW) ratio pages:
        - Coin/All and Singles/All vs x
        - Coin/All and Singles/All vs y
    """
    c = ROOT.TCanvas("c_xy", "c_xy", 1000, 800)
    c.Print(outpdf + "[")

    # 1) x overlay
    c.Clear()
    setup_pad(mright=0.04)
    hx_list[0].GetXaxis().SetTitle(f"{xlab} (m)")
    hx_list[0].GetYaxis().SetTitle("Normalized counts" if norm1d else "Counts")
    hx_list[0].Draw("hist")
    for h in hx_list[1:]:
        h.Draw("hist same")

    leg = ROOT.TLegend(0.55, 0.62, 0.94, 0.92)
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    for h, lab, n in zip(hx_list, labels, counts):
        leg.AddEntry(h, f"{lab} (N={n})", "l")
    leg.Draw()
    label_ndc(f"1D overlay: {xlab} (All + triggers)", 0.14, 0.93, 0.045)
    c.Print(outpdf)

    # 2) y overlay
    c.Clear()
    setup_pad(mright=0.04)
    hy_list[0].GetXaxis().SetTitle(f"{ylab} (m)")
    hy_list[0].GetYaxis().SetTitle("Normalized counts" if norm1d else "Counts")
    hy_list[0].Draw("hist")
    for h in hy_list[1:]:
        h.Draw("hist same")

    leg = ROOT.TLegend(0.55, 0.62, 0.94, 0.92)
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    for h, lab, n in zip(hy_list, labels, counts):
        leg.AddEntry(h, f"{lab} (N={n})", "l")
    leg.Draw()
    label_ndc(f"1D overlay: {ylab} (All + triggers)", 0.14, 0.93, 0.045)
    c.Print(outpdf)

    # 3.. ) 2D y vs x per category
    for i, lab in enumerate(labels):
        c.Clear()
        setup_pad(mright=0.14)
        hxy_list[i].GetXaxis().SetTitle(f"{xlab} (m)")
        hxy_list[i].GetYaxis().SetTitle(f"{ylab} (m)")
        hxy_list[i].Draw("COLZ")
        label_ndc(f"2D: {ylab} vs {xlab} | {lab}  (N={counts[i]})", 0.14, 0.93, 0.040)
        c.Print(outpdf)

    # ---- NEW: ratio pages (Coin/All and Singles/All) ----
    i_all = find_index_by_key(labels, all_key, prefer_exact=True)
    if i_all is None:
        i_all = find_index_by_key(labels, all_key, prefer_exact=False)

    i_coin = find_index_by_key(labels, coin_key, prefer_exact=False)
    i_sing = find_index_by_key(labels, singles_key, prefer_exact=False)

    # Only do ratio pages if we found All and at least one of Coin/Singles
    if i_all is not None and (i_coin is not None or i_sing is not None):
        # helper to make a ratio hist with binomial errors
        def make_ratio(num, den, name):
            h = num.Clone(name)
            h.Reset("ICESM")
            # h = num/den with binomial errors ("B") is meaningful if num is subset of den
            h.Divide(num, den, 1.0, 1.0, "B")
            return h

        # Page: ratios vs x
        c.Clear()
        setup_pad(mright=0.04)
        frame = hx_list[i_all].Clone("ratio_frame_x")
        frame.Reset("ICESM")
        frame.GetXaxis().SetTitle(f"{xlab} (m)")
        frame.GetYaxis().SetTitle("Fraction")
        frame.SetMinimum(0.0)
        frame.SetMaximum(1.2)
        frame.Draw("AXIS")

        leg = ROOT.TLegend(0.55, 0.72, 0.94, 0.92)
        leg.SetBorderSize(0); leg.SetFillStyle(0)

        drawn_any = False
        if i_coin is not None:
            hcx = make_ratio(hx_list[i_coin], hx_list[i_all], "h_coin_over_all_x")
            style_ratio_hist(hcx, colors[i_coin], mstyle=20, msize=1.2)
            hcx.Draw("SAME PE1")
            leg.AddEntry(hcx, f"{labels[i_coin]} / {labels[i_all]}", "pe")
            drawn_any = True

        if i_sing is not None:
            hsx = make_ratio(hx_list[i_sing], hx_list[i_all], "h_sing_over_all_x")
            style_ratio_hist(hsx, colors[i_sing], mstyle=21, msize=1.2)
            hsx.Draw("SAME PE1")
            leg.AddEntry(hsx, f"{labels[i_sing]} / {labels[i_all]}", "pe")
            drawn_any = True

        if drawn_any:
            leg.Draw()
            label_ndc(f"Ratios vs {xlab}: Coin/All and Singles/All", 0.14, 0.93, 0.045)
            c.Print(outpdf)

        # Page: ratios vs y
        c.Clear()
        setup_pad(mright=0.04)
        frame = hy_list[i_all].Clone("ratio_frame_y")
        frame.Reset("ICESM")
        frame.GetXaxis().SetTitle(f"{ylab} (m)")
        frame.GetYaxis().SetTitle("Fraction")
        frame.SetMinimum(0.0)
        frame.SetMaximum(1.2)
        frame.Draw("AXIS")

        leg = ROOT.TLegend(0.55, 0.72, 0.94, 0.92)
        leg.SetBorderSize(0); leg.SetFillStyle(0)

        drawn_any = False
        if i_coin is not None:
            hcy = make_ratio(hy_list[i_coin], hy_list[i_all], "h_coin_over_all_y")
            style_ratio_hist(hcy, colors[i_coin], mstyle=20, msize=1.2)
            hcy.Draw("SAME PE1")
            leg.AddEntry(hcy, f"{labels[i_coin]} / {labels[i_all]}", "pe")
            drawn_any = True

        if i_sing is not None:
            hsy = make_ratio(hy_list[i_sing], hy_list[i_all], "h_sing_over_all_y")
            style_ratio_hist(hsy, colors[i_sing], mstyle=21, msize=1.2)
            hsy.Draw("SAME PE1")
            leg.AddEntry(hsy, f"{labels[i_sing]} / {labels[i_all]}", "pe")
            drawn_any = True

        if drawn_any:
            leg.Draw()
            label_ndc(f"Ratios vs {ylab}: Coin/All and Singles/All", 0.14, 0.93, 0.045)
            c.Print(outpdf)

    else:
        # still produce the PDF, but tell the user why ratios didn't appear
        print("[NOTE] Ratio pages (Coin/All, Singles/All) not added to xy PDF.")
        print(f"       Could not find required labels: all='{all_key}', coin='{coin_key}', singles='{singles_key}'.")
        print(f"       Labels present: {labels}")

    c.Print(outpdf + "]")

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

    # y/x ranges for 2D (also used for trigger plots)
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

    ap.add_argument("--print_every", type=int, default=20000)
    ap.add_argument("--profile_color", default="black", help="black/white/red/blue")

    # Trigger plots options
    ap.add_argument("--make_trig_plots", action="store_true",
                    help="Also make xHCAL/yHCAL and xHCAL vs yHCAL plots for All + trigger categories")
    ap.add_argument("--out_xy", default="xyHCAL_byTrig.pdf",
                    help="Output PDF for trigger x/y plots")
    ap.add_argument("--trigvar", default="trigbits",
                    help="Trigger variable/branch (default: trigbits)")
    ap.add_argument("--trigs", default="",
                    help='Comma list of trigger values; optionally include "all". '
                         'e.g. "all,1,4,8" or "1,4,8". If empty, defaults to "all" when --make_trig_plots is set.')
    ap.add_argument("--triglabels", default="",
                    help='Comma labels matching --trigs-expanded order. Example: "All,Singles,Coin,Other"')
    ap.add_argument("--trig_mode", default="eq", choices=["eq", "bitand"],
                    help='eq -> trigvar==val, bitand -> (trigvar & val)!=0')
    ap.add_argument("--trig_basecut", default="__USE_CUT1__",
                    help='Basecut for trigger plots; default uses cut1.')

    # NEW: ratio matching keys
    ap.add_argument("--all_key", default="all", help='Substring used to identify the "All" category label (default: all)')
    ap.add_argument("--coin_key", default="coin", help='Substring used to identify the "Coin" category label (default: coin)')
    ap.add_argument("--singles_key", default="single", help='Substring used to identify the "Singles" category label (default: single)')

    args = ap.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetNumberContours(80)

    prof_color = {
        "black": ROOT.kBlack,
        "white": ROOT.kWhite,
        "red":   ROOT.kRed+1,
        "blue":  ROOT.kBlue+1,
    }.get(args.profile_color.lower(), ROOT.kBlack)

    f = ROOT.TFile.Open(args.input)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open file: {args.input}")
    t = f.Get(args.tree)
    if not t:
        raise RuntimeError(f"Could not find tree '{args.tree}' in {args.input}")

    cuts   = [args.cut1, args.cut2, args.cut3]
    labels = [args.label1, args.label2, args.label3]
    colors_cut = [ROOT.kBlack, ROOT.kRed+1, ROOT.kBlue+1]

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
        he = ROOT.TH1D(f"h_e_{i+1}", "", args.nbinsE, args.emin, args.emax)
        he.Sumw2()
        style_hist1d(he, colors_cut[i])
        h1_e.append(he)

        hy_e = ROOT.TH2D(f"h2_e_vs_y_{i+1}", "", args.nbinsY, args.ymin, args.ymax,
                         args.nbinsE, args.emin, args.emax)
        hx_e = ROOT.TH2D(f"h2_e_vs_x_{i+1}", "", args.nbinsX, args.xmin, args.xmax,
                         args.nbinsE, args.emin, args.emax)
        hy_e.Sumw2(); hx_e.Sumw2()
        h2y_e.append(hy_e); h2x_e.append(hx_e)

        hr = ROOT.TH1D(f"h_r_{i+1}", "", args.nbinsR, args.rmin, args.rmax)
        hr.Sumw2()
        style_hist1d(hr, colors_cut[i])
        h1_r.append(hr)

        hy_r = ROOT.TH2D(f"h2_r_vs_y_{i+1}", "", args.nbinsY, args.ymin, args.ymax,
                         args.nbinsR, args.rmin, args.rmax)
        hx_r = ROOT.TH2D(f"h2_r_vs_x_{i+1}", "", args.nbinsX, args.xmin, args.xmax,
                         args.nbinsR, args.rmin, args.rmax)
        hy_r.Sumw2(); hx_r.Sumw2()
        h2y_r.append(hy_r); h2x_r.append(hx_r)

    # -------- Trigger plots setup --------
    make_trig = bool(args.make_trig_plots)

    trig_labels = []
    trig_colors = []
    fTrigCuts = []
    hx_trig, hy_trig, hxy_trig = [], [], []
    trig_counts = []

    if make_trig:
        trig_basecut = args.cut1 if args.trig_basecut == "__USE_CUT1__" else args.trig_basecut

        trig_items = parse_trigger_items(
            trigs_str=args.trigs,
            trigvar=args.trigvar,
            trig_mode=args.trig_mode,
            trig_basecut=trig_basecut
        )

        default_tlabels = [lab for (lab, _) in trig_items]
        trig_labels = apply_user_labels(default_tlabels, args.triglabels)

        base_colors = [ROOT.kBlack, ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1,
                       ROOT.kOrange+7, ROOT.kCyan+2, ROOT.kViolet+1]

        for i, (_, expr) in enumerate(trig_items):
            fTrigCuts.append(ROOT.TTreeFormula(f"fTrig_{i}", expr, t))

            hx = ROOT.TH1D(f"hx_trig_{i}", "", args.nbinsX, args.xmin, args.xmax); hx.Sumw2()
            hy = ROOT.TH1D(f"hy_trig_{i}", "", args.nbinsY, args.ymin, args.ymax); hy.Sumw2()
            h2 = ROOT.TH2D(f"hxy_trig_{i}", "", args.nbinsX, args.xmin, args.xmax,
                           args.nbinsY, args.ymin, args.ymax); h2.Sumw2()

            col = base_colors[i % len(base_colors)]
            trig_colors.append(col)

            style_hist1d(hx, col, width=3)
            style_hist1d(hy, col, width=3)

            hx_trig.append(hx)
            hy_trig.append(hy)
            hxy_trig.append(h2)

        trig_counts = [0] * len(fTrigCuts)

    # --- Event loop (single pass) ---
    nentries = t.GetEntries()
    acc = [0, 0, 0]
    t0 = time.time()

    if nentries > 0:
        t.GetEntry(0)
        fE.GetNdata(); fX.GetNdata(); fY.GetNdata(); fQ2.GetNdata()
        for fc in fCuts:
            fc.GetNdata()
        if make_trig:
            for ft in fTrigCuts:
                ft.GetNdata()

    for i in range(nentries):
        t.GetEntry(i)

        e  = float(fE.EvalInstance())
        x  = float(fX.EvalInstance())
        y  = float(fY.EvalInstance())
        q2 = float(fQ2.EvalInstance())

        has_ratio = (math.isfinite(q2) and q2 != 0.0 and math.isfinite(e))
        r = (e * 2.0 * MP / q2) if has_ratio else float("nan")

        # Fill eHCAL/ratio cut histos
        for ic, fc in enumerate(fCuts):
            if fc.EvalInstance():
                acc[ic] += 1
                h1_e[ic].Fill(e)
                h2y_e[ic].Fill(y, e)
                h2x_e[ic].Fill(x, e)

                if has_ratio and math.isfinite(r):
                    h1_r[ic].Fill(r)
                    h2y_r[ic].Fill(y, r)
                    h2x_r[ic].Fill(x, r)

        # Fill trigger x/y histos
        if make_trig:
            for it, ft in enumerate(fTrigCuts):
                if ft.EvalInstance():
                    trig_counts[it] += 1
                    hx_trig[it].Fill(x)
                    hy_trig[it].Fill(y)
                    hxy_trig[it].Fill(x, y)

        if (i % args.print_every) == 0 or (i == nentries - 1):
            frac = (i + 1) / float(nentries) if nentries > 0 else 1.0
            dt = time.time() - t0
            rate = (i + 1) / dt if dt > 0 else 0.0
            sys.stdout.write(f"\r{render_bar(frac)}  {i+1}/{nentries}  ({rate:,.0f} ev/s)")
            sys.stdout.flush()

    print("\nDone filling.\n")

    # Normalize 1D if requested
    if args.norm1d:
        for hh in h1_e:
            if hh.Integral() > 0:
                hh.Scale(1.0 / hh.Integral())
        for hh in h1_r:
            if hh.Integral() > 0:
                hh.Scale(1.0 / hh.Integral())
        if make_trig:
            for hh in hx_trig:
                if hh.Integral() > 0:
                    hh.Scale(1.0 / hh.Integral())
            for hh in hy_trig:
                if hh.Integral() > 0:
                    hh.Scale(1.0 / hh.Integral())

    # Consistent overlay maxima for trigger overlays
    if make_trig:
        set_overlay_max(hx_trig)
        set_overlay_max(hy_trig)

    # --- Write eHCAL PDF ---
    make_multipage_pdf(
        outpdf=args.out_e,
        prefix="eHCAL",
        var_label="eHCAL",
        xlab_1d="eHCAL",
        ylab_1d=("Normalized counts" if args.norm1d else "Counts"),
        xlab_y2d=f"{args.yvar}  [{args.ymin:g},{args.ymax:g}]",
        xlab_x2d=f"{args.xvar}  [{args.xmin:g},{args.xmax:g}]",
        ylab_2d=f"eHCAL  [{args.emin:g},{args.emax:g}]",
        labels=labels,
        counts=acc,
        h1_list=h1_e, h2y_list=h2y_e, h2x_list=h2x_e,
        prof_color=prof_color
    )

    # --- Write ratio PDF ---
    make_multipage_pdf(
        outpdf=args.out_r,
        prefix="ratio",
        var_label="R = eHCAL*2*0.938 / Q2",
        xlab_1d="R = eHCAL*2*0.938 / Q2",
        ylab_1d=("Normalized counts" if args.norm1d else "Counts"),
        xlab_y2d=f"{args.yvar}  [{args.ymin:g},{args.ymax:g}]",
        xlab_x2d=f"{args.xvar}  [{args.xmin:g},{args.xmax:g}]",
        ylab_2d=f"R  [{args.rmin:g},{args.rmax:g}]",
        labels=labels,
        counts=acc,
        h1_list=h1_r, h2y_list=h2y_r, h2x_list=h2x_r,
        prof_color=prof_color
    )

    # --- Write trigger PDF ---
    if make_trig:
        make_xy_trigger_pdf(
            outpdf=args.out_xy,
            xlab=args.xvar,
            ylab=args.yvar,
            labels=trig_labels,
            counts=trig_counts,
            colors=trig_colors,
            hx_list=hx_trig,
            hy_list=hy_trig,
            hxy_list=hxy_trig,
            norm1d=args.norm1d,
            all_key=args.all_key,
            coin_key=args.coin_key,
            singles_key=args.singles_key
        )

    # --- Summary ---
    print("=== Accepted counts (cuts) ===")
    for lab, n in zip(labels, acc):
        print(f"{lab:>10s}: N = {n}")
    print(f"Saved: {args.out_e}")
    print(f"Saved: {args.out_r}")

    if make_trig:
        print("=== Trigger counts (trig categories) ===")
        for lab, n in zip(trig_labels, trig_counts):
            print(f"{lab:>18s}: N = {n}")
        print(f"Saved: {args.out_xy}")

if __name__ == "__main__":
    main()
