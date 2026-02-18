#!/usr/bin/env python3
# overlay_dx.py
#
# Example:
#   python3 overlay_dx.py input.root --tree TOut --branch dx \
#     --cut1 'ntrack>0&&abs(vz)<0.27&&ePS>0.2&&eHCAL>0.2&&abs((eSH+ePS)/trP-1)<0.2&&W2<1.6&&abs(dy)<0.4&&abs(coin_time-185)<5' \
#     --cut2 'ntrack>0&&abs(vz)<0.27&&ePS>0.2&&eHCAL>0.2&&abs((eSH+ePS)/trP-1)<0.2&&W2<1.6&&abs(dy)<0.4&&abs(coin_time-185)<5&&trigbits==4' \
#     --cut3 'ntrack>0&&abs(vz)<0.27&&ePS>0.2&&eHCAL>0.2&&abs((eSH+ePS)/trP-1)<0.2&&W2<1.6&&abs(dy)<0.4&&abs(coin_time-185)<5&&trigbits==1' \
#     --label1 "All" --label2 "trigbits==4" --label3 "trigbits==1" \
#     --out dx_overlay.pdf
#
import argparse
import ROOT

def integral_in_window(h, xmin, xmax):
    # integrate including bins whose centers lie in [xmin, xmax]
    b1 = h.FindBin(xmin)
    b2 = h.FindBin(xmax) 
    return h.Integral(b1, b2)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="ROOT file")
    ap.add_argument("--tree", default="Tout", help="Tree name (default: TOut)")
    ap.add_argument("--branch", default="dx", help="Branch/expression to histogram (default: dx)")
    ap.add_argument("--nbins", type=int, default=100)
    ap.add_argument("--xmin", type=float, default=-4.0)
    ap.add_argument("--xmax", type=float, default=3.0)
    ap.add_argument("--window", type=float, default=0.4, help="Count window half-width (default: 0.4 => -0.4..0.4)")
    ap.add_argument("--title", default="", help="Plot title (optional)")
    ap.add_argument("--norm", action="store_true", help="Normalize each hist to area=1")
    ap.add_argument("--out", default="dx_overlay.pdf")

    ap.add_argument("--cut1", required=True)
    ap.add_argument("--cut2", required=True)
    ap.add_argument("--cut3", required=True)
    ap.add_argument("--label1", default="all")
    ap.add_argument("--label2", default="coin")
    ap.add_argument("--label3", default="singles")
    args = ap.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleFontSize(0.04)

    f = ROOT.TFile.Open(args.input)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open file: {args.input}")

    t = f.Get(args.tree)
    if not t:
        raise RuntimeError(f"Could not find tree '{args.tree}' in {args.input}")

    cuts   = [args.cut1, args.cut2, args.cut3]
    labels = [args.label1, args.label2, args.label3]
    colors = [ROOT.kBlack, ROOT.kRed+1, ROOT.kBlue+1]

    hists = []
    counts = []

    # make hists + fill
    for i, (cut, lab, col) in enumerate(zip(cuts, labels, colors), start=1):
        hname = f"h_dx_{i}"
        h = ROOT.TH1D(hname, "", args.nbins, args.xmin, args.xmax)
        h.Sumw2()
        h.SetLineColor(col)
        h.SetLineWidth(3)
        h.SetFillStyle(0)

        # Fill with TTree::Draw
        draw_expr = f"{args.branch}>>{hname}"
        n = t.Draw(draw_expr, cut, "goff")


        # optional normalization
        if args.norm and h.Integral() > 0:
            h.Scale(1.0 / h.Integral())

        # count in -window..+window (sum of weights; for unweighted this is entries)
        cwin = integral_in_window(h, -args.window, +args.window)

        hists.append(h)
        counts.append((n, cwin))

    # canvas
    c = ROOT.TCanvas("c", "c", 900, 700)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.07)

    # set axes using first hist
    h0 = hists[0]
    h0.GetXaxis().SetTitle("dx (m)")
    h0.GetYaxis().SetTitle("Normalized counts" if args.norm else "Counts")
    h0.GetXaxis().SetTitleSize(0.045)
    h0.GetYaxis().SetTitleSize(0.045)
    h0.GetXaxis().SetLabelSize(0.04)
    h0.GetYaxis().SetLabelSize(0.04)
    h0.SetTitle(args.title)

    # draw
    h0.Draw("hist")
    for h in hists[1:]:
        h.Draw("hist same")

    # window lines
    lineL = ROOT.TLine(-args.window, 0, -args.window, h0.GetMaximum()*1.05)
    lineR = ROOT.TLine(+args.window, 0, +args.window, h0.GetMaximum()*1.05)
    for ln in (lineL, lineR):
        ln.SetLineStyle(2)
        ln.SetLineWidth(2)
        ln.SetLineColor(ROOT.kGray+2)
        ln.Draw("same")

    # legend
    leg = ROOT.TLegend(0.55, 0.68, 0.94, 0.92)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for h, lab, col, (ntot, cwin) in zip(hists, labels, colors, counts):
        # legend includes counts in window
        if args.norm:
            entry = f"{lab}  |  N[-{args.window:.1f},+{args.window:.1f}]={cwin:.4g}"
        else:
            entry = f"{lab}  |  N[-{args.window:.1f},+{args.window:.1f}]={cwin:.0f}"
        leg.AddEntry(h, entry, "l")
    leg.Draw()

    # text box with totals + window counts
    tex = ROOT.TLatex()
    tex.SetNDC(True)
    tex.SetTextSize(0.035)
    tex.SetTextFont(42)

    y0 = 0.64
    tex.DrawLatex(0.55, y0, f"Window:  {args.window:.1f} < dx < +{args.window:.1f}")
    for i, (lab, col, (ntot, cwin)) in enumerate(zip(labels, colors, counts), start=1):
        tex.SetTextColor(col)
        if args.norm:
            s = f"{lab}:  filled={ntot:d},  in-window={cwin:.4g}"
        else:
            s = f"{lab}:  filled={ntot:d},  in-window={cwin:.0f}"
        tex.DrawLatex(0.55, y0 - 0.05*i, s)
    tex.SetTextColor(ROOT.kBlack)

    c.Modified()
    c.Update()
    c.SaveAs(args.out)

    # also print to terminal
    print("\n=== dx overlay summary ===")
    for lab, (ntot, cwin) in zip(labels, counts):
        if args.norm:
            print(f"{lab:>12s} : filled={ntot:8d}   in [-{args.window},{args.window}] = {cwin:.6g}")
        else:
            print(f"{lab:>12s} : filled={ntot:8d}   in [-{args.window},{args.window}] = {cwin:.0f}")
    print(f"Saved: {args.out}\n")

if __name__ == "__main__":
    main()

