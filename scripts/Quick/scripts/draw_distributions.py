#!/usr/bin/env python3
"""
draw_distributions.py

Reads a ROOT TTree (default: Tout) and draws 1D (and optional 2D) distributions
for selected branches, with up to THREE ROOT-ish cut strings overlaid.

- Uses uproot + awkward (no PyROOT needed)
- Works with awkward v1 and v2
- Handles scalar and jagged branches:
    * scalar branches: event-level mask applied directly
    * jagged branches: event-level mask applied BEFORE flattening

Notes:
  - For jagged branches, overlay counts are element-counts after flattening.
  - For 2D, we use ONLY cut1 (common use-case). You can extend similarly if needed.

Install:
  pip install uproot awkward matplotlib numexpr
"""

import os
import re
import argparse
import numpy as np
import uproot
import awkward as ak
import matplotlib.pyplot as plt

try:
    import numexpr as ne
except ImportError:
    ne = None


# -----------------------------
# Helpers
# -----------------------------
def rootish_to_numexpr(expr: str) -> str:
    """
    Convert ROOT-ish boolean ops to numexpr-friendly ops.
    Supports: &&, ||, ! (but preserves !=)

    IMPORTANT:
      numexpr '&' and '|' can have precedence surprises with comparisons.
      So provide parentheses around each term, e.g.:
        (a>0)&&(abs(b)<1)&&(c==2)
      Your defaults do this.
    """
    expr = expr.replace("&&", "&").replace("||", "|")
    expr = expr.replace("!=", "<>")   # protect
    expr = expr.replace("!", "~")     # not
    expr = expr.replace("<>", "!=")   # restore
    return expr


def is_jagged(a) -> bool:
    try:
        return ("var *" in str(ak.type(a)))
    except Exception:
        return False


def flatten_to_numpy(a) -> np.ndarray:
    if isinstance(a, np.ndarray):
        return a
    if is_jagged(a):
        return ak.to_numpy(ak.flatten(a, axis=None))
    return ak.to_numpy(a)


def evaluate_cut(arrs: dict, cut: str):
    """
    Evaluate event-level boolean mask from a cut string using numexpr.
    Only scalar branches can be used in the cut (jagged ignored).
    """
    if cut is None or cut.strip() == "":
        return None
    if ne is None:
        raise RuntimeError("numexpr not installed. Install with: pip install numexpr")

    expr = rootish_to_numexpr(cut)

    local = {}
    for k, v in arrs.items():
        if is_jagged(v):
            continue
        local[k] = ak.to_numpy(v)

    mask = ne.evaluate(expr, local_dict=local)
    return np.asarray(mask, dtype=bool)


def select_1d_values(a, mask):
    """
    Apply event-level mask to awkward array 'a' correctly, then return numpy 1D:
      - jagged: apply mask to outer axis then flatten
      - scalar: apply mask directly
    """
    if mask is None:
        return flatten_to_numpy(a)

    if is_jagged(a):
        return flatten_to_numpy(a[mask])
    else:
        x = flatten_to_numpy(a)
        if x.shape[0] == mask.shape[0]:
            return x[mask]
        return x


# -----------------------------
# Plotting
# -----------------------------
def plot_1d_overlay_multi(xs, labels, name: str, outdir: str, bins=100, range_=None, logy=False):
    clean = []
    for x in xs:
        if x is None:
            clean.append(None)
        else:
            x = np.asarray(x)
            x = x[np.isfinite(x)]
            clean.append(x)

    if all((x is None or x.size == 0) for x in clean):
        print(f"[warn] {name}: no finite entries in any selection")
        return

    colors = ["0.75", "0.45", "0.15"]

    # ---- NEW: common bin edges ----
    allx = np.concatenate([x for x in clean if x is not None and x.size > 0])
    if range_ is None:
        xmin, xmax = float(allx.min()), float(allx.max())
    else:
        xmin, xmax = range_
    edges = np.linspace(xmin, xmax, bins + 1)
    # ------------------------------

    plt.figure()
    for x, lab, c in zip(clean, labels, colors):
        if x is None or x.size == 0:
            continue
        plt.hist(
            x, bins=edges,              # <-- use same edges for all
            histtype="stepfilled", linewidth=2,
            color=c,
            label=f"{lab} (N={x.size})"
        )

    plt.xlabel(name, fontsize=18)
    plt.ylabel("Counts", fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if logy:
        plt.yscale("log")
    plt.legend(frameon=False, fontsize=15)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{name}.png"), dpi=200)
    plt.savefig(os.path.join(outdir, f"{name}.pdf"))
    plt.close()



def plot_2d(x: np.ndarray, y: np.ndarray, name: str, outdir: str, bins=200, range_=None):
    x = np.asarray(x)
    y = np.asarray(y)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size == 0:
        print(f"[warn] {name}: no finite entries")
        return

    plt.figure()
    plt.hist2d(x, y, bins=bins, range=range_)
    plt.xlabel(name.split("_vs_")[0], fontsize=14)
    plt.ylabel(name.split("_vs_")[1], fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    cbar = plt.colorbar(label="Counts")
    cbar.ax.tick_params(labelsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{name}.png"), dpi=200)
    plt.savefig(os.path.join(outdir, f"{name}.pdf"))
    plt.close()


# -----------------------------
# Main
# -----------------------------

#gmn cuts "(ntrack>0)&&(abs(vz)<0.075)&&(eHCAL>0.025)&&(ePS>0.2)&&(abs((ePS+eSH)/trP-1)<0.2)"
#                "&&(abs(dy-0.2)<0.4)&&(abs(dx+1)<0.8)&&(abs(coin_time-0)<50)&&(W2>-1)&&(W2<3)"

#gen cuts "(ntrack>0)&&(abs(vz)<0.27)&&(eHCAL>0.025)&&(ePS>0.2)&&(abs((ePS+eSH)/trP-1)<0.2)"
#                "&&(abs(dy+0.15)<0.4)&&(abs(dx+1.53)<0.4)&&(abs(coin_time-121)<5)&&(W2>-1)&&(W2<3)"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rootfile")
    ap.add_argument("--tree", default="Tout")
    ap.add_argument("--outdir", default="plots")

    # ap.add_argument(
    #     "--cut1",
    #     default="(ntrack>0)&&(abs(vz)<0.27)&&(eHCAL>0.025)&&(ePS>0.2)&&(abs((ePS+eSH)/trP-1)<0.2)"
    #             "&&((((dx+1.1)/0.4)**2+((dy+0.0)/0.4)**2)<=1)&&(abs(coin_time-185)<5)&&(W2>-1)&&(W2<1.6)&&(eHCAL<2)",
    #     help="cut1 (base)"
    # )
    # ap.add_argument("--cut2", default="(ntrack>0)&&(abs(vz)<0.27)&&(eHCAL>0.025)&&(ePS>0.2)&&(abs((ePS+eSH)/trP-1)<0.2)"
    #             "&&((((dx+1.1)/0.4)**2+((dy+0.0)/0.4)**2)<=1)&&(abs(coin_time-185)<5)&&(W2>-1)&&(W2<1.6)&&(eHCAL<2)&&(hcal_best_idx==0)", help="cut2 (optional)")
    # ap.add_argument("--cut3", default="(ntrack>0)&&(abs(vz)<0.27)&&(eHCAL>0.025)&&(ePS>0.2)&&(abs((ePS+eSH)/trP-1)<0.2)"
    #             "&&((((dx+1.1)/0.4)**2+((dy+0.0)/0.4)**2)<=1)&&(abs(coin_time-185)<5)&&(W2>-1)&&(W2<1.6)&&(eHCAL<2)&&(hcal_best_idx!=0)", help="cut3 (optional)")

    ap.add_argument(
        "--cut1",
        default="(ntrack>0)&&(abs(vz)<0.27)&&(eHCAL>0.025)&&(ePS>0.2)&&(abs((ePS+eSH)/trP-1)<0.2)"
                "&&((((dx+1.1)/0.4)**2+((dy+0.15)/0.4)**2)<=1)&&(abs(coin_time-185)<5)&&(W2>-1)&&(W2<1.4)&&(eHCAL<2)",
        help="cut1 (base)"
    )
    ap.add_argument("--cut2", default="(ntrack>0)&&(abs(vz)<0.27)&&(eHCAL>0.025)&&(ePS>0.2)&&(abs((ePS+eSH)/trP-1)<0.2)"
                "&&((((dx+1.1)/0.4)**2+((dy+0.15)/0.4)**2)<=1)&&(abs(coin_time-185)<5)&&(W2>-1)&&(W2<1.4)&&(eHCAL<2)&&(hcal_best_idx==0)", help="cut2 (optional)")
    ap.add_argument("--cut3", default="(ntrack>0)&&(abs(vz)<0.27)&&(eHCAL>0.025)&&(ePS>0.2)&&(abs((ePS+eSH)/trP-1)<0.2)"
                "&&((((dx+1.1)/0.4)**2+((dy+0.15)/0.4)**2)<=1)&&(abs(coin_time-185)<5)&&(W2>-1)&&(W2<1.4)&&(eHCAL<2)&&(hcal_best_idx!=0)", help="cut3 (optional)")

    # ap.add_argument(
    #     "--cut1",
    #     default="(ntrack>0)&&(abs(vz)<0.075)&&(eHCAL>0.025)&&(ePS>0.2)&&(abs((ePS+eSH)/trP-1)<0.2)"
    #             "&&((((dx+1.0)/0.7)**2+((dy-0.2)/0.4)**2)<=1)&&(abs(coin_time-0)<5)&&(W2>-1)&&(W2<1.4)&&(eHCAL<2)",
    #     help="cut1 (base)"
    # )
    # ap.add_argument("--cut2", default="(ntrack>0)&&(abs(vz)<0.075)&&(eHCAL>0.025)&&(ePS>0.2)&&(abs((ePS+eSH)/trP-1)<0.2)"
    #             "&&((((dx+1.0)/0.7)**2+((dy-0.2)/0.4)**2)<=1)&&(abs(coin_time-0)<5)&&(W2>-1)&&(W2<1.4)&&(eHCAL<2)&&(hcal_best_idx==0)", help="cut2 (optional)")
    # ap.add_argument("--cut3", default="(ntrack>0)&&(abs(vz)<0.075)&&(eHCAL>0.025)&&(ePS>0.2)&&(abs((ePS+eSH)/trP-1)<0.2)"
    #             "&&((((dx+1.0)/0.7)**2+((dy-0.2)/0.4)**2)<=1)&&(abs(coin_time-0)<5)&&(W2>-1)&&(W2<1.4)&&(eHCAL<2)&&(hcal_best_idx!=0)", help="cut3 (optional)")




    ap.add_argument("--label1", default="cut1")
    ap.add_argument("--label2", default="cut2")
    ap.add_argument("--label3", default="cut3")

    ap.add_argument("--max-events", type=int, default=None, help="optional entry_stop")
    ap.add_argument(
        "--vars",
        default="Q2,W2,vz,dx,dy,eHCAL,xHCAL,yHCAL,coin_time,trigbits,nclus_HCAL",
        help="comma-separated 1D variables",
    )
    ap.add_argument(
        "--vars2d",
        default="",
        help="comma-separated 2D as x:y (e.g. dx:dy,xHCAL:yHCAL). 2D uses cut1 only.",
    )
    ap.add_argument("--bins", type=int, default=200)
    ap.add_argument("--logy", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    f = uproot.open(args.rootfile)
    tree = f[args.tree]

    vars_1d = [v.strip() for v in args.vars.split(",") if v.strip()]
    vars_2d = [v.strip() for v in args.vars2d.split(",") if v.strip()]

    branches = set(vars_1d)
    for v in vars_2d:
        x, y = v.split(":")
        branches.add(x.strip())
        branches.add(y.strip())

    # Add cut tokens that exist in tree
    tree_keys = set(tree.keys())
    for cut in [args.cut1, args.cut2, args.cut3]:
        if cut and cut.strip():
            tokens = re.findall(r"[A-Za-z_]\w*", cut)
            for t in tokens:
                if t in tree_keys:
                    branches.add(t)

    branches = sorted(branches)
    arrs = tree.arrays(branches, how=dict, entry_stop=args.max_events)

    mask1 = evaluate_cut(arrs, args.cut1) if (args.cut1 and args.cut1.strip()) else None
    mask2 = evaluate_cut(arrs, args.cut2) if (args.cut2 and args.cut2.strip()) else None
    mask3 = evaluate_cut(arrs, args.cut3) if (args.cut3 and args.cut3.strip()) else None

    # Print event yields (event-level counts)
    def _yield(mask, name):
        if mask is None:
            return
        print(f"[info] Events passing {name}: {mask.sum()} / {mask.size}")

    if mask1 is None:
        print("[info] cut1 empty -> using all events")
    else:
        _yield(mask1, "cut1")
    if mask2 is not None:
        _yield(mask2, "cut2")
    if mask3 is not None:
        _yield(mask3, "cut3")

    labels = [args.label1, args.label2, args.label3]

    # 1D overlay plots
    for name in vars_1d:
        a = arrs[name]
        x1 = select_1d_values(a, mask1)
        x2 = select_1d_values(a, mask2) if mask2 is not None else None
        x3 = select_1d_values(a, mask3) if mask3 is not None else None

        plot_1d_overlay_multi(
            [x1, x2, x3],
            labels,
            name,
            args.outdir,
            bins=args.bins,
            logy=args.logy
        )

    # 2D plots: cut1 only, scalar only
    for spec in vars_2d:
        xname, yname = [s.strip() for s in spec.split(":")]
        ax = arrs[xname]
        ay = arrs[yname]

        if is_jagged(ax) or is_jagged(ay):
            print(f"[skip] {xname}:{yname} is jagged; implement jagged 2D if needed.")
            continue

        x = ak.to_numpy(ax)
        y = ak.to_numpy(ay)

        if mask1 is not None:
            x = x[mask1]
            y = y[mask1]

        plot_2d(x, y, f"{xname}_vs_{yname}", args.outdir, bins=args.bins)

    print(f"Saved plots to: {args.outdir}/")


if __name__ == "__main__":
    main()
