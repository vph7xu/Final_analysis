#!/usr/bin/env python3
"""
Plot G_E^n/G_D world data with method-colored markers and overlay a model
read from a lookup table providing (Q2, GEn/GD, d(GEn/GD), dGEn_Par/GD).

What it does
------------
- y-axis is **GEn/GD** (not GEn).
- Transforms data errors:  sigma(GEn/GD) = sigma(GEn) / G_D(Q^2).
- Uses the lookup’s GEn/GD directly for the model line and bands.
- Colors by method: polHe3 (blue), polH2/ND3 (red), recoil (black), other (gray).
- Distinct marker shape per paper and a paper legend “Paper (year)” colored by method,
  placed directly under the method legend.

Usage
-----
python plot_GEn_over_GD.py \
  --points "Previous GEn measurements - Sheet1.csv" \
  --lookup "neutron_lookup.dat" \
  --out-prefix "GEn_over_GD" \
  --clip \
  --dpi 450
"""

import re
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------- Helpers ----------------

def _norm_cols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [re.sub(r"\s+", " ", str(c).strip().lower()) for c in df.columns]
    return df

def _find_col(df: pd.DataFrame, patterns, required=True):
    for c in df.columns:
        for pat in patterns:
            if re.search(pat, c, flags=re.I):
                return c
    if required:
        raise KeyError(f"Missing a column matching: {patterns}")
    return None

def _method_category(m: str) -> str:
    m = (str(m) if m is not None else "").lower()
    if "polhe3" in m:
        return "polHe3"
    if ("polh2" in m) or ("nd3" in m):
        return "polH2/ND3"
    if "recoil" in m:
        return "recoil"
    return "other"

def G_D(q2, lam=0.71):
    return 1.0 / (1.0 + q2/lam)**2

# ---------------- Core ----------------

def load_points(points_csv: Path):
    """Return DataFrame with: Paper, Q2, GEn/GD (GEn_over_GD), its total error,
       method, year, and plotting helpers.
    """
    df = pd.read_csv(points_csv)
    df = _norm_cols(df)

    paper_col = _find_col(df, [r"\bpaper\b", r"\bexperiment\b", r"\bref\b"])
    method_col = _find_col(df, [r"\bmethod\b"], required=False)
    q2_col     = _find_col(df, [r"\bq\^?2\b", r"\bq2\b"])
    gen_col    = _find_col(df, [r"\bgen\b", r"\bg e n\b"])
    year_col   = _find_col(df, [r"\byear\b", r"\bpub.*year\b", r"\bdate\b"], required=False)

    # Prefer explicit 'gen stat' / 'gen sys'; else nearest to the right of GEn
    gen_idx = list(df.columns).index(gen_col)
    stat_explicit = [c for c in df.columns if re.search(r"\bgen.*stat\b|\bstat.*gen\b", c, flags=re.I)]
    sys_explicit  = [c for c in df.columns if re.search(r"\bgen.*sys\b|\bsys.*gen\b",  c, flags=re.I)]
    if stat_explicit and sys_explicit:
        gen_stat_col, gen_sys_col = stat_explicit[0], sys_explicit[0]
    else:
        stat_candidates = [(i, c) for i, c in enumerate(df.columns) if re.search(r"\bstat\b", c, flags=re.I)]
        sys_candidates  = [(i, c) for i, c in enumerate(df.columns) if re.search(r"\bsys\b",  c, flags=re.I)]
        stat_right = [(i - gen_idx, c) for i, c in stat_candidates if i > gen_idx]
        sys_right  = [(i - gen_idx, c) for i, c in sys_candidates  if i > gen_idx]
        gen_stat_col = min(stat_right, default=(1e9, None), key=lambda t: t[0])[1] or (stat_candidates[0][1] if stat_candidates else None)
        gen_sys_col  = min(sys_right,  default=(1e9, None), key=lambda t: t[0])[1] or (sys_candidates[0][1]  if sys_candidates  else None)

    out = df[[paper_col, q2_col, gen_col]].copy()
    out.columns = ["Paper", "Q2", "GEn"]
    out["GEn_stat"] = pd.to_numeric(df[gen_stat_col], errors="coerce")
    out["GEn_sys"]  = pd.to_numeric(df[gen_sys_col],  errors="coerce")
    out["Q2"]       = pd.to_numeric(out["Q2"], errors="coerce")
    out["GEn"]      = pd.to_numeric(out["GEn"], errors="coerce")
    out["method"]   = df[method_col].astype(str) if method_col else "other"

    if year_col:
        out["year"] = pd.to_numeric(df[year_col], errors="coerce").astype("Int64")
    else:
        out["year"] = pd.Series([pd.NA]*len(out), dtype="Int64")

    # Valid rows + total error
    out = out.dropna(subset=["Q2","GEn","GEn_stat","GEn_sys"]).copy()
    out["GEn_tot_err"] = np.sqrt(out["GEn_stat"]**2 + out["GEn_sys"]**2)

    # Transform to GEn/GD and its error
    gd = G_D(out["Q2"].values)
    out["GEn_over_GD"]     = out["GEn"].values / gd
    out["GEn_over_GD_err"] = out["GEn_tot_err"].values / gd

    # Method category for coloring
    out["method_cat"] = out["method"].map(_method_category)
    return out

def load_model_lookup(lookup_path: Path):
    """Read whitespace-delimited file with at least: Q2, GEn/GD, d(GEn/GD), dGEn_Par/GD.
       Return DataFrame with those and a total error column.
    """
    dfm = pd.read_csv(lookup_path, comment="#", delim_whitespace=True, header=None)
    dfm = dfm.dropna(axis=1, how="all")
    for c in dfm.columns:
        dfm[c] = pd.to_numeric(dfm[c], errors="coerce")
    dfm = dfm.dropna(how="any", subset=[dfm.columns[0], dfm.columns[1]]).copy()

    dfm = dfm.iloc[:, :4]
    dfm.columns = ["Q2", "GEn_over_GD", "dGEn_over_GD", "dGEnPar_over_GD"]
    dfm["dGEn_total"] = np.sqrt(dfm["dGEn_over_GD"]**2 + dfm["dGEnPar_over_GD"]**2)
    return dfm

def plot(points_csv: Path, lookup_path: Path, out_prefix: Path, clip_to_data: bool=True, dpi: int=300):
    pts = load_points(points_csv)
    model = load_model_lookup(lookup_path)

    # Clip model to points' Q2 range (for readability)
    if clip_to_data and len(pts) > 0:
        q2_max = pts["Q2"].max()
        model = model[(model["Q2"] >= 0) & (model["Q2"] <= q2_max*1.05)].copy()

    # Colors and markers
    color_map = {"polHe3":"blue", "polH2/ND3":"red", "recoil":"black", "other":"gray"}
    markers = ['o','s','^','D','v','>','<','p','P','X','*','h','H','+','x','d']
    unique_papers = list(pts["Paper"].astype(str).unique())
    marker_map = {paper: markers[i % len(markers)] for i, paper in enumerate(unique_papers)}

    # Representative year per paper (first non-null)
    paper_year = {}
    for p in unique_papers:
        y = pts.loc[pts["Paper"] == p, "year"].dropna()
        paper_year[p] = int(y.iloc[0]) if len(y) else None

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(9.8,6.8))

    # Data points: y = GEn/GD with transformed errors
    for paper in unique_papers:
        sub = pts[pts["Paper"] == paper]
        for _, r in sub.iterrows():
            ax.errorbar(
                [r["Q2"]], [r["GEn_over_GD"]], yerr=[r["GEn_over_GD_err"]],
                fmt=marker_map[paper], linestyle="None", markersize=7, capsize=3,
                color=color_map.get(r["method_cat"], "gray")
            )

    # Bands + line in GEn/GD
    band_total = ax.fill_between(
        model["Q2"].values,
        (model["GEn_over_GD"] - model["dGEn_total"]).values,
        (model["GEn_over_GD"] + model["dGEn_total"]).values,
        alpha=0.20, label="Model ± total"
    )
    band_stat = ax.fill_between(
        model["Q2"].values,
        (model["GEn_over_GD"] - model["dGEn_over_GD"]).values,
        (model["GEn_over_GD"] + model["dGEn_over_GD"]).values,
        alpha=0.35, label="Model ± stat"
    )
    (line_model,) = ax.plot(model["Q2"].values, model["GEn_over_GD"].values, linewidth=2, label="Model (GEn/GD)")

    ax.set_xlabel(r"$Q^2\ \mathrm{(GeV^2)}$")
    ax.set_ylabel(r"$G_E^n/G_D$")
    ax.grid(True, alpha=0.35)

    # ---- Legends ----
    # 1) Method legend (top-right)
    method_handles = [
        Line2D([0],[0], marker='o', color='blue',  linestyle='None', label='polHe3'),
        Line2D([0],[0], marker='o', color='red',   linestyle='None', label='polH2 / ND3'),
        Line2D([0],[0], marker='o', color='black', linestyle='None', label='recoil'),
        line_model,  # model line
        band_stat,   # stat band
        band_total   # total band
    ]
    leg1 = ax.legend(handles=method_handles, title="Legend",
                     loc='upper right', bbox_to_anchor=(1, 1),
                     fontsize=9, title_fontsize=10)

    # 2) Paper legend (marker shapes) – colored by each paper's method, with YEAR
    shape_handles = []
    for p in unique_papers:
        meth = pts.loc[pts['Paper'] == p, 'method_cat'].mode().iloc[0]
        col  = color_map.get(meth, 'gray')
        mk   = marker_map[p]
        yr   = paper_year.get(p, None)
        label = f"{p} ({yr})" if yr is not None else p
        shape_handles.append(Line2D([0],[0],
                                    marker=mk, linestyle='None',
                                    markerfacecolor=col, markeredgecolor=col, color=col,
                                    label=label))

    # Place the paper legend just under the method legend (same right edge)
    fig.canvas.draw()  # so leg1 has a size
    bbox_pix = leg1.get_window_extent(fig.canvas.get_renderer())
    _, y0 = ax.transAxes.inverted().transform((0, bbox_pix.y0))
    y_under = max(0, y0 - 0.02)

    leg2 = ax.legend(handles=shape_handles, title="Paper (marker shape, year)",
                     loc='upper left', bbox_to_anchor=(0, 1),
                     fontsize=10, title_fontsize=9, ncol=2,
                     handletextpad=0.6, borderpad=0.4, labelspacing=0.4)

    # Keep both legends
    #ax.add_artist(leg1)

    # ---- Save figure ----
    out_png = Path(f"{out_prefix}.png")
    out_pdf = Path(f"{out_prefix}.pdf")
    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")

    # Save the model data used (in GEn/GD)
    model_out = Path(f"{out_prefix}_model.csv")
    model[["Q2","GEn_over_GD","dGEn_over_GD","dGEnPar_over_GD","dGEn_total"]].to_csv(model_out, index=False)

    print(f"Wrote: {out_png}\n       {out_pdf}\n       {model_out}")

def main():
    p = argparse.ArgumentParser(description="Plot G_E^n/G_D with model overlay and uncertainty bands")
    p.add_argument("--points", type=Path, default=Path("Previous GEn measurements - Sheet1.csv"),
                   help="CSV with columns: Paper, Method, YEAR, Q2, GEn, (stat), (sys)")
    p.add_argument("--lookup", type=Path, default=Path("neutron_lookup.dat"),
                   help="Whitespace file with at least: Q2, GEn/GD, d(GEn/GD), dGEn_Par/GD")
    p.add_argument("--out-prefix", type=Path, default=Path("GEn_over_GD"),
                   help="Output file prefix for PNG/PDF/CSV")
    p.add_argument("--clip", action="store_true", help="Clip the model line/bands to the data Q2-range")
    p.add_argument("--dpi", type=int, default=300, help="Image DPI")
    args = p.parse_args()

    plot(points_csv=args.points, lookup_path=args.lookup,
         out_prefix=args.out_prefix, clip_to_data=args.clip, dpi=args.dpi)

if __name__ == "__main__":
    main()

