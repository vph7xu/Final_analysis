#!/usr/bin/env python3
"""
plot_GMn_over_muGD.py

World data: plot GMn in the standard normalization:

    R(Q2) = G_M^n(Q2) / (mu_n * G_D(Q2))

Inputs
------
1) CSV world points:
   Columns (case/space insensitive):
     Paper (or Author/Ref/Experiment), Method (optional), Year (optional), Q2
     and either:
       GMn, GMn_stat, GMn_sys   (raw GMn in nuclear magnetons)
     or:
       R, R_stat, R_sys         (already GMn/(mu_n*G_D))

2) Lookup file (Ye-style):
   Your neutron_lookup.dat format (commented header starting with ##):
     Q2  GEn/GD  dGEn/GD  dGEn_Par/GD  GMn/muGD  dGMn/muGD  dGMn_Par/muGD

   Or generic 4-column whitespace:
     Q2  R  dR_stat  dR_par

Plot settings
-------------
- Axis range to match your figure: Q2 in [0,12], R in [0.2,1.25]

Usage
-----
python3 plot_GMn_over_muGD.py --no-hepdata --points GMn_world_points.csv --lookup neutron_lookup.dat
"""

import re
import argparse
from pathlib import Path
from typing import Optional, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# requests is only needed if you enable HEPData ingestion (not used here by default)
try:
    import requests
except Exception:
    requests = None


# ---------------- Physics helpers ----------------

def G_D(q2, lam=0.71):
    """Standard dipole form factor: G_D(Q2) = 1/(1+Q2/0.71)^2, Q2 in GeV^2."""
    q2 = np.asarray(q2, dtype=float)
    return 1.0 / (1.0 + q2 / lam) ** 2

MU_N = -1.9130427  # neutron magnetic moment (nuclear magnetons), sign included


# ---------------- DataFrame helpers ----------------

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

def _fallback_first_col(df: pd.DataFrame) -> str:
    # If the sheet export gave "Unnamed: 0" or similar, prefer that, else first column
    for c in df.columns:
        if re.search(r"^unnamed", c, flags=re.I):
            return c
    return df.columns[0]

def _method_category(m: str) -> str:
    m = (str(m) if m is not None else "").lower()
    if "lt" in m:
        return "lt"
    if "he3" in m or "3he" in m or "pol" in m:
        return "polHe3"
    if "ratio" in m or "en/ep" in m or "d(e,e'n)/d(e,e'p)" in m:
        return "ratio"
    if "inclusive" in m:
        return "inclusive"

    return "other"


# ---------------- Local CSV ingestion ----------------

def load_points_csv(points_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(points_csv)
    df = _norm_cols(df)

    paper_col  = _find_col(df, [r"\bpaper\b", r"\bauthor\b", r"\bexperiment\b", r"\bref\b"], required=False)
    if paper_col is None:
        paper_col = _fallback_first_col(df)

    method_col = _find_col(df, [r"\bmethod\b"], required=False)
    year_col   = _find_col(df, [r"\byear\b", r"\bpub.*year\b", r"\bdate\b"], required=False)
    q2_col     = _find_col(df, [r"\bq\^?2\b", r"\bq2\b"])

    gmn_col = _find_col(df, [r"\bgmn\b", r"\bg m n\b", r"\bg_mn\b"], required=False)
    r_col   = _find_col(df, [r"^\s*r\s*$", r"\br\b", r"mu.*gd", r"over.*mu", r"gmn.*/.*mu"], required=False)

    def err_cols(prefix):
        stat = [c for c in df.columns if re.search(fr"\b{prefix}.*stat\b|\bstat.*{prefix}\b", c, flags=re.I)]
        sysc = [c for c in df.columns if re.search(fr"\b{prefix}.*sys\b|\bsys.*{prefix}\b",  c, flags=re.I)]
        return (stat[0] if stat else None, sysc[0] if sysc else None)

    gmn_stat_col, gmn_sys_col = err_cols("gmn")
    r_stat_col,   r_sys_col   = err_cols("r")

    out = pd.DataFrame()
    out["Paper"]  = df[paper_col].astype(str).str.strip()
    out["method"] = df[method_col].astype(str).str.strip() if method_col else "other"
    out["year"]   = pd.to_numeric(df[year_col], errors="coerce").astype("Int64") if year_col else pd.Series([pd.NA]*len(df), dtype="Int64")
    out["Q2"]     = pd.to_numeric(df[q2_col], errors="coerce")

    # drop accidental embedded header rows, if any
    out = out[~out["Paper"].str.fullmatch("paper|author|ref", case=False, na=False)].copy()

    if r_col is not None and df[r_col].astype(str).str.strip().replace("nan", "").ne("").any():
        out["y_kind"] = "R"
        out["y"]      = pd.to_numeric(df[r_col], errors="coerce")
        out["y_stat"] = pd.to_numeric(df[r_stat_col], errors="coerce") if r_stat_col else 0.0
        out["y_sys"]  = pd.to_numeric(df[r_sys_col],  errors="coerce") if r_sys_col  else 0.0
    elif gmn_col is not None and df[gmn_col].astype(str).str.strip().replace("nan", "").ne("").any():
        out["y_kind"] = "GMn"
        out["y"]      = pd.to_numeric(df[gmn_col], errors="coerce")
        out["y_stat"] = pd.to_numeric(df[gmn_stat_col], errors="coerce") if gmn_stat_col else 0.0
        out["y_sys"]  = pd.to_numeric(df[gmn_sys_col],  errors="coerce") if gmn_sys_col  else 0.0
    else:
        raise KeyError("CSV must contain either GMn (and errors) or R (and errors).")

    out["method_cat"] = out["method"].map(_method_category)
    out = out.dropna(subset=["Paper", "Q2", "y"]).copy()
    return out


# ---------------- Transform to R ----------------

def to_R_table(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    q2 = df["Q2"].to_numpy(dtype=float)
    gd = G_D(q2)
    mu = MU_N

    is_gmn = df["y_kind"].str.upper().eq("GMN")
    is_r   = df["y_kind"].str.upper().eq("R")

    R = np.full(len(df), np.nan, dtype=float)
    R_stat = np.full(len(df), np.nan, dtype=float)
    R_sys  = np.full(len(df), np.nan, dtype=float)

    denom = np.abs(mu * gd)

    # raw GMn -> R
    if is_gmn.any():
        R[is_gmn]      = df.loc[is_gmn, "y"].to_numpy(dtype=float) / (mu * gd[is_gmn])
        R_stat[is_gmn] = df.loc[is_gmn, "y_stat"].to_numpy(dtype=float) / denom[is_gmn]
        R_sys[is_gmn]  = df.loc[is_gmn, "y_sys"].to_numpy(dtype=float)  / denom[is_gmn]

    # already R
    if is_r.any():
        R[is_r]      = df.loc[is_r, "y"].to_numpy(dtype=float)
        R_stat[is_r] = df.loc[is_r, "y_stat"].to_numpy(dtype=float)
        R_sys[is_r]  = df.loc[is_r, "y_sys"].to_numpy(dtype=float)

    df["R"] = R
    df["R_stat"] = R_stat
    df["R_sys"]  = R_sys
    df["R_tot"]  = np.sqrt(df["R_stat"]**2 + df["R_sys"]**2)

    df = df.dropna(subset=["R", "R_tot"]).copy()
    return df


# ---------------- Lookup loader (supports your neutron_lookup.dat) ----------------

def load_lookup_any(lookup_path: Path) -> Tuple[pd.DataFrame, str]:
    """
    Returns (df, label_kind)
      df columns: Q2, R, dR_stat, dR_par, dR_tot
      label_kind: "ye2018" or "generic"
    """
    # comment="#" will skip the "## ..." header line in your file
    dfm = pd.read_csv(lookup_path, comment="#", delim_whitespace=True, header=None)
    dfm = dfm.dropna(axis=1, how="all").copy()

    # numeric conversion
    for c in dfm.columns:
        dfm[c] = pd.to_numeric(dfm[c], errors="coerce")
    dfm = dfm.dropna(how="any").copy()

    # Detect your 7-column neutron_lookup format:
    # Q2, GEn/GD, dGEn/GD, dGEn_Par/GD, GMn/muGD, dGMn/muGD, dGMn_Par/muGD
    if dfm.shape[1] >= 7:
        out = pd.DataFrame({
            "Q2": dfm.iloc[:, 0].astype(float).values,
            "R":  dfm.iloc[:, 4].astype(float).values,
            "dR_stat": dfm.iloc[:, 5].astype(float).values,
            "dR_par":  dfm.iloc[:, 6].astype(float).values,
        })
        out["dR_tot"] = np.sqrt(out["dR_stat"]**2 + out["dR_par"]**2)
        out = out.sort_values("Q2").reset_index(drop=True)
        return out, "ye2018"

    # Generic 4-column: Q2 R dR_stat dR_par
    if dfm.shape[1] >= 4:
        out = pd.DataFrame({
            "Q2": dfm.iloc[:, 0].astype(float).values,
            "R":  dfm.iloc[:, 1].astype(float).values,
            "dR_stat": dfm.iloc[:, 2].astype(float).values,
            "dR_par":  dfm.iloc[:, 3].astype(float).values,
        })
        out["dR_tot"] = np.sqrt(out["dR_stat"]**2 + out["dR_par"]**2)
        out = out.sort_values("Q2").reset_index(drop=True)
        return out, "generic"

    raise ValueError(f"Lookup file {lookup_path} needs >=4 numeric columns; found {dfm.shape[1]}.")


# ---------------- Plotting ----------------

def plot(points_csv: Optional[Path],
         lookup_path: Optional[Path],
         out_prefix: Path,
         dpi: int = 300,
         no_fit: bool = False):

    frames: List[pd.DataFrame] = []

    if points_csv is not None and points_csv.exists():
        frames.append(load_points_csv(points_csv))

    if not frames:
        raise RuntimeError("No data loaded. Provide --points <world_data.csv>.")

    raw = pd.concat(frames, ignore_index=True)
    pts = to_R_table(raw)

    fit = None
    fit_kind = None
    if (not no_fit) and (lookup_path is not None) and lookup_path.exists():
        fit, fit_kind = load_lookup_any(lookup_path)

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(9.8, 6.8))

    # world data points (color by method category, marker by paper)
    color_map = {"polHe3":"blue", "ratio":"red", "lt":"black", "inclusive":"purple", "other":"gray"}
    unique_papers = list(pts["Paper"].astype(str).unique())
    markers = ['o','s','^','D','v','>','<','p','P','X','*','h','H','+','x','d']
    marker_map = {paper: markers[i % len(markers)] for i, paper in enumerate(unique_papers)}

    # Plot each paper in one call so legend handles are clean
    for paper in unique_papers:
        sub = pts[pts["Paper"] == paper].copy()
        if sub.empty:
            continue

        # choose a representative method_cat for color (mode)
        method_cat = sub["method_cat"].mode().iloc[0] if "method_cat" in sub.columns else "other"
        col = color_map.get(method_cat, "gray")

        ax.errorbar(
            sub["Q2"].values,
            sub["R"].values,
            yerr=sub["R_tot"].values,
            fmt=marker_map[paper],
            linestyle="None",
            markersize=6.5,
            capsize=3,
            color=col,
            alpha=0.95,
            label=paper,   # <-- this is the per-paper legend label
        )

    # Ye (or lookup) band
    band = line = None
    if fit is not None:
        # clip to displayed range (avoid Q2<=0 for log axis)
        fitp = fit[(fit["Q2"] > 0.0) & (fit["Q2"] <= 12.0)].copy()

        band = ax.fill_between(
            fitp["Q2"].values,
            (fitp["R"] - fitp["dR_tot"]).values,
            (fitp["R"] + fitp["dR_tot"]).values,
            alpha=0.25,
            label="Global fit (Ye 2018) ± total" if fit_kind == "ye2018" else "Fit ± total"
        )
        (line,) = ax.plot(
            fitp["Q2"].values,
            fitp["R"].values,
            linewidth=2.2,
            color="black",
            label="Global fit (Ye 2018)" if fit_kind == "ye2018" else "Fit"
        )

    # ---- axes ----
    ax.set_xscale("log")                 # <-- log x-axis
    ax.set_xlim(0.01, 12.0)               # log scale cannot include 0
    ax.set_ylim(0.2, 1.25)

    ax.set_xlabel(r"$Q^2\ \mathrm{(GeV/c)^2}$",fontsize=18)
    ax.set_ylabel(r"$G_M^n / (\mu_n\,G_D)$",fontsize=18)
    ax.grid(True, which="both", alpha=0.35)

    # ---- Legends ----
    # (1) Fit legend (bottom-left), like your thesis figure
    handles_fit = []
    handles_fit.append(Line2D([0],[0], marker='o', color='black', linestyle='None', label='World data'))
    if line is not None:
        handles_fit.append(Line2D([0],[0], color='black', linewidth=2.2,
                                  label="Global fit (Ye 2018)" if fit_kind == "ye2018" else "Fit"))
    if band is not None:
        handles_fit.append(band)

    #leg_fit = ax.legend(handles=handles_fit, loc="lower left", fontsize=10, frameon=True)
    #ax.add_artist(leg_fit)

    # (2) Per-paper legend (data points). Put it on the right; split into columns.
    # If there are many papers, ncol=2 keeps it compact.
    leg_papers = ax.legend(
        loc="lower left",
        fontsize=15,
        frameon=True,
        ncol=2,
        title="",
        title_fontsize=9,
        handletextpad=0.4,
        borderpad=0.4,
        labelspacing=0.35
    )

    # ---- Save ----
    out_png = Path(f"{out_prefix}.png")
    out_pdf = Path(f"{out_prefix}.pdf")
    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")

    # Save plotted points (R-space)
    out_points = Path(f"{out_prefix}_points.csv")
    pts_out = pts[["Paper","year","method","Q2","R","R_stat","R_sys","R_tot","y_kind"]].copy()
    pts_out.to_csv(out_points, index=False)

    print(f"Wrote:\n  {out_png}\n  {out_pdf}\n  {out_points}")


def main():
    p = argparse.ArgumentParser(description="Plot G_M^n/(mu_n G_D) world data + Ye(2018) band from neutron_lookup.dat")
    p.add_argument("--points", type=Path, required=True,
                   help="CSV with world data points")
    p.add_argument("--lookup", type=Path, default=None,
                   help="Lookup file (supports neutron_lookup.dat format)")
    p.add_argument("--out-prefix", type=Path, default=Path("GMn_over_muGD"),
                   help="Output prefix")
    p.add_argument("--dpi", type=int, default=300, help="DPI for PNG")
    p.add_argument("--no-fit", action="store_true", help="Disable lookup/fit overlay")
    args = p.parse_args()

    plot(points_csv=args.points,
         lookup_path=args.lookup,
         out_prefix=args.out_prefix,
         dpi=args.dpi,
         no_fit=args.no_fit)

if __name__ == "__main__":
    main()
