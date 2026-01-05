import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

pol_path = "He3_pol.csv"
runs_path = "Copy of GEn Good Runs - All Good Runs.csv"

kin_ranges = [('Kin2', '2022-10-17', '2022-10-30'), ('Kin3', '2022-11-11', '2022-12-14'), ('Kin4a', '2023-01-15', '2023-03-13'), ('Kin4b', '2023-09-16', '2023-10-30')]

table_vals = {'Kin2': (34.91, 1.73), 'Kin3': (41.76, 1.09), 'Kin4a': (53.18, 1.92), 'Kin4b': (36.52, 2.19)}

kin_color = {'Kin2': '#4C6FFF', 'Kin3': '#FFD400', 'Kin4a': '#5CFF40', 'Kin4b': '#FF66FF'}
pt_color = "#FF00CC"

# Choose one: "start", "mid", "end"
sample_mode = "mid"

dfp = pd.read_csv(pol_path)
dfp["Interpolated Dates"] = pd.to_datetime(dfp["Interpolated Dates"], errors="coerce")
dfp["Interpolated Data"] = pd.to_numeric(dfp["Interpolated Data"], errors="coerce")
dfp["Interpolated Errors"] = pd.to_numeric(dfp["Interpolated Errors"], errors="coerce")
dfp = dfp.dropna(subset=["Interpolated Dates","Interpolated Data","Interpolated Errors"]).copy()
dfp = dfp.sort_values("Interpolated Dates").set_index("Interpolated Dates")

t_pol = dfp.index.view("int64") / 1e9
p_pol = dfp["Interpolated Data"].values.astype(float)
e_pol = dfp["Interpolated Errors"].values.astype(float)

raw = pd.read_csv(runs_path, header=None).rename(columns={0:"Run",1:"Date",2:"Start",3:"End"})
raw["Run_num"] = pd.to_numeric(raw["Run"], errors="coerce")
runs = raw.dropna(subset=["Run_num","Date","Start","End"]).copy()
runs["Run_num"] = runs["Run_num"].astype(int)

runs["Start_dt"] = pd.to_datetime(runs["Date"].astype(str) + " " + runs["Start"].astype(str), errors="coerce")
runs["End_dt"]   = pd.to_datetime(runs["Date"].astype(str) + " " + runs["End"].astype(str), errors="coerce")
runs = runs.dropna(subset=["Start_dt","End_dt"]).copy()
roll = runs["End_dt"] < runs["Start_dt"]
runs.loc[roll, "End_dt"] = runs.loc[roll, "End_dt"] + pd.Timedelta(days=1)
runs = runs.sort_values("Start_dt").reset_index(drop=True)

def assign_kin(start_dt):
    for kin, s, e in kin_ranges:
        sdt = pd.Timestamp(s)
        edt = pd.Timestamp(e) + pd.Timedelta(days=1) - pd.Timedelta(seconds=1)
        if (start_dt >= sdt) and (start_dt <= edt):
            return kin
    return None

runs["Kinematic"] = runs["Start_dt"].apply(assign_kin)
runs = runs.dropna(subset=["Kinematic"]).copy()

if sample_mode == "start":
    runs["Sample_dt"] = runs["Start_dt"]
elif sample_mode == "end":
    runs["Sample_dt"] = runs["End_dt"]
else:
    runs["Sample_dt"] = runs["Start_dt"] + (runs["End_dt"] - runs["Start_dt"]) / 2

t_samp = runs["Sample_dt"].view("int64") / 1e9
tmin, tmax = t_pol.min(), t_pol.max()
runs = runs.loc[(t_samp >= tmin) & (t_samp <= tmax)].copy()
t_samp = runs["Sample_dt"].view("int64") / 1e9

runs["P"] = np.interp(t_samp, t_pol, p_pol)
runs["dP"] = np.interp(t_samp, t_pol, e_pol)

runs[["Run_num","Kinematic","Start_dt","End_dt","Sample_dt","P","dP"]].to_csv(
    f"He3_pol_sampled_per_run_{sample_mode}.csv", index=False
)

left_df = runs[runs["Kinematic"].isin(["Kin2","Kin3","Kin4a"])].copy()
right_df = runs[runs["Kinematic"].isin(["Kin4b"])].copy()

kin_run_spans = {}
for kin, _, _ in kin_ranges:
    sub = runs[runs["Kinematic"] == kin]
    kin_run_spans[kin] = (sub["Run_num"].min(), sub["Run_num"].max()) if len(sub) else (np.nan, np.nan)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.9), sharey=True,
                               gridspec_kw={"width_ratios":[3.8, 1.4]})

def shade_kin(ax, kin_list, y_top=69.0):
    for kin in kin_list:
        x0, x1 = kin_run_spans.get(kin, (np.nan, np.nan))
        if not (np.isfinite(x0) and np.isfinite(x1)):
            continue
        ax.axvspan(x0, x1, alpha=0.20, linewidth=0, facecolor=kin_color.get(kin, "#eeeeee"), zorder=0)
        ax.axvline(x0, linewidth=2.2, color="k", zorder=1)
        ax.axvline(x1, linewidth=2.2, color="k", zorder=1)
        ax.text((x0+x1)/2, y_top, kin, ha="center", va="top", fontsize=9, fontweight="bold")

shade_kin(ax1, ["Kin2","Kin3","Kin4a"])
shade_kin(ax2, ["Kin4b"])

def draw_points(ax, df):
    ax.errorbar(df["Run_num"].values, df["P"].values, yerr=df["dP"].values,
                fmt="o", linestyle="none",
                markersize=3.0, elinewidth=1.2, capsize=3.0, capthick=1.2,
                color=pt_color, ecolor=pt_color, zorder=3)

for ax, sub in [(ax1, left_df), (ax2, right_df)]:
    draw_points(ax, sub)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
    ax.set_ylim(0, 70)

def add_table_band(ax, kin):
    mean, unc = table_vals[kin]
    x0, x1 = kin_run_spans.get(kin, (np.nan, np.nan))
    if not (np.isfinite(x0) and np.isfinite(x1)):
        return
    ax.fill_between([x0, x1], [mean-unc, mean-unc], [mean+unc, mean+unc], alpha=0.15, zorder=2)
    ax.plot([x0, x1], [mean, mean], linewidth=2.2, color="k", zorder=2.5)

for kin in ["Kin2","Kin3","Kin4a"]:
    add_table_band(ax1, kin)
add_table_band(ax2, "Kin4b")

ax1.set_ylabel("Density Corrected Polarization (%)")
ax1.set_title("GEn-IIa", pad=6)
ax2.set_title("GEn-IIb", pad=6)
ax2.tick_params(labelleft=False)
ax1.set_xlabel("Run number")
ax2.set_xlabel("Run number")
ax1.xaxis.set_major_locator(mticker.MaxNLocator(nbins=8, integer=True))
ax2.xaxis.set_major_locator(mticker.MaxNLocator(nbins=5, integer=True))

ax1.spines["right"].set_visible(False)
ax2.spines["left"].set_visible(False)
d = 0.012
kwargs = dict(transform=ax1.transAxes, color="k", clip_on=False, linewidth=1.6)
ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
kwargs.update(transform=ax2.transAxes)
ax2.plot((-d, +d), (-d, +d), **kwargs)
ax2.plot((-d, +d), (1-d, 1+d), **kwargs)

kin_patches = [Patch(facecolor=kin_color[k], edgecolor="k", label=k) for k in ["Kin2","Kin3","Kin4a","Kin4b"]]
mean_handle = Line2D([0],[0], color="k", linewidth=2.2, label="⟨P⟩ (table)")
band_handle = Patch(facecolor="k", alpha=0.15, edgecolor="none", label="±ΔP (table)")
err_handle = ax1.errorbar([], [], yerr=[1], fmt="none", ecolor="k", elinewidth=1.2, capsize=3.0, capthick=1.2, label="Error")
ax1.legend(handles=kin_patches + [mean_handle, band_handle, err_handle],
           loc="lower left", fontsize=8, framealpha=0.95, handlelength=1.8)

fig.text(0.5, 0.02, f"Polarization sampled at run {sample_mode}point (no averaging)", ha="center", va="bottom", fontsize=12)
fig.tight_layout(rect=(0, 0.05, 1, 1))
fig.savefig(f"He3_pol_vs_run_sampled_{sample_mode}.png", dpi=300)
fig.savefig(f"He3_pol_vs_run_sampled_{sample_mode}.pdf")
print("Wrote plots + CSV for sample_mode =", sample_mode)
