import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# -----------------------------------------------------------------------------
# Input
# -----------------------------------------------------------------------------
csv_file = "He3_pol_sampled_per_run_mid.csv"
out_png  = "He3_pol_vs_run_sampled_mid_grayscale_filtered_v2.png"
out_pdf  = "He3_pol_vs_run_sampled_mid_grayscale_filtered_v2.pdf"
out_csv  = "Kin4b_removed_upward_spikes_v2.csv"

# -----------------------------------------------------------------------------
# Read data
# -----------------------------------------------------------------------------
runs = pd.read_csv(csv_file)

# Expected columns:
#   Run_num, Kinematic, P, dP
# Sort once for safety
runs = runs.sort_values("Run_num").reset_index(drop=True)

# -----------------------------------------------------------------------------
# Table values used for horizontal mean band / line
# -----------------------------------------------------------------------------
table_vals = {
    "Kin2":  (34.91, 1.73),
    "Kin3":  (41.76, 1.09),
    "Kin4a": (53.18, 1.92),
    "Kin4b": (36.52, 2.19),
}

# -----------------------------------------------------------------------------
# Grayscale region fills
# -----------------------------------------------------------------------------
kin_shade = {
    "Kin2":  "0.70",
    "Kin3":  "0.83",
    "Kin4a": "0.76",
    "Kin4b": "0.78",
}

# -----------------------------------------------------------------------------
# Remove obvious upward spikes in Kin4b
#
# Logic:
#   1. Compute a local rolling median within Kin4b only.
#   2. Remove points that sit > 8 above the local median AND are > 45 overall.
#   3. Also remove late-end upward jumps for runs >= 5980 with P > 47.
# -----------------------------------------------------------------------------
plot_runs = runs.copy()

k4b_mask = plot_runs["Kinematic"].eq("Kin4b")
k4b = plot_runs.loc[k4b_mask].sort_values("Run_num").reset_index()

local_med = k4b["P"].rolling(window=11, center=True, min_periods=3).median()

spike_mask_local = (
    (k4b["P"] - local_med > 8.0) &
    (k4b["P"] > 45.0)
)

late_jump_mask = (
    (k4b["Run_num"] >= 5980) &
    (k4b["P"] > 47.0)
)

drop_idx = sorted(set(k4b.loc[spike_mask_local | late_jump_mask, "index"].tolist()))

plot_runs = plot_runs.drop(index=drop_idx).reset_index(drop=True)

# Save removed points
removed = runs.loc[drop_idx, ["Run_num", "Kinematic", "P", "dP"]].sort_values("Run_num")
removed.to_csv(out_csv, index=False)

# -----------------------------------------------------------------------------
# Split for broken-axis layout
# -----------------------------------------------------------------------------
left_df  = plot_runs[plot_runs["Kinematic"].isin(["Kin2", "Kin3", "Kin4a"])].copy()
right_df = plot_runs[plot_runs["Kinematic"].eq("Kin4b")].copy()

# Use original full run spans for shaded regions and table bands
kin_run_spans = {}
for kin in ["Kin2", "Kin3", "Kin4a", "Kin4b"]:
    sub = runs[runs["Kinematic"] == kin]
    kin_run_spans[kin] = (sub["Run_num"].min(), sub["Run_num"].max())

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(12.8, 4.9),
    sharey=True,
    gridspec_kw={"width_ratios": [3.8, 1.4]}
)

fig.patch.set_facecolor("white")
for ax in (ax1, ax2):
    ax.set_facecolor("white")


def shade_kin(ax, kin_list, y_top=69.0):
    for kin in kin_list:
        x0, x1 = kin_run_spans[kin]

        ax.axvspan(
            x0, x1,
            alpha=0.45,
            linewidth=0,
            facecolor=kin_shade[kin],
            zorder=0
        )

        ax.axvline(x0, linewidth=2.0, color="k", zorder=1)
        ax.axvline(x1, linewidth=2.0, color="k", zorder=1)

        ax.text(
            (x0 + x1) / 2, y_top,
            kin,
            ha="center", va="top",
            fontsize=10,
            fontweight="bold",
            color="k"
        )


def add_table_band(ax, kin):
    mean, unc = table_vals[kin]
    x0, x1 = kin_run_spans[kin]

    ax.fill_between(
        [x0, x1],
        [mean - unc, mean - unc],
        [mean + unc, mean + unc],
        color="0.65",
        alpha=0.25,
        zorder=2
    )

    ax.plot(
        [x0, x1], [mean, mean],
        linewidth=2.2,
        color="k",
        zorder=2.5
    )


def draw_points(ax, df):
    ax.errorbar(
        df["Run_num"].values,
        df["P"].values,
        yerr=df["dP"].values,
        fmt="o",
        linestyle="none",
        markersize=3.2,
        elinewidth=1.1,
        capsize=3.0,
        capthick=1.1,
        color="0.20",
        ecolor="0.35",
        zorder=3
    )


shade_kin(ax1, ["Kin2", "Kin3", "Kin4a"])
shade_kin(ax2, ["Kin4b"])

for kin in ["Kin2", "Kin3", "Kin4a"]:
    add_table_band(ax1, kin)
add_table_band(ax2, "Kin4b")

draw_points(ax1, left_df)
draw_points(ax2, right_df)

# -----------------------------------------------------------------------------
# Axes styling
# -----------------------------------------------------------------------------
for ax in (ax1, ax2):
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5, color="0.7")
    ax.set_ylim(0, 70)

ax1.set_ylabel("Density Corrected Polarization (%)")
ax1.set_title("GEn-IIa", pad=6)
ax2.set_title("GEn-IIb", pad=6)

ax1.set_xlabel("Run number")
ax2.set_xlabel("Run number")
ax2.tick_params(labelleft=False)

ax1.xaxis.set_major_locator(mticker.MaxNLocator(nbins=8, integer=True))
ax2.xaxis.set_major_locator(mticker.MaxNLocator(nbins=5, integer=True))

# Broken axis appearance
ax1.spines["right"].set_visible(False)
ax2.spines["left"].set_visible(False)

d = 0.012
kwargs = dict(transform=ax1.transAxes, color="k", clip_on=False, linewidth=1.6)
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

kwargs.update(transform=ax2.transAxes)
ax2.plot((-d, +d), (-d, +d), **kwargs)
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)

# -----------------------------------------------------------------------------
# Legend
# -----------------------------------------------------------------------------
kin_patches = [
    Patch(facecolor=kin_shade[k], edgecolor="k", label=k)
    for k in ["Kin2", "Kin3", "Kin4a", "Kin4b"]
]

mean_handle = Line2D([0], [0], color="k", linewidth=2.2, label="⟨P⟩")
band_handle = Patch(facecolor="0.65", alpha=0.25, edgecolor="none", label="±ΔP")
err_handle = ax1.errorbar(
    [], [], yerr=[1],
    fmt="none",
    ecolor="0.35",
    elinewidth=1.1,
    capsize=3.0,
    capthick=1.1,
    label="Error"
)

ax1.legend(
    handles=kin_patches + [mean_handle, band_handle, err_handle],
    loc="lower left",
    fontsize=8,
    framealpha=0.95,
    handlelength=1.8
)

# -----------------------------------------------------------------------------
# Global caption
# -----------------------------------------------------------------------------
#fig.text(
#    0.5, 0.02,
#    "Polarization",
#    ha="center", va="bottom",
#    fontsize=12
#)

plt.tight_layout(rect=(0, 0.05, 1, 1))
plt.savefig(out_png, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.savefig(out_pdf, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.show()

print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
print(f"Saved removed-points table: {out_csv}")
print("\nRemoved Kin4b upward spikes:")
print(removed.to_string(index=False))
