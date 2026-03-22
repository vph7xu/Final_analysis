import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("02010954.dat", sep=r"\s+", header=None, names=["x", "freq", "signal"])
ticks = np.arange(len(df))

# Detect sharp transitions
jump_idx = np.where(df["freq"].diff().abs().fillna(0) > 20)[0]

fig, ax = plt.subplots(figsize=(10, 4.75), dpi=200)

# Background colors in grayscale
#fig.patch.set_facecolor("#d9d9d9")
#ax.set_facecolor("#efefef")

# Highlight transition regions
band_halfwidth = 5
for j in jump_idx:
    ax.axvspan(j - band_halfwidth, j + band_halfwidth,
               color="0.75", alpha=0.7, zorder=0)

# Main curve and points
ax.plot(ticks, df["freq"], color="0.15", linewidth=1.6, zorder=2)
ax.scatter(ticks, df["freq"], s=10, color="0.35",
           edgecolors="none", zorder=3, label="Frequency Data")

# Labels and title
ax.set_title(r"$^{39}$K EPR AFP Spectrum", fontsize=15, pad=10)
ax.set_xlabel("Ticks", fontsize=12)
ax.set_ylabel("Frequency (kHz)", fontsize=12)

# Text annotations
ax.text(14, 18909.8, r"$^{3}$He anti-aligned", fontsize=11, color="0.1")
ax.text(76, 18867.2, r"$^{3}$He aligned", fontsize=11, color="0.1")
ax.text(51, 18884.3, r"$\Delta \nu_1$", fontsize=15, color="0.1")
ax.text(211, 18884.0, r"$\Delta \nu_2$", fontsize=15, color="0.1")

# Limits and legend
ax.set_ylim(18852, 18921.5)
ax.legend(loc="lower left", frameon=False, fontsize=11)

plt.tight_layout()
plt.savefig("02010954_grayscale_plot.pdf", dpi=300,
            facecolor=fig.get_facecolor(), bbox_inches="tight")
plt.show()
