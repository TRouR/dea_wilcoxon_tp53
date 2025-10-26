import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Resolve the directory
try:
    SCRIPT_DIR = Path(__file__).resolve().parent
except NameError:
    SCRIPT_DIR = Path.cwd()

# Parameters
alpha = 0.05          # FDR threshold
log2fc_min = 0.5     # effect-size guard (median log2 diff threshold)
group_a = "WT"        # reference
group_b = "TP53-mut"  # case
csv_path = SCRIPT_DIR / "DEA_results_tails.csv"
out_path = SCRIPT_DIR / f"volcano_{group_b}_vs_{group_a}.png"
out_path_svg = SCRIPT_DIR / f"volcano_{group_b}_vs_{group_a}.svg"

# Load DE results
df = pd.read_csv(csv_path)

# Check required columns
required = {"gene", "med_diff", "p_two", "padj_two"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Input must contain columns: {', '.join(sorted(missing))}")

plot_df = df.copy()

# Prepare plotting columns
eps = np.finfo(float).tiny  # avoid -log10(0)
plot_df["neglog10_p"] = -np.log10(np.maximum(plot_df["p_two"].values, eps))  # raw p on y
plot_df["log2FC"] = plot_df["med_diff"]

# Masks: significance (FDR) + side
sig = (plot_df["padj_two"] < alpha) & (plot_df["log2FC"].abs() >= log2fc_min)
left_sig  = sig & (plot_df["log2FC"] <= -log2fc_min)   # under-expressed in TP53-mut
right_sig = sig & (plot_df["log2FC"] >=  log2fc_min)   # over-expressed  in TP53-mut
ns_mask   = ~(left_sig | right_sig)

n_left  = int(left_sig.sum())
n_right = int(right_sig.sum())
n_ns    = int(ns_mask.sum())

# Colors (more saturated)
COL_NS    = "#6f6f6f"   # grey
COL_LEFT  = "#2b8cbe"   # blue   (under)
COL_RIGHT = "#e45757"   # red    (over)
ALPHA_NS  = 0.6
ALPHA_SIG = 0.9
SIZE_NS   = 12
SIZE_SIG  = 14

# Plot
fig, ax = plt.subplots(figsize=(8, 6))

# Non-significant (grey)
ax.scatter(plot_df.loc[ns_mask, "log2FC"],
           plot_df.loc[ns_mask, "neglog10_p"],
           c=COL_NS, alpha=ALPHA_NS, s=SIZE_NS, edgecolors="none",
           label=f"Not significant (n={n_ns})")

# Significant left (Under-expressed in TP53-mut)
ax.scatter(plot_df.loc[left_sig, "log2FC"],
           plot_df.loc[left_sig, "neglog10_p"],
           c=COL_LEFT, alpha=ALPHA_SIG, s=SIZE_SIG, edgecolors="none",
           label=f"Under-expressed in {group_b} (n={n_left})")

# Significant right (Over-expressed in TP53-mut)
ax.scatter(plot_df.loc[right_sig, "log2FC"],
           plot_df.loc[right_sig, "neglog10_p"],
           c=COL_RIGHT, alpha=ALPHA_SIG, s=SIZE_SIG, edgecolors="none",
           label=f"Over-expressed in {group_b} (n={n_right})")

# Threshold guide lines (no legend entry for p-line)
guide_color = "black"
p_y = -np.log10(alpha)
ax.axhline(p_y, linestyle="--", linewidth=0.8, color=guide_color)  # no label
ax.axvline(log2fc_min,  linestyle="--", linewidth=0.8, color=guide_color)
ax.axvline(-log2fc_min, linestyle="--", linewidth=0.8, color=guide_color)

# Place "p = α" text just above the horizontal line near the right
x_min, x_max = ax.get_xlim()
y_min, y_max = ax.get_ylim()
y_offset = 0.02 * (y_max - y_min)
ax.text(x_max, p_y + y_offset, f"p = {alpha}",
        ha="right", va="bottom", fontsize=10, color=guide_color, clip_on=False)

# Labels
ax.set_xlabel(f"log2 Fold Change ({group_b} - {group_a})")
ax.set_ylabel("-log10(two-sided p-value) [Wilcoxon]")
ax.set_title(f"{group_b} vs {group_a}")

# Compact legend (only the three point groups)
leg = ax.legend(frameon=True, loc="upper left", markerscale=1.0)
for lh in leg.legend_handles:
    lh.set_alpha(1.0)

plt.tight_layout()
fig.savefig(out_path, dpi=600)
fig.savefig(out_path_svg, dpi=600)
