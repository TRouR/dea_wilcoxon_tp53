import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch


# Horizontal bar plot of 30 enrichment terms grouped by GO category

# Load data 
csv_path = "top10_terms_each.csv"
df = pd.read_csv(csv_path)

# Basic validation 
required_cols = {"source", "term_name", "negative_log10_of_adjusted_p_value"}
missing = required_cols - set(df.columns)
if missing:
    raise SystemExit(f"Missing required columns in CSV: {missing}. Found: {df.columns.tolist()}")

# Ensure numeric
valcol = "negative_log10_of_adjusted_p_value"
df[valcol] = pd.to_numeric(df[valcol], errors="coerce")

# Drop completely invalid rows (if any)
df = df.dropna(subset=[valcol, "term_name", "source"]).copy()

# Order by category, then value 
cat_order = ["GO:BP", "GO:MF", "GO:CC"]
df["source"] = pd.Categorical(df["source"], categories=cat_order, ordered=True)
df_sorted = df.sort_values(["source", valcol], ascending=[True, False]).reset_index(drop=True)

# Color mapping: 3 base colors, intensity scales with value
base_colors = {
    "GO:BP": "tab:blue",
    "GO:MF": "tab:green",
    "GO:CC": "tab:orange",
}

# Normalize values globally (avoid zero division)
vmin = float(df_sorted[valcol].min())
vmax = float(df_sorted[valcol].max())
rng = vmax - vmin if vmax > vmin else 1.0

def tinted(color_name, value, vmin=vmin, vmax=vmax):
    # strength in [0.4, 1.0] to keep colors readable
    norm = (value - vmin) / (vmax - vmin) if vmax > vmin else 1.0
    strength = 0.4 + 0.6 * norm
    base_rgb = np.array(mcolors.to_rgb(color_name))
    white = np.array([1.0, 1.0, 1.0])
    rgb = white * (1 - strength) + base_rgb * strength
    return tuple(rgb)

colors = [
    tinted(base_colors.get(cat, "tab:gray"), val)
    for cat, val in zip(df_sorted["source"], df_sorted[valcol])
]

# Plot 
fig, ax = plt.subplots(figsize=(11, 12))

ypos = np.arange(len(df_sorted))
ax.barh(ypos, df_sorted[valcol].to_numpy(), color=colors, edgecolor="none")

# y labels (term names)
ax.set_yticks(ypos)
ax.set_yticklabels(df_sorted["term_name"], fontsize=7)

# x label and title
ax.set_xlabel("-log10(adjusted p-value)")
ax.set_title("Functional Enrichment (Gene Ontology)", pad=10)

# Group separators between categories
sources = df_sorted["source"].astype(str).to_list()
break_idxs = [i for i in range(1, len(sources)) if sources[i] != sources[i-1]]
for i in break_idxs:
    ax.axhline(i - 0.5, color="0.85", linewidth=0.8)

# Legend with base colors
legend_handles = [
    Patch(facecolor=base_colors[c], label=c) for c in cat_order if c in set(sources)
]
if legend_handles:
    ax.legend(handles=legend_handles, title="Source", loc="lower right", frameon=False)

# Grid for readability
ax.grid(axis="x", linestyle=":", linewidth=0.8, alpha=0.4)
ax.invert_yaxis()  # top item at top
plt.tight_layout()

# Save a copy
dpi = 600
out_path = f"functional_enrichment_barplot_{dpi}dpi.png"
svg_out = f"functional_enrichment_barplot_{dpi}dpi.svg"
fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
fig.savefig(svg_out, bbox_inches="tight")