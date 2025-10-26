from pathlib import Path
import pandas as pd
from sklearn.impute import KNNImputer 
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Resolve the directory
try:
    SCRIPT_DIR = Path(__file__).resolve().parent
except NameError:
    SCRIPT_DIR = Path.cwd()
 
in_path  = SCRIPT_DIR / "DEA_results_tails.csv"
out_filtered = SCRIPT_DIR / "filtered_DEA_genes_.csv"
out_ranked = SCRIPT_DIR / "filtered_ranked_DEA_genes.csv"

PADJ =  0.05
LFC  =  0.5

df = pd.read_csv(in_path)

# basic sanity check
if not {"med_diff", "padj_two"}.issubset(df.columns):
    raise SystemExit("ERROR: required columns 'med_diff' and 'padj_two' not found in input.")

# Use med_diff as log2FC
df["log2FC"] = df["med_diff"]

# Filter
filt = (df["padj_two"] < PADJ) & (df["log2FC"].abs() >= LFC)
out = df.loc[filt].copy()

# Rank
out["abs_log2FC"] = out["log2FC"].abs()
out.sort_values(["abs_log2FC"], ascending=[False], inplace=True)
out.to_csv(out_filtered, index=False)

# SUBSET: keep only the top 100 by biggest abs_log2FC
top100 = out.nlargest(100, "abs_log2FC").copy()

# Present them sorted by abs_log2FC descending
top100.sort_values("abs_log2FC", ascending=False, inplace=True)
top100.to_csv(out_ranked, index=False)

# Read expression matrix and subset TOP 100 genes
expr_path = SCRIPT_DIR / "exprs.csv"
expr = pd.read_csv(expr_path)

# Use a gene identifier column if present, else first column as index
gene_cols = ["gene","Gene","symbol","Symbol","gene_id","GeneID","ENSEMBL","ensembl","id","ID","Unnamed: 0"]
expr_gene_col = next((c for c in gene_cols if c in expr.columns), None)
if expr_gene_col is not None:
    expr = expr.set_index(expr_gene_col)
else:
    expr = expr.set_index(expr.columns[0])

# Your DE file has a 'gene' column, if not, adjust de_gene_col accordingly
de_gene_col = "gene" if "gene" in top100.columns else next((c for c in gene_cols if c in top100.columns), None)

if de_gene_col is None:
    raise SystemExit("Could not find a gene ID column in the DE table (e.g. 'gene').")

top100_genes = top100[de_gene_col].astype(str).tolist()

# Subset and keep the same order as in top100
top100_exprs = expr.loc[expr.index.intersection(top100_genes)].copy()
top100_exprs = top100_exprs.reindex(top100_genes)

# kNN IMPUTATION
# Ensure numeric for imputation
top100_exprs = top100_exprs.apply(pd.to_numeric, errors='coerce')

# kNN impute across samples (treat rows=samples for sklearn, so transpose)
imputer = KNNImputer(n_neighbors=5, weights='distance', metric='nan_euclidean')
X_T = top100_exprs.T.values  # samples x genes

X_imp_T = imputer.fit_transform(X_T)
top100_exprs_imputed = pd.DataFrame(X_imp_T, index=top100_exprs.T.index, columns=top100_exprs.T.columns).T

# Save imputed matrix for heatmap
(top100_exprs_imputed).to_csv(SCRIPT_DIR / "top100_expression_imputed.csv")

# Save original (with NAs) as well if you still want it
top100_exprs.to_csv(SCRIPT_DIR / "top100_expression.csv")

# HEATMAP (log2 values) 
X = top100_exprs_imputed  # imputed log2-normalized expression

# Detect groups from column names
cols = X.columns.astype(str)
is_mut = cols.str.match(r'(?i)^mut')
is_wt  = cols.str.match(r'(?i)^(wt|norm|control|ctrl)')

group = np.where(is_mut, "Mut", np.where(is_wt, "WT", "Unknown"))
group_colors = {"WT": "#6BA3D6", "Mut": "#D66B6B"}

# reorder columns: WT/Norm first, then Mut
order_idx = np.r_[np.where(group == "WT")[0],
                  np.where(group == "Mut")[0]]
X = X.iloc[:, order_idx]
group = group[order_idx]
col_colors = [group_colors[g] for g in group]

out_png = SCRIPT_DIR / "heatmap_top100_log2_1.png"
out_svg = SCRIPT_DIR / "heatmap_top100_log2_1.svg"

import seaborn as sns
sns.set(context="notebook")
g = sns.clustermap(
    X,
    cmap="RdBu_r",
    center=0.0,
    col_colors=col_colors,
    xticklabels=False,           
    yticklabels=False,            
    figsize=(12, 10),
    metric="euclidean",
    method="average",
    col_cluster=False,
    row_cluster=True
)

# legend for the top color bar
for lab, color in group_colors.items():
    g.ax_col_dendrogram.bar(0, 0, color=color, label=lab, linewidth=0)
g.ax_col_dendrogram.legend(title="Group", loc="center", ncols=len(group_colors))

# vertical separator between WT and Mut
n_wt = (group == "WT").sum()
n_mut = (group == "Mut").sum()
if 0 < n_wt < X.shape[1]:
    g.ax_heatmap.axvline(n_wt - 0.5, linestyle="--", linewidth=1, color="k", alpha=0.6)
if 0 < (n_wt + n_mut) < X.shape[1]:  # separator before Unknowns (if any)
    g.ax_heatmap.axvline(n_wt + n_mut - 0.5, linestyle="--", linewidth=1, color="k", alpha=0.4)

g.fig.suptitle("Top 100 Differentially Expressed Genes", y=1.02, fontsize=20)
g.savefig(out_png, dpi=150, bbox_inches="tight")
g.savefig(out_svg, dpi=150, bbox_inches="tight")

plt.close(g.fig)
print(f"Saved heatmap to {out_png}")