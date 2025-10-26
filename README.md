# TP53-mut vs WT in LUAD — Differential Gene Expression + Enrichment Analyses (Gitools-replicated pipeline)

A compact, reproducible pipeline to:

1) Run per-gene **Mann–Whitney–Wilcoxon** differential expression (DGE) in **R** with **BH FDR** (mirrors Gitools “Group Comparison”), using a log2 expression matrix from the **Gitools LUAD** source heatmap dataset,  
2) Visualize results (volcano plot, top-genes heatmap) and **GO enrichment**.

---

## What’s inside

### R / Differential Expression
- **`DiffExprAnalysis.R`** – Wilcoxon rank-sum per gene, left/right/two-sided p, **BH over tested genes**, filters: ≥3 non-NA per group & non-constant, writes `DEA_results_tails.csv`.  
  *Also includes a commented block showing how to run BH over the **full universe** by padding untestable genes with p=1.*

### Python / Data prep & figures
- **`data_processing.py`** – converts the Gitools `.tdm` + `Target.csv` TP53 status into a **genes×samples** matrix (`exprs.csv`) with columns renamed `Mut*` and `Norm*`.  
- **`volcano.py`** – volcano plot using **two-sided p** (y) and **median log2 fold change** (TP53-mut − WT) (x), highlights FDR<0.05 & |log2FC|≥0.5.  
- **`heatmap.py`** – selects top-100 DE genes by |log2FC| within FDR<0.05, **kNN-imputes** missing log2 values *for visualization only*, and draws a clustered heatmap with columns ordered **WT → Mut**. Saves image and the imputed matrix.  
- **`barplot.py`** – horizontal bar chart for **top 10 GO terms per ontology (BP/MF/CC)** using −log10(FDR).

### Data 

- Gitools LUAD heatmap export (log2 expression layer): `expression_data.tdm`.  
- Sample status table: **`Target.csv`** (rows = sample IDs as in TDM, one column `Mut` with 1 = TP53-mut, 0 = WT).

> The pipeline expects **log2-normalized** expression (as provided by the Gitools heatmap).

---

## Quick start

### 0) Environment

- **R** ≥ 4.0 (uses `stats::wilcox.test`, `p.adjust`, and `matrixStats`).  
- **Python** ≥ 3.9 with: `pandas`, `numpy`, `matplotlib`, `seaborn`, `scikit-learn`.

```bash
pip install pandas numpy matplotlib seaborn scikit-learn
```

### 1) Build the expression matrix

```bash
python data_processing.py
# → writes exprs.csv with columns: Mut1..MutN, Norm1..NormM
```

### 2) Differential expression in R

```r
# edit setwd() inside the script or run from the repo root
source("DiffExprAnalysis.R")
# → writes DEA_results_tails.csv with:
#   gene, n_wt, n_mut, med_wt, med_mut, med_diff, auc,
#   p_left, padj_left, p_right, padj_right, p_two, padj_two
```

**Method details (implemented in R script):**

- **Test**: two-sample Mann–Whitney–Wilcoxon.  
- **Per-gene NA handling**: NAs **dropped within each group**, a gene is tested only if **n_wt ≥ 3** and **n_mut ≥ 3**, and both groups are **non-constant**.  
- **Effect size**:  
  - `med_diff` = median(Mut) − median(WT) (log2 scale).  
  - `auc` = U / (n₁·n₂) ≈ P(Mut > WT).  
- **FDR**: default is **BH over tested genes**. Optional block demonstrates BH over **all genes** by padding untestable rows with p=1.

### 3) Figures

```bash
# Volcano plot (two-sided p, FDR<0.05 & |log2FC|≥0.5)
python volcano.py
# → volcano_TP53-mut_vs_WT.png/.svg

# Heatmap of top-100 DE genes (FDR<0.05, |log2FC|≥0.5)
python heatmap.py
# → heatmap_top100_log2.png/.svg, top100_expression_imputed.csv

# GO bar plot (expects your prepared enrichment CSV with top 10 BP/MF/CC)
python barplot.py
# → functional_enrichment_barplot.png/.svg
```

---

## File specs

- **`exprs.csv`** – genes×samples matrix (log2) produced from the Gitools TDM by `data_processing.py`. Mutated samples = `Mut*`, non-mutated = `Norm*`.  
- **`DEA_results_tails.csv`** – per-gene DE results from R (left/right/two-sided p and BH-FDR).  
- **`volcano_*.png/.svg`** – volcano plot showing over/under-expressed genes in TP53-mut.  
- **`top100_expression_imputed.csv`** – imputed log2 matrix for the selected top genes (**visualization only**).  
- **`heatmap_top100_log2*.png/.svg`** – clustered heatmap with columns ordered WT → Mut.  
- **`functional_enrichment_barplot*.png/.svg`** – grouped bar chart for GO:BP/MF/CC.

---

## Reproducing Gitools behavior

- **Same statistical test**: two-sample MWW, per-gene NAs are discarded.  
- **Same eligibility rule**: gene tested only if **≥3** non-missing values in **each** group.  
- **Outputs**: left (`Mut < WT`), right (`Mut > WT`), and two-sided p-values.  
- **R vs Gitools differences**: expect **minor numeric jitter** from implementation details (NA/constant-row filtering, tie/continuity corrections in Wilcoxon, two-sided construction, and BH ordering).  
- **Why “BH over tested genes” (default)?** It controls FDR on the set of **actual decisions** and avoids over-penalizing by including untestable genes.

---

## Parameter notes

- **Significance**: default **FDR < 0.05** (two-sided) and **|log2FC| ≥ 0.5** (median-based: ≈1.41× change).  
- **Heatmap imputation**: kNN (k=5, nan-euclidean) on **log2** values **after** DE filtering. Statistics are computed **without** imputation.  
- **Volcano**: y = −log10(two-sided p), x = log2FC (median difference), with FDR and effect cutoffs annotated.

---

