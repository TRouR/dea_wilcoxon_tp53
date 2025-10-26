# ===============================================================
# Mann–Whitney–Wilcoxon DGE
# - Input:  exprs.csv  (first column = gene IDs, rest = samples)
# - Groups: inferred from column names (Mut/WT).
# - Filters: per gene, per group: >=3 non-NA and non-constant
# - P-values: left (Mut < WT), right (Mut > WT), two-sided
# - BH FDR: only tested genes
# - Output: DEA_results_tails.csv
# ===============================================================

library(matrixStats)

# Load expression matrix (first column = gene IDs)
setwd("...")
expr_df <- read.csv("exprs.csv", check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_df) <- expr_df[[1]]
expr_df <- expr_df[-1]
expr <- as.matrix(apply(expr_df, 2, function(x) suppressWarnings(as.numeric(x))))
rownames(expr) <- rownames(expr_df)

# Define groups
cn <- colnames(expr)
is_mut <- grepl("(?i)mut|tp53mut|tp53_mut", cn, perl = TRUE)
is_wt  <- grepl("(?i)^(wt|norm|control|ctrl)|wt|norm|control|ctrl", cn, perl = TRUE)
group <- ifelse(is_mut, "Mut", ifelse(is_wt, "WT", NA))
keep_cols <- !is.na(group)
expr  <- expr[, keep_cols, drop = FALSE]

group <- factor(group[keep_cols], levels = c("WT","Mut"))

ixW <- which(group == "WT")
ixM <- which(group == "Mut")
stopifnot(length(ixW) >= 3, length(ixM) >= 3)

# Filters (per gene, per group)
min_n <- 3
nW <- rowSums(!is.na(expr[, ixW, drop = FALSE]))
nM <- rowSums(!is.na(expr[, ixM, drop = FALSE]))
not_const <- function(v) length(unique(stats::na.omit(v))) > 1
varW <- apply(expr[, ixW, drop = FALSE], 1, not_const)
varM <- apply(expr[, ixM, drop = FALSE], 1, not_const)
keep_genes <- (nW >= min_n) & (nM >= min_n) & varW & varM
expr_f <- expr[keep_genes, , drop = FALSE]
sum(vapply(expr_df, function(x) all(is.na(x)), logical(1)))


# Wilcoxon per gene: left/right/two-sided
wilcox_three <- function(v) {
  a <- v[ixW]; b <- v[ixM]
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  # two-sided
  p_two <- suppressWarnings(wilcox.test(a, b, alternative="two.sided", exact=FALSE)$p.value)
  # Right tail: Mut > WT  (equivalently a<b)
  p_right <- suppressWarnings(wilcox.test(b, a, alternative="greater", exact=FALSE)$p.value)
  # Left tail: Mut < WT   (equivalently a>b)
  p_left  <- suppressWarnings(wilcox.test(b, a, alternative="less",    exact=FALSE)$p.value)
  # effects
  n1 <- length(a); n2 <- length(b)
  W  <- suppressWarnings(wilcox.test(a, b, alternative="two.sided", exact=FALSE)$statistic)
  U  <- as.numeric(W) - n2*(n2+1)/2
  auc <- U / (n1*n2)
  med_wt  <- median(a); med_mut <- median(b)
  c(p_left=p_left, p_right=p_right, p_two=p_two,
    med_wt=med_wt, med_mut=med_mut,
    med_diff=med_mut - med_wt,
    auc=auc)
}

res <- t(apply(expr_f, 1, wilcox_three))
res <- as.data.frame(res)

# BH per tail (only across tested genes)
res$padj_left  <- p.adjust(res$p_left,  method = "BH")
res$padj_right <- p.adjust(res$p_right, method = "BH")
res$padj_two   <- p.adjust(res$p_two,   method = "BH")

# # BH with FULL UNIVERSE padding
# # p-vectors for the genes you tested
# p_left_named  <- setNames(res$p_left,  rownames(expr_f))
# p_right_named <- setNames(res$p_right, rownames(expr_f))
# p_two_named   <- setNames(res$p_two,   rownames(expr_f))
# 
# # helper: fill missing genes with p=1 so BH runs over ALL genes
# pad_vec <- function(p_named, all_ids) {
#   v <- rep(1, length(all_ids))
#   names(v) <- all_ids
#   overlap <- intersect(names(p_named), all_ids)
#   v[overlap] <- p_named[overlap]
#   v
# }
# 
# pL_all <- pad_vec(p_left_named,  all_genes)
# pR_all <- pad_vec(p_right_named, all_genes)
# pT_all <- pad_vec(p_two_named,   all_genes)
# 
# padj_left_all  <- p.adjust(pL_all, method = "BH")
# padj_right_all <- p.adjust(pR_all, method = "BH")
# padj_two_all   <- p.adjust(pT_all, method = "BH")
# 
# # map adjusted p's back to tested genes
# res$padj_left  <- padj_left_all [rownames(expr_f)]
# res$padj_right <- padj_right_all[rownames(expr_f)]
# res$padj_two   <- padj_two_all  [rownames(expr_f)]

# Tidy output
out <- data.frame(
  gene = rownames(expr_f),
  n_wt = nW[keep_genes],
  n_mut = nM[keep_genes],
  med_wt = res$med_wt, med_mut = res$med_mut, med_diff = res$med_diff,
  auc = res$auc,
  p_left = res$p_left,  padj_left  = res$padj_left,
  p_right = res$p_right, padj_right = res$padj_right,
  p_two = res$p_two,     padj_two   = res$padj_two,
  check.names = FALSE
)
out <- out[order(out$padj_two, out$p_two), ]

write.csv(out, "DEA_results_tails.csv", row.names = FALSE)
