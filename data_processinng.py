import pandas as pd
from pathlib import Path


# path to your wide TDM
tdm_path = Path("expression_data.tdm")   # change if needed

# 1) read as strings so we can clean safely; treat empty fields as NaN
df = pd.read_csv(
    tdm_path,
    sep="\t",
    index_col=0,           
    dtype=str,             
    keep_default_na=True,
    na_values=[""]         
)

# 2) tidy headers (strip quotes/whitespace)
df.columns = df.columns.str.strip().str.strip('"')
df.index   = df.index.str.strip().str.strip('"')

# 3) convert all data cells to numeric (non-numeric -> NaN)
df = df.apply(pd.to_numeric, errors="coerce")
df.dropna(axis=1, how='all', inplace=True)

# --- read target (assumes a column 'Mut' with 1/0 and the sample IDs in the first column or a column named 'Sample') ---
t = pd.read_csv("Target.csv", dtype=str)

# sample column: first column if unnamed, else a column named like 'Sample'
sample_col = t.columns[0] if t.columns[0].startswith(("Unnamed", "")) or t.columns[0] == "" else (
    "Sample" if "Sample" in t.columns else t.columns[0]
)

mut_col = "Mut"  # per your format

# normalize types
t[mut_col] = t[mut_col].astype(str).str.strip().map({"1": 1, "0": 0}).astype(int)
t[sample_col] = t[sample_col].astype(str).str.strip()

# keep only samples present in df (and de-duplicate, preserving order)
samples_in_df = set(df.columns.astype(str))
t = t[t[sample_col].isin(samples_in_df)].drop_duplicates(subset=[sample_col], keep="first")

# lists in the order they appear in Target.csv
mut_samples  = t.loc[t[mut_col] == 1, sample_col].tolist()
norm_samples = t.loc[t[mut_col] == 0, sample_col].tolist()

# build rename map (only for samples present in df)
rename_map = {s: f"Mut{i+1}" for i, s in enumerate(mut_samples)}
rename_map.update({s: f"Norm{i+1}" for i, s in enumerate(norm_samples)})

# apply (reorder: Mut* first, then Norm*)
ordered_cols = [s for s in mut_samples + norm_samples if s in df.columns]
df_renamed = df.reindex(columns=ordered_cols).rename(columns=rename_map)

# Save
df_renamed.to_csv("exprs.csv")
print(f"Mutated: {len(mut_samples)} | Not mutated: {len(norm_samples)}")
print(df_renamed.iloc[:5, :8])  # preview first 5 genes × first 8 samples


