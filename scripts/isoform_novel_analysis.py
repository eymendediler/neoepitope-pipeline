import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path("/Users/eymen/Desktop/v105")

MERGED = BASE / "mhcflurry_vs_netmhcpan_8mer_merged.csv"
NOVEL = BASE / "novel_peptides_enriched.csv"
NOVEL_TISSUE = BASE / "novel_peptides_enriched_with_tissue.csv"
DOMINANT = BASE / "gtex_most_dominant_per_gene.tsv"

print("Loading files...")
df = pd.read_csv(MERGED)
novel = pd.read_csv(NOVEL)
novel_tissue = pd.read_csv(NOVEL_TISSUE)
dom = pd.read_csv(DOMINANT, sep="\t")

print("merged:", df.shape)
print("novel:", novel.shape)
print("novel_tissue:", novel_tissue.shape)
print("dominant:", dom.shape)

# -----------------------------
# Basic cleaning
# -----------------------------
for x in [df, novel, novel_tissue]:
    for col in ["peptide", "transcript_id", "gene_id", "gene_name"]:
        if col in x.columns:
            x[col] = x[col].astype(str).str.strip()

    if "peptide" in x.columns:
        x["peptide"] = x["peptide"].str.upper()

for col in ["mhcflurry_rank", "ba_rank"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df["rank_gap"] = df["ba_rank"] - df["mhcflurry_rank"]
df["abs_rank_gap"] = df["rank_gap"].abs()

dom.columns = ["gene_id", "most_tx"]
dom["gene_id"] = dom["gene_id"].astype(str).str.strip()
dom["most_tx"] = dom["most_tx"].astype(str).str.strip()

# -----------------------------
# 1. Novel catalog labels
# -----------------------------
novel_peptides = set(novel["peptide"].dropna().unique())
df["is_novel_catalog"] = df["peptide"].isin(novel_peptides)

# transcript count inside novel catalog
novel_pep_tx = (
    novel.groupby("peptide")["transcript_id"]
    .nunique()
    .reset_index(name="novel_transcript_count")
)
df = df.merge(novel_pep_tx, on="peptide", how="left")
df["novel_transcript_count"] = df["novel_transcript_count"].fillna(0).astype(int)

# tissue breadth
if "tissue" in novel_tissue.columns:
    novel_tissue_breadth = (
        novel_tissue.groupby("peptide")["tissue"]
        .nunique()
        .reset_index(name="novel_tissue_count")
    )
    df = df.merge(novel_tissue_breadth, on="peptide", how="left")
    df["novel_tissue_count"] = df["novel_tissue_count"].fillna(0).astype(int)
else:
    df["novel_tissue_count"] = 0

# -----------------------------
# 2. Dominant transcript label
# -----------------------------
df = df.merge(dom, on="gene_id", how="left")
df["is_major_transcript"] = df["transcript_id"] == df["most_tx"]

# -----------------------------
# 3. Isoform multiplicity
# -----------------------------
pep_tx = (
    df.groupby("peptide")["transcript_id"]
    .nunique()
    .reset_index(name="peptide_transcript_count")
)

pep_prot = (
    df.groupby("peptide")["protein_id"]
    .nunique()
    .reset_index(name="peptide_protein_count")
)

df = df.merge(pep_tx, on="peptide", how="left")
df = df.merge(pep_prot, on="peptide", how="left")

df["is_isoform_specific"] = df["peptide_transcript_count"] == 1
df["is_multi_isoform"] = df["peptide_transcript_count"] > 1

# -----------------------------
# 4. Summary table
# -----------------------------
summary_rows = []

def add_summary(label, mask):
    sub = df.loc[mask].copy()
    if len(sub) == 0:
        return
    summary_rows.append({
        "group": label,
        "n": len(sub),
        "mean_abs_rank_gap": sub["abs_rank_gap"].mean(),
        "median_abs_rank_gap": sub["abs_rank_gap"].median(),
        "mean_mhcflurry_rank": sub["mhcflurry_rank"].mean(),
        "mean_ba_rank": sub["ba_rank"].mean(),
    })

add_summary("novel_catalog", df["is_novel_catalog"])
add_summary("non_novel_catalog", ~df["is_novel_catalog"])
add_summary("major_transcript", df["is_major_transcript"] == True)
add_summary("non_major_transcript", df["is_major_transcript"] == False)
add_summary("isoform_specific", df["is_isoform_specific"])
add_summary("multi_isoform", df["is_multi_isoform"])

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(BASE / "isoform_novel_summary.csv", index=False)
print("Saved: isoform_novel_summary.csv")
print(summary_df)

# -----------------------------
# 5. Gene-level table
# -----------------------------
gene_df = (
    df.groupby("gene_name")
    .agg(
        n=("peptide", "size"),
        unique_peptides=("peptide", "nunique"),
        unique_transcripts=("transcript_id", "nunique"),
        mean_abs_rank_gap=("abs_rank_gap", "mean"),
        mean_rank_gap=("rank_gap", "mean"),
        novel_fraction=("is_novel_catalog", "mean"),
    )
    .reset_index()
    .sort_values("mean_abs_rank_gap", ascending=False)
)

gene_df.to_csv(BASE / "isoform_gene_level_bias.csv", index=False)
print("Saved: isoform_gene_level_bias.csv")

# -----------------------------
# 6. Peptide-level summary
# -----------------------------
peptide_df = (
    df.groupby("peptide")
    .agg(
        gene_name=("gene_name", lambda x: ",".join(sorted(set(map(str, x.dropna()))))[:300]),
        peptide_transcript_count=("peptide_transcript_count", "max"),
        peptide_protein_count=("peptide_protein_count", "max"),
        is_novel_catalog=("is_novel_catalog", "max"),
        novel_tissue_count=("novel_tissue_count", "max"),
        mean_abs_rank_gap=("abs_rank_gap", "mean"),
        max_abs_rank_gap=("abs_rank_gap", "max"),
    )
    .reset_index()
    .sort_values("max_abs_rank_gap", ascending=False)
)

peptide_df.to_csv(BASE / "isoform_peptide_level_summary.csv", index=False)
print("Saved: isoform_peptide_level_summary.csv")

# -----------------------------
# 7. Interesting candidate tables
# -----------------------------
novel_extreme = df[df["is_novel_catalog"]].sort_values("abs_rank_gap", ascending=False)
novel_extreme.head(500).to_csv(BASE / "novel_extreme_disagreement_top500.csv", index=False)

nonmajor_extreme = df[df["is_major_transcript"] == False].sort_values("abs_rank_gap", ascending=False)
nonmajor_extreme.head(500).to_csv(BASE / "nonmajor_transcript_disagreement_top500.csv", index=False)

multi_isoform_extreme = df[df["is_multi_isoform"]].sort_values("abs_rank_gap", ascending=False)
multi_isoform_extreme.head(500).to_csv(BASE / "multi_isoform_disagreement_top500.csv", index=False)

print("Saved top500 candidate tables")

# -----------------------------
# 8. Figures
# -----------------------------
# Novel vs non-novel
plot1 = pd.DataFrame({
    "group": np.where(df["is_novel_catalog"], "Novel", "Non-novel"),
    "abs_rank_gap": df["abs_rank_gap"]
}).dropna()

plt.figure(figsize=(6, 4))
plot1.boxplot(by="group", column="abs_rank_gap")
plt.title("Novel vs non-novel: absolute rank disagreement")
plt.suptitle("")
plt.ylabel("abs(ba_rank - mhcflurry_rank)")
plt.tight_layout()
plt.savefig(BASE / "fig_novel_vs_non_novel_rankgap.png", dpi=300)
plt.close()

# Major vs non-major
plot2 = df[df["is_major_transcript"].notna()].copy()
plot2["group"] = np.where(plot2["is_major_transcript"], "Major transcript", "Non-major transcript")

plt.figure(figsize=(6, 4))
plot2.boxplot(by="group", column="abs_rank_gap")
plt.title("Major vs non-major transcript: disagreement")
plt.suptitle("")
plt.ylabel("abs(ba_rank - mhcflurry_rank)")
plt.tight_layout()
plt.savefig(BASE / "fig_major_vs_nonmajor_rankgap.png", dpi=300)
plt.close()

# Transcript count effect
tx_effect = (
    df.groupby("peptide_transcript_count")["abs_rank_gap"]
    .mean()
    .reset_index()
    .sort_values("peptide_transcript_count")
)

plt.figure(figsize=(6, 4))
plt.plot(tx_effect["peptide_transcript_count"], tx_effect["abs_rank_gap"], marker="o")
plt.xlabel("Number of transcripts per peptide")
plt.ylabel("Mean absolute rank gap")
plt.title("Isoform complexity vs model disagreement")
plt.tight_layout()
plt.savefig(BASE / "fig_isoform_complexity_rankgap.png", dpi=300)
plt.close()

print("Saved figures")

# -----------------------------
# 9. Console outputs
# -----------------------------
print("\n=== TOP GENES BY DISAGREEMENT ===")
print(gene_df.head(20).to_string(index=False))

print("\n=== TOP NOVEL EXTREME CASES ===")
cols = [
    "peptide", "allele", "gene_name", "transcript_id",
    "is_novel_catalog", "is_major_transcript",
    "peptide_transcript_count", "mhcflurry_rank", "ba_rank", "abs_rank_gap"
]
print(novel_extreme[cols].head(20).to_string(index=False))
