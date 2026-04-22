import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import spearmanr

BASE = Path("/Users/eymen/Desktop/v105")

MERGED = BASE / "mhcflurry_vs_netmhcpan_8mer_merged.csv"
GENE = BASE / "isoform_gene_level_bias.csv"
PEPTIDE = BASE / "isoform_peptide_level_summary.csv"

outdir = BASE / "figures_scientific_v2"
outdir.mkdir(exist_ok=True)

# -----------------------------
# Load
# -----------------------------
df = pd.read_csv(MERGED)
gene_df = pd.read_csv(GENE)
pep_df = pd.read_csv(PEPTIDE)

# Clean
for col in ["mhcflurry_rank", "ba_rank"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.dropna(subset=["mhcflurry_rank", "ba_rank"]).copy()
df = df[(df["mhcflurry_rank"] > 0) & (df["ba_rank"] > 0)].copy()

df["rank_gap_signed"] = df["ba_rank"] - df["mhcflurry_rank"]
df["abs_rank_gap"] = df["rank_gap_signed"].abs()
df["mean_rank"] = (df["mhcflurry_rank"] + df["ba_rank"]) / 2

# -----------------------------
# Figure 1: density-like scatter
# -----------------------------
plt.figure(figsize=(6.5, 6))
plt.scatter(
    np.log10(df["mhcflurry_rank"]),
    np.log10(df["ba_rank"]),
    alpha=0.08,
    s=6
)
mn = min(np.log10(df["mhcflurry_rank"]).min(), np.log10(df["ba_rank"]).min())
mx = max(np.log10(df["mhcflurry_rank"]).max(), np.log10(df["ba_rank"]).max())
plt.plot([mn, mx], [mn, mx], linestyle="--", linewidth=1)
plt.xlabel("log10 MHCflurry rank")
plt.ylabel("log10 NetMHCpan BA rank")
plt.title("Global agreement across 8-mer peptide–HLA pairs")
plt.tight_layout()
plt.savefig(outdir / "fig1_global_rank_agreement.png", dpi=300)
plt.close()

# -----------------------------
# Figure 2: rank gap vs mean rank
# -----------------------------
plot_df = df.sample(min(15000, len(df)), random_state=42).copy()

plt.figure(figsize=(7, 5))
plt.scatter(
    np.log10(plot_df["mean_rank"]),
    plot_df["rank_gap_signed"],
    alpha=0.10,
    s=8
)
plt.axhline(0, linestyle="--", linewidth=1)
plt.xlabel("log10 mean rank")
plt.ylabel("NetMHCpan rank - MHCflurry rank")
plt.title("Model disagreement as a function of ranking strength")
plt.tight_layout()
plt.savefig(outdir / "fig2_rankgap_vs_meanrank.png", dpi=300)
plt.close()

# -----------------------------
# Figure 3: allele-level effect size
# -----------------------------
allele_rows = []
for allele, g in df.groupby("allele"):
    if len(g) < 20:
        continue
    rho, _ = spearmanr(g["mhcflurry_rank"], g["ba_rank"])
    allele_rows.append({
        "allele": allele,
        "n": len(g),
        "spearman_rho": rho,
        "mean_abs_rank_gap": g["abs_rank_gap"].mean()
    })

allele_df = pd.DataFrame(allele_rows).sort_values("mean_abs_rank_gap", ascending=False)
allele_df.to_csv(outdir / "table_allele_effects.csv", index=False)

plt.figure(figsize=(7, 5))
sizes = 20 + 0.03 * allele_df["n"]
plt.scatter(
    allele_df["mean_abs_rank_gap"],
    allele_df["spearman_rho"],
    s=sizes,
    alpha=0.8
)
for _, r in allele_df.iterrows():
    plt.text(r["mean_abs_rank_gap"], r["spearman_rho"], r["allele"], fontsize=8)
plt.xlabel("Mean absolute rank gap")
plt.ylabel("Spearman rho")
plt.title("Allele-specific agreement and disagreement")
plt.tight_layout()
plt.savefig(outdir / "fig3_allele_effect_size.png", dpi=300)
plt.close()

# -----------------------------
# Figure 4: gene-level landscape
# -----------------------------
gene_plot = gene_df.copy()
gene_plot = gene_plot.dropna(subset=["unique_transcripts", "mean_abs_rank_gap"])
gene_plot = gene_plot[gene_plot["n"] >= 5].copy()

plt.figure(figsize=(7, 5))
plt.scatter(
    gene_plot["unique_transcripts"],
    gene_plot["mean_abs_rank_gap"],
    alpha=0.5,
    s=20
)

top_gene = gene_plot.sort_values("mean_abs_rank_gap", ascending=False).head(15)
for _, r in top_gene.iterrows():
    plt.text(r["unique_transcripts"], r["mean_abs_rank_gap"], str(r["gene_name"]), fontsize=8)

plt.xlabel("Number of unique transcripts per gene")
plt.ylabel("Mean absolute rank gap")
plt.title("Gene-level disagreement vs transcript complexity")
plt.tight_layout()
plt.savefig(outdir / "fig4_gene_landscape.png", dpi=300)
plt.close()

# -----------------------------
# Figure 5: top peptide discordance
# -----------------------------
pep_plot = pep_df.copy()
pep_plot = pep_plot.sort_values("max_abs_rank_gap", ascending=False).head(20).copy()
pep_plot = pep_plot.sort_values("max_abs_rank_gap", ascending=True)

plt.figure(figsize=(8, 6))
sc = plt.scatter(
    pep_plot["max_abs_rank_gap"],
    pep_plot["peptide"],
    s=40 + 4 * pep_plot["peptide_transcript_count"],
    c=pep_plot["peptide_transcript_count"]
)
plt.xlabel("Maximum absolute rank gap")
plt.ylabel("Peptide")
plt.title("Top discordant peptides and transcript multiplicity")
cbar = plt.colorbar(sc)
cbar.set_label("Transcript count")
plt.tight_layout()
plt.savefig(outdir / "fig5_top_peptide_discordance.png", dpi=300)
plt.close()

# -----------------------------
# Figure 6: isoform-specific vs multi-isoform
# -----------------------------
iso_df = pep_df.copy()
iso_df["group"] = np.where(
    iso_df["peptide_transcript_count"] == 1,
    "Isoform-specific",
    "Multi-isoform"
)

vals1 = iso_df.loc[iso_df["group"] == "Isoform-specific", "mean_abs_rank_gap"].dropna()
vals2 = iso_df.loc[iso_df["group"] == "Multi-isoform", "mean_abs_rank_gap"].dropna()

plt.figure(figsize=(6, 5))
plt.boxplot([vals1, vals2], tick_labels=["Isoform-specific", "Multi-isoform"])
plt.ylabel("Mean absolute rank gap")
plt.title("Disagreement by peptide isoform context")
plt.tight_layout()
plt.savefig(outdir / "fig6_isoform_context_boxplot.png", dpi=300)
plt.close()

# -----------------------------
# Helpful summary files
# -----------------------------
summary = pd.DataFrame([{
    "n_pairs": len(df),
    "n_unique_peptides": df["peptide"].nunique(),
    "n_alleles": df["allele"].nunique(),
    "global_mean_abs_rank_gap": df["abs_rank_gap"].mean(),
    "global_median_abs_rank_gap": df["abs_rank_gap"].median()
}])
summary.to_csv(outdir / "table_figure_summary.csv", index=False)

print("Done. Files saved in:", outdir)
print(sorted([p.name for p in outdir.iterdir()]))
