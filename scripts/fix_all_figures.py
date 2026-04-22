import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Optional: better label placement
try:
    from adjustText import adjust_text
    HAS_ADJUST = True
except ImportError:
    HAS_ADJUST = False

plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 11

df = pd.read_csv("mhcflurry_vs_netmhcpan_8mer_merged.csv")

# Clean numeric columns
for col in ["mhcflurry_rank", "ba_rank"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.dropna(subset=["mhcflurry_rank", "ba_rank"]).copy()
df = df[(df["mhcflurry_rank"] > 0) & (df["ba_rank"] > 0)].copy()

df["rank_gap"] = df["ba_rank"] - df["mhcflurry_rank"]
df["abs_rank_gap"] = df["rank_gap"].abs()
df["mean_rank"] = (df["mhcflurry_rank"] + df["ba_rank"]) / 2

# -----------------------------
# FIG 1 — Global agreement
# -----------------------------
plot_df = df.sample(min(12000, len(df)), random_state=42)

plt.figure(figsize=(6.5, 6))
plt.scatter(
    np.log10(plot_df["mhcflurry_rank"]),
    np.log10(plot_df["ba_rank"]),
    alpha=0.12,
    s=8
)

mn = min(np.log10(df["mhcflurry_rank"]).min(), np.log10(df["ba_rank"]).min())
mx = max(np.log10(df["mhcflurry_rank"]).max(), np.log10(df["ba_rank"]).max())
plt.plot([mn, mx], [mn, mx], linestyle="--", linewidth=1)

plt.xlabel("log10 MHCflurry rank")
plt.ylabel("log10 NetMHCpan BA rank")
plt.title("Global agreement across 8-mer peptide–HLA pairs")
plt.tight_layout()
plt.savefig("fig1_global_rank_agreement_clean.png")
plt.close()

# -----------------------------
# FIG 2 — Disagreement vs ranking strength
# -----------------------------
plot_df2 = df.sample(min(15000, len(df)), random_state=42)

plt.figure(figsize=(7, 5))
plt.scatter(
    np.log10(plot_df2["mean_rank"]),
    plot_df2["rank_gap"],
    alpha=0.10,
    s=8
)
plt.axhline(0, linestyle="--", linewidth=1)
plt.xlabel("log10 mean rank")
plt.ylabel("NetMHCpan rank - MHCflurry rank")
plt.title("Model disagreement as a function of ranking strength")
plt.tight_layout()
plt.savefig("fig2_rankgap_vs_meanrank_clean.png")
plt.close()

# -----------------------------
# FIG 3 — Allele-specific effect size
# -----------------------------
from scipy.stats import spearmanr

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

plt.figure(figsize=(8, 5.8))
sizes = 20 + 0.03 * allele_df["n"]
plt.scatter(
    allele_df["mean_abs_rank_gap"],
    allele_df["spearman_rho"],
    s=sizes,
    alpha=0.8
)

texts = []
for _, r in allele_df.iterrows():
    texts.append(
        plt.text(r["mean_abs_rank_gap"], r["spearman_rho"], r["allele"], fontsize=8)
    )

if HAS_ADJUST:
    adjust_text(texts, arrowprops=dict(arrowstyle="-", lw=0.5))

plt.xlabel("Mean absolute rank gap")
plt.ylabel("Spearman rho")
plt.title("Allele-specific agreement and disagreement")
plt.tight_layout()
plt.savefig("fig3_allele_effect_size_clean.png")
plt.close()

# -----------------------------
# FIG 4 — Gene-level landscape
# -----------------------------
gene_df = (
    df.groupby("gene_name")
    .agg(
        n=("peptide", "size"),
        unique_transcripts=("transcript_id", "nunique"),
        mean_abs_rank_gap=("abs_rank_gap", "mean")
    )
    .reset_index()
)

gene_plot = gene_df[gene_df["n"] >= 5].copy()

plt.figure(figsize=(8, 6))
plt.scatter(
    gene_plot["unique_transcripts"],
    gene_plot["mean_abs_rank_gap"],
    alpha=0.5,
    s=22
)

top_gene = gene_plot.sort_values("mean_abs_rank_gap", ascending=False).head(12)
texts = []
for _, r in top_gene.iterrows():
    texts.append(
        plt.text(r["unique_transcripts"], r["mean_abs_rank_gap"], str(r["gene_name"]), fontsize=8)
    )

if HAS_ADJUST:
    adjust_text(texts, arrowprops=dict(arrowstyle="-", lw=0.5))

plt.xlabel("Number of unique transcripts per gene")
plt.ylabel("Mean absolute rank gap")
plt.title("Gene-level disagreement vs transcript complexity")
plt.tight_layout()
plt.savefig("fig4_gene_landscape_clean.png")
plt.close()

# -----------------------------
# FIG 5 — Top discordant peptides
# -----------------------------
pep_df = (
    df.groupby("peptide")
    .agg(
        peptide_transcript_count=("transcript_id", "nunique"),
        max_abs_rank_gap=("abs_rank_gap", "max")
    )
    .reset_index()
    .sort_values("max_abs_rank_gap", ascending=False)
    .head(20)
    .sort_values("max_abs_rank_gap", ascending=True)
)

plt.figure(figsize=(8, 6))
sc = plt.scatter(
    pep_df["max_abs_rank_gap"],
    pep_df["peptide"],
    s=40 + 6 * pep_df["peptide_transcript_count"],
    c=pep_df["peptide_transcript_count"]
)
plt.xlabel("Maximum absolute rank gap")
plt.ylabel("Peptide")
plt.title("Top discordant peptides and transcript multiplicity")
cbar = plt.colorbar(sc)
cbar.set_label("Transcript count")
plt.tight_layout()
plt.savefig("fig5_top_peptide_discordance_clean.png")
plt.close()

# -----------------------------
# FIG 6 — Isoform context boxplot
# -----------------------------
iso_df = (
    df.groupby("peptide")
    .agg(
        peptide_transcript_count=("transcript_id", "nunique"),
        mean_abs_rank_gap=("abs_rank_gap", "mean")
    )
    .reset_index()
)

vals1 = iso_df.loc[iso_df["peptide_transcript_count"] == 1, "mean_abs_rank_gap"].dropna()
vals2 = iso_df.loc[iso_df["peptide_transcript_count"] > 1, "mean_abs_rank_gap"].dropna()

plt.figure(figsize=(6, 5))
plt.boxplot([vals1, vals2], tick_labels=["Isoform-specific", "Multi-isoform"], showfliers=False)
plt.ylabel("Mean absolute rank gap")
plt.title("Disagreement by peptide isoform context")
plt.tight_layout()
plt.savefig("fig6_isoform_context_boxplot_clean.png")
plt.close()

print("All clean figures saved.")
