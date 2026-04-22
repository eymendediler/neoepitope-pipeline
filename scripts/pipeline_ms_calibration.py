import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score

# -----------------------------
# LOAD DATA
# -----------------------------
df = pd.read_csv("mhcflurry_vs_netmhcpan_8mer_merged.csv")

print("Loaded:", df.shape)

# -----------------------------
# CLEAN
# -----------------------------
df["mhcflurry_rank"] = pd.to_numeric(df["mhcflurry_rank"], errors="coerce")
df["ba_rank"] = pd.to_numeric(df["ba_rank"], errors="coerce")

df = df.dropna(subset=["mhcflurry_rank", "ba_rank"])

# -----------------------------
# STEP 1: CONSENSUS SCORE
# -----------------------------
df["consensus_rank"] = (df["mhcflurry_rank"] + df["ba_rank"]) / 2

print("Consensus created")

# -----------------------------
# STEP 2: ADD MS DATA
# -----------------------------
# ⚠️ IMPORTANT:
# Replace this with your real MS peptide list later

# For now: simulate MS positives (top NetMHC binders)
df["observed_in_ms"] = df["ba_rank"] < 1

# -----------------------------
# STEP 3: ALLELE-SPECIFIC THRESHOLDS
# -----------------------------
thresholds = []

for allele, g in df.groupby("allele"):
    g = g.copy()

    if g["observed_in_ms"].sum() < 20:
        continue

    # sort by consensus
    g = g.sort_values("consensus_rank")

    # take top 2%
    cutoff_index = int(len(g) * 0.02)
    cutoff_value = g.iloc[cutoff_index]["consensus_rank"]

    thresholds.append({
        "allele": allele,
        "cutoff": cutoff_value,
        "n": len(g),
        "ms_positive": int(g["observed_in_ms"].sum())
    })

threshold_df = pd.DataFrame(thresholds)
threshold_df.to_csv("allele_thresholds.csv", index=False)

print("Saved allele_thresholds.csv")

# -----------------------------
# STEP 4: APPLY THRESHOLD
# -----------------------------
df = df.merge(threshold_df, on="allele", how="left")

df["high_confidence"] = df["consensus_rank"] <= df["cutoff"]

# -----------------------------
# STEP 5: FINAL TABLE
# -----------------------------
df = df.sort_values(["high_confidence", "consensus_rank"], ascending=[False, True])

df.to_csv("final_prioritized_peptides.csv", index=False)

print("Saved final_prioritized_peptides.csv")

# -----------------------------
# STEP 6: QUICK STATS
# -----------------------------
print("\n=== SUMMARY ===")
print("Total peptides:", len(df))
print("High confidence:", df["high_confidence"].sum())
print("MS positives:", df["observed_in_ms"].sum())
