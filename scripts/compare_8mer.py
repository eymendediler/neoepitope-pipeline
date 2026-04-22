import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import cohen_kappa_score, confusion_matrix

# --- Files ---
mhc_file = "/Users/eymen/Desktop/v105/mhcflurry_multi_length_results/mhcflurry_results_len8.csv"
net_file = "/Users/eymen/netmhcpan_v105_compact_run/results/chunk_results_8mer/combined_8mer.xls"

# --- Load MHCflurry ---
mhc = pd.read_csv(mhc_file)
mhc = mhc.rename(columns={
    "affinity": "mhcflurry_affinity",
    "rank": "mhcflurry_rank",
    "presentation_score": "mhcflurry_presentation_score"
})

print("MHCflurry:", mhc.shape)
print(mhc.head())

# --- Load NetMHCpan headers ---
h1 = pd.read_csv(net_file, sep="\t", header=None, nrows=1)
net = pd.read_csv(net_file, sep="\t", header=None, skiprows=2)

alleles = [h1.iloc[0, i] for i in range(3, 123, 6)]
print("n_alleles:", len(alleles))
print("alleles:", alleles)

# --- Convert NetMHCpan wide -> long ---
records = []

for _, row in net.iterrows():
    peptide = row[1]
    for j, allele in enumerate(alleles):
        start = 3 + j * 6
        records.append({
            "peptide": peptide,
            "allele": allele,
            "score": row[start + 2],
            "rank": row[start + 3],
            "ba_score": row[start + 4],
            "ba_rank": row[start + 5],
        })

net_long = pd.DataFrame(records)

for c in ["score", "rank", "ba_score", "ba_rank"]:
    net_long[c] = pd.to_numeric(net_long[c], errors="coerce")

print("NetMHCpan long:", net_long.shape)
print(net_long.head())

# --- Merge ---
df = pd.merge(
    mhc,
    net_long,
    on=["peptide", "allele"],
    how="inner"
)

print("Merged:", df.shape)
print(df.head())

# --- Rank-based comparison ---
sub = df[["mhcflurry_rank", "ba_rank"]].dropna()
rho, pval = spearmanr(sub["mhcflurry_rank"], sub["ba_rank"])
print("Spearman (mhcflurry_rank vs netmhcpan_ba_rank):", rho, pval)

# --- Binder agreement ---
df["mhcflurry_binder"] = df["mhcflurry_rank"] < 2
df["netmhcpan_binder"] = df["ba_rank"] < 2

sub2 = df[["mhcflurry_binder", "netmhcpan_binder"]].dropna()

print("Confusion matrix:")
print(confusion_matrix(sub2["mhcflurry_binder"], sub2["netmhcpan_binder"]))
print("Cohen kappa:", cohen_kappa_score(sub2["mhcflurry_binder"], sub2["netmhcpan_binder"]))

# --- Save ---
df.to_csv("mhcflurry_vs_netmhcpan_8mer_merged.csv", index=False)
print("Saved: mhcflurry_vs_netmhcpan_8mer_merged.csv")
