import pandas as pd
from pathlib import Path

# Paths
MHC_FILE = Path("/Users/eymen/Desktop/v105/mhcflurry_multi_length_results/mhcflurry_results_len8.csv")
NET_FILE = Path("/Users/eymen/netmhcpan_v105_compact_run/results/chunk_results_8mer/netmhcpan_8mer_long_sample.csv")
OUT_FILE = Path("/Users/eymen/Desktop/v105/mhcflurry_vs_netmhcpan_8mer_merged.csv")

# Load NetMHCpan long once
print("Loading NetMHCpan long...")
net_long = pd.read_csv(NET_FILE)
net_long["allele"] = net_long["allele"].astype(str).str.replace("*", "", regex=False)
net_long["peptide"] = net_long["peptide"].astype(str).str.strip().str.upper()

for c in ["score", "rank", "ba_score", "ba_rank"]:
    net_long[c] = pd.to_numeric(net_long[c], errors="coerce")

print("NetMHCpan:", net_long.shape)
print(net_long.head())

# Remove old output if exists
if OUT_FILE.exists():
    OUT_FILE.unlink()

chunk_size = 200000
first_write = True
total_merged = 0

print("Streaming MHCflurry and merging in chunks...")

for i, chunk in enumerate(pd.read_csv(MHC_FILE, chunksize=chunk_size)):
    print(f"\nChunk {i+1} raw shape:", chunk.shape)

    chunk = chunk.rename(columns={
        "affinity": "mhcflurry_affinity",
        "presentation_score": "mhcflurry_presentation_score",
        "rank": "mhcflurry_rank",
    })

    chunk["allele"] = chunk["allele"].astype(str).str.replace("*", "", regex=False)
    chunk["peptide"] = chunk["peptide"].astype(str).str.strip().str.upper()

    merged = pd.merge(
        chunk,
        net_long,
        on=["peptide", "allele"],
        how="inner"
    )

    print(f"Chunk {i+1} merged shape:", merged.shape)

    total_merged += len(merged)

    if len(merged) > 0:
        merged.to_csv(
            OUT_FILE,
            mode="w" if first_write else "a",
            header=first_write,
            index=False
        )
        first_write = False

print("\nDone.")
print("Total merged rows:", total_merged)
print("Saved to:", OUT_FILE)
