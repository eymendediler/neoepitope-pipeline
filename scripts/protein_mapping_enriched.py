import pandas as pd

MHC_FILE = "mhcflurry_results_merged.csv"
MAP_FILE = "ensg_enst_ensp_geneName_v105.tsv.gz"
OUTPUT = "mhcflurry_results_enriched_mapping.csv"

print("[INFO] Loading MHCflurry merged data...")
df = pd.read_csv(MHC_FILE)

print("[INFO] Loading Ensembl v105 mapping table...")
map_df = pd.read_csv(MAP_FILE, sep="\t")

# Kolon isimlerini görelim (debug)
print("[INFO] Mapping file columns:", list(map_df.columns))

# Ensembl mapping tablosundaki gerekli kolonlar:
# 'Protein stable ID'
# 'Transcript stable ID'
# 'Gene stable ID'
# 'Gene name'

# Bu kolonları daha kolay eşleşme için yeniden adlandıralım
map_df = map_df.rename(columns={
    "Protein stable ID": "protein_id_map",
    "Transcript stable ID": "transcript_id_map",
    "Gene stable ID": "gene_id_map",
    "Gene name": "gene_name_map"
})

print("[INFO] Merging on protein_id...")
merged = df.merge(
    map_df,
    how="left",
    left_on="protein_id",
    right_on="protein_id_map"
)

print("[INFO] Merge complete.")

# Eksik eşleşen proteinler
missing = merged[merged["protein_id_map"].isna()]
missing_count = len(missing)
print(f"[INFO] Missing mappings: {missing_count}")

if missing_count > 0:
    missing.to_csv("missing_proteins.csv", index=False)
    print("[INFO] Missing proteins saved to missing_proteins.csv")

print("[INFO] Saving enriched mapping file...")
merged.to_csv(OUTPUT, index=False)

print(f"[DONE] Enriched mapping saved as: {OUTPUT}")

