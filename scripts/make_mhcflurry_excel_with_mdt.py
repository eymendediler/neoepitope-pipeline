#!/usr/bin/env python3
import os, re, glob
import pandas as pd

MHC_DIR = "results/mhcflurry_by_tissue"
MDT_DIR = "results/mdt_by_tissue"
OUT_XLSX = "results/mhcflurry_by_tissue_WITH_MDT.xlsx"

IMPORTANT_COLS = [
    "tissue",
    "peptide", "allele",
    "presentation_score", "affinity", "rank",
    "gene_name", "gene_id", "transcript_id", "protein_id", "start_position",
    # mapping columns (if present)
    "gene_id_map", "transcript_id_map", "protein_id_map", "gene_name_map",
    # MDT columns we will add
    "mdt_enst", "is_tissue_MDT",
    "top1_median_tpm", "top1_mean_tpm", "top2_enst", "top2_median_tpm", "delta_median_tpm",
]

def safe_sheet_name(name: str) -> str:
    # Excel sheet name max 31 chars, no : \ / ? * [ ]
    name = re.sub(r"[:\\/?*\[\]]", "_", str(name))
    name = re.sub(r"\s+", " ", name).strip()
    if len(name) > 31:
        name = name[:31]
    return name or "Sheet"

def load_mdt_map():
    # Build dict: tissue -> MDT dataframe keyed by ensg
    m = {}
    for path in glob.glob(os.path.join(MDT_DIR, "*.tsv")):
        df = pd.read_csv(path, sep="\t")
        # expect columns: tissue, ensg, mdt_enst, top1_median_tpm, ...
        if "tissue" not in df.columns:
            # tissue name can be inferred from filename
            tissue = os.path.splitext(os.path.basename(path))[0]
            df.insert(0, "tissue", tissue)
        tissue = str(df["tissue"].iloc[0])
        m[tissue] = df
    return m

def strip_version(x):
    if pd.isna(x): 
        return x
    return str(x).split(".")[0]

def main():
    os.makedirs("results", exist_ok=True)

    if not os.path.isdir(MHC_DIR):
        raise SystemExit(f"Missing folder: {MHC_DIR}")
    if not os.path.isdir(MDT_DIR):
        raise SystemExit(f"Missing folder: {MDT_DIR}")

    mdt_by_tissue = load_mdt_map()
    mhc_files = sorted(glob.glob(os.path.join(MHC_DIR, "*.csv")))
    if not mhc_files:
        raise SystemExit(f"No CSV files found in {MHC_DIR}")

    all_rows = []
    summary = []

    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as xw:
        for f in mhc_files:
            tissue = os.path.splitext(os.path.basename(f))[0]
            df = pd.read_csv(f)

            # ensure tissue column exists and is consistent
            if "tissue" not in df.columns:
                df["tissue"] = tissue

            # normalize IDs (versionless compare)
            df["_gene_nov"] = df["gene_id"].map(strip_version) if "gene_id" in df.columns else None
            df["_tx_nov"] = df["transcript_id"].map(strip_version) if "transcript_id" in df.columns else None

            # Find matching MDT table:
            # MDT files have a 'tissue' column that is the original GTEx tissue name; filenames are sanitized.
            # We'll match by filename tissue first; if not found, try to match by tissue column values.
            mdt = None
            if tissue in mdt_by_tissue:
                mdt = mdt_by_tissue[tissue]
            else:
                # try match by tissue name inside MDT tables
                for k, v in mdt_by_tissue.items():
                    if safe_sheet_name(k).lower() == tissue.lower():
                        mdt = v
                        break

            if mdt is not None:
                mdt = mdt.copy()
                mdt["_gene_nov"] = mdt["ensg"].map(strip_version) if "ensg" in mdt.columns else None
                mdt["_mdt_tx_nov"] = mdt["mdt_enst"].map(strip_version) if "mdt_enst" in mdt.columns else None

                keep_cols = [c for c in ["_gene_nov","mdt_enst","top1_median_tpm","top1_mean_tpm","top2_enst","top2_median_tpm","delta_median_tpm"] if c in mdt.columns]
                mdt_small = mdt[keep_cols].drop_duplicates("_gene_nov")

                df = df.merge(mdt_small, on="_gene_nov", how="left")
                df["is_tissue_MDT"] = (df["_tx_nov"] == df.get("mdt_enst").map(strip_version))
            else:
                df["mdt_enst"] = pd.NA
                df["top1_median_tpm"] = pd.NA
                df["top1_mean_tpm"] = pd.NA
                df["top2_enst"] = pd.NA
                df["top2_median_tpm"] = pd.NA
                df["delta_median_tpm"] = pd.NA
                df["is_tissue_MDT"] = pd.NA

            # pick important columns that exist
            cols = [c for c in IMPORTANT_COLS if c in df.columns]
            out = df[cols].copy()

            # Write tissue sheet
            sheet = safe_sheet_name(tissue)
            out.to_excel(xw, sheet_name=sheet, index=False)

            all_rows.append(out)

            # summary stats
            n = len(out)
            mdt_match = out["is_tissue_MDT"].dropna().mean() if "is_tissue_MDT" in out.columns else float("nan")
            summary.append({"tissue": tissue, "rows": n, "mdt_match_fraction": mdt_match})

        # ALL sheet
        all_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
        all_df.to_excel(xw, sheet_name="ALL", index=False)

        # SUMMARY sheet
        summary_df = pd.DataFrame(summary).sort_values("rows", ascending=False)
        summary_df.to_excel(xw, sheet_name="SUMMARY", index=False)

    print("Wrote Excel:", OUT_XLSX)
    print("Sheets:", len(mhc_files) + 2)

if __name__ == "__main__":
    main()
