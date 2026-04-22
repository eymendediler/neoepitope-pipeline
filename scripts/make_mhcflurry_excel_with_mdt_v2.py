#!/usr/bin/env python3
import os, re, glob
import pandas as pd
import numpy as np

MHC_DIR = "results/mhcflurry_by_tissue"
MDT_DIR = "results/mdt_by_tissue"
OUT_XLSX = "results/mhcflurry_by_tissue_WITH_MDT_v2.xlsx"

RANK_STRONG = 2.0
TOP_N_PER_TISSUE = 200

BASE_COLS = [
    "tissue",
    "peptide","allele",
    "presentation_score","affinity","rank",
    "gene_name","gene_id","transcript_id","protein_id","start_position",
    "gene_id_map","transcript_id_map","protein_id_map","gene_name_map",
    "mdt_enst","is_tissue_MDT",
    "top1_median_tpm","top1_mean_tpm","top2_enst","top2_median_tpm","delta_median_tpm",
]

def safe_sheet_name(name: str) -> str:
    name = re.sub(r"[:\\/?*\[\]]", "_", str(name))
    name = re.sub(r"\s+", " ", name).strip()
    return (name[:31] if len(name) > 31 else name) or "Sheet"

def strip_version(x):
    if pd.isna(x): return x
    return str(x).split(".")[0]

def load_mdt_tables():
    m = {}
    for path in glob.glob(os.path.join(MDT_DIR, "*.tsv")):
        df = pd.read_csv(path, sep="\t")
        if "tissue" not in df.columns:
            df.insert(0, "tissue", os.path.splitext(os.path.basename(path))[0])
        tissue = str(df["tissue"].iloc[0])
        df["_gene_nov"] = df["ensg"].map(strip_version) if "ensg" in df.columns else None
        df["_mdt_tx_nov"] = df["mdt_enst"].map(strip_version) if "mdt_enst" in df.columns else None
        keep = [c for c in ["_gene_nov","mdt_enst","top1_median_tpm","top1_mean_tpm","top2_enst","top2_median_tpm","delta_median_tpm"] if c in df.columns]
        m[tissue] = df[keep].drop_duplicates("_gene_nov")
    return m

def main():
    os.makedirs("results", exist_ok=True)
    mdt_by_tissue = load_mdt_tables()

    mhc_files = sorted(glob.glob(os.path.join(MHC_DIR, "*.csv")))
    if not mhc_files:
        raise SystemExit(f"No CSV files found in {MHC_DIR}")

    all_rows = []
    strong_rows = []
    top_rows = []
    summary_rows = []

    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as xw:
        for f in mhc_files:
            tissue_file = os.path.splitext(os.path.basename(f))[0]
            df = pd.read_csv(f)
            if "tissue" not in df.columns:
                df["tissue"] = tissue_file

            df["_gene_nov"] = df["gene_id"].map(strip_version) if "gene_id" in df.columns else None
            df["_tx_nov"] = df["transcript_id"].map(strip_version) if "transcript_id" in df.columns else None

            # try match MDT table by exact filename tissue first
            mdt = mdt_by_tissue.get(tissue_file)
            if mdt is None:
                # fallback: match by sanitized tissue name
                for k, v in mdt_by_tissue.items():
                    if safe_sheet_name(k).lower() == tissue_file.lower():
                        mdt = v
                        break

            if mdt is not None:
                df = df.merge(mdt, on="_gene_nov", how="left")
                df["is_tissue_MDT"] = (df["_tx_nov"] == df["mdt_enst"].map(strip_version))
            else:
                for c in ["mdt_enst","top1_median_tpm","top1_mean_tpm","top2_enst","top2_median_tpm","delta_median_tpm"]:
                    df[c] = pd.NA
                df["is_tissue_MDT"] = pd.NA

            # select columns that exist
            cols = [c for c in BASE_COLS if c in df.columns]
            out = df[cols].copy()

            # write tissue sheet
            out.to_excel(xw, sheet_name=safe_sheet_name(tissue_file), index=False)

            # accumulate ALL
            all_rows.append(out)

            # STRONG binders
            if "rank" in out.columns:
                r = pd.to_numeric(out["rank"], errors="coerce")
                strong = out[r <= RANK_STRONG].copy()
                if not strong.empty:
                    strong_rows.append(strong)

            # TOP N per tissue by rank (lowest)
            if "rank" in out.columns:
                tmp = out.copy()
                tmp["_rank_num"] = pd.to_numeric(tmp["rank"], errors="coerce")
                tmp = tmp.sort_values(["_rank_num","presentation_score"], ascending=[True, False])
                tmp = tmp.drop(columns=["_rank_num"])
                top = tmp.head(TOP_N_PER_TISSUE)
                top_rows.append(top)

            # summary per tissue
            n = len(out)
            uniq_pep = out["peptide"].nunique() if "peptide" in out.columns else np.nan
            uniq_gene = out["gene_id"].nunique() if "gene_id" in out.columns else np.nan
            mdt_avail = int(out["mdt_enst"].notna().sum()) if "mdt_enst" in out.columns else 0
            mdt_frac = (mdt_avail / n) if n else np.nan
            if "is_tissue_MDT" in out.columns:
                x = out["is_tissue_MDT"].dropna()
                mdt_match = float((x==True).mean()) if len(x) else np.nan
            else:
                mdt_match = np.nan
            summary_rows.append({
                "tissue": tissue_file,
                "rows": n,
                "unique_peptides": uniq_pep,
                "unique_genes": uniq_gene,
                "mdt_available_rows": mdt_avail,
                "mdt_available_fraction": mdt_frac,
                "mdt_match_fraction": mdt_match,
            })

        all_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
        all_df.to_excel(xw, sheet_name="ALL", index=False)

        if strong_rows:
            strong_df = pd.concat(strong_rows, ignore_index=True)
            strong_df.to_excel(xw, sheet_name=f"STRONG_rankLE{int(RANK_STRONG)}", index=False)

        if top_rows:
            top_df = pd.concat(top_rows, ignore_index=True)
            top_df.to_excel(xw, sheet_name=f"TOP{TOP_N_PER_TISSUE}", index=False)

        summary_df = pd.DataFrame(summary_rows).sort_values("rows", ascending=False)
        summary_df.to_excel(xw, sheet_name="SUMMARY", index=False)

    print("Wrote Excel:", OUT_XLSX)

if __name__ == "__main__":
    main()
