#!/usr/bin/env python3

import os
import csv
from typing import Dict, Any, Optional, List

import pandas as pd
from Bio import SeqIO
from mhcflurry import Class1PresentationPredictor


# ====================================================
# USER SETTINGS
# ====================================================

FASTA_PATH = "/Users/eymen/Desktop/v105/Homo_sapiens.GRCh38.pep.all.fa"
MAPPING_PATH = "/Users/eymen/Desktop/v105/ensg_enst_ensp_geneName_v105.tsv.gz"
OUTPUT_CSV = "/Users/eymen/Desktop/v105/mhcflurry_results_ensembl_v105.csv"

PEPTIDE_LENGTH = 9
BATCH_SIZE = 20000

# Filter by presentation rank. Set None to keep all predictions.
RANK_THRESHOLD: Optional[float] = 2.0

# 20 alleles you provided:
ALLELES = [
    "HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*24:02",
    "HLA-A*26:01", "HLA-A*30:01", "HLA-A*23:01", "HLA-A*68:01", "HLA-A*33:01",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*27:05", "HLA-B*35:01",
    "HLA-B*38:01", "HLA-B*40:01", "HLA-B*44:02", "HLA-B*51:01", "HLA-B*18:01",
]

# Only allow standard amino acids
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


# ====================================================
# MAPPING FILE LOADER (v105-specific)
# ====================================================

def load_mapping(path: str) -> Dict[str, Dict[str, Any]]:
    """
    Load Ensembl v105 mapping file:
    Protein stable ID -> Transcript stable ID -> Gene stable ID -> Gene name
    """

    print(f"[INFO] Loading mapping file: {path}")
    df = pd.read_csv(path, sep="\t")

    protein_col = "Protein stable ID"
    transcript_col = "Transcript stable ID"
    gene_col = "Gene stable ID"
    gene_name_col = "Gene name"

    mapping: Dict[str, Dict[str, Any]] = {}

    for _, row in df.iterrows():
        ensp = str(row[protein_col])          # e.g. ENSP00000354587.5
        ensp_novers = ensp.split(".")[0]      # ENSP00000354587

        info = {
            "protein_id": ensp,
            "transcript_id": str(row[transcript_col]),
            "gene_id": str(row[gene_col]),
            "gene_name": str(row[gene_name_col]),
        }

        # Store both versioned and non-versioned keys
        mapping[ensp] = info
        mapping[ensp_novers] = info

    print(f"[INFO] Mapping loaded. Proteins mapped (keys): {len(mapping)}")
    return mapping


# ====================================================
# PEPTIDE GENERATOR
# ====================================================

def generate_peptide_batches(
    fasta_path: str,
    mapping: Dict[str, Dict[str, Any]],
    peptide_length: int,
    batch_size: int
):
    """
    Yield batches of (peptides, metas) from the FASTA file.

    peptides: List[str]
    metas: List[dict] with:
        protein_id, transcript_id, gene_id, gene_name, start_position, peptide
    """
    print(f"[INFO] Reading FASTA: {fasta_path}")

    batch_peptides: List[str] = []
    batch_meta: List[Dict[str, Any]] = []
    total = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq)
        if len(seq) < peptide_length:
            continue

        ensp_id = record.id
        ensp_novers = ensp_id.split(".")[0]

        info = mapping.get(ensp_id) or mapping.get(ensp_novers) or {
            "protein_id": ensp_id,
            "transcript_id": None,
            "gene_id": None,
            "gene_name": None,
        }

        for i in range(len(seq) - peptide_length + 1):
            peptide = seq[i:i + peptide_length]

            # Skip any peptide that contains non-standard amino acids
            if any(aa not in VALID_AA for aa in peptide):
                continue

            batch_peptides.append(peptide)
            batch_meta.append({
                "protein_id": info["protein_id"],
                "transcript_id": info["transcript_id"],
                "gene_id": info["gene_id"],
                "gene_name": info["gene_name"],
                "start_position": i,
                "peptide": peptide,
            })

            total += 1
            if len(batch_peptides) >= batch_size:
                yield batch_peptides, batch_meta
                batch_peptides, batch_meta = [], []

    if batch_peptides:
        yield batch_peptides, batch_meta

    print(f"[INFO] Total peptides generated (after filtering): {total}")


# ====================================================
# MAIN PIPELINE
# ====================================================

def main():
    print("[INFO] Loading MHCflurry predictor...")
    predictor = Class1PresentationPredictor.load()
    print("[INFO] Predictor loaded.")

    mapping = load_mapping(MAPPING_PATH)

    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    out = open(OUTPUT_CSV, "w", newline="")
    writer = csv.DictWriter(out, fieldnames=[
        "protein_id", "transcript_id", "gene_id", "gene_name",
        "start_position", "peptide", "allele",
        "affinity", "presentation_score", "rank",
    ])
    writer.writeheader()

    total_written = 0
    batch_index = 0

    for peptides, metas in generate_peptide_batches(
        FASTA_PATH, mapping, PEPTIDE_LENGTH, BATCH_SIZE
    ):
        batch_index += 1
        print(f"[INFO] Processing batch {batch_index}, size={len(peptides)}")

        for allele in ALLELES:
            print(f"  [INFO] Predicting allele: {allele}")

            # Modern MHCflurry API: alleles=[...]
            pred = predictor.predict(
                peptides=peptides,
                alleles=[allele],
            )

            # Robust column identification
            cols = {c.lower(): c for c in pred.columns}

            affinity_col = next(
                (cols[c] for c in cols if "affinity" in c or "ic50" in c),
                None
            )
            score_col = next(
                (cols[c] for c in cols if "presentation" in c or "score" in c),
                None
            )
            rank_col = next(
                (cols[c] for c in cols if "rank" in c or "percentile" in c),
                None
            )

            for i, meta in enumerate(metas):
                affinity = pred.iloc[i][affinity_col] if affinity_col else None
                score = pred.iloc[i][score_col] if score_col else None
                rank = pred.iloc[i][rank_col] if rank_col else None

                # Optional rank filtering
                if RANK_THRESHOLD is not None and rank is not None:
                    try:
                        if float(rank) > float(RANK_THRESHOLD):
                            continue
                    except Exception:
                        # If rank cannot be parsed, don't filter
                        pass

                writer.writerow({
                    **meta,
                    "allele": allele,
                    "affinity": affinity,
                    "presentation_score": score,
                    "rank": rank,
                })
                total_written += 1

        print(f"[INFO] Total rows written so far: {total_written}")

    out.close()
    print(f"[DONE] Final output: {OUTPUT_CSV}")
    print(f"[DONE] Total rows written: {total_written}")


if __name__ == "__main__":
    main()


