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
OUTPUT_DIR = "/Users/eymen/Desktop/v105/mhcflurry_multi_length_results"

PEPTIDE_LENGTHS = [8, 9, 10, 11, 12]
BATCH_SIZE = 20000

# Set None to keep all predictions
RANK_THRESHOLD: Optional[float] = 2.0

ALLELES = [
    "HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*24:02",
    "HLA-A*26:01", "HLA-A*30:01", "HLA-A*23:01", "HLA-A*68:01", "HLA-A*33:01",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*27:05", "HLA-B*35:01",
    "HLA-B*38:01", "HLA-B*40:01", "HLA-B*44:02", "HLA-B*51:01", "HLA-B*18:01",
]

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


# ====================================================
# MAPPING FILE LOADER
# ====================================================

def load_mapping(path: str) -> Dict[str, Dict[str, Any]]:
    print(f"[INFO] Loading mapping file: {path}")
    df = pd.read_csv(path, sep="\t", low_memory=False)

    protein_col = "Protein stable ID"
    transcript_col = "Transcript stable ID"
    gene_col = "Gene stable ID"
    gene_name_col = "Gene name"

    mapping: Dict[str, Dict[str, Any]] = {}

    for _, row in df.iterrows():
        ensp = str(row[protein_col])                    # versioned
        ensp_novers = ensp.split(".")[0]                # non-versioned

        transcript_id = str(row[transcript_col]).split(".")[0]
        gene_id = str(row[gene_col]).split(".")[0]
        gene_name = str(row[gene_name_col])

        info = {
            "protein_id": ensp_novers,
            "transcript_id": transcript_id,
            "gene_id": gene_id,
            "gene_name": gene_name,
        }

        mapping[ensp] = info
        mapping[ensp_novers] = info

    print(f"[INFO] Mapping loaded. Protein keys: {len(mapping)}")
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
    print(f"[INFO] Reading FASTA for peptide length {peptide_length}: {fasta_path}")

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
            "protein_id": ensp_novers,
            "transcript_id": None,
            "gene_id": None,
            "gene_name": None,
        }

        for i in range(len(seq) - peptide_length + 1):
            peptide = seq[i:i + peptide_length]

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
                "peptide_length": peptide_length,
            })

            total += 1
            if len(batch_peptides) >= batch_size:
                yield batch_peptides, batch_meta
                batch_peptides, batch_meta = [], []

    if batch_peptides:
        yield batch_peptides, batch_meta

    print(f"[INFO] Total peptides generated for length {peptide_length}: {total}")


# ====================================================
# RUN FOR A SINGLE PEPTIDE LENGTH
# ====================================================

def run_for_length(
    predictor: Class1PresentationPredictor,
    mapping: Dict[str, Dict[str, Any]],
    peptide_length: int
) -> str:
    output_csv = os.path.join(OUTPUT_DIR, f"mhcflurry_results_len{peptide_length}.csv")
    print(f"[INFO] Output file for {peptide_length}-mer: {output_csv}")

    out = open(output_csv, "w", newline="")
    writer = csv.DictWriter(out, fieldnames=[
        "protein_id", "transcript_id", "gene_id", "gene_name",
        "start_position", "peptide", "peptide_length", "allele",
        "affinity", "presentation_score", "rank",
    ])
    writer.writeheader()

    total_written = 0
    batch_index = 0

    for peptides, metas in generate_peptide_batches(
        FASTA_PATH, mapping, peptide_length, BATCH_SIZE
    ):
        batch_index += 1
        print(f"[INFO] Processing {peptide_length}-mer batch {batch_index}, size={len(peptides)}")

        for allele in ALLELES:
            print(f"  [INFO] Predicting allele: {allele}")

            pred = predictor.predict(
                peptides=peptides,
                alleles=[allele],
            )

            cols = {c.lower(): c for c in pred.columns}

            affinity_col = next(
                (cols[c] for c in cols if "affinity" in c or "ic50" in c),
                None
            )
            score_col = next(
                (cols[c] for c in cols if "presentation" in c or c == "score" or "presentation_score" in c),
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

                if RANK_THRESHOLD is not None and rank is not None:
                    try:
                        if float(rank) > float(RANK_THRESHOLD):
                            continue
                    except Exception:
                        pass

                writer.writerow({
                    **meta,
                    "allele": allele,
                    "affinity": affinity,
                    "presentation_score": score,
                    "rank": rank,
                })
                total_written += 1

        print(f"[INFO] Total rows written so far for {peptide_length}-mer: {total_written}")

    out.close()
    print(f"[DONE] Final output for {peptide_length}-mer: {output_csv}")
    print(f"[DONE] Total rows written for {peptide_length}-mer: {total_written}")
    return output_csv


# ====================================================
# MERGE OUTPUTS
# ====================================================

def merge_outputs(files: List[str], merged_path: str) -> None:
    print("[INFO] Merging output files...")
    first = True

    with open(merged_path, "w", newline="") as fout:
        writer = None

        for f in files:
            print(f"  [INFO] Reading {f}")
            with open(f, "r", newline="") as fin:
                reader = csv.DictReader(fin)
                if first:
                    writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
                    writer.writeheader()
                    first = False
                for row in reader:
                    writer.writerow(row)

    print(f"[DONE] Merged output written to: {merged_path}")


# ====================================================
# MAIN PIPELINE
# ====================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("[INFO] Loading MHCflurry predictor...")
    predictor = Class1PresentationPredictor.load()
    print("[INFO] Predictor loaded.")

    mapping = load_mapping(MAPPING_PATH)

    output_files = []
    for peptide_length in PEPTIDE_LENGTHS:
        print("=" * 70)
        print(f"[INFO] Running peptide length: {peptide_length}")
        print("=" * 70)
        out_csv = run_for_length(
            predictor=predictor,
            mapping=mapping,
            peptide_length=peptide_length,
        )
        output_files.append(out_csv)

    merged_path = os.path.join(OUTPUT_DIR, "mhcflurry_results_ensembl_v105_8to12mer.csv")
    merge_outputs(output_files, merged_path)

    print("[DONE] All peptide lengths completed.")
    print(f"[DONE] Per-length files are in: {OUTPUT_DIR}")
    print(f"[DONE] Combined file: {merged_path}")


if __name__ == "__main__":
    main()
