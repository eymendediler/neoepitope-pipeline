import argparse
import pandas as pd
import sys

PATH="/Users/eymen/Desktop/v105/v105_second/unique_peptides_with_gtex_plus.tsv"

def main():
    ap = argparse.ArgumentParser(description="Query unique peptides table by peptide and/or ENSG.")
    ap.add_argument("--peptide", help="Peptide sequence (e.g., RSHYERIIY)")
    ap.add_argument("--ensg", help="Gene ID ENSG... (e.g., ENSG00000012817)")
    ap.add_argument("--limit", type=int, default=50, help="Max rows to show (default 50)")
    args = ap.parse_args()

    if not args.peptide and not args.ensg:
        print("Provide --peptide and/or --ensg", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(PATH, sep="\t")

    q = df
    if args.peptide:
        q = q[q["peptide"] == args.peptide.strip()]
    if args.ensg:
        # tabloda hem gene_id hem gtex_gene_id var; ikisine de bak
        ensg = args.ensg.strip()
        if "gtex_gene_id" in q.columns:
            q = q[(q["gene_id"] == ensg) | (q["gtex_gene_id"] == ensg)]
        else:
            q = q[q["gene_id"] == ensg]

    if len(q) == 0:
        print("No matches.")
        sys.exit(0)

    # okunur kolon sırası
    preferred = [
        "peptide",
        "gtex_gene_id",
        "gtex_most_dominant_tx",
        "gene_id",
        "transcript_id",
        "n_minor_transcripts",
        "minor_specificity",
        "unique_reason",
        "isoforms",
    ]
    cols = [c for c in preferred if c in q.columns]
    q2 = q[cols].drop_duplicates()

    print("Matches:", len(q2))
    print(q2.head(args.limit).to_string(index=False))

    if len(q2) > args.limit:
        print(f"\n... truncated, showing first {args.limit} rows. Use --limit to increase.")

if __name__ == "__main__":
    main()
