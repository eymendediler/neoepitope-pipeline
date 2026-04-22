#!/usr/bin/env python3
import os, re, gzip
import pandas as pd
import numpy as np

EXPR_GZ = "GTEx_Analysis_v10_RSEMv1.3.3_transcripts_tpm.txt.gz"
SAMPLE_TISSUE = "sample_tissue.tsv"
MAP_GZ = "ensg_enst_ensp_geneName_v105.tsv.gz"
OUTDIR = "results/mdt_by_tissue"
MIN_TOP1_MEDIAN = 0.0

def safe_name(s: str) -> str:
    s = re.sub(r"[^\w]+", "_", s.strip())
    return re.sub(r"_+", "_", s).strip("_")

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    st = pd.read_csv(SAMPLE_TISSUE, sep="\t")
    tissue_to_samples = st.groupby("SMTSD")["SAMPID"].apply(list).to_dict()

    # expression header
    with gzip.open(EXPR_GZ, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    all_samples = header[1:]
    all_samples_set = set(all_samples)

    tissue_samples_idx = {}
    tissue_nsamples = {}
    for tissue, samples in tissue_to_samples.items():
        s = [x for x in samples if x in all_samples_set]
        if not s:
            continue
        tissue_nsamples[tissue] = len(s)
        tissue_samples_idx[tissue] = [header.index(x) for x in s]

    if not tissue_samples_idx:
        raise SystemExit("No overlapping samples between expression matrix and SampleAttributes mapping.")

    # ENST -> ENSG
    m = pd.read_csv(MAP_GZ, sep="\t", compression="gzip", low_memory=False)
    m = m[["Transcript stable ID","Gene stable ID"]].dropna().drop_duplicates()
    m.columns = ["enst","ensg"]
    m["enst"] = m["enst"].astype(str).str.split(".").str[0]
    enst2ensg = dict(zip(m["enst"], m["ensg"]))

    best = {t:{} for t in tissue_samples_idx.keys()}

    def update_gene(tissue, ensg, enst, med, mean):
        d = best[tissue]
        if ensg not in d:
            d[ensg] = [enst, med, mean, None, -1.0, -1.0]
            return
        t1, m1, mu1, t2, m2, mu2 = d[ensg]
        if (med > m1) or (med == m1 and mean > mu1):
            d[ensg] = [enst, med, mean, t1, m1, mu1]
        elif (med > m2) or (med == m2 and mean > mu2):
            d[ensg][3:] = [enst, med, mean]

    with gzip.open(EXPR_GZ, "rt") as f:
        _ = f.readline()
        n = 0
        for line in f:
            n += 1
            parts = line.rstrip("\n").split("\t")
            enst = parts[0].split(".")[0]
            ensg = enst2ensg.get(enst)
            if not ensg:
                continue

            for tissue, idxs in tissue_samples_idx.items():
                vals = []
                for j in idxs:
                    v = parts[j]
                    if v == "NA" or v == "":
                        continue
                    x = float(v)
                    if x < 0.001:
                        x = 0.0
                    vals.append(x)
                if not vals:
                    continue
                med = float(np.median(vals))
                mean = float(np.mean(vals))
                update_gene(tissue, ensg, enst, med, mean)

            if n % 200000 == 0:
                print(f"[INFO] processed {n:,} transcripts")

    for tissue, d in best.items():
        rows = []
        for ensg, (t1, m1, mu1, t2, m2, mu2) in d.items():
            if m1 < MIN_TOP1_MEDIAN:
                continue
            delta = m1 - (m2 if t2 else 0.0)
            rows.append([tissue, tissue_nsamples[tissue], ensg, t1, m1, mu1, t2, (m2 if t2 else np.nan), (mu2 if t2 else np.nan), delta])

        out = pd.DataFrame(rows, columns=[
            "tissue","n_samples","ensg","mdt_enst","top1_median_tpm","top1_mean_tpm",
            "top2_enst","top2_median_tpm","top2_mean_tpm","delta_median_tpm"
        ]).sort_values("ensg")

        fn = os.path.join(OUTDIR, f"{safe_name(tissue)}.tsv")
        out.to_csv(fn, sep="\t", index=False)

    print("[DONE] Wrote per-tissue MDT files to", OUTDIR)

if __name__ == "__main__":
    main()
