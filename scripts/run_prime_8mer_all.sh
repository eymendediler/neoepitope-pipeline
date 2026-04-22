#!/bin/bash

PRIME_BIN="/Users/eymen/Desktop/prime_pipeline/PRIME2.1/PRIME"
MIX_BIN="/Users/eymen/Desktop/prime_pipeline/MixMHCpred/MixMHCpred"
IN_DIR="/Users/eymen/Desktop/v105/prime_input_chunks"
OUT_DIR="/Users/eymen/Desktop/v105/prime_output_chunks_rerun"
LOG_DIR="/Users/eymen/Desktop/v105/prime_logs_rerun"

ALLELES="A0201,A0301,A2402,A2601,A3001,B0801,B1501,B2705,B3501,B3801,B4001,B4402,B5101"

mkdir -p "$OUT_DIR" "$LOG_DIR"

for f in "$IN_DIR"/peptides_8mer_chunk_*.txt; do
    [ -e "$f" ] || continue
    base=$(basename "$f" .txt)

    echo "Running $base"

    "$PRIME_BIN" \
        -i "$f" \
        -o "$OUT_DIR/${base}_prime.txt" \
        -a "$ALLELES" \
        -mix "$MIX_BIN" \
        > "$LOG_DIR/${base}.log" 2>&1
done
