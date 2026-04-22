#!/usr/bin/env bash
set -euo pipefail

# Example NetMHCpan wrapper script
# Edit these paths and alleles to match your local setup.

INPUT="data/example_input/v105_8mer_from_mdt_clean.txt"
OUTDIR="results/netmhcpan"
OUTFILE="${OUTDIR}/netmhcpan_results.txt"

# Use NetMHCpan allele formatting that matches your installation
ALLELES="HLA-A02:01,HLA-A11:01,HLA-B07:02"

mkdir -p "${OUTDIR}"

echo "Running NetMHCpan..."
echo "Input: ${INPUT}"
echo "Alleles: ${ALLELES}"
echo "Output: ${OUTFILE}"

# Example command; adjust flags to your installation/version
netMHCpan -p "${INPUT}" -a "${ALLELES}" -BA > "${OUTFILE}"

echo "Done."
