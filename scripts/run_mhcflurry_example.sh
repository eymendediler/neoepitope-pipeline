#!/usr/bin/env bash
set -euo pipefail

# Example wrapper for MHCflurry pipeline

SCRIPT="scripts/mhcflurry_pipeline_v105_full.py"
LOG="results/mhcflurry_8to12.log"

mkdir -p results

echo "Running MHCflurry pipeline..."
echo "Script: ${SCRIPT}"
echo "Log: ${LOG}"

python "${SCRIPT}" > "${LOG}" 2>&1

echo "Done."
