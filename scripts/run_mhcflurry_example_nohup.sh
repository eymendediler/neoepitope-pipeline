#!/usr/bin/env bash
set -euo pipefail

SCRIPT="scripts/mhcflurry_pipeline_v105_full.py"
LOG="results/mhcflurry_8to12.log"

mkdir -p results

nohup python "${SCRIPT}" > "${LOG}" 2>&1 &
echo "MHCflurry started in background. Log: ${LOG}"
