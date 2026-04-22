# Neoepitope Prediction Pipeline

This repository contains a neoepitope prioritization workflow based on:
- Most Dominant Transcript (MDT) selection from GTEx
- peptide generation from MDT-derived transcripts/proteins
- peptide-HLA binding prediction with NetMHCpan and MHCflurry
- immunogenicity estimation with PRIME
- integration of final candidate tables

## Repository structure
- `scripts/` executable and helper scripts
- `docs/` workflow documentation
- `metadata/` reference metadata files
- `data/final_tables/` selected final result tables
- `database/` SQLite schema and build scripts

## Main workflow
1. Select Most Dominant Transcripts (MDT)
2. Generate MDT-derived peptides
3. Run NetMHCpan
4. Run MHCflurry
5. Run PRIME
6. Merge and prioritize candidates

## Documentation
- `docs/mdt_workflow.md`
- `docs/netmhcpan_workflow.md`
- `docs/mhcflurry_workflow.md`
- `docs/prime_workflow.md`
