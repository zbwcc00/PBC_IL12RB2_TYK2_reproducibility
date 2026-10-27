# PBC IL12RB2–TYK2 reproducibility repository

This repository accompanies the manuscript **“Human genetics and liver single-cell analyses prioritize the IL12RB2–TYK2 axis for therapeutic testing in primary biliary cholangitis: a multi-cohort translational genomics study.”** It contains code, compact derived source data, figure-generation scripts, and audit records needed to inspect and reproduce the reported computational evidence.

Release: **v0.9.1**. Repository file and directory names use portable English ASCII paths. Analysis scripts derive the repository root from their own location and use repository-relative inputs and outputs.

## What is included

- Portable plotting scripts for Figures 1–7 and Supplementary Figures S1–S13.
- Compact, deidentified or aggregate figure-source data.
- Selected analysis scripts for fine-mapping, signal-level colocalization, external disease replication, donor-level single-cell analysis, bulk meta-analysis, proteomic robustness and TYK2 LD sensitivity.
- Data-source registry, figure manifests, software snapshots, random-seed registry and SHA-256 manifests.
- A source-project code audit that records selection decisions without requiring the original workstation layout.

## What is not included

Raw human sequencing files, complete GWAS summary statistics, full single-cell objects, complete external-QTL releases, third-party source-code mirrors, software environments, temporary downloads and QA rendering folders are intentionally excluded. Obtain these assets from their official accessions and follow the applicable data-use terms; see `DATA_AVAILABILITY.md`.

## Reproduce the figures

Python 3.12 is recommended.

```bash
python -m venv .venv
python -m pip install -r environment/requirements-public.txt
python code/00_release_audit/verify_repository.py
python code/04_figures/run_all_figures.py
python code/00_release_audit/verify_repository.py
```

Generated figures are written to `results/figures/`. The committed reference renders support visual and numerical comparison. PNG files are the primary pixel-audit format because PDF and SVG metadata can vary across backends.

## Reproduction levels

- **Level A:** figures can be redrawn directly from compact derived tables in this repository.
- **Level B:** full upstream analyses require downloading public or controlled-access source data from the accessions listed in `data/metadata/data_source_registry.tsv` and placing them under `data/raw/` as described in `REPRODUCIBILITY.md`.

## Citation

Use `CITATION.cff`. The repository is prepared for GitHub release and subsequent Zenodo archiving. The final GitHub URL and Zenodo DOI must be inserted only after those records exist.
