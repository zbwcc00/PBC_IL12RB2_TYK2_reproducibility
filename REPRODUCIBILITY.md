# Reproducibility guide

## Scope

The public package separates figure-level reproduction from full upstream reanalysis. This distinction prevents large or access-controlled data from being silently redistributed while retaining an auditable path from source accession to reported result.

## Level A: compact figure reproduction

Run `python code/04_figures/run_all_figures.py`. The scripts read only repository-relative paths. Figure 5 uses a barcode-free cell-level visualization table; inferential panels remain donor-level. Supplementary Figure S2 uses compact credible-set LD blocks exported from the exact signed matrices used in the locked analysis.

## Level B: full upstream analysis

1. Download the source data listed in `data/metadata/data_source_registry.tsv`.
2. Place downloaded inputs under the English relative layout in `data/raw/`; `config/paths.example.env` documents the expected top-level locations.
3. Run the selected scripts indexed in `manifests/analysis_code_index.tsv` in numerical order within each module.
4. Compare outputs against the compact derived tables and the locked evidence statements.

The analysis scripts in modules 01–03 derive the repository root from the script location. Raw inputs are read from `data/raw/`, intermediate files are written to `work/`, and compact release assets are stored under `data/derived/analysis/` and `data/derived/figures/`. No workstation-specific drive letter or path substitution is required.

The full upstream workflows remain accession-dependent: some inputs are controlled-access, some public releases must be downloaded separately, and selected external software must be installed by the user. The compact Level A figure workflow does not require these raw assets.

## Determinism

Random seeds are recorded in `environment/05_random_seed_registry.tsv`. File integrity is recorded in `manifests/SHA256SUMS.tsv`. Numerical equality should be checked for TSV outputs. PNG pixel identity is recorded where the rendering environment is unchanged; visual or semantic equivalence is appropriate for PDF and SVG outputs because backend metadata can vary.
