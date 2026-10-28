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

Install the upstream-analysis dependencies with:

```bash
python -m pip install -r environment/requirements-upstream.txt
```

The following workflows close the main public-data gaps in the prior archive:

```bash
python code/02_single_cell/07_cardon_antigen_context.py
python code/02_single_cell/08_andrews_external_donor_pseudobulk.py --download
python code/03_bulk_proteomics_pharmacology/07_mesa_cross_platform_pqtl.py
```

- `07_cardon_antigen_context.py` expects the GSE270739 processed metadata, expression matrix and clinical supplement under `data/raw/single_cell/gse270739/` and reproduces the five-patient paired antigen-context tables.
- `08_andrews_external_donor_pseudobulk.py` can download the two version-locked CELLxGENE H5AD assets, verifies their byte sizes, aggregates raw counts by donor and reproduces the scRNA and snRNA IL12RB2 pseudobulk tables.
- `07_mesa_cross_platform_pqtl.py` expects the Zenodo 17642645 Olink and SomaScan cis/trans summary files under `data/raw/molecular_qtl/mesa/` and reproduces the 8,647-variant allele-aligned comparison and single-causal ABF sensitivity analysis.

The analysis scripts in modules 01–03 derive the repository root from the script location. Raw inputs are read from `data/raw/`, intermediate files are written to `work/`, and compact release assets are stored under `data/derived/analysis/` and `data/derived/figures/`. No workstation-specific drive letter or path substitution is required.

The full upstream workflows remain accession-dependent: some inputs are controlled-access, some public releases must be downloaded separately, and selected external software must be installed by the user. The compact Level A figure workflow does not require these raw assets. Full Finnish multi-signal `coloc.susie` remains out of scope because complete study-specific signed FinnGen LD was not available; the repository retains the stated single-causal regional compatibility analysis instead.

## Determinism

Random seeds are recorded in `environment/05_random_seed_registry.tsv`. File integrity is recorded in `manifests/SHA256SUMS.tsv`. Numerical equality should be checked for TSV outputs. PNG pixel identity is recorded where the rendering environment is unchanged; visual or semantic equivalence is appropriate for PDF and SVG outputs because backend metadata can vary.
