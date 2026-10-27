# Public-release asset audit

## Decision

The 306-GB working project is not suitable for direct publication. A whitelist release was constructed instead. The release includes original analysis code selected for the manuscript, portable figure scripts, compact derived data and final renders. It excludes raw/controlled human data, complete GWAS and QTL releases, software installations, third-party repositories, caches and QA artifacts.

## Key risk controls

1. No BAM, CRAM, VCF, BCF, RDS, QS, H5AD or Matrix Market count object is included.
2. HRA008003 cell barcodes were removed from the compact UMAP table; the annotation table uses synthetic row identifiers.
3. The 373-MB PBC GWAS is replaced by the deterministic subset used by the final Manhattan plot.
4. Ten full signed-LD matrices are replaced by compact credible-set LD blocks; the exact PBC04 lead-r² vector is embedded in the Figure 2 regional table.
5. The 148-MB Figure S12 TIFF and all files larger than 100 MB are excluded.
6. Third-party code mirrors are not relicensed or redistributed.
7. Analysis and figure scripts derive the repository root from their own location and use English repository-relative paths.

## Remaining publication action

Before a remote GitHub release, the authors should choose a repository URL, add it to `CITATION.cff`, create a tagged release and archive that tag in Zenodo. No DOI is fabricated in this package.
