# Data availability

The repository distributes only compact derived data required for figure reproduction and audit. Raw human sequencing data and large external datasets are not redistributed.

The versioned reproducibility archive is available from Zenodo under concept DOI [10.5281/zenodo.23690864](https://doi.org/10.5281/zenodo.23690864), with source control and release assets hosted at [GitHub](https://github.com/zbwcc00/PBC_IL12RB2_TYK2_reproducibility).

| Resource | Accession/version | Repository use | Redistribution decision |
|---|---|---|---|
| PBC GWAS | GCST90061440 | Discovery association and fine-mapping | Full file excluded; deterministic plotting subset included |
| OneK1K/eQTL Catalogue | eQTL Catalogue r8 beta | Single-cell eQTL fine-mapping and colocalization | Official releases must be downloaded by the user |
| 1000 Genomes | Phase 3 EUR, n=503 | External signed LD | Full genotype/LD panels excluded; compact derived LD blocks included |
| FinnGen | R12 CHIRBIL_PRIM | Independent disease replication | Regional derived summaries included; official source remains authoritative |
| HRA008003 | GSA-Human HRA008003 | PBC/control liver single-cell analysis | BAM and full objects excluded; barcodes removed from compact visualization data |
| HRA002347/GSE270739/CELLxGENE | See source registry | External single-cell context | Aggregate donor-level results included; public upstream workflows provided |
| GEO bulk datasets | GSE119600, GSE159676, GSE79850, GSE138888 | Bulk validation and pathway meta-analysis | Derived effect tables included |
| PRIDE | PXD083842 | Proteomic robustness and confounding audit | Derived results included; processed matrices remain at PRIDE |
| MESA Olink/SomaScan | Zenodo 17642645 | Cross-platform pQTL comparison | Allele-aligned derived cis table and upstream analysis code included; source files remain at Zenodo |

All users must comply with the original repositories' consent, access and redistribution requirements. The MIT license in this repository applies to original code only, not to third-party data.
