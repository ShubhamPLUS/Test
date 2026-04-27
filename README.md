# ProteoForge

**Reproducible, publication-grade LFQ proteomics and phosphoproteomics analysis.**

[![R CMD check](https://github.com/shubhamplus/test/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/shubhamplus/test/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## What it does

ProteoForge is an R package and Shiny application for end-to-end label-free
quantitative (LFQ) proteomics and phosphoproteomics analysis. It supports
every major upstream search engine and provides a complete, reproducible
pipeline from raw quantification matrices to publication-ready reports.

## Supported inputs

| Software | Files |
|----------|-------|
| DIA-NN | `report.pg_matrix.tsv`, `report.tsv`, `report.parquet` |
| Spectronaut | Long-format TSV/CSV/parquet |
| FragPipe / MSFragger / IonQuant | `combined_protein.tsv` |
| MaxQuant | `proteinGroups.txt`, `Phospho (STY)Sites.txt` |
| Proteome Discoverer | Proteins/Peptide Groups TSV export |
| Skyline | Custom report export |
| Generic matrix | Feature × sample matrix + metadata |

## Installation

```r
# Install from GitHub (requires remotes or pak)
pak::pkg_install("shubhamplus/test")

# Or install all dependencies with renv
Rscript scripts/install.R
```

## Quick start — CLI

```bash
Rscript scripts/run_demo.R
# or with your own data:
Rscript scripts/run_pipeline.R \
  --config inst/config/default_lfq.yml \
  --quant  path/to/protein_matrix.tsv \
  --metadata path/to/metadata.tsv \
  --output results \
  --organism human
```

## Quick start — Shiny

```bash
Rscript scripts/launch_shiny.R
```

Then open http://localhost:3838 in your browser.

## Output structure

```
results/<project>_<timestamp>/
├── 00_logs/          run.log, warnings.tsv, session_info.txt
├── 02_qc/            QC tables and plots
├── 03_preprocessing/ filtering, normalization, imputation
├── 04_statistics/    differential abundance results
├── 05_overlap/       Venn / UpSet diagrams
├── 06_enrichment/    ORA and GSEA results
├── 07_networks/      STRING / OmniPath networks
├── 08_phospho/       phosphosite results (if applicable)
├── 09_reports/       final_report.docx, .html, .pptx
├── 10_excel/         complete_analysis_results.xlsx
└── 11_serialized/    R objects for downstream analysis
```

## Statistical methods

**Default engine:** limma (moderated t-statistic) + DEqMS (peptide-count
weighted variance correction). Alternative engines: msqrob2, MSstats,
Welch t-test.

**Default imputation:** mixed strategy — MinProb (downshift) for
missing-not-at-random (MNAR) proteins, kNN for missing-completely-at-random
(MCAR) proteins, classified per protein per condition.

**Default normalization:** median centering (configurable: quantile, VSN,
cyclic loess).

## Phosphoproteomics

When `analysis_level: phosphosite` is set:
- Class-1 site filter (localization probability ≥ 0.75)
- Proper MaxQuant multiplicity expansion
- Parent-protein abundance correction
- Kinase activity inference (KSEA + KEA3 + decoupleR consensus)
- Motif enrichment (rmotifx)
- PhosphoSitePlus annotation overlay (download script bundled; data not included)

## Reproducibility

Every run writes:
- `parameter_hash.txt` — sha256 of all analysis parameters
- `session_info.txt` — full R session information
- `config_used.yml` — exact configuration used
- All database versions and access dates in result tables

## Citation

If you use ProteoForge in your research, please cite:
```
ProteoForge: Reproducible LFQ proteomics analysis platform.
https://github.com/shubhamplus/test
```

## Known limitations

- TMT/iTRAQ/SILAC quantification is out of scope for v1.
- SaaS frontend, billing, and multi-tenant authentication are out of scope for v1.
- AI narrative generation is stubbed; the insertion point is clean for future integration.
- PhosphoSitePlus data requires registration and must be downloaded separately.
- KEA3 kinase activity requires internet access.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Please follow the
[Code of Conduct](CODE_OF_CONDUCT.md).
