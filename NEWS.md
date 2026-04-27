# proteoforge 0.1.0

## Initial Release

* Package skeleton with core utilities.
* Standardized `ProteoForgeData` S4 class.
* Generic LFQ matrix importer with auto-detection.
* Sample metadata validation (batch confounding, pairing, replicate checks).
* Preprocessing pipeline: filter → log2 → normalize → impute (mixed MNAR/MCAR strategy).
* QC module: missingness, ID counts, PCA, correlation, CV.
* Differential abundance via limma + DEqMS.
* Volcano, heatmap, Venn/UpSet plots with `theme_proteoforge()`.
* Excel, Word, PPTX report export.
* Full `targets` pipeline.
* Production Shiny app.
* Docker images.
* GitHub Actions CI.
