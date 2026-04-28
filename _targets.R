## ProteoForge targets pipeline
## Run with: targets::tar_make()
## Visualise: targets::tar_visnetwork()

library(targets)
library(tarchetypes)

# Source all package functions (dev mode without install)
if (!requireNamespace("proteoforge", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(proteoforge)
}

# ── Pipeline configuration ─────────────────────────────────────────────────
tar_option_set(
  packages   = c("proteoforge", "data.table"),
  format     = "rds",
  error      = "continue",   # don't halt on single-target failure
  memory     = "transient",  # free memory after each target
  garbage_collection = TRUE
)

# Load user config (override via Sys.setenv(PF_CONFIG=...))
config_path <- Sys.getenv("PF_CONFIG", unset = "config.yml")
config <- if (file.exists(config_path)) {
  yaml::read_yaml(config_path)
} else {
  yaml::read_yaml(
    system.file("config", "default_lfq.yml", package = "proteoforge")
  )
}

# ── Targets list ──────────────────────────────────────────────────────────
list(

  # ── 1. Configuration hash ────────────────────────────────────────────────
  tar_target(config_hash, hash_params(config)),

  # ── 2. Import data ───────────────────────────────────────────────────────
  tar_target(
    quant_file,
    config$input$quant_file,
    format = "file"
  ),
  tar_target(
    metadata_file,
    config$input$metadata_file,
    format = "file"
  ),
  tar_target(
    pfd_raw,
    {
      sw <- config$input$software %||% "auto"
      if (sw == "auto") sw <- detect_software(quant_file)
      switch(sw,
        diann        = import_diann(quant_file, metadata_file, config = config),
        spectronaut  = import_spectronaut(quant_file, metadata_file, config = config),
        fragpipe     = import_fragpipe(quant_file, metadata_file, config = config),
        maxquant     = import_maxquant(quant_file, metadata_file, config = config),
        pd           = import_pd(quant_file, metadata_file, config = config),
        skyline      = import_skyline(quant_file, metadata_file, config = config),
        import_generic(quant_file, metadata_file, config = config)
      )
    }
  ),

  # ── 3. QC metrics ────────────────────────────────────────────────────────
  tar_target(qc_metrics, compute_qc_metrics(pfd_raw)),
  tar_target(qc_plots_out, {
    dirs <- make_output_dirs(config)
    generate_qc_plots(pfd_raw, qc_metrics, dirs, config)
  }),

  # ── 4. Preprocessing ─────────────────────────────────────────────────────
  tar_target(pfd_filtered,  filter_features(pfd_raw, config)),
  tar_target(pfd_log2,      log2_transform(pfd_filtered)),
  tar_target(pfd_norm,      normalise_matrix(pfd_log2,
                              method = config$preprocessing$normalisation %||% "median")),
  tar_target(missingness,   classify_missingness(pfd_norm, config)),
  tar_target(pfd_imputed,   impute_missing(pfd_norm,
                              method = config$preprocessing$imputation %||% "mixed",
                              config = config)),

  # ── 5. Differential abundance ─────────────────────────────────────────────
  tar_target(design_mat,   build_design_matrix(pfd_imputed, config)),
  tar_target(contrast_mat, build_contrast_matrix(design_mat, config)),
  tar_target(
    stat_results,
    run_limma_deqms(pfd_imputed, design_mat, contrast_mat, config)
  ),
  tar_target(contrast_tables, {
    dirs <- make_output_dirs(config)
    write_contrast_tables(stat_results, dirs, config)
  }),

  # ── 6. Plots ─────────────────────────────────────────────────────────────
  tar_target(volcano_plots, {
    dirs <- make_output_dirs(config)
    lapply(names(stat_results), function(cname) {
      plot_volcano(stat_results[[cname]], contrast = cname, dirs = dirs)
    })
  }),
  tar_target(heatmap_plot, {
    dirs <- make_output_dirs(config)
    plot_heatmap(pfd_imputed, stat_results[[1]], dirs = dirs)
  }),

  # ── 7. Enrichment ────────────────────────────────────────────────────────
  tar_target(
    enrichment_results,
    {
      if (isFALSE(config$enrichment$enabled)) return(NULL)
      lapply(stat_results, function(dt) {
        list(
          ora  = run_ora(dt, pfd_imputed, config),
          gsea = run_gsea(dt, config)
        )
      })
    }
  ),

  # ── 8. Networks ──────────────────────────────────────────────────────────
  tar_target(
    network_results,
    {
      if (isFALSE(config$networks$enabled)) return(NULL)
      dirs <- make_output_dirs(config)
      list(
        string   = tryCatch(run_string_network(stat_results[[1]], dirs, config), error = function(e) NULL),
        omnipath = tryCatch(run_omnipath(stat_results[[1]], config), error = function(e) NULL)
      )
    }
  ),

  # ── 9. Excel report ──────────────────────────────────────────────────────
  tar_target(excel_report, {
    dirs <- make_output_dirs(config)
    export_excel(stat_results, qc_metrics, enrichment_results, dirs, config)
  }),

  # ── 10. Word / PPTX ──────────────────────────────────────────────────────
  tar_target(word_report, {
    dirs <- make_output_dirs(config)
    export_word(stat_results, qc_metrics, dirs, config)
  }),
  tar_target(pptx_report, {
    dirs <- make_output_dirs(config)
    export_pptx(stat_results, dirs, config)
  }),

  # ── 11. Quarto HTML ──────────────────────────────────────────────────────
  tar_target(quarto_report, {
    dirs <- make_output_dirs(config)
    render_quarto_report(
      stat_results      = stat_results,
      qc_metrics        = qc_metrics,
      enrichment_results= enrichment_results,
      dirs              = dirs,
      config            = config
    )
  }),

  # ── 12. AI summary JSON ──────────────────────────────────────────────────
  tar_target(ai_summary, {
    dirs <- make_output_dirs(config)
    write_ai_summary_json(
      stat_results       = stat_results,
      qc_metrics         = qc_metrics,
      enrichment_results = enrichment_results,
      config_hash        = config_hash,
      dirs               = dirs
    )
  }),

  # ── 13. Param hash (reproducibility) ─────────────────────────────────────
  tar_target(param_hash_file, {
    dirs <- make_output_dirs(config)
    write_param_hash(config, dirs$results_dir)
  })
)
