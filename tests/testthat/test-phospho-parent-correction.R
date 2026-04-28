test_that("parent_protein_correction returns corrected_pfd and correction_stats", {
  set.seed(42)
  n_sites   <- 20L
  n_samples <- 8L
  snames    <- paste0("S", seq_len(n_samples))

  # Build a minimal phospho ProteoForgeData
  ph_mat <- matrix(rnorm(n_sites * n_samples, mean = 20, sd = 1),
                   nrow = n_sites, ncol = n_samples,
                   dimnames = list(paste0("P1_S", seq_len(n_sites)), snames))

  ph_fmd <- data.table::data.table(
    feature_id = rownames(ph_mat),
    ProteinID  = "P1",
    GeneSymbol = "GENE1"
  )
  ph_smd <- data.table::data.table(
    sample_id = snames,
    condition = rep(c("A", "B"), each = n_samples / 2)
  )

  ph_pfd <- ProteoForgeData(
    raw_matrix       = ph_mat,
    feature_metadata = ph_fmd,
    sample_metadata  = ph_smd,
    analysis_level   = "phosphosite",
    quant_type       = "LFQ",
    source_software  = "generic"
  )

  # Build a paired protein ProteoForgeData with same parent protein
  pr_mat <- matrix(rnorm(n_samples, mean = 20, sd = 0.5), nrow = 1,
                   dimnames = list("P1", snames))
  pr_fmd <- data.table::data.table(feature_id = "P1", GeneSymbol = "GENE1")
  pr_pfd <- ProteoForgeData(
    raw_matrix       = pr_mat,
    feature_metadata = pr_fmd,
    sample_metadata  = ph_smd,
    analysis_level   = "protein",
    quant_type       = "LFQ",
    source_software  = "generic"
  )

  result <- parent_protein_correction(ph_pfd, pr_pfd)

  expect_named(result, c("corrected_pfd", "correction_stats"))
  expect_s4_class(result$corrected_pfd, "ProteoForgeData")
  expect_true(result$corrected_pfd@parameters$parent_corrected)

  stats_dt <- result$correction_stats
  expect_true(data.table::is.data.table(stats_dt))
  expect_true(all(c("site_id", "prot_id", "r_squared", "n_obs") %in% names(stats_dt)))
  expect_true(nrow(stats_dt) > 0)
  expect_true(all(stats_dt$r_squared >= 0 & stats_dt$r_squared <= 1))
})

test_that("parent_protein_correction stops with < 3 shared samples", {
  n_samples <- 2L
  snames    <- paste0("S", seq_len(n_samples))

  ph_mat <- matrix(rnorm(4), nrow = 2,
                   dimnames = list(c("P1_S1", "P1_S2"), snames))
  ph_fmd <- data.table::data.table(feature_id = rownames(ph_mat), ProteinID = "P1")
  ph_smd <- data.table::data.table(sample_id = snames, condition = c("A", "B"))

  ph_pfd <- ProteoForgeData(
    raw_matrix = ph_mat, feature_metadata = ph_fmd,
    sample_metadata = ph_smd, analysis_level = "phosphosite",
    quant_type = "LFQ", source_software = "generic"
  )
  pr_mat <- matrix(rnorm(2), nrow = 1, dimnames = list("P1", snames))
  pr_fmd <- data.table::data.table(feature_id = "P1")
  pr_pfd <- ProteoForgeData(
    raw_matrix = pr_mat, feature_metadata = pr_fmd,
    sample_metadata = ph_smd, analysis_level = "protein",
    quant_type = "LFQ", source_software = "generic"
  )

  expect_error(parent_protein_correction(ph_pfd, pr_pfd), "< 3 shared samples")
})

test_that("parent_protein_correction handles missing parent protein gracefully", {
  set.seed(1)
  n_samples <- 6L
  snames    <- paste0("S", seq_len(n_samples))

  ph_mat <- matrix(rnorm(12), nrow = 2,
                   dimnames = list(c("MISSING_S1", "MISSING_S2"), snames))
  ph_fmd <- data.table::data.table(
    feature_id = rownames(ph_mat),
    ProteinID  = "NOTEXIST"
  )
  ph_smd <- data.table::data.table(sample_id = snames,
                                    condition = rep(c("A","B"), 3))

  ph_pfd <- ProteoForgeData(
    raw_matrix = ph_mat, feature_metadata = ph_fmd,
    sample_metadata = ph_smd, analysis_level = "phosphosite",
    quant_type = "LFQ", source_software = "generic"
  )
  pr_mat <- matrix(rnorm(6), nrow = 1, dimnames = list("P_OTHER", snames))
  pr_fmd <- data.table::data.table(feature_id = "P_OTHER")
  pr_pfd <- ProteoForgeData(
    raw_matrix = pr_mat, feature_metadata = pr_fmd,
    sample_metadata = ph_smd, analysis_level = "protein",
    quant_type = "LFQ", source_software = "generic"
  )

  # No matching parent → correction_stats is empty but no error
  result <- parent_protein_correction(ph_pfd, pr_pfd)
  expect_equal(nrow(result$correction_stats), 0L)
})
