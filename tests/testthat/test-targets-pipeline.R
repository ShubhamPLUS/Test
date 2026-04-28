test_that("pf_targets_exists returns FALSE for non-existent store", {
  expect_false(pf_targets_exists(store = tempfile("_notexist")))
})

test_that("pf_targets_exists returns TRUE when store directory exists", {
  td <- tempfile("_targets_test")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE))
  file.create(file.path(td, "dummy.rds"))
  expect_true(pf_targets_exists(store = td))
})

test_that("pf_tar_read returns NULL when targets not installed or store missing", {
  # This should not error — just return NULL gracefully
  result <- suppressWarnings(pf_tar_read("nonexistent_target",
                                          store = tempfile("_fake")))
  expect_null(result)
})

test_that("pf_tar_status returns a data.table (possibly empty)", {
  dt <- suppressWarnings(pf_tar_status(store = tempfile("_fake")))
  expect_true(data.table::is.data.table(dt))
})

test_that("render_quarto_report returns NULL when quarto is not installed", {
  skip_if(requireNamespace("quarto", quietly = TRUE),
          "quarto installed; skipping absent-package test")
  dirs <- list(results_dir = tempdir(), reports_dir = tempdir())
  res <- suppressWarnings(
    render_quarto_report(list(), list(), dirs = dirs, config = list())
  )
  expect_null(res)
})

test_that("write_ai_summary_json writes valid JSON", {
  set.seed(11)
  n <- 10L
  stat_results <- list(
    CondA_vs_B = data.table::data.table(
      feature_id  = paste0("P", seq_len(n)),
      GeneSymbol  = paste0("GENE", seq_len(n)),
      log2FC      = rnorm(n),
      adj.P.Val   = runif(n),
      significant = runif(n) > 0.5,
      direction   = sample(c("up","down","ns"), n, replace = TRUE)
    )
  )
  qc_metrics <- list(
    n_samples   = 6L,
    n_features  = 100L,
    median_cv   = 0.15,
    outlier_samples = character(0),
    pca_variance    = c(0.45, 0.20)
  )
  tmp_dir <- tempdir()
  dirs    <- list(results_dir = tmp_dir)

  out <- write_ai_summary_json(
    stat_results = stat_results,
    qc_metrics   = qc_metrics,
    dirs         = dirs,
    config       = list(project_name = "TestProject")
  )

  expect_true(file.exists(out))
  parsed <- jsonlite::fromJSON(out)
  expect_equal(parsed$schema_version, "1.0")
  expect_equal(parsed$project_name, "TestProject")
  expect_equal(length(parsed$contrast_summaries), 1L)
  expect_equal(parsed$qc_summary$n_samples, 6L)

  unlink(out)
})
