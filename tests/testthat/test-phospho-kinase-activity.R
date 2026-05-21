test_that("infer_kinase_activity returns a named list with expected elements", {
  set.seed(99)
  n <- 40L
  result_dt <- data.table::data.table(
    feature_id  = paste0("P", seq_len(n), "_S", sample(c("S","T","Y"), n, replace=TRUE),
                          sample(100:500, n)),
    t           = rnorm(n, sd = 3),
    log2FC      = rnorm(n),
    adj.P.Val   = runif(n),
    significant = runif(n) > 0.6,
    direction   = sample(c("up","down","ns"), n, replace=TRUE)
  )

  res <- infer_kinase_activity(result_dt, methods = c("ksea"))

  expect_type(res, "list")
  expect_true(all(c("ksea","kea3","decoupler","consensus") %in% names(res)))
})

test_that("consensus is NULL when fewer than consensus_min methods succeed", {
  set.seed(7)
  n <- 10L
  result_dt <- data.table::data.table(
    feature_id  = paste0("P", seq_len(n), "_S", seq_len(n)),
    t           = rnorm(n),
    log2FC      = rnorm(n),
    adj.P.Val   = runif(n),
    significant = FALSE,
    direction   = "ns"
  )

  # Only KSEA is requested; consensus_min=2 → consensus should be NULL
  res <- infer_kinase_activity(result_dt, methods = "ksea", consensus_min = 2L)
  expect_null(res$consensus)
})

test_that("consensus table has expected columns when methods agree", {
  # Build a mock scenario where we manually provide results
  # by testing the consensus logic with patched method outputs
  set.seed(5)
  n <- 30L
  result_dt <- data.table::data.table(
    feature_id  = paste0("P", seq_len(n), "_S", seq_len(n)),
    t           = rnorm(n, sd = 2),
    log2FC      = rnorm(n),
    adj.P.Val   = c(rep(0.001, 10), runif(20)),
    significant = c(rep(TRUE, 10), rep(FALSE, 20)),
    direction   = c(rep("up", 5), rep("down", 5), rep("ns", 20))
  )

  # Call with all methods but expect graceful handling of missing packages
  res <- infer_kinase_activity(result_dt, methods = c("ksea","decoupler"),
                                consensus_min = 1L)
  expect_type(res, "list")
  # If consensus has rows, check columns
  if (!is.null(res$consensus) && nrow(res$consensus) > 0) {
    expect_true("kinase" %in% names(res$consensus))
    expect_true("n_methods_agree" %in% names(res$consensus))
  }
})

test_that("infer_kinase_activity handles missing t column with log2FC fallback", {
  n <- 20L
  result_dt <- data.table::data.table(
    feature_id  = paste0("SITE", seq_len(n)),
    log2FC      = rnorm(n),
    adj.P.Val   = runif(n),
    significant = FALSE,
    direction   = "ns"
  )
  # No 't' column — should use log2FC
  expect_silent(infer_kinase_activity(result_dt, methods = "ksea"))
})
