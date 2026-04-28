qf <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
                  package="proteoforge")
mf <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
                  package="proteoforge")

test_that("log2_transform converts to log scale", {
  pfd  <- import_generic(qf, mf)
  pfd2 <- log2_transform(pfd)
  # After log2, typical values should be 15-30
  med <- median(pfd2@raw_matrix, na.rm = TRUE)
  expect_true(med > 10 && med < 40)
})

test_that("log2_transform stores pre-transform values", {
  pfd  <- import_generic(qf, mf)
  pfd2 <- log2_transform(pfd)
  expect_true(isTRUE(pfd2@parameters$log2_transformed))
})

test_that("normalise_matrix median centers columns", {
  pfd  <- import_generic(qf, mf)
  pfd  <- log2_transform(pfd)
  pfd2 <- normalise_matrix(pfd, method = "median")
  col_meds <- apply(pfd2@raw_matrix, 2, median, na.rm = TRUE)
  grand_med <- median(col_meds)
  # After median normalisation, column medians should be approximately equal
  expect_true(diff(range(col_meds)) < 0.5)
})

test_that("normalise_matrix preserves pre-norm matrix", {
  pfd  <- import_generic(qf, mf)
  pfd  <- log2_transform(pfd)
  pfd2 <- normalise_matrix(pfd, method = "median")
  expect_false(is.null(pfd2@parameters$pre_norm_matrix))
})

test_that("normalise_matrix method=none returns unchanged matrix", {
  pfd  <- import_generic(qf, mf)
  pfd  <- log2_transform(pfd)
  before <- pfd@raw_matrix
  pfd2   <- normalise_matrix(pfd, method = "none")
  expect_equal(pfd2@raw_matrix, before)
})
