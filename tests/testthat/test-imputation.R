qf <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
                  package="proteoforge")
mf <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
                  package="proteoforge")

.prep_pfd <- function() {
  pfd <- import_generic(qf, mf)
  pfd <- log2_transform(pfd)
  normalise_matrix(pfd, method = "median")
}

test_that("impute_missing method=none leaves NAs unchanged", {
  pfd  <- .prep_pfd()
  n_na_before <- sum(is.na(pfd@raw_matrix))
  pfd2 <- impute_missing(pfd, method = "none")
  expect_equal(sum(is.na(pfd2@raw_matrix)), n_na_before)
})

test_that("impute_missing method=minprob fills all NAs", {
  pfd  <- .prep_pfd()
  pfd2 <- impute_missing(pfd, method = "minprob", seed = 42)
  expect_equal(sum(is.na(pfd2@raw_matrix)), 0L)
})

test_that("impute_missing method=knn fills all NAs", {
  pfd  <- .prep_pfd()
  pfd2 <- impute_missing(pfd, method = "knn")
  expect_equal(sum(is.na(pfd2@raw_matrix)), 0L)
})

test_that("impute_missing method=mixed fills all NAs", {
  pfd  <- .prep_pfd()
  pfd2 <- impute_missing(pfd, method = "mixed", seed = 42)
  expect_equal(sum(is.na(pfd2@raw_matrix)), 0L)
})

test_that("impute_missing preserves pre-impute matrix", {
  pfd  <- .prep_pfd()
  pfd2 <- impute_missing(pfd, method = "minprob", seed = 42)
  expect_false(is.null(pfd2@parameters$pre_impute_matrix))
  expect_true(sum(is.na(pfd2@parameters$pre_impute_matrix)) > 0)
})

test_that("impute_missing records imputation statistics", {
  pfd  <- .prep_pfd()
  pfd2 <- impute_missing(pfd, method = "mixed", seed = 42)
  stats <- pfd2@parameters$imputation_stats
  expect_true(!is.null(stats))
  expect_true(stats$n_imputed > 0)
})

test_that("impute_missing is reproducible with same seed", {
  pfd   <- .prep_pfd()
  pfd2a <- impute_missing(pfd, method = "minprob", seed = 999)
  pfd2b <- impute_missing(pfd, method = "minprob", seed = 999)
  expect_equal(pfd2a@raw_matrix, pfd2b@raw_matrix)
})
