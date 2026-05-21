qf <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
                  package="proteoforge")
mf <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
                  package="proteoforge")

test_that("validate_sample_metadata passes clean data", {
  pfd <- import_generic(qf, mf)
  pfd2 <- validate_sample_metadata(pfd)
  expect_true(pfd2@validation$passed)
})

test_that("validate_sample_metadata populates validation slot", {
  pfd  <- import_generic(qf, mf)
  pfd2 <- validate_sample_metadata(pfd)
  expect_true(is.list(pfd2@validation))
  expect_true("checks" %in% names(pfd2@validation))
})

test_that("validate_sample_metadata errors on wrong sample IDs", {
  pfd <- import_generic(qf, mf)
  pfd@sample_metadata[1, sample_id := "WRONG_ID"]
  expect_error(validate_sample_metadata(pfd))
})

test_that("validate_sample_metadata errors on single condition", {
  pfd <- import_generic(qf, mf)
  pfd@sample_metadata[, condition := "OnlyOne"]
  expect_error(validate_sample_metadata(pfd))
})

test_that("validate_sample_metadata writes JSON when output_dir given", {
  pfd  <- import_generic(qf, mf)
  tmp  <- tempdir()
  pfd2 <- validate_sample_metadata(pfd, output_dir = tmp)
  expect_true(file.exists(file.path(tmp, "validation_report.json")))
})
