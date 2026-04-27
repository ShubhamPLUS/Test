test_that("assert_columns passes when all columns present", {
  dt <- data.table::data.table(a = 1, b = 2, c = 3)
  expect_silent(assert_columns(dt, c("a", "b"), "test"))
})

test_that("assert_columns errors on missing columns", {
  dt <- data.table::data.table(a = 1)
  expect_error(assert_columns(dt, c("a", "b", "c"), "test"),
               "Missing required columns")
})

test_that("validate_sample_id_match passes on exact match", {
  expect_true(validate_sample_id_match(c("S1","S2"), c("S1","S2")))
})

test_that("validate_sample_id_match errors on mismatch", {
  expect_error(validate_sample_id_match(c("S1","S2","S3"), c("S1","S2")))
  expect_error(validate_sample_id_match(c("S1"), c("S1","S2")))
})

test_that("check_batch_confounding errors on fully confounded design", {
  smd <- data.table::data.table(
    condition = c("A","A","B","B"),
    batch     = c("b1","b1","b2","b2")
  )
  expect_error(check_batch_confounding(smd), "fully confounded")
})

test_that("check_batch_confounding passes on non-confounded design", {
  smd <- data.table::data.table(
    condition = c("A","A","B","B"),
    batch     = c("b1","b2","b1","b2")
  )
  expect_silent(check_batch_confounding(smd))
})

test_that("check_batch_confounding returns TRUE when no batch column", {
  smd <- data.table::data.table(condition = c("A","B"))
  expect_true(check_batch_confounding(smd))
})

test_that("validate_replicate_counts errors below min_reps", {
  smd <- data.table::data.table(condition = c("A","B","B"))
  expect_error(validate_replicate_counts(smd, min_reps = 2),
               "fewer than 2 replicates")
})

test_that("validate_replicate_counts passes with sufficient replicates", {
  smd <- data.table::data.table(condition = rep(c("A","B"), each = 3))
  counts <- validate_replicate_counts(smd)
  expect_equal(nrow(counts), 2L)
  expect_true(all(counts$N == 3L))
})
