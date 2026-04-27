test_that("hash_params returns 64-char hex string", {
  h <- hash_params(list(seed = 42, method = "median"))
  expect_type(h, "character")
  expect_equal(nchar(h), 64L)
  expect_match(h, "^[0-9a-f]+$")
})

test_that("hash_params is deterministic", {
  params <- list(seed = 1234, alpha = 0.05, engine = "limma_deqms")
  expect_equal(hash_params(params), hash_params(params))
})

test_that("hash_params differs for different inputs", {
  h1 <- hash_params(list(a = 1))
  h2 <- hash_params(list(a = 2))
  expect_false(h1 == h2)
})

test_that("write_param_hash writes a file", {
  tmp <- tempdir()
  p   <- write_param_hash("abc123def456", tmp)
  expect_true(fs::file_exists(p))
  expect_equal(readLines(p), "abc123def456")
})
