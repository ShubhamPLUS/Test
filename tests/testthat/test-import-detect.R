test_that("detect_software identifies DIA-NN from filename", {
  f <- system.file("extdata","diann_example","report.pg_matrix.tsv",
                   package="proteoforge")
  expect_equal(detect_software(f), "diann")
})

test_that("detect_software identifies FragPipe from filename", {
  f <- system.file("extdata","fragpipe_example","combined_protein.tsv",
                   package="proteoforge")
  expect_equal(detect_software(f), "fragpipe")
})

test_that("detect_software identifies MaxQuant from filename", {
  f <- system.file("extdata","maxquant_example","proteinGroups.txt",
                   package="proteoforge")
  expect_equal(detect_software(f), "maxquant")
})

test_that("detect_software uses hint when provided", {
  f <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
                   package="proteoforge")
  expect_equal(detect_software(f, hint="spectronaut"), "spectronaut")
})

test_that("detect_software falls back to generic for plain matrix", {
  f <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
                   package="proteoforge")
  expect_equal(detect_software(f, hint="auto"), "generic")
})
