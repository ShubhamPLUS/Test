qf_pg <- system.file("extdata","diann_example","report.pg_matrix.tsv",
                     package="proteoforge")
qf_rep <- system.file("extdata","diann_example","report.tsv",
                      package="proteoforge")
mf     <- system.file("extdata","diann_example","sample_metadata.tsv",
                      package="proteoforge")

test_that("import_diann pg_matrix returns ProteoForgeData", {
  pfd <- import_diann(qf_pg, mf)
  expect_s4_class(pfd, "ProteoForgeData")
  expect_equal(pfd@source_software, "diann")
})

test_that("import_diann pg_matrix has correct sample count", {
  pfd <- import_diann(qf_pg, mf)
  expect_equal(ncol(pfd@raw_matrix), 8L)
})

test_that("import_diann pg_matrix zeros become NA", {
  pfd <- import_diann(qf_pg, mf)
  expect_false(any(pfd@raw_matrix == 0, na.rm=TRUE))
})

test_that("import_diann report.tsv returns ProteoForgeData", {
  pfd <- import_diann(qf_rep, mf)
  expect_s4_class(pfd, "ProteoForgeData")
})
