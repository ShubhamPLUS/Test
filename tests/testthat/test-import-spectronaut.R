qf <- system.file("extdata","spectronaut_example","spectronaut_report.tsv",
                  package="proteoforge")
mf <- system.file("extdata","spectronaut_example","sample_metadata.tsv",
                  package="proteoforge")

test_that("import_spectronaut returns ProteoForgeData", {
  pfd <- import_spectronaut(qf, mf)
  expect_s4_class(pfd, "ProteoForgeData")
  expect_equal(pfd@source_software, "spectronaut")
})

test_that("import_spectronaut has correct sample count", {
  pfd <- import_spectronaut(qf, mf)
  expect_equal(ncol(pfd@raw_matrix), 8L)
})
