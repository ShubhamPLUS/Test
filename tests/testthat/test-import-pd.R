qf <- system.file("extdata","pd_example","proteins.tsv",
                  package="proteoforge")
mf <- system.file("extdata","pd_example","sample_metadata.tsv",
                  package="proteoforge")

test_that("import_pd returns ProteoForgeData", {
  pfd <- import_pd(qf, mf)
  expect_s4_class(pfd, "ProteoForgeData")
  expect_equal(pfd@source_software, "pd")
})

test_that("import_pd filters contaminants", {
  pfd <- import_pd(qf, mf)
  expect_lt(nrow(pfd@raw_matrix), 110L)
})

test_that("import_pd correct sample count", {
  pfd <- import_pd(qf, mf)
  expect_equal(ncol(pfd@raw_matrix), 8L)
})
