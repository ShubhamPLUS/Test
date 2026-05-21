qf <- system.file("extdata","fragpipe_example","combined_protein.tsv",
                  package="proteoforge")
mf <- system.file("extdata","fragpipe_example","sample_metadata.tsv",
                  package="proteoforge")

test_that("import_fragpipe returns ProteoForgeData", {
  pfd <- import_fragpipe(qf, mf)
  expect_s4_class(pfd, "ProteoForgeData")
  expect_equal(pfd@source_software, "fragpipe")
})

test_that("import_fragpipe correct dimensions", {
  pfd <- import_fragpipe(qf, mf)
  expect_equal(ncol(pfd@raw_matrix), 8L)
  expect_equal(nrow(pfd@raw_matrix), 130L)
})
