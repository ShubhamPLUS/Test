qf <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
                  package="proteoforge")
mf <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
                  package="proteoforge")

test_that("import_generic returns ProteoForgeData", {
  pfd <- import_generic(qf, mf)
  expect_s4_class(pfd, "ProteoForgeData")
})

test_that("import_generic reads correct dimensions", {
  pfd <- import_generic(qf, mf)
  expect_equal(nrow(pfd@raw_matrix), 300L)
  expect_equal(ncol(pfd@raw_matrix), 12L)
})

test_that("import_generic sample IDs match metadata", {
  pfd <- import_generic(qf, mf)
  expect_equal(sort(colnames(pfd@raw_matrix)),
               sort(pfd@sample_metadata$sample_id))
})

test_that("import_generic converts zeros to NA when configured", {
  pfd <- import_generic(qf, mf, config = list(preprocessing = list(zero_to_na = TRUE)))
  # Should have no zeros in matrix (all become NA)
  expect_false(any(pfd@raw_matrix == 0, na.rm = TRUE))
})

test_that("import_generic feature_metadata has feature_id column", {
  pfd <- import_generic(qf, mf)
  expect_true("feature_id" %in% names(pfd@feature_metadata))
})

test_that("import_generic raw_long is populated", {
  pfd <- import_generic(qf, mf)
  expect_true(nrow(pfd@raw_long) > 0)
  expect_true(all(c("feature_id","sample_id","intensity") %in%
                    names(pfd@raw_long)))
})
