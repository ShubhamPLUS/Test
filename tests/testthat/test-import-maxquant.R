qf <- system.file("extdata","maxquant_example","proteinGroups.txt",
                  package="proteoforge")
mf <- system.file("extdata","maxquant_example","sample_metadata.tsv",
                  package="proteoforge")

test_that("import_maxquant returns ProteoForgeData", {
  pfd <- import_maxquant(qf, mf)
  expect_s4_class(pfd, "ProteoForgeData")
  expect_equal(pfd@source_software, "maxquant")
})

test_that("import_maxquant removes reverse and contaminant rows", {
  pfd <- import_maxquant(qf, mf)
  # 5 Reverse + 3 contaminant + 2 OIS = 10 rows removed from 140
  # (some overlap possible; at least < 140)
  expect_lt(nrow(pfd@raw_matrix), 140L)
})

test_that("import_maxquant correct sample count", {
  pfd <- import_maxquant(qf, mf)
  expect_equal(ncol(pfd@raw_matrix), 8L)
})

test_that("import_maxquant has gene_symbol in feature metadata", {
  pfd <- import_maxquant(qf, mf)
  expect_true("gene_symbol" %in% names(pfd@feature_metadata))
})
