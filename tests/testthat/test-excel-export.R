test_that("export_excel creates an xlsx file", {
  skip_if_not_installed("openxlsx2")
  qf  <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
                     package="proteoforge")
  mf  <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
                     package="proteoforge")
  pfd <- import_generic(qf, mf)
  pfd <- log2_transform(pfd)
  pfd <- normalise_matrix(pfd)
  pfd <- impute_missing(pfd, method="minprob", seed=42)

  out <- tempfile(fileext=".xlsx")
  export_excel(pfd, output_path=out)
  expect_true(file.exists(out))
  expect_gt(file.size(out), 1000L)
})
