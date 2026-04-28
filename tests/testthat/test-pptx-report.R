test_that("export_pptx creates a pptx file", {
  skip_if_not_installed("officer")
  qf  <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
                     package="proteoforge")
  mf  <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
                     package="proteoforge")
  pfd <- import_generic(qf, mf)
  pfd <- log2_transform(pfd)
  out <- tempfile(fileext=".pptx")
  export_pptx(pfd, output_path=out)
  expect_true(file.exists(out))
  expect_gt(file.size(out), 1000L)
})
