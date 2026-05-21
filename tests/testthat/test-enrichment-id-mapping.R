test_that("map_ids_to_entrez returns named character vector", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")
  ids <- c("TP53","BRCA1","EGFR","MYC")
  res <- map_ids_to_entrez(ids, organism="human", id_type="symbol")
  expect_type(res, "character")
  expect_named(res)
  expect_true(any(!is.na(res)))
})

test_that("detect_id_type recognises UniProt", {
  expect_equal(proteoforge:::.detect_id_type(c("P38398","O15350","Q9Y6K9")),
               "uniprot")
})

test_that("detect_id_type recognises Ensembl", {
  expect_equal(proteoforge:::.detect_id_type(
    c("ENSG00000141510","ENSG00000012048")), "ensembl")
})

test_that("detect_id_type defaults to symbol", {
  expect_equal(proteoforge:::.detect_id_type(c("TP53","BRCA1","EGFR")),
               "symbol")
})

test_that("build_gsea_ranking produces descending sorted vector", {
  dt <- data.table::data.table(
    feature_id  = paste0("P",1:20),
    gene_symbol = paste0("G",1:20),
    log2FC      = rnorm(20),
    P.Value     = runif(20,0.001,0.5)
  )
  ranked <- build_gsea_ranking(dt)
  expect_false(is.unsorted(rev(ranked)))
  expect_named(ranked)
})
