qf <- system.file("extdata","generic_lfq_small","protein_matrix.tsv",
                  package="proteoforge")
mf <- system.file("extdata","generic_lfq_small","sample_metadata.tsv",
                  package="proteoforge")

.prep_stat <- function() {
  pfd <- import_generic(qf, mf)
  pfd <- log2_transform(pfd)
  pfd <- normalise_matrix(pfd, method="median")
  pfd <- impute_missing(pfd, method="minprob", seed=42)
  cfg <- list(statistics = list(
    design_formula  = "~ 0 + condition",
    contrasts = list(list(name="Treatment_vs_Control",
                          numerator="Treatment", denominator="Control")),
    p_adjust_method = "BH", alpha = 0.05, lfc_cutoff = 1
  ))
  list(pfd=pfd, cfg=cfg)
}

test_that("build_design_matrix produces correct dimensions", {
  smd <- data.table::data.table(
    sample_id = paste0("S",1:6),
    condition = rep(c("Control","Treatment"), each=3)
  )
  dm <- build_design_matrix(smd, "~ 0 + condition")
  expect_equal(nrow(dm), 6L)
  expect_equal(ncol(dm), 2L)
})

test_that("run_limma_deqms returns named list", {
  skip_if_not_installed("limma")
  x <- .prep_stat()
  res <- run_limma_deqms(x$pfd, x$cfg)
  expect_type(res, "list")
  expect_named(res, "Treatment_vs_Control")
})

test_that("run_limma_deqms result has required columns", {
  skip_if_not_installed("limma")
  x <- .prep_stat()
  res <- run_limma_deqms(x$pfd, x$cfg)
  dt  <- res[[1]]
  required <- c("feature_id","log2FC","P.Value","adj.P.Val",
                "significant","direction")
  expect_true(all(required %in% names(dt)))
})

test_that("run_limma_deqms detects known DE proteins", {
  skip_if_not_installed("limma")
  x      <- .prep_stat()
  res    <- run_limma_deqms(x$pfd, x$cfg)
  dt     <- res[[1]]
  truth  <- data.table::fread(system.file("extdata","generic_lfq_small",
                                           "ground_truth.tsv",
                                           package="proteoforge"))
  sig_ids  <- dt[significant==TRUE, feature_id]
  true_ids <- truth[true_de==TRUE,  ProteinID]
  # Expect at least 50% of true DE proteins to be recovered
  recovery <- length(intersect(sig_ids, true_ids)) / length(true_ids)
  expect_gt(recovery, 0.5)
})

test_that("build_gsea_ranking returns sorted named numeric vector", {
  dt <- data.table::data.table(
    feature_id  = paste0("P",1:10),
    gene_symbol = paste0("G",1:10),
    log2FC      = c(2,-2,1,-1,rep(0,6)),
    P.Value     = c(0.001,0.01,0.05,0.1,rep(0.5,6))
  )
  ranked <- build_gsea_ranking(dt)
  expect_type(ranked, "double")
  expect_true(!is.unsorted(rev(ranked)))
})
