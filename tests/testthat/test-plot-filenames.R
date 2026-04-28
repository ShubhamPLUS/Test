test_that("volcano plot is a ggplot object", {
  dt <- data.table::data.table(
    feature_id  = paste0("P",1:50),
    gene_symbol = paste0("G",1:50),
    log2FC      = rnorm(50,0,2),
    adj.P.Val   = runif(50,0,1),
    P.Value     = runif(50,0,1),
    significant = FALSE,
    direction   = "ns"
  )
  p <- plot_volcano(dt, title="Test")
  expect_s3_class(p, "gg")
})

test_that("theme_proteoforge returns a theme object", {
  th <- theme_proteoforge()
  expect_s3_class(th, "theme")
})

test_that("pf_palette returns correct number of colours", {
  expect_length(pf_palette(4), 4)
  expect_length(pf_palette(8), 8)
})

test_that("get_intersection_table returns correct columns", {
  dt1 <- data.table::data.table(feature_id=c("P1","P2","P3"),
    gene_symbol=c("G1","G2","G3"), significant=TRUE, direction="up")
  dt2 <- data.table::data.table(feature_id=c("P2","P3","P4"),
    gene_symbol=c("G2","G3","G4"), significant=TRUE, direction="up")
  tbl <- get_intersection_table(list(A=dt1,B=dt2))
  expect_true(all(c("intersection_id","set_membership","count") %in% names(tbl)))
  expect_true(any(tbl$count > 1))
})

test_that("save_plot creates files", {
  skip_if_not_installed("ggplot2")
  p <- ggplot2::ggplot(data.frame(x=1:5,y=1:5),
                        ggplot2::aes(x,y)) + ggplot2::geom_point()
  tmp <- tempdir()
  saved <- save_plot(p, tmp, "test_plot", formats="png", width=4, height=3)
  expect_true(any(fs::file_exists(saved)))
})
