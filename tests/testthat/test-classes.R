test_that("ProteoForgeData constructor returns valid object", {
  mat <- matrix(c(1,2,3,4), nrow=2,
                dimnames=list(c("P1","P2"), c("S1","S2")))
  smd <- data.table::data.table(sample_id = c("S1","S2"),
                                condition = c("A","B"))
  pfd <- ProteoForgeData(raw_matrix = mat, sample_metadata = smd)
  expect_s4_class(pfd, "ProteoForgeData")
  expect_equal(pfd@analysis_level, "protein")
  expect_equal(pfd@source_software, "generic")
})

test_that("ProteoForgeData show method runs without error", {
  pfd <- ProteoForgeData()
  expect_output(show(pfd), "ProteoForgeData")
})

test_that("ProteoForgeData validity rejects bad analysis_level", {
  expect_error(
    methods::new("ProteoForgeData", analysis_level = "bad_level"),
    "analysis_level"
  )
})

test_that("ProteoForgeData validity rejects bad source_software", {
  expect_error(
    methods::new("ProteoForgeData", source_software = "unknown_tool"),
    "source_software"
  )
})

test_that("ProteoForgeData auto-computes parameter_hash", {
  pfd <- ProteoForgeData(parameters = list(seed = 42, alpha = 0.05))
  expect_true(nchar(pfd@parameter_hash) == 64L)
})

test_that("ProteoForgeData dimensions are correct", {
  mat <- matrix(runif(60), nrow=10, ncol=6,
                dimnames=list(paste0("P",1:10), paste0("S",1:6)))
  smd <- data.table::data.table(
    sample_id = paste0("S",1:6),
    condition = rep(c("A","B"), each=3)
  )
  pfd <- ProteoForgeData(raw_matrix = mat, sample_metadata = smd)
  expect_equal(nrow(pfd@raw_matrix), 10L)
  expect_equal(ncol(pfd@raw_matrix), 6L)
})
