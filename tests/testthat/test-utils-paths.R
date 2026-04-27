test_that("sanitise_name replaces bad characters", {
  expect_equal(sanitise_name("hello world!"), "hello_world")
  expect_equal(sanitise_name("A vs B"), "A_vs_B")
  expect_equal(sanitise_name("OK-name_123"), "OK-name_123")
  expect_equal(sanitise_name("__leading__"), "leading")
})

test_that("safe_filename builds correct format", {
  ts <- as.POSIXct("2026-04-27 14:30:00", tz = "UTC")
  fn <- safe_filename("MyStudy", "STAT", "protein", "Tx_vs_Ctrl", "volcano", "pdf",
                      timestamp = ts)
  expect_match(fn, "^MyStudy__STAT__protein__Tx_vs_Ctrl__volcano__20260427_143000\\.pdf$")
})

test_that("safe_filename sanitises project name", {
  ts <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  fn <- safe_filename("My Study!", "QC", "protein", "all", "pca", "png", timestamp = ts)
  expect_false(grepl(" ", fn))
  expect_false(grepl("!", fn))
})

test_that("make_output_dirs creates all subdirectories", {
  tmp <- tempfile()
  dirs <- make_output_dirs(tmp, "TestProject",
                           timestamp = as.POSIXct("2026-01-01", tz="UTC"))
  expect_true(fs::dir_exists(dirs["logs"]))
  expect_true(fs::dir_exists(dirs["stat_tables"]))
  expect_true(fs::dir_exists(dirs["serialized"]))
  fs::dir_delete(tmp)
})
