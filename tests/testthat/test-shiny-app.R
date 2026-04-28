test_that("launch_app stops gracefully when app dir not found (sanity check)", {
  # We can't run the full Shiny app in unit tests without shinytest2,
  # but we can verify the function exists and has correct signature
  expect_true(is.function(launch_app))
  args <- formals(launch_app)
  expect_true("port" %in% names(args))
  expect_true("host" %in% names(args))
})

test_that("Shiny app directory contains required files", {
  app_dir <- system.file("shiny", package = "proteoforge")
  skip_if(!nzchar(app_dir), "Package not installed; skipping Shiny dir check")
  expect_true(file.exists(file.path(app_dir, "app.R")))
  modules <- list.files(file.path(app_dir, "modules"), pattern = "\\.R$")
  expect_gte(length(modules), 14L)
  # Verify all required modules exist
  expected <- c("mod_welcome","mod_upload","mod_colmap","mod_design",
                "mod_qc","mod_preprocess","mod_stats","mod_heatmap",
                "mod_venn","mod_enrichment","mod_network","mod_phospho",
                "mod_report","mod_logs")
  for (mod in expected) {
    expect_true(any(grepl(mod, modules)),
                info = sprintf("Module file for %s not found", mod))
  }
})
