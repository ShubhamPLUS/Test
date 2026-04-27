#!/usr/bin/env Rscript
## Regenerate all synthetic test datasets.
## Usage: Rscript scripts/make_test_data.R

message("Regenerating synthetic test datasets...")
source(file.path("data-raw", "generate_synthetic_datasets.R"))
message("Done.")
