# Contributing to ProteoForge

Thank you for your interest in contributing!

## Getting started

1. Fork the repository and create a feature branch.
2. Install dependencies: `Rscript scripts/install.R`
3. Make your changes on a branch named `feat/<description>` or `fix/<description>`.
4. Write tests for all new exported functions.
5. Run `R CMD check` locally — 0 errors, 0 warnings.
6. Run `lintr::lint_package()` and fix any issues.
7. Open a pull request with a clear description.

## Commit style

Use [Conventional Commits](https://www.conventionalcommits.org/):
- `feat:` new feature
- `fix:` bug fix
- `test:` tests only
- `docs:` documentation
- `chore:` CI, deps, build
- `refactor:` code reorganisation without behaviour change

## Code style

- Follow the [tidyverse style guide](https://style.tidyverse.org/).
- Run `styler::style_pkg()` before committing.
- Document all exported functions with roxygen2, including at least one `@examples` block.

## Testing

All exported functions require a `testthat` test in `tests/testthat/`.
Target ≥ 80% line coverage (`covr::package_coverage()`).

## Statistical decisions

Before adding or changing a default normalization, imputation, or batch
correction strategy, open an issue to discuss. Defaults materially affect
downstream biology and must be documented with a rationale.
