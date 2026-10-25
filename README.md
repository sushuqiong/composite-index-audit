# composite-index-audit

R package **CompositeIndexAudit** (version 0.34.0, tag v34): pre-fitting
algebraic classification of composite indices entered alongside their own
components in a regression model.

Many composite indices (TyG, HOMA-IR, LAP, CVAI, NLR, SII, GLM7, ...) are
built from the same variables that enter the fitted model. Whether the index
then carries information beyond those components is a question about the
index *and* the model, and it can be answered from the formulas before any
data are used. This package answers it.

Given a registered index expression and a model specification, the audit
returns one of four classes:

| class | meaning |
|---|---|
| `rank_preserving_exact` | the index lies in the span of the base model; the added coefficient is not identifiable and the increment is zero |
| `rank_adding_exact` | the decomposition is exact but the index carries a direction outside the base model; its coefficient is the residual's, not an independent contribution |
| `reducible_to_interaction` | a native-scale product that expands into interactions of columns already in the model |
| `outside_supported_class` | the symbolic layer cannot decide, and says which terms blocked it |

It also returns the closed-form coefficient map `gamma = K beta` where one
exists, the named interaction bundle, the intercept remapping
`beta_0* = beta_0 - k beta_C` for a unit conversion `C* = C + k`, and an
enumerated failure state (`approximate_collinearity`, `fit_failed`,
`not_decomposable`) that is never silently downgraded.

## Term identity

Terms are matched on a five-field identity - variable, transformation, log
base, units, encoding - so `log x` and `1/x`, or `ln x` and `log10 x`, are
different terms even when the column names coincide. A term whose units are
not declared never matches a declared one, so unknown provenance fails
closed.

## Install and test

```r
# from a source checkout
# install.packages(c("survey", "haven", "digest"))  # the only dependencies
# R CMD INSTALL CompositeIndexAudit_0.34.0.tar.gz   # or install the source dir
library(CompositeIndexAudit)

# GLM7 = metab + log10(age) + log10(BMI); base = sex + log10(age) + log10(BMI)
a <- audit_index_in_model(c(metab = 1, log_age = 1, log_bmi = 1),
                          c(log_age = 1, log_bmi = 1))
a$classification                     # "rank_adding_exact" (no data needed)

# same index, base now containing metab: the class follows the model
b <- audit_index_in_model(c(metab = 1, log_age = 1, log_bmi = 1),
                          c(metab = 1, log_age = 1, log_bmi = 1))
b$classification                     # "rank_preserving_exact"
```

`testthat::test_dir("tests/testthat")` runs the distributed suite (46
expectations across four files); `R CMD check` returns 0 errors, 0 warnings
and 0 notes under R 4.5.1.

## Layout

- `R/` - the package sources: the symbolic and numerical layers
  (`algebra_audit.R`), the frozen 25-expression registry
  (`index_registry.R`), the five-field term identity (`term_identity.R`),
  the public API (`public_api.R`), `audit_fit()` for fitted-model objects
  (`audit_fit.R`), and the independent registry re-derivation
  (`registry_audit.R`)
- `man/` - roxygen-generated help pages; `inst/CITATION`; `NEWS.md`
- `tests/testthat/` - the test suite
- `analysis/` - the NHANES demonstration pipeline: `build_results.R`,
  `simulation.R` (the 12-cell ADEMP study), `reanalysis_cases.R`
  (reproduction of published model specifications), `benchmarks.R`
  (runtime and the measured comparison with `alias()`, QR rank, `kappa()`,
  `car::vif()` and `deriv()`), `flow_counts.R`, `verify_unit_factor.R`,
  `reclassify.R`, `dump_registry.R`, the standalone driver `test_audit.R`
  and the runnable example `vignette.R` (Box 1 of the manuscript). The
  analysis scripts are self-contained script versions that `source()` each
  other; they do not require the package to be installed.

## Reproducibility

The analysis scripts take their data root from the `CIA_DATA_ROOT`
environment variable or an explicit argument; no absolute path is
hard-coded. They write to `CIA_OUT` (default `../results`).

## Requirements

R >= 4.5; CRAN packages `survey`, `haven`, `digest`. No compiled code.

## License

MIT
