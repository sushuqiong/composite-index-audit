# composite-index-audit

R package **CompositeIndexAudit** (version 0.39.0): pre-fitting
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
| `rank_adding_exact` | under the declared symbolic terms, the index carries a direction outside the base model; its coefficient is the residual's, not an independent contribution; observed sample rank is verified separately |
| `reducible_to_interaction` | a native-scale product that expands into interactions of columns already in the model |
| `outside_supported_class` | the symbolic layer cannot decide, and says which terms blocked it |

It also returns the closed-form coefficient map `gamma = K beta` where one
exists, the named interaction bundle, the intercept remapping
`beta_0* = beta_0 - k beta_C` for a unit conversion `C* = C + k`, and an
enumerated non-classification state (`approximate_collinearity`, `fit_failed`,
`not_decomposable`, `not_certified`) that is never silently downgraded.

## Term identity

Terms are matched on a five-field identity - variable, transformation, log
base, units, encoding - so `log x` and `1/x`, or `ln x` and `log10 x`, are
different terms even when the column names coincide. A term whose units are
not declared never matches a declared one, so unknown provenance fails
closed.

## Install and use

```r
# The runtime package imports only stats; tests use the included base R runner.
# R CMD INSTALL /path/to/composite-index-audit
library(CompositeIndexAudit)

# GLM7 = metab + log10(age) + log10(BMI); base = log10(age) + log10(BMI)
a <- audit_index_in_model(c(metab = 1, log_age = 1, log_bmi = 1),
                          c(log_age = 1, log_bmi = 1))
a$classification                     # "rank_adding_exact" (no data needed)

# The same index against a base containing metab
b <- audit_index_in_model(c(metab = 1, log_age = 1, log_bmi = 1),
                          c(metab = 1, log_age = 1, log_bmi = 1))
b$classification                     # "rank_preserving_exact"
```

`audit_registered_index()` and `audit_index_in_model()` are the pre-fitting
entry points. `audit_fit()` is a convenience interface for an already-fitted
object. Exact certification through the public registered/expression interfaces
requires a named vector of five-field term keys and declared units. Incomplete
metadata return `not_certified`; `strict = FALSE` is reserved for explicitly
exploratory legacy comparisons.

Version 0.39.0 separates the declared symbolic class, certification of the supplied
column, observed sample rank and numerical fit state. `verify_index_column()`
returns `input_certification`, `verification_status` and `failure_states`; a
mismatched supplied column therefore need not change the declared class or prevent
model fitting. Boundary tests cover composite base expressions, no-intercept
models, categorical encoding, existing interactions, centering, unit mismatch and
coefficient mismatch. The analysis drivers also provide a common-domain alias/QR
comparison, a covariance-referenced CVAI decomposition and survey-weighted
threshold/net-benefit sensitivity analyses.

The earlier DOI 10.5281/zenodo.24217746 applies only to archived version 0.36.0.
Do not assign it to the 0.39.0 source or results. Validation claims for this
submission are limited to the logs produced from the final 0.39.0 source archive.

## Layout

- `R/` - the package sources: the symbolic and numerical layers
  (`algebra_audit.R`), the frozen 25-expression registry
  (`index_registry.R`), the five-field term identity (`term_identity.R`),
  the public API (`public_api.R`), `audit_fit()` for fitted-model objects
  (`audit_fit.R`), and the registry concordance checks based on hand-written
  algebraic expansions (`registry_audit.R`)
- `man/` - roxygen-generated help pages; `inst/CITATION`; `NEWS.md`
- `tests/testthat/` - the test suite
- `analysis/` - the NHANES demonstration pipeline: `build_results.R`,
  `simulation.R` (the 12-cell ADEMP study), `reanalysis_cases.R`
  (reproduction of published model specifications), `benchmarks.R`
  (runtime and the measured comparison with `alias()`, QR rank, `kappa()`,
  `car::vif()` and `deriv()`), `flow_counts.R`, `verify_unit_factor.R`,
  `reclassify.R`, `dump_registry.R`, `hba1c_outcome.R`, the standalone drivers `test_audit.R` and `test_analysis_helpers.R`
  and the runnable example `vignette.R` (Box 1 of the manuscript). The
  analysis scripts are self-contained script versions that `source()` each
  other; they do not require the package to be installed.

## Reproducibility

The analysis scripts take their data root from the `CIA_DATA_ROOT`
environment variable or an explicit argument; no absolute path is
hard-coded. They write to `CIA_OUT` (default `../results`).

## Requirements

The runtime package imports base R `stats`. The formal check executes the
testthat-compatible files with the included base R runner; the optional standard
entry point uses `testthat`. Model-class tests additionally use `survey` and
`survival`. The NHANES analysis scripts use the base/recommended `foreign`
XPT reader, `survey` and `digest`; they do not load `haven`.
No compiled code is included.

## License

MIT
