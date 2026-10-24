# composite-index-audit

Pre-fitting algebraic classification of composite indices entered alongside
their own components in a regression model.

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
enumerated failure state that is never silently downgraded.

## Term identity

Terms are matched on a five-field identity - variable, transformation, log
base, units, encoding - so `log x` and `1/x`, or `ln x` and `log10 x`, are
different terms even when the column names coincide. A term whose units are
not declared never matches a declared one, so unknown provenance fails
closed.

## Quick start (no NHANES download)

```r
source("algebra_audit_v33.R"); source("term_identity_v33.R")
source("index_registry_v33.R"); source("public_api_v33.R")

# GLM7 = metab + log10(age) + log10(BMI); base = sex + log10(age) + log10(BMI)
a <- audit_index_in_model(c(metab = 1, log_age = 1, log_bmi = 1),
                          c(log_age = 1, log_bmi = 1))
a$classification                     # "rank_adding_exact" (no data needed)

# same index, base now containing metab: the class follows the model
b <- audit_index_in_model(c(metab = 1, log_age = 1, log_bmi = 1),
                          c(metab = 1, log_age = 1, log_bmi = 1))
b$classification                     # "rank_preserving_exact"

# a native-scale product: LAP = (waist - 65|58) x TG
e <- audit_product_in_model("LAP", c("sex", "waist_cm", "tg_mmol_L"),
                            product_expansions()$LAP, data = d)
e$classification                     # "reducible_to_interaction"
```

`vignette_v33.R` is the full self-verifying example; `test_audit_v33.R`
runs the unit tests, including the v33 regression cases (wrong log base,
missing intercept, extra base column, zero-variance column, pre-existing
interaction, deliberately failing fit, and the no-data C = x + y case that
the v32 `length(residual) == 0` test mis-classified).

## Files

- `algebra_audit_v33.R` - symbolic and numerical layers
  (`audit_index_in_model`, `audit_product_in_model`)
- `index_registry_v33.R` - the frozen 25-expression registry and product
  expansions
- `term_identity_v33.R` - the five-field term identity and its three-state
  comparison (match / unkeyed / mismatch)
- `public_api_v33.R` - `parse_index_expression`, `register_index`,
  `audit_registered_index`, `verify_index_column`, `normalize_log_base`
- `registry_audit_v33.R` - independent re-derivation of all 26 registry
  classifications (this is what makes the class counts a measurement)
- `model_evaluation_v33.R`, `data_build_v33.R`, `build_results_v33.R` -
  the design-aware NHANES pipeline
- `reanalysis_cases_v33.R` - reproduction of published model specifications
  on public NHANES data
- `simulation_v33.R` - the 12-cell ADEMP study
- `benchmarks_v33.R` - runtime and the measured comparison with
  `alias()`, QR rank, `kappa()`, `car::vif()` and `deriv()`
- `vignette_v33.R`, `test_audit_v33.R` - Box 1 and the unit tests

## Reproducibility

`data_build_v33.R` takes its data root from the `CIA_DATA_ROOT` environment
variable or an explicit argument; no absolute path is hard-coded. The
scripts write to `CIA_OUT` (default `../results`).

## Requirements

R >= 4.5; CRAN packages `survey`, `haven`, `digest`. No compiled code.

## License

MIT
