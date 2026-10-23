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

## Quick start (15 lines, no NHANES download)

```r
source("algebra_audit_v32.R"); source("public_api_v32.R")

# GLM7 = metab + log10(age) + log10(BMI); base = sex + log10(age) + log10(BMI)
a <- audit_index_in_model(c(metab = 1, log_age = 1, log_bmi = 1),
                          c(log_age = 1, log_bmi = 1), data = d)
a$classification                     # "rank_adding_exact"

# same index, base now containing metab: the classification follows the model
b <- audit_index_in_model(c(metab = 1, log_age = 1, log_bmi = 1),
                          c(metab = 1, log_age = 1, log_bmi = 1), data = d)
b$classification                     # "rank_preserving_exact"
fit <- glm(dx ~ metab + log_age + log_bmi + GLM7, data = d, family = binomial())
coef(fit)["GLM7"]                    # NA - aliased, where the audit predicted
```

`vignette_v32.R` is the full, self-verifying example; `test_audit_v32.R`
runs 43 assertions.

## Files

- `algebra_audit_v32.R` - symbolic and numerical layers
  (`audit_index_in_model`, `audit_product_in_model`)
- `index_registry_v32.R` - the frozen 25-expression registry and product
  expansions
- `public_api_v32.R` - `parse_index_expression`, `register_index`,
  `audit_registered_index`, `verify_index_column`, `normalize_log_base`
- `model_evaluation_v32.R`, `data_build_v32.R`, `build_results_v32.R` -
  the design-aware NHANES pipeline
- `reanalysis_cases_v32.R` - reproduction of two published model
  specifications on public NHANES data
- `vignette_v32.R`, `test_audit_v32.R` - Box 1 and the unit tests

## Requirements

R >= 4.5; CRAN packages `survey`, `haven`, `digest`. No compiled code.

## License

MIT
