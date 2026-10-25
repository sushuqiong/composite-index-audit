# CompositeIndexAudit 0.34.0

## New entry point

* `audit_fit(fit, index)` audits an index against an **already-fitted**
  model object (`glm`, `lm`, `svyglm`, `coxph`, `survreg`). It reads the
  design matrix, the terms, the aliased columns, the family, the offset and
  the weights off the fit, so the user no longer has to transcribe the base
  model columns by hand. Offsets and aliased columns are reported as explicit
  failure states and are never downgraded.

## Clarifications

* The four classes and the non-classification states are now reported
  separately (`state_type`, `is_classification`), so "one of four classes, or
  an explicit non-classification state" is a statement the code enforces
  rather than a claim in the text.
* `outside_supported_class` is now distinguishable into
  `outside_grammar` (the parser cannot express the registered expression, for
  any base), `outside_no_registered_expression` (the registry stores a label
  rather than a formula) and `component_not_in_model` (the decomposition is
  known but this base lacks a column). The coarse label is unchanged.
* `glm7_full_expression()` writes out the seven routine measures of GLM7,
  which the registered string had hidden behind a sub-expression name.
* The registry now carries operational scale columns (`scale_of_C`,
  `scale_of_components_in_W`) plus `provenance` and `endorsement`.

## Fixes

* Benjamini-Hochberg values are no longer rounded to six decimals, which had
  printed every corrected value below 5e-7 as 0 and made the corrected value
  look smaller than the raw one.
* Base model terms are reported as model *terms*, not as design-matrix dummy
  column names, so a factor term is matched as `sex` and not as `sexMale`.
