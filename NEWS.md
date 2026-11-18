# CompositeIndexAudit 0.39.0

## Specification and reporting decisions

* `verify_index_column()` now reports supplied-column certification separately from the declared symbolic class and fit state. A discordant column returns `input_certification = 'withheld'` and `verification_status = 'input_column_mismatch'`; it does not erase a valid symbolic decomposition or imply a failed fit.
* Added boundary tests for composite base expressions, no-intercept models, categorical encodings, pre-existing interactions, centering, unit mismatch and coefficient mismatch. Certification applies to explicitly declared terms rather than arbitrary implicit column-space identities.
* The ADEMP driver reports symbolic-class, input-certification, input-mismatch, structural-alias and failure-state rates separately. Its target coefficient and coverage denominator are explicit; no zero target is assigned to an aliased or discordant full-rank composite coefficient.
* Added the covariance-referenced CVAI decomposition, common-domain alias/QR checks and survey-weighted threshold/net-benefit sensitivity outputs. Unit conversion and intercept mapping are reported separately from discrimination.
* Submission metadata are aligned to version 0.39.0 and R 4.5.3. The earlier DOI 10.5281/zenodo.24217746 identifies only version 0.36.0. A new DOI must not be inferred before a public release is archived.

# CompositeIndexAudit 0.36.0

## Decision integrity and strict examples

* Replaced `haven` with the base/recommended `foreign` XPT reader, eliminating the observed R exit crash. All 24 official files passed an independent parser cross-check for dimensions, SEQN order, missingness and numerical agreement within 1e-12.
* Reran the complete NHANES pipeline with neutral `nhanes_` output names. Corrected the former AUC labelling error: `P_apparent_auc_*` now contains development-period apparent AUC and `L_auc_*` contains frozen-coefficient later-period transported AUC.
* Added the explicitly post hoc HbA1c >=5.7% outcome sensitivity driver and complete 24-index output. All contrasts have 1,000 valid paired Rao-Wu replicates; HbA1c is not an index component, and the outcome combines prediabetes- and diabetes-range values.
* Restored synthetic regression tests for survey reparameterization and unit-scaling helpers.

* Added 270 deterministic expectations covering 100 log-base linear-predictor equivalences, 50 round trips, 40 five-field identity comparisons, 30 scale-error detections, 30 correct-column checks and 20 symbolic classifications on degenerate samples.
* The runnable `vignette.R` now includes the manuscript's strict `DemoGLM7` public-API example with complete five-field keys for both the minimal and embedded bases.
* Documentation now distinguishes a symbolic `rank_adding_exact` direction from observed sample rank and reports the exact bundled-versus-unrestricted interpretation for fixed-ratio interaction bundles.
* Added deterministic recomputation of unweighted referral yield as `(368 - cases_missed) / n_flagged` and algebraic sensitivity for insulin conversion factor 6.945. The current release subsequently reran the complete NHANES analysis; the 6.945 calculation remains an algebraic sensitivity without a new threshold table.
* Validation uses R 4.5.2 with `--vanilla` to isolate user-profile and user-library state. The check command is `R CMD check --no-manual`; PDF-manual rendering is outside the available LaTeX environment.
* The base R runner completed 30 test blocks and 329 expectations with no failures. Source installation, `vignette.R`, both standalone audit drivers and all 12 simulation cells (16,000 repetitions) exited normally. `R CMD check --no-manual` returned `Status: OK`. Unit-test coverage was measured as 54.90%; the covr process then exited abnormally, so this is not reported as a successful coverage gate or as meeting an 85% target.

# CompositeIndexAudit 0.35.0

## Certification and classification fixes

* Symbolic classes no longer depend on observed QR rank; sample degeneration is reported separately as `verification_status` and `symbolic_numeric_rank_mismatch`.
* `strict = TRUE` is now the default for registered and expression-based public audits. Incomplete five-field metadata return `not_certified` without a coefficient map. Legacy bare-name exploration requires explicit `strict = FALSE`.
* `audit_fit(expression = ...)` now dispatches product expressions to the product audit and shares metadata matching with `audit_registered_index()`.
* Product audits subtract interaction columns already present in the base; a wholly embedded product is `rank_preserving_exact`, not a new interaction bundle.
* NLR and PLR are labelled as published standard forms; the constructed-illustration count is eight.
* The threshold-shift analysis now emits its domain size, case count, unweighted PPV and both unweighted and survey-weighted net benefit directly, preventing denominator or weighting ambiguity.
* The v35 simulation no longer assigns a zero target to a full-rank discordant composite; its bias and coverage are undefined rather than compared with the residual coefficient as if they shared an estimand.

## Validation status of this source revision

With R 4.5.2, the included base R runner executed all 24 test blocks and 59 expectations with no failures, and `R CMD check --no-manual` returned `Status: OK`. The original testthat-compatible files and optional standard entry point are retained. Source installation, Box 1, `vignette.R`, `test_audit.R` and `test_audit_fit.R` completed successfully. The complete simulation was separately rerun with R 4.6.2: all 16,000 replications completed and every assertion passed. NHANES outputs were not regenerated and remain archived v34 results. The maintainer email is `sushuqiong@sr.gxmu.edu.cn`; a standalone DESCRIPTION comment was removed for R 4.5 DCF compatibility, and the simulation aggregation now assigns the discordant composite target in the correct scope.


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
