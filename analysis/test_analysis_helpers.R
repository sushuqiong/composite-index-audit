# Regression test for analysis-only helpers used by build_results.R.
# These tests use synthetic data and require only survey plus recommended/base R.
CODE_DIR <- Sys.getenv("CIA_CODE_DIR", unset = "")
if (!nzchar(CODE_DIR)) stop("set CIA_CODE_DIR before running test_analysis_helpers.R")
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")
source(file.path(CODE_DIR, "model_evaluation.R"), local = TRUE)
source(file.path(CODE_DIR, "algebra_audit.R"), local = TRUE)
stopifnot(is.function(check_survey_reparameterization),
          is.function(run_unit_scaling_sensitivity))

make_side <- function(seed, n = 240L) {
  set.seed(seed)
  stra <- factor(rep(sprintf("S%02d", 1:12), each = n / 12))
  psu <- paste(stra, rep(rep(1:2, each = 10), 12), sep = "_")
  sex <- factor(sample(c("Male", "Female"), n, replace = TRUE), levels = c("Male", "Female"))
  log_age <- log10(runif(n, 20, 80))
  log_bmi <- log10(runif(n, 18, 38))
  glucose_mg_dL <- runif(n, 70, 180)
  insulin_pmol_L <- runif(n, 20, 180)
  tg_mmol_L <- runif(n, 0.5, 4.5)
  LDL <- runif(n, 1.2, 5.5)
  hdl_mmol_L <- runif(n, 0.6, 2.2)
  GLM7 <- log10(glucose_mg_dL) + log10(insulin_pmol_L) + log10(tg_mmol_L) +
    log10(LDL) - log10(hdl_mmol_L) + log_age + log_bmi
  metab <- GLM7 - log_age - log_bmi
  eta <- -2.2 + 0.25 * GLM7 + 0.3 * log_age + 0.15 * (sex == "Female")
  dx <- rbinom(n, 1L, plogis(eta))
  data.frame(dx, weight = runif(n, 0.5, 2), stra, psu, sex, log_age, log_bmi,
             GLM7, metab, glucose_mg_dL, insulin_pmol_L, tg_mmol_L, LDL, hdl_mmol_L)
}
p <- make_side(3601L)
l <- make_side(3602L)
p_des <- make_survey_design(p)
l_des <- make_survey_design(l)
survey_check <- check_survey_reparameterization(p, design = p_des)
stopifnot(survey_check$status == "survey_weighted_equivalence_confirmed",
          survey_check$max_linear_predictor_difference < 1e-8,
          survey_check$max_coefficient_mapping_difference < 1e-8,
          survey_check$max_covariance_mapping_difference < 1e-8)
unit_check <- run_unit_scaling_sensitivity(
  p, l, development_design = p_des, later_design = l_des)
stopifnot(is.finite(unit_check$unit_shift),
          abs(unit_check$unit_shift - unit_check$expected_unit_shift) < 1e-12,
          unit_check$intercept_mapping_difference < 1e-8,
          unit_check$development_refit_max_lp_difference < 1e-8,
          unit_check$later_refit_max_lp_difference < 1e-8,
          abs(unit_check$coefficient_difference) < 1e-8,
          abs(unit_check$original_later_metrics$auc -
              unit_check$incorrectly_transferred_later_metrics$auc) < 1e-12,
          length(unit_check$lp_l_correct) == nrow(l),
          length(unit_check$lp_l_incorrect_transfer) == nrow(l))
cat("ANALYSIS_HELPER_REGRESSION_OK: function existence; synthetic survey equivalence; unit scaling; prediction mapping\n")
quit(save = "no", status = 0L, runLast = FALSE)
