# Covariance-referenced CVAI Wald decomposition; analysis-only current artifact.
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")
CODE_DIR <- "E:/KimiData/kimi/Workspaces/论文写作/课题2_v38_BMC_MIDM投稿包/04_代码与GitHub仓库/composite-index-audit/analysis"
OUT <- "E:/KimiData/kimi/Workspaces/论文写作/work_revision_20261117_1831/topic2_results"
source(file.path(CODE_DIR, "index_registry.R"))
source(file.path(CODE_DIR, "model_evaluation.R"))
source(file.path(CODE_DIR, "data_build.R"))
trms <- c("age_f", "bmi_f", "w_f", "ltg_f", "hdl_f")
fixed <- c(1.03, 4.20, -2.88, 17.76, 4.66)
rows <- list()
for (period in c("P", "L")) {
  d <- build_period(period)
  vars <- c("age", "sex", "BMI", "waist_cm", "log10_tg", "hdl_mmol_L")
  dd <- d[complete.cases(d[, c(vars, "dx_diab", "weight", "stra", "psu")]), ]
  dd$sexFemale <- as.numeric(dd$sex == "Female")
  dd$age_f <- dd$age * dd$sexFemale
  dd$bmi_f <- dd$BMI * dd$sexFemale
  dd$w_f <- dd$waist_cm * dd$sexFemale
  dd$ltg_f <- dd$log10_tg * dd$sexFemale
  dd$hdl_f <- dd$hdl_mmol_L * dd$sexFemale
  des <- make_survey_design(dd)
  form <- dx_diab ~ age + sex + BMI + waist_cm + log10_tg + hdl_mmol_L +
    age_f + bmi_f + w_f + ltg_f + hdl_f
  fit <- svyglm(form, design = des, family = quasibinomial())
  b <- coef(fit)[trms]
  V <- vcov(fit)[trms, trms, drop = FALSE]
  ee <- eigen(V, symmetric = TRUE)
  stopifnot(min(ee$values) > 0)
  Vhalf <- ee$vectors %*% (t(ee$vectors) * sqrt(ee$values))
  Vinvhalf <- ee$vectors %*% (t(ee$vectors) / sqrt(ee$values))
  ct <- as.numeric(Vinvhalf %*% fixed)
  Q <- qr.Q(qr(t(matrix(ct, nrow = 1))), complete = TRUE)
  Dw <- t(Q[, 2:5, drop = FALSE]) %*% Vinvhalf
  Vi <- solve(V)
  precision_c <- as.numeric(crossprod(fixed, Vi %*% fixed))
  direction <- as.numeric(Vi %*% fixed / sqrt(precision_c))
  theta <- as.numeric(crossprod(fixed, Vi %*% b) / precision_c)
  z <- as.numeric(Dw %*% b)
  q_direction <- as.numeric(crossprod(direction, b))^2
  q_constraint <- sum(z^2)
  q_total <- as.numeric(crossprod(b, Vi %*% b))
  unrestricted <- regTermTest(fit, ~ age_f + bmi_f + w_f + ltg_f + hdl_f)
  denominator_df <- unname(unrestricted$ddf)
  stopifnot(is.finite(denominator_df), denominator_df > 0)
  bundled_fit <- svyglm(update(form, . ~ . - age_f - bmi_f - w_f - ltg_f - hdl_f + CVAI),
                         design = des, family = quasibinomial())
  bundled <- regTermTest(bundled_fit, ~CVAI)
  max_orthogonality_error <- max(abs(Dw %*% V %*% t(Dw) - diag(4)))
  max_null_error <- max(abs(Dw %*% fixed))
  decomposition_error <- abs(q_total - q_constraint - q_direction)
  stopifnot(max_orthogonality_error < 1e-8, max_null_error < 1e-10,
            decomposition_error < 1e-8)
  rows[[period]] <- data.frame(
    period = period, n = nrow(dd), design_df = degf(des),
    numerator_df_constraint = 4L, numerator_df_direction = 1L,
    numerator_df_total = 5L, denominator_df = denominator_df,
    theta_covariance_direction = theta,
    wald_chisq_4df = q_constraint, wald_chisq_1df = q_direction,
    wald_chisq_5df = q_total,
    F_4df = q_constraint / 4, F_1df = q_direction,
    F_5df = q_total / 5,
    constraint_p = pf(q_constraint / 4, 4, denominator_df, lower.tail = FALSE),
    covariance_direction_p = pf(q_direction, 1, denominator_df, lower.tail = FALSE),
    unrestricted_5df_p = pf(q_total / 5, 5, denominator_df, lower.tail = FALSE),
    bundled_composite_1df_p = unname(bundled$p),
    max_abs_D_times_fixed = max_null_error,
    max_abs_DVDt_minus_I = max_orthogonality_error,
    absolute_wald_decomposition_error = decomposition_error,
    note = paste("Covariance-referenced decomposition in the unrestricted interaction model;",
                 "the bundled composite-model coefficient test is a separate estimand."))
}
out <- do.call(rbind, rows)
write.csv(out, file.path(OUT, "cvai_covariance_wald_decomposition.csv"), row.names = FALSE)
print(out, row.names = FALSE)
cat("TOPIC2_CVAI_COVARIANCE_DECOMPOSITION_COMPLETE\n")
