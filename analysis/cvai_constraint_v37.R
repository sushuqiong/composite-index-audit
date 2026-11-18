# Topic 2 v37: CVAI fixed-direction constraints; no evidence-of-no-effect claim.
# Four independent contrasts D gamma = 0, D c = 0, tested in the unrestricted model.
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")
CODE_DIR <- Sys.getenv("CIA_CODE_DIR", unset = "")
OUT <- Sys.getenv("CIA_OUT", unset = "")
if (!nzchar(CODE_DIR) || !nzchar(OUT)) stop("set CIA_CODE_DIR and CIA_OUT")
source(file.path(CODE_DIR, "index_registry.R"))
source(file.path(CODE_DIR, "model_evaluation.R"))
source(file.path(CODE_DIR, "data_build.R"))
trms <- c("age_f", "bmi_f", "w_f", "ltg_f", "hdl_f")
fixed <- c(1.03, 4.20, -2.88, 17.76, 4.66)
# A deterministic orthonormal basis for the null space of t(fixed).
Q <- qr.Q(qr(t(matrix(fixed, nrow = 1))), complete = TRUE)
D <- t(Q[, 2:5, drop = FALSE])
stopifnot(all(dim(D) == c(4L, 5L)), max(abs(D %*% fixed)) < 1e-10,
          qr(D)$rank == 4L)
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
  z <- as.numeric(D %*% b)
  VD <- D %*% V %*% t(D)
  eig <- eigen(VD, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(min(eig) > 0)
  wald_chisq <- as.numeric(crossprod(z, solve(VD, z)))
  unrestricted <- regTermTest(fit, ~ age_f + bmi_f + w_f + ltg_f + hdl_f)
  denominator_df <- unname(unrestricted$ddf)
  stopifnot(is.finite(denominator_df), denominator_df > 0)
  F_constraint <- wald_chisq / 4
  p_constraint <- pf(F_constraint, 4, denominator_df, lower.tail = FALSE)
  bundled_form <- dx_diab ~ age + sex + BMI + waist_cm + log10_tg + hdl_mmol_L + CVAI
  bundled_fit <- svyglm(bundled_form, design = des, family = quasibinomial())
  bundled <- regTermTest(bundled_fit, ~CVAI)
  rows[[period]] <- data.frame(
    period = period, n = nrow(dd), rank_D = qr(D)$rank,
    max_abs_D_times_fixed = max(abs(D %*% fixed)), min_eigenvalue_DVDt = min(eig),
    wald_chisq_4df = wald_chisq, F_4df = F_constraint, numerator_df = 4L,
    denominator_df = denominator_df, design_df = degf(des),
    constraint_p = p_constraint, unrestricted_5df_p = unname(unrestricted$p),
    bundled_1df_p = unname(bundled$p),
    note = "Non-rejection of fixed-ratio constraints is not proof that they are correct.")
}
out <- do.call(rbind, rows)
write.csv(out, file.path(OUT, "nhanes_cvai_constraint_v37.csv"), row.names = FALSE)
print(out, row.names = FALSE)
cat("CVAI_CONSTRAINT_V37_COMPLETE\n")
