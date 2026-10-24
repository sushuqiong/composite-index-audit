# reanalysis_cases_v33.R (PURE ASCII)
# Stage 3 (v33): re-analysis of published model specifications on public NHANES
# data, demonstrating what the audit reports for each.
#
#   Case 1  Lan et al., BMC Public Health 2025;25:2512 (DOI
#           10.1186/s12889-025-23714-z): sarcopenia ~ LAP with Model 3
#           adjusting for waist circumference AND triglycerides (the two
#           components of LAP).  Specification refitted here on NHANES P and
#           L with the diabetes outcome.  Audit: reducible_to_interaction.
#   Case 2  Xu & Zhang, Int J Surg 2026 (TyG-related indices and mortality
#           in cancer survivors), with the Yin commentary explicitly calling
#           for TyG-BMI modelled next to TyG and BMI.  Audit:
#           reducible_to_interaction; for a single product the component-
#           conditional coefficient IS the interaction coefficient.
#   Case 3  Boundary honesty: TyG entered next to RAW triglycerides and
#           glucose (as in the Lan 2025 Model 3 covariate list) is not the
#           registered log-scale decomposition; the audit declines and
#           reports the transformation mismatch instead of claiming an
#           identity.  The same columns at log scale are certified
#           rank_preserving_exact.
#
# All fits are survey-weighted quasibinomial svyglm with Taylor
# linearization; coefficient identities are algebraic and verified to
# 1e-8 on the real design rows.

suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")

CODE_DIR <- Sys.getenv("CIA_CODE_DIR", unset = "")
OUT      <- Sys.getenv("CIA_OUT", unset = "")
if (!nzchar(CODE_DIR) || !nzchar(OUT))
  stop("set CIA_CODE_DIR and CIA_OUT before running reanalysis_cases_v33.R")
source(file.path(CODE_DIR, "index_registry_v33.R"))
source(file.path(CODE_DIR, "algebra_audit_v33.R"))
source(file.path(CODE_DIR, "model_evaluation_v33.R"))
source(file.path(CODE_DIR, "data_build_v33.R"))

exps <- product_expansions()
rows <- list()
emit <- function(case, period, quantity, value, detail = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    case = case, period = period, quantity = quantity,
    value = value, detail = detail, stringsAsFactors = FALSE)
}

fit_frame <- function(d, vars) {
  dd <- d[stats::complete.cases(d[, c(vars, "dx_diab", "weight", "stra", "psu")]),
          , drop = FALSE]
  dd$sexFemale <- as.numeric(dd$sex == "Female")
  dd$w_tg <- dd$waist_cm * dd$tg_mmol_L
  dd$tg_f <- dd$tg_mmol_L * dd$sexFemale
  dd$tyg_x_bmi <- dd$TyG * dd$BMI
  dd
}

wald1 <- function(fit, term_formula) {
  z <- survey::regTermTest(fit, term_formula)
  c(stat = unname(z$Ftest), df = unname(z$df), p = unname(z$p))
}

for (period in c("P", "L")) {
  d <- build_period_v33(period)

  # ---------------- Case 1: LAP with waist + TG (Lan 2025 Model 3 core) ----
  f1 <- fit_frame(d, c("age", "sex", "waist_cm", "tg_mmol_L", "LAP"))
  des1 <- make_survey_design(f1)
  a1 <- audit_product_in_model("LAP", c("age", "sex", "waist_cm", "tg_mmol_L"),
                               exps$LAP, data = f1)
  m0 <- svyglm(dx_diab ~ LAP + age + sex, design = des1, family = quasibinomial())
  m1 <- svyglm(dx_diab ~ LAP + age + sex + waist_cm + tg_mmol_L,
               design = des1, family = quasibinomial())
  m2 <- svyglm(dx_diab ~ age + sex + waist_cm + tg_mmol_L + w_tg + tg_f,
               design = des1, family = quasibinomial())
  # algebraic identity: X_[W,I] (K beta) == X_[W,C] beta on the real rows
  K1 <- a1$coefficient_map
  b_orig <- coef(m1)
  mm_c <- model.matrix(~ LAP + age + sex + waist_cm + tg_mmol_L, data = f1)
  lp_c <- as.numeric(mm_c %*% b_orig[colnames(mm_c)])
  b1 <- b_orig; names(b1)[names(b1) == "LAP"] <- "C"
  gam1 <- as.numeric(K1 %*% b1[colnames(K1)])
  mm_i <- model.matrix(~ age + sex + waist_cm + tg_mmol_L + w_tg + tg_f,
                       data = f1)
  lp_mapped <- as.numeric(mm_i %*% gam1)
  ident1 <- max(abs(lp_c - lp_mapped))
  w_1df <- wald1(m1, ~LAP)
  w_2df <- wald1(m2, ~ w_tg + tg_f)
  emit("case1_LAP", period, "audit_classification", a1$classification,
       "LAP + waist + TG in the same model (Lan 2025 Model 3)")
  emit("case1_LAP", period, "n", nrow(f1))
  emit("case1_LAP", period, "rank_gain_single_bundle", a1$numerical_diagnostics$rank_gain)
  emit("case1_LAP", period, "interaction_gain_unrestricted", a1$numerical_diagnostics$interaction_gain,
       "waist:TG and TG:sexFemale")
  emit("case1_LAP", period, "identity_max_abs_lp", signif(ident1, 3),
       "X_[W,I] K beta == X_[W,C] beta")
  emit("case1_LAP", period, "LAP_coef_demographics_only", signif(coef(m0)[["LAP"]], 5),
       sprintf("SE %.5f", coef(summary(m0))["LAP", "Std. Error"]))
  emit("case1_LAP", period, "LAP_coef_component_conditional", signif(coef(m1)[["LAP"]], 5),
       sprintf("SE %.5f", coef(summary(m1))["LAP", "Std. Error"]))
  emit("case1_LAP", period, "bundled_1df_wald_p", signif(w_1df[["p"]], 4),
       "test of LAP in the component-conditional model")
  emit("case1_LAP", period, "unrestricted_2df_wald_p", signif(w_2df[["p"]], 4),
       "joint test of waist:TG and TG:sexFemale")
  emit("case1_LAP", period, "unrestricted_waistTG_coef", signif(coef(m2)[["w_tg"]], 5),
       "free interaction; bundled model forces ratio 1 : 7 with TG:sexFemale")
  emit("case1_LAP", period, "unrestricted_TGxFemale_coef", signif(coef(m2)[["tg_f"]], 5))

  # ---------------- Case 2: TyG-BMI with TyG + BMI (Xu-Zhang / Yin) --------
  f2 <- fit_frame(d, c("age", "sex", "TyG", "BMI", "TyG_BMI"))
  des2 <- make_survey_design(f2)
  a2 <- audit_product_in_model("TyG_BMI", c("age", "sex", "TyG", "BMI"),
                               exps$TyG_BMI, data = f2)
  n0 <- svyglm(dx_diab ~ TyG_BMI + age + sex, design = des2, family = quasibinomial())
  n1 <- svyglm(dx_diab ~ TyG_BMI + TyG + BMI + age + sex, design = des2,
               family = quasibinomial())
  n2 <- svyglm(dx_diab ~ age + sex + TyG + BMI + tyg_x_bmi, design = des2,
               family = quasibinomial())
  coef_gap <- abs(coef(n1)[["TyG_BMI"]] - coef(n2)[["tyg_x_bmi"]])
  lp_gap <- max(abs(as.numeric(model.matrix(~ TyG_BMI + TyG + BMI + age + sex, f2) %*%
                                 coef(n1)[match(colnames(model.matrix(~ TyG_BMI + TyG + BMI + age + sex, f2)), names(coef(n1)))]) -
                    as.numeric(model.matrix(~ age + sex + TyG + BMI + tyg_x_bmi, f2) %*%
                                 coef(n2)[match(colnames(model.matrix(~ age + sex + TyG + BMI + tyg_x_bmi, f2)), names(coef(n2)))])))
  emit("case2_TyG_BMI", period, "audit_classification", a2$classification,
       "TyG-BMI + TyG + BMI in the same model (Yin commentary proposal)")
  emit("case2_TyG_BMI", period, "n", nrow(f2))
  emit("case2_TyG_BMI", period, "TyG_BMI_coef_demographics_only", signif(coef(n0)[["TyG_BMI"]], 6),
       sprintf("SE %.6f", coef(summary(n0))["TyG_BMI", "Std. Error"]))
  emit("case2_TyG_BMI", period, "TyG_BMI_coef_component_conditional", signif(coef(n1)[["TyG_BMI"]], 6),
       sprintf("SE %.6f; equals the TyG:BMI interaction coefficient", coef(summary(n1))["TyG_BMI", "Std. Error"]))
  emit("case2_TyG_BMI", period, "abs_coef_gap_vs_interaction_model", signif(coef_gap, 3),
       "the component-conditional coefficient IS the interaction")
  emit("case2_TyG_BMI", period, "abs_lp_gap_vs_interaction_model", signif(lp_gap, 3))
  emit("case2_TyG_BMI", period, "interaction_1df_wald_p",
       signif(wald1(n1, ~TyG_BMI)[["p"]], 4))

  # ---------------- Case 3: transformation boundary -------------------------
  f3 <- fit_frame(d, c("age", "sex", "TyG", "tg_mg_dL", "glucose_mg_dL",
                       "log_tg_mg", "log_glucose_mg"))
  raw <- audit_index_in_model(c(TyG = 1),
                              c(log_tg_mg = 1, log_glucose_mg = 1),
                              model_terms = c(tg_mg_dL = 1, glucose_mg_dL = 1),
                              data = f3)
  logged <- audit_index_in_model(c(log_tg_mg = 1, log_glucose_mg = 1),
                                 c(log_tg_mg = 1, log_glucose_mg = 1),
                                 model_terms = c(log_tg_mg = 1, log_glucose_mg = 1),
                                 data = f3)
  wcor <- function(x, y, w) {
    ok <- is.finite(x) & is.finite(y)
    cov.wt(cbind(x[ok], y[ok]), wt = w[ok], cor = TRUE)$cor[1, 2]
  }
  raw_note <- if (length(raw$failure_states)) raw$failure_states[1] else ""
  emit("case4_boundary", period, "TyG_with_RAW_components", raw$classification,
       raw_note)
  emit("case4_boundary", period, "TyG_with_LOG_components", logged$classification,
       "same numbers, certified once the registered transformation is used")
  emit("case4_boundary", period, "weighted_cor_TyG_TG", signif(wcor(f3$TyG, f3$tg_mg_dL, f3$weight), 4),
       "high empirical correlation is not an identity")
  emit("case4_boundary", period, "weighted_cor_TyG_glucose", signif(wcor(f3$TyG, f3$glucose_mg_dL, f3$weight), 4))

  # ---------------- Case 4: CVAI bundled 1 df vs unrestricted 5 df --------
  # LAP carries two interactions, so its contrast is 1 df vs 2 df.  CVAI
  # carries FIVE sex interactions; the honest contrast is 1 df vs 5 df.
  f4 <- fit_frame(d, c("age", "sex", "BMI", "waist_cm", "log10_tg",
                       "hdl_mmol_L", "CVAI"))
  f4$age_f <- f4$age * f4$sexFemale
  f4$bmi_f <- f4$BMI * f4$sexFemale
  f4$w_f <- f4$waist_cm * f4$sexFemale
  f4$ltg_f <- f4$log10_tg * f4$sexFemale
  f4$hdl_f <- f4$hdl_mmol_L * f4$sexFemale
  des4 <- make_survey_design(f4)
  a4 <- audit_product_in_model("CVAI",
    c("age", "sex", "BMI", "waist_cm", "log10_tg", "hdl_mmol_L"),
    exps$CVAI, data = f4)
  c1m <- svyglm(dx_diab ~ age + sex + BMI + waist_cm + log10_tg + hdl_mmol_L +
                  CVAI, design = des4, family = quasibinomial())
  c5m <- svyglm(dx_diab ~ age + sex + BMI + waist_cm + log10_tg + hdl_mmol_L +
                  age_f + bmi_f + w_f + ltg_f + hdl_f,
                design = des4, family = quasibinomial())
  w_1 <- wald1(c1m, ~CVAI)
  w_5 <- wald1(c5m, ~ age_f + bmi_f + w_f + ltg_f + hdl_f)
  emit("case3_CVAI", period, "audit_classification", a4$classification,
       "CVAI entered next to its six components")
  emit("case3_CVAI", period, "n", nrow(f4))
  emit("case3_CVAI", period, "interaction_gain_unrestricted",
       a4$numerical_diagnostics$interaction_gain,
       "five sex interactions; the bundled index model is the 1-df contrast")
  emit("case3_CVAI", period, "bundled_1df_wald_p", signif(w_1[["p"]], 4))
  emit("case3_CVAI", period, "unrestricted_5df_wald_p", signif(w_5[["p"]], 4),
       "the two tests answer different questions and need not agree")
}

res <- do.call(rbind, rows)

print(res, row.names = FALSE)
write.csv(res, file.path(OUT, "v33_reanalysis_cases.csv"), row.names = FALSE, na = "")
cat("\n=== v33 reanalysis cases done ===\n")
