# test_audit_v32.R -- synthetic-data locks for the v32 symbolic layer.
# Covers: four-class registry counts, reducible_to_interaction classification,
# LAP main-effect shift (-65/+7), CVAI sex-segmented expansion, coefficient
# maps K, rank bookkeeping, and regression of the v31 GLM7 behaviour.
# Pure base R; ASCII only.

source("algebra_audit_v32.R")
source("index_registry_v32.R")

tol <- 1e-8
pass <- function(msg) cat("[PASS]", msg, "\n")
check <- function(cond, msg) {
  if (!isTRUE(cond)) stop("[FAIL] ", msg)
  pass(msg)
}

# --- synthetic cohort -------------------------------------------------------
set.seed(20261023)
n <- 3000
d <- data.frame(
  sex = factor(sample(c("Male", "Female"), n, replace = TRUE), levels = c("Male", "Female")),
  age = runif(n, 20, 80),
  BMI = rnorm(n, 27, 5),
  waist_cm = rnorm(n, 95, 14),
  height_cm = rnorm(n, 168, 10),
  TG_mmol_L = rlnorm(n, log(1.3), 0.5),
  glucose_mmol_L = rnorm(n, 5.6, 1.2),
  insulin_uIU_mL = rlnorm(n, log(9), 0.5),
  hdl_mmol_L = rnorm(n, 1.3, 0.35)
)
d$TyG <- log((d$TG_mmol_L * 88.57) * (d$glucose_mmol_L * 18) / 2)
d$tg_mmol_L <- d$TG_mmol_L  # model-column name used by the registry term map
d$log10_tg <- log10(d$TG_mmol_L)
d$sexFemale <- as.numeric(d$sex == "Female")

# published formulas, exactly as compute_candidate_indices()
d$TyG_BMI <- d$TyG * d$BMI
d$TyG_WC <- d$TyG * d$waist_cm
d$LAP <- ifelse(d$sex == "Male", d$waist_cm - 65, d$waist_cm - 58) * d$TG_mmol_L
d$CVAI <- ifelse(d$sex == "Male",
                 -267.93 + 0.68 * d$age + 0.03 * d$BMI + 4.00 * d$waist_cm +
                   22.00 * d$log10_tg - 16.32 * d$hdl_mmol_L,
                 -187.32 + 1.71 * d$age + 4.23 * d$BMI + 1.12 * d$waist_cm +
                   39.76 * d$log10_tg - 11.66 * d$hdl_mmol_L)
d$HOMA_IR <- d$glucose_mmol_L * d$insulin_uIU_mL / 22.5
d$TyG_WHtR <- d$TyG * (d$waist_cm / d$height_cm)
d$log_mets_sum <- log(2 * (d$glucose_mmol_L * 18) + d$TG_mmol_L * 88.57)
d$log_hdl_mg <- log(d$hdl_mmol_L * 38.67)
d$METS_IR <- d$log_mets_sum * d$BMI / d$log_hdl_mg

# --- 1. registry four-class counts ------------------------------------------
reg <- v32_registry()
check(nrow(reg) == 25L, "registry has 25 rows")
check(all(c("audit_class_v32", "added_columns_v32") %in% names(reg)),
      "registry carries v32 audit columns")
# GLM7 dual registration adds one rank_preserving_exact record on top.
cls <- table(reg$audit_class_v32)
records <- c(rank_preserving_exact = unname(cls[["rank_preserving_exact"]]) + 1,
             rank_adding_exact = unname(cls[["rank_adding_exact"]]),
             reducible_to_interaction = unname(cls[["reducible_to_interaction"]]),
             outside_supported_class = unname(cls[["outside_supported_class"]]))
check(identical(as.integer(records), c(15L, 1L, 5L, 5L)),
      "26 audit records split 15 / 1 / 5 / 5")
check(identical(sort(names(v32_product_expansions())),
                sort(c("TyG_BMI", "TyG_WC", "LAP", "CVAI", "HOMA_IR"))),
      "expansions registered for exactly the five product-form indices")

# --- 2. classification + identity + ranks for the five reducible indices ----
exps <- v32_product_expansions()
model_terms_of <- function(idx) v28_component_model_terms(idx, reg)
for (idx in names(exps)) {
  a <- audit_product_in_model(idx, model_terms_of(idx), exps[[idx]], data = d, tol = tol)
  check(a$classification == "reducible_to_interaction",
        paste0(idx, " classified reducible_to_interaction"))
  check(is.finite(a$numerical_diagnostics$identity_error) &&
          a$numerical_diagnostics$identity_error < 1e-8,
        paste0(idx, " expansion reproduces the stored column"))
  k_new <- length(exps[[idx]]$products)
  check(a$numerical_diagnostics$rank_gain == 1L,
        paste0(idx, " augmenting adds exactly one bundled dimension"))
  check(a$numerical_diagnostics$interaction_gain == k_new,
        paste0(idx, " unrestricted interaction space adds ", k_new, " column(s)"))
  check(is.finite(a$numerical_diagnostics$c_span_residual) &&
          a$numerical_diagnostics$c_span_residual < 1e-6,
        paste0(idx, " C lies in span[W, I] to tolerance"))
  check(!is.null(a$coefficient_map), paste0(idx, " coefficient map K returned"))
}

# --- 3. outside_supported_class honesty --------------------------------------
a_wh <- audit_product_in_model("TyG_WHtR", model_terms_of("TyG_WHtR"),
                               list(products = c("TyG:waist_per_height" = 1)),
                               data = d, tol = tol)
check(a_wh$classification == "outside_supported_class" &&
        grepl("TyG:waist_per_height", a_wh$blocking_reason),
      "TyG_WHtR (ratio second factor) stays outside the supported class")
a_mets <- audit_product_in_model("METS_IR", model_terms_of("METS_IR"),
                                 list(products = c("log_mets_sum:BMI" = 1,
                                                   "log_hdl_mg_reciprocal" = 1)),
                                 data = d, tol = tol)
check(a_mets$classification == "outside_supported_class",
      "METS_IR (denominator log term) stays outside the supported class")

# --- 4. LAP: main-effect shift -65 and +7 interaction, algebraic identity ---
d$lp <- -3 + 0.02 * d$waist_cm + 0.4 * d$tg_mmol_L + 0.02 * d$LAP +
  0.3 * d$sexFemale
d$dx <- rbinom(n, 1, plogis(d$lp))
m_aug <- glm(dx ~ sex + waist_cm + tg_mmol_L + LAP, data = d, family = binomial())
d$w_tg <- d$waist_cm * d$tg_mmol_L
d$tg_f <- d$tg_mmol_L * d$sexFemale
K <- audit_product_in_model("LAP", model_terms_of("LAP"), exps$LAP)$coefficient_map
check(K["tg_mmol_L", "C"] == -65 && K["waist_cm:tg_mmol_L", "C"] == 1 &&
        K["tg_mmol_L:sexFemale", "C"] == 7,
      "LAP: K encodes the -65 main-effect shift and the 1 / 7 interaction ratios")
b <- coef(m_aug); names(b)[names(b) == "LAP"] <- "C"
gam_mapped <- as.numeric(K %*% b[colnames(K)])
names(gam_mapped) <- rownames(K)
mm_wi <- model.matrix(~ sex + waist_cm + tg_mmol_L + w_tg + tg_f, data = d)
stopifnot(ncol(mm_wi) == nrow(K))  # same order by construction
lp_mapped <- as.numeric(mm_wi %*% gam_mapped)
lp_diff <- max(abs(predict(m_aug, type = "link") - lp_mapped))
check(lp_diff < 1e-8,
      "LAP: X_[W,I] K beta reproduces the augmented linear predictor exactly")

# --- 5. CVAI: intercept and sexFemale main shifts -----------------------------
d$lp2 <- -2 + 0.01 * d$age + 0.05 * d$BMI + 0.02 * d$CVAI
d$dx2 <- rbinom(n, 1, plogis(d$lp2))
m_aug2 <- glm(dx2 ~ sex + age + BMI + waist_cm + log10_tg + hdl_mmol_L + CVAI,
              data = d, family = binomial())
d$age_f <- d$age * d$sexFemale; d$bmi_f <- d$BMI * d$sexFemale
d$w_f <- d$waist_cm * d$sexFemale; d$ltg_f <- d$log10_tg * d$sexFemale
d$hdl_f <- d$hdl_mmol_L * d$sexFemale
K2 <- audit_product_in_model("CVAI", model_terms_of("CVAI"), exps$CVAI)$coefficient_map
check(K2["(Intercept)", "C"] == -267.93 && K2["sexFemale", "C"] == 80.61 &&
        K2["age:sexFemale", "C"] == 1.03 && K2["hdl_mmol_L:sexFemale", "C"] == 4.66,
      "CVAI: K encodes the intercept, sexFemale and interaction published constants")
b2 <- coef(m_aug2); names(b2)[names(b2) == "CVAI"] <- "C"
gam2 <- as.numeric(K2 %*% b2[colnames(K2)])
names(gam2) <- rownames(K2)
mm_wi2 <- model.matrix(~ sex + age + BMI + waist_cm + log10_tg + hdl_mmol_L +
                         age_f + bmi_f + w_f + ltg_f + hdl_f, data = d)
stopifnot(ncol(mm_wi2) == nrow(K2))  # same order by construction
lp_mapped2 <- as.numeric(mm_wi2 %*% gam2)
check(max(abs(predict(m_aug2, type = "link") - lp_mapped2)) < 1e-8,
      "CVAI: X_[W,I] K beta reproduces the augmented linear predictor exactly")

# --- 6. v31 GLM7 behaviour is unchanged ---------------------------------------
# The additive audit keys on the index/embedded terms; other model columns Z
# (here sex) are outside its scope, matching how the pipeline calls it.
z <- audit_index_in_model(
  c(log_age = 1, log_bmi = 1, metab = 1),
  c(log_age = 1, log_bmi = 1))
check(z$classification == "rank_adding_exact",
      "GLM7 against the minimal base stays rank_adding_exact")
check(!is.null(z$coefficient_map) && z$coefficient_map["log_age", "C"] == 1,
      "GLM7 coefficient map still embeds log_age with coefficient 1")
# data = NULL leaves rank_gain NA, which the v31 branch maps to
# rank_adding_exact; supply rows so the rank check confirms gain 0.
d_tyg <- data.frame(log_tg_mg = rnorm(200), log_glucose_mg = rnorm(200))
z0 <- audit_index_in_model(
  c(log_tg_mg = 1, log_glucose_mg = 1),
  c(log_tg_mg = 1, log_glucose_mg = 1),
  data = d_tyg)
check(z0$classification == "rank_preserving_exact",
      "TyG against its own components stays rank_preserving_exact")

cat("=== all v32 audit tests passed ===\n")
