# test_audit_fit.R (PURE ASCII)
# v34: locks for audit_fit(), the entry point that reads a fitted model object.
#
# Covers: extraction from glm / lm / svyglm / coxph / survreg, the aliased
# column drop, the offset and weights reports, ad-hoc expressions, refusal of
# an unsupported model class, and -- the one that matters -- record-by-record
# equivalence with the registry audit over every fittable registry record.
#
# Pure base R + survey + survival for the model objects; no NHANES data.

suppressPackageStartupMessages({
  library(survey); library(survival)
})
CD <- Sys.getenv("CIA_CODE_DIR", ".")
for (f in c("index_registry.R", "term_identity.R", "algebra_audit.R",
            "public_api.R", "registry_audit.R", "audit_fit.R"))
  source(file.path(CD, f))

PASS <- 0L; FAILS <- character()
pass <- function(msg) { PASS <<- PASS + 1L; cat("[PASS]", msg, "\n") }
check <- function(cond, msg) {
  if (!isTRUE(cond)) { FAILS <<- c(FAILS, msg); cat("[FAIL]", msg, "\n") }
  else pass(msg)
}

set.seed(20261025)
n <- 1200
mk <- function() rnorm(n)
d <- data.frame(
  sex = factor(sample(c("Male", "Female"), n, TRUE)),
  log_age = mk() + 1.65, log_bmi = mk() + 1.45,
  log_tg_mg = mk() + 4.7, log_glucose_mg = mk() + 4.4,
  tg_mmol_L = rlnorm(n, 0.1, 0.3), hdl_mmol_L = rlnorm(n, 0.2, 0.3),
  waist_cm = rnorm(n, 98, 14), height_cm = rnorm(n, 168, 9),
  BMI = rnorm(n, 28, 5), age = rnorm(n, 48, 15),
  log10_tg = mk(), log_tg = mk(), log_hdl = mk(),
  log_tg_mg_2 = mk(), log_mets_sum = mk(), log_hdl_mg = mk(),
  log_hscrp = mk(), log_albumin = mk(), log_monocytes = mk(),
  log_TyG = mk(), ln_bmi = mk(), log_waist = mk(), log_SII = mk(),
  log_neutrophils = mk(), log_lymphocytes = mk(), log_platelets = mk(),
  albumin_g_dL = rnorm(n, 4.2, 0.4),
  lymphocytes_cells_mm3 = rnorm(n, 2100, 600),
  glucose_mmol_L = rnorm(n, 5.4, 1.1), insulin_uIU_mL = rnorm(n, 11, 6),
  log_glucose_mmol = mk(), log_insulin = mk(),
  stringsAsFactors = FALSE)
d$TyG <- d$log_tg_mg + d$log_glucose_mg - log(2)
d$metab <- mk()
d$log_glucose_mol <- d$log_glucose_mmol
lp <- -6 + 0.9 * d$log_age + 0.7 * d$log_bmi + 0.4 * d$log_tg_mg
d$dx <- rbinom(n, 1, plogis(lp))
d$time <- rexp(n, exp(lp - 6)); d$status <- rbinom(n, 1, 0.8)
d$pw <- runif(n, 0.5, 3); d$off <- rnorm(n, 0, 0.2)

# ---- 1. glm, two different bases ------------------------------------------
f_min  <- glm(dx ~ sex + log_age + log_bmi, data = d, family = binomial())
f_comp <- glm(dx ~ sex + log_tg_mg + log_glucose_mg, data = d, family = binomial())

a1 <- suppressWarnings(audit_fit(f_min, "TyG"))
a2 <- suppressWarnings(audit_fit(f_comp, "TyG"))
check(identical(a1$classification, "rank_adding_exact"),
      "glm: TyG against the minimal base is rank_adding_exact")
check(identical(a2$classification, "rank_preserving_exact"),
      "glm: TyG against its registered component base is rank_preserving_exact")
check(isTRUE(a1$is_classification) && isTRUE(a2$is_classification),
      "glm: both classifications are one of the four classes")
check(!is.null(a2$coefficient_map), "glm: a coefficient map exists when decomposable")
check(identical(a1$model_object$family, "binomial"), "glm: family is read off the fit")
check(identical(a1$model_object$link, "logit"), "glm: link is read off the fit")
check(isFALSE(a1$model_object$index_is_in_fit), "glm: index is not in the minimal-base fit")
check(identical(a1$model_object$base_model_terms, c("sex", "log_age", "log_bmi")),
      "glm: base is reported as model TERMS, not dummy column names")

# ---- 2. lm, svyglm, coxph, survreg -----------------------------------------
f_lm <- lm(log_bmi ~ sex + log_tg_mg + log_glucose_mg, data = d)
check(identical(suppressWarnings(audit_fit(f_lm, "TyG"))$classification,
                "rank_preserving_exact"),
      "lm: same class as glm against the same base")

des <- svydesign(ids = ~1, strata = NULL, weights = ~pw, data = d)
f_svy <- svyglm(dx ~ sex + log_tg_mg + log_glucose_mg, design = des,
                family = quasibinomial())
a_svy <- suppressWarnings(audit_fit(f_svy, "TyG"))
check(identical(a_svy$classification, "rank_preserving_exact"),
      "svyglm: TyG against its component base is rank_preserving_exact")
check(grepl("survey design weights", a_svy$model_object$weights$weights_source),
      "svyglm: the weights report names survey design weights")

f_cox <- coxph(Surv(time, status) ~ sex + log_tg_mg + log_glucose_mg,
               data = d, weights = pw)
a_cox <- suppressWarnings(audit_fit(f_cox, "TyG"))
check(identical(a_cox$classification, "rank_preserving_exact"),
      "coxph: TyG against its component base is rank_preserving_exact")
check(grepl("coxph", a_cox$model_object$weights$weights_source),
      "coxph: the weights report names survival case weights")

f_sr <- survreg(Surv(time, status) ~ sex + log_tg_mg + log_glucose_mg, data = d)
check(identical(suppressWarnings(audit_fit(f_sr, "TyG"))$classification,
                "rank_preserving_exact"),
      "survreg: TyG against its component base is rank_preserving_exact")

# ---- 3. aliased columns are reported, not silently matched -----------------
d2 <- d; d2$log_tg_mg_dup <- d2$log_tg_mg
f_alias <- glm(dx ~ sex + log_tg_mg + log_glucose_mg + log_tg_mg_dup,
               data = d2, family = binomial())
a_alias <- suppressWarnings(audit_fit(f_alias, "TyG"))
check(length(a_alias$model_object$aliased_columns) > 0L,
      "aliasing: the fit reports an aliased column")
check("aliased_columns_present" %in% a_alias$failure_states,
      "aliasing: reported as an explicit failure state, never downgraded")

# ---- 4. offset --------------------------------------------------------------
f_off <- glm(dx ~ sex + log_tg_mg + log_glucose_mg, data = d,
             family = binomial(), offset = off)
a_off <- suppressWarnings(audit_fit(f_off, "TyG"))
check(isTRUE(a_off$model_object$offset$has_offset), "offset: detected")
check("offset_present" %in% a_off$failure_states,
      "offset: reported as an explicit failure state, never downgraded")

# ---- 5. ad-hoc expressions, no registry entry needed ------------------------
f_ad <- glm(dx ~ sex + tg_mmol_L + waist_cm, data = d, family = binomial())
a_ad1 <- suppressWarnings(audit_fit(f_ad, "my_sum",
                                    expression = "tg_mmol_L + waist_cm"))
check(identical(a_ad1$classification, "rank_preserving_exact"),
      "ad-hoc expression: fully embedded in this base is rank_preserving_exact")
f_ad2 <- glm(dx ~ sex + tg_mmol_L, data = d, family = binomial())
a_ad2 <- suppressWarnings(audit_fit(f_ad2, "my_sum",
                                    expression = "tg_mmol_L + waist_cm"))
check(identical(a_ad2$classification, "rank_adding_exact"),
      "ad-hoc expression: partly embedded in this base is rank_adding_exact")
err <- tryCatch(audit_fit(structure(list(), class = "nls"), "TyG"),
                error = function(e) conditionMessage(e))
check(grepl("unsupported model class", err),
      "unsupported model class is refused with a named error")

# ---- 6. THE ONE THAT MATTERS ------------------------------------------------
# Every fittable registry record: fit the registered component base as a real
# glm, read the class off the fitted object with audit_fit(), and require it
# to equal the class the registry audit derives from first principles.
ref <- reference_audit_classes()
exps <- product_expansions()
mismatch <- character(); covered <- 0L
for (id in names(ref)) {
  base_id <- sub("-b$", "", id)
  if (base_id %in% names(exps)) {
    mt <- tryCatch(v28_component_model_terms(base_id, v28_registry()),
                   error = function(e) NULL)
  } else {
    t <- derive_additive_terms(id)
    if (is.null(t)) next                      # outside class: nothing to fit
    mt <- tryCatch(v28_component_model_terms(base_id, v28_registry()),
                   error = function(e) NULL)
  }
  if (is.null(mt) || !all(mt %in% names(d))) next
  fit <- suppressWarnings(tryCatch(
    glm(dx ~ ., data = d[, unique(c(mt, "dx")), drop = FALSE],
        family = binomial()),
    error = function(e) NULL))
  if (is.null(fit)) next
  got <- suppressWarnings(tryCatch(
    audit_fit(fit, base_id)$classification, error = function(e) NA_character_))
  want <- audit_one_record(id)$actual_class
  covered <- covered + 1L
  if (!identical(got, want))
    mismatch <- c(mismatch, sprintf("%s: audit_fit=%s registry=%s", id, got, want))
}
check(length(mismatch) == 0L,
      paste0("audit_fit equals the registry audit on every fittable record (n=",
             covered, ")",
             if (length(mismatch)) paste0(" -- ", paste(mismatch, collapse = "; "))
             else ""))
check(covered >= 18L,
      paste0("the equivalence loop actually covered the registry (n=", covered, ")"))

cat("\nPASS", PASS, " FAIL", length(FAILS), "\n")
if (length(FAILS)) quit(status = 1)
