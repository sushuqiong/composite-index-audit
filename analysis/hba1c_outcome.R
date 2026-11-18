# HbA1c >= 5.7% post hoc outcome sensitivity analysis for CompositeIndexAudit.
# The binary outcome means HbA1c in the prediabetes-or-diabetes range; it is
# NOT a pure prediabetes outcome that excludes diabetes-range values.
# All 24 primary indices use one shared complete-case domain per period.

suppressPackageStartupMessages({library(foreign); library(survey); library(digest)})
options(survey.lonely.psu = "adjust")
need_env <- function(nm) {
  v <- Sys.getenv(nm, unset = "")
  if (!nzchar(v)) stop("set ", nm, " before running hba1c_outcome.R")
  v
}
CODE_DIR <- need_env("CIA_CODE_DIR")
DATA_ROOT <- need_env("CIA_DATA_ROOT")
OUT <- need_env("CIA_OUT")
BREP <- as.integer(Sys.getenv("CIA_BOOTSTRAP", unset = "1000"))
SEED <- 20261023L
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
for (f in c("index_registry.R", "model_evaluation.R", "algebra_audit.R",
            "term_identity.R", "public_api.R", "registry_audit.R", "data_build.R"))
  source(file.path(CODE_DIR, f), local = TRUE)

registry <- v28_registry()
primary_registry <- registry[registry$analysis_role == "primary_candidate", , drop = FALSE]
stopifnot(nrow(primary_registry) == 24L)
minimal_base_terms <- c("sex", "log_age", "log_bmi")

read_ghb <- function(period) {
  file <- if (period == "P") "P_GHB.XPT" else "GHB_L.XPT"
  x <- foreign::read.xport(file.path(DATA_ROOT, file))
  stopifnot(all(c("SEQN", "LBXGH") %in% names(x)), !anyDuplicated(x$SEQN))
  x[, c("SEQN", "LBXGH")]
}
Pd <- build_period("P"); Ld <- build_period("L")
Pd <- merge(Pd, read_ghb("P"), by = "SEQN", all.x = TRUE, sort = FALSE)
Ld <- merge(Ld, read_ghb("L"), by = "SEQN", all.x = TRUE, sort = FALSE)
add_outcome <- function(d) {
  d$dx_hba1c <- ifelse(is.finite(d$LBXGH), as.integer(d$LBXGH >= 5.7), NA_integer_)
  d
}
Pd <- add_outcome(Pd); Ld <- add_outcome(Ld)
stopifnot(all(Pd$dx_hba1c[is.finite(Pd$LBXGH)] == as.integer(Pd$LBXGH[is.finite(Pd$LBXGH)] >= 5.7)),
          all(Ld$dx_hba1c[is.finite(Ld$LBXGH)] == as.integer(Ld$LBXGH[is.finite(Ld$LBXGH)] >= 5.7)))
P_full <- make_full_survey_design(transform(Pd, design_support = TRUE), "design_support")
L_full <- make_full_survey_design(transform(Ld, design_support = TRUE), "design_support")

v28_score_is_valid <- function(score, index) {
  if (!is.numeric(score)) score <- as.numeric(score)
  score[!is.finite(score)] <- NA_real_
  if (index %in% c("CAR", "MHR", "TyG_hsCRP", "hsCRP_BMI", "hsCRP_WC",
                   "SII_BMI", "NLR", "PLR", "SII", "SIRI")) score[score <= 0] <- NA_real_
  score
}
score_vector <- function(d, index) {
  score <- v28_score_is_valid(d[[index]], index)
  scale <- registry$model_scale[match(index, registry$index)]
  if (identical(scale, "log") && index != "HOMA_IR_log")
    score <- ifelse(is.finite(score) & score > 0, log(score), NA_real_)
  score
}
audit_class_for_base <- function(index, base_terms) {
  exps <- product_expansions()
  if (index %in% names(exps))
    return(audit_product_in_model(index, base_terms, exps[[index]])$classification)
  t <- derive_additive_terms(index)
  if (is.null(t)) return("outside_supported_class")
  emb <- t$index_terms[intersect(names(t$index_terms), base_terms)]
  audit_index_in_model(t$index_terms, emb, emb, intercept = t$intercept)$classification
}
outcome_design <- function(des) {
  out <- des; out$variables$dx <- out$variables$dx_hba1c; out
}
repdes_cache <- new.env()
repdes_for <- function(des) {
  key <- digest::digest(paste(sort(des$variables$SEQN), collapse = ","))
  if (!exists(key, envir = repdes_cache))
    assign(key, make_raowu_replicates(des, replicates = BREP, seed = SEED), envir = repdes_cache)
  get(key, envir = repdes_cache)
}

build_side <- function(des_full) {
  des <- outcome_design(des_full); v <- des$variables
  S <- as.data.frame(lapply(primary_registry$index, function(idx) score_vector(v, idx)))
  names(S) <- primary_registry$index
  req <- c("dx", "weight", "stra", "psu", minimal_base_terms)
  mask <- complete.cases(S) & complete.cases(v[, req, drop = FALSE]) &
    is.finite(v$metab) & is.finite(v$GLM7)
  vv <- v[mask, , drop = FALSE]
  dd <- des; dd$variables$.shared_domain <- mask
  sub <- subset_analysis_domain(dd, ".shared_domain")
  list(df = vv, des = sub, scores = S[mask, , drop = FALSE])
}
dom <- list(P = build_side(P_full), L = build_side(L_full))
domain <- do.call(rbind, lapply(c("P", "L"), function(period) {
  d <- dom[[period]]$df
  data.frame(period = period, n = nrow(d), cases = sum(d$dx == 1L),
             noncases = sum(d$dx == 0L), weighted_case_pct = 100 * weighted.mean(d$dx, d$weight),
             design_df = survey::degf(dom[[period]]$des), stringsAsFactors = FALSE)
}))
print(domain)
stopifnot(all(domain$n >= 100), all(domain$cases >= 20), all(domain$noncases >= 20))
write.csv(domain, file.path(OUT, "hba1c_domain.csv"), row.names = FALSE)

eval_candidate <- function(p, l, p_des, l_des, base_terms) {
  bt <- paste(base_terms, collapse = " + ")
  f0 <- stats::as.formula(paste("dx ~", bt)); f1 <- stats::as.formula(paste("dx ~", bt, "+ score"))
  fit0 <- fit_survey_logistic(p, f0, design = p_des)
  fit1 <- fit_survey_logistic(p, f1, design = p_des)
  c0 <- stats::coef(fit0); c1 <- stats::coef(fit1)
  linear_predictor <- function(fit, formula, data, coefs) {
    mm <- stats::model.matrix(stats::delete.response(stats::terms(formula)), data = data)
    estimable <- intersect(colnames(mm), names(coefs)[is.finite(coefs)])
    as.numeric(mm[, estimable, drop = FALSE] %*% coefs[estimable])
  }
  p_lp0 <- linear_predictor(fit0, f0, p, c0)
  p_lp1 <- linear_predictor(fit1, f1, p, c1)
  lp0 <- linear_predictor(fit0, f0, l, c0)
  lp1 <- linear_predictor(fit1, f1, l, c1)
  p_y <- as.integer(p$dx); p_w <- as.numeric(p$weight)
  y <- as.integer(l$dx); w <- as.numeric(l$weight)
  p_a0 <- weighted_auc(p_lp0, p_y, p_w)
  p_a1 <- weighted_auc(p_lp1, p_y, p_w)
  a0 <- weighted_auc(lp0, y, w); a1 <- weighted_auc(lp1, y, w); d_auc <- a1 - a0
  repdes <- repdes_for(l_des); rw <- weights(repdes, type = "analysis")
  boot <- vapply(seq_len(ncol(rw)), function(b) {
    wb <- as.numeric(rw[, b]); weighted_auc(lp1, y, wb) - weighted_auc(lp0, y, wb)
  }, numeric(1))
  okb <- is.finite(boot); ci <- as.numeric(stats::quantile(boot[okb], c(.025, .975)))
  wt <- rao_wu_wald_test(d_auc, boot, repdes$scale, repdes$rscales, mse = repdes$mse, minimum_valid = 1L)
  pvt <- if (is.finite(wt$p_value)) 2 * stats::pt(-abs(d_auc / wt$se), df = survey::degf(l_des)) else NA_real_
  list(P_apparent_auc_base = p_a0, P_apparent_auc_aug = p_a1,
       L_auc_base = a0, L_auc_aug = a1, L_delta_auc = d_auc,
       L_lo = ci[1], L_hi = ci[2], wald_p = wt$p_value, wald_p_t15 = pvt,
       valid = sum(okb), aliased = !("score" %in% names(c1)) || !is.finite(c1[["score"]]),
       max_abs_lp_difference = max(abs(lp1 - lp0)))
}

# TyG_hsCRP contains TyG and therefore glucose information. This membership is
# fixed from the registered formulas, independently of the clinical outcome.
glucose_indices <- c("TyG", "TyG_WHtR", "TyG_BMI", "TyG_WC", "METS_IR", "HOMA_IR", "GLM7", "TyG_hsCRP")
rows <- vector("list", nrow(primary_registry))
for (i in seq_len(nrow(primary_registry))) {
  index <- primary_registry$index[i]
  p <- dom$P$df; l <- dom$L$df
  p$score <- dom$P$scores[[index]]; l$score <- dom$L$scores[[index]]
  p_sub <- dom$P$des; l_sub <- dom$L$des
  p_sub$variables$score <- p$score; l_sub$variables$score <- l$score
  z <- eval_candidate(p, l, p_sub, l_sub, minimal_base_terms)
  rows[[i]] <- data.frame(
    arm = "H_HbA1c_5.7_minimalbase", index = index,
    audit_class_registry = reference_audit_classes()[[index]],
    audit_class_fitted_base = audit_class_for_base(index, minimal_base_terms),
    model_scale = registry$model_scale[match(index, registry$index)],
    P_n = nrow(p), P_cases = sum(p$dx == 1), L_n = nrow(l), L_cases = sum(l$dx == 1),
    P_apparent_auc_base = z$P_apparent_auc_base, P_apparent_auc_aug = z$P_apparent_auc_aug,
    L_auc_base = z$L_auc_base, L_auc_aug = z$L_auc_aug,
    L_delta_auc = z$L_delta_auc, L_lo = z$L_lo, L_hi = z$L_hi,
    wald_p = z$wald_p, wald_p_t15 = z$wald_p_t15, valid = z$valid,
    glucose_component = index %in% glucose_indices,
    hba1c_outcome_component = FALSE, aliased = z$aliased,
    max_abs_lp_difference = z$max_abs_lp_difference,
    status = if (z$aliased) "aliased_review_required" else "estimated", stringsAsFactors = FALSE)
  cat(sprintf("[HbA1c] %-12s dAUC=%+.4f (%+.4f, %+.4f) p_t=%.6g\n",
              index, z$L_delta_auc, z$L_lo, z$L_hi, z$wald_p_t15), flush = TRUE)
}
arm <- do.call(rbind, rows)
bh <- function(p) {
  ok <- is.finite(p); out <- rep(NA_real_, length(p)); n <- sum(ok)
  if (!n) return(out)
  o <- order(p[ok]); ranked <- p[ok][o] * n / seq_len(n)
  ranked <- rev(cummin(rev(ranked))); out[ok][o] <- pmin(ranked, 1); out
}
arm$bh_q_24 <- bh(arm$wald_p_t15)
arm$bh_q_non_glucose_16 <- NA_real_
ng <- !arm$glucose_component
arm$bh_q_non_glucose_16[ng] <- bh(arm$wald_p_t15[ng])
arm$family_m_24 <- nrow(arm); arm$family_m_non_glucose_16 <- sum(ng)
stopifnot(length(unique(arm$P_n)) == 1L, length(unique(arm$L_n)) == 1L,
          nrow(arm) == 24L, sum(ng) == 16L, all(!arm$hba1c_outcome_component),
          all(arm$status == "estimated"), all(arm$valid >= 0.9 * BREP),
          all(is.finite(arm$L_delta_auc)), all(is.finite(arm$bh_q_24)),
          all(is.finite(arm$P_apparent_auc_base)), all(is.finite(arm$L_auc_base)),
          all(arm$bh_q_24 >= arm$wald_p_t15 - 1e-12),
          all(arm$bh_q_non_glucose_16[ng] >= arm$wald_p_t15[ng] - 1e-12),
          all(abs(arm$L_auc_aug - arm$L_auc_base - arm$L_delta_auc) < 1e-12))
arm <- arm[order(arm$wald_p_t15), , drop = FALSE]
write.csv(arm, file.path(OUT, "hba1c_armB.csv"), row.names = FALSE)

summary <- data.frame(
  metric = c("P_n", "L_n", "P_cases", "L_cases", "P_apparent_base_AUC",
             "L_base_AUC", "positive_delta_n", "negative_delta_n", "zero_delta_n",
             "q24_lt_0.05_n", "q16_lt_0.05_n", "minimum_valid_replicates"),
  value = c(domain$n[domain$period == "P"], domain$n[domain$period == "L"],
            domain$cases[domain$period == "P"], domain$cases[domain$period == "L"],
            arm$P_apparent_auc_base[1], arm$L_auc_base[1],
            sum(arm$L_delta_auc > 0), sum(arm$L_delta_auc < 0), sum(arm$L_delta_auc == 0),
            sum(arm$bh_q_24 < 0.05), sum(arm$bh_q_non_glucose_16 < 0.05, na.rm = TRUE), min(arm$valid)),
  stringsAsFactors = FALSE)
stopifnot(length(unique(arm$P_apparent_auc_base)) == 1L,
          length(unique(arm$L_auc_base)) == 1L)
write.csv(summary, file.path(OUT, "hba1c_summary.csv"), row.names = FALSE)
saveRDS(list(domain = domain, arm = arm, summary = summary,
             outcome = "LBXGH >= 5.7%; prediabetes-or-diabetes range, not pure prediabetes",
             evaluation = "P apparent AUC and L transported AUC are separately labelled; inference uses L delta AUC",
             replicates = BREP, seed = SEED), file.path(OUT, "hba1c_results.rds"))
print(summary)
cat("\nHBA1C_RUN_COMPLETE: 24 indices, 8 containing glucose and 16 not containing glucose; shared domains and paired Rao-Wu replicates.\n")
