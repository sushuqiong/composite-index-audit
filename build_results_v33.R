# build_results_v33.R  (PURE ASCII)
# v33 single source of truth for every NHANES number in the manuscript.
#
# WHAT CHANGED FROM v32 (see the manuscript's revision notes)
#   1. A failed augmented fit is NO LONGER folded into "aliased".  In v32,
#      eval_candidate() set lp1 = lp0 whenever the augmented fit errored OR
#      the score column was dropped as aliased, which made every structural
#      zero deterministic by construction.  In v33 a failed fit returns
#      status = "fit_failed" and is excluded from the structural-zero count;
#      a converged fit with an aliased score column builds lp1 from the
#      ESTIMABLE columns of the augmented fit itself, and the zero is
#      verified by measurement (max |lp1 - lp0| and the replicate deltas).
#   2. Arm B now runs on ONE shared complete-case domain across all 24
#      primary indices (plus the single-component benchmarks and the metab
#      control), so "the same participants" is literally true: P_n and L_n
#      are identical in every row, and every contrast shares one set of
#      Rao-Wu replicate weights per period.
#   3. The audit class is computed against the base the row was actually
#      fitted with (audit_class_fitted_base), alongside the registered-base
#      reference label (audit_class_registry).
#   4. New exports: single-component benchmarks (Arm B0), testable
#      predictions (GLM7 == metab; SII == SII*BMI on the log scale), the
#      record-flow table, the multiplicity-family table, the registry
#      concordance audit, and the decision-threshold shift table.
#
# Paths come from the environment and fail closed:
#   CIA_CODE_DIR   directory holding the v33 .R sources
#   CIA_DATA_ROOT  directory holding the official NHANES XPT files
#   CIA_OUT        output directory for CSV / RDS
#   CIA_BOOTSTRAP  replicate count (default 1000; set 50 for a debug run)
#
# Seeds: one fixed seed (20261023) for every domain; replicate designs are
# built ONCE per analysis domain (keyed by the sorted SEQN list) and shared
# by every contrast on that domain.

suppressPackageStartupMessages({
  library(haven); library(survey); library(digest)
})
options(survey.lonely.psu = "adjust")

need_env <- function(nm) {
  v <- Sys.getenv(nm, unset = "")
  if (!nzchar(v)) stop("set ", nm, " before running build_results_v33.R")
  v
}
CODE_DIR  <- need_env("CIA_CODE_DIR")
DATA_ROOT <- need_env("CIA_DATA_ROOT")
OUT       <- need_env("CIA_OUT")
BREP      <- as.integer(Sys.getenv("CIA_BOOTSTRAP", unset = "1000"))
SEED      <- 20261023L
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(DATA_ROOT)) stop("CIA_DATA_ROOT does not exist: ", DATA_ROOT)

for (f in c("index_registry_v33.R", "model_evaluation_v33.R",
            "algebra_audit_v33.R", "term_identity_v33.R",
            "public_api_v33.R", "registry_audit_v33.R"))
  source(file.path(CODE_DIR, f), local = TRUE)

registry <- v28_registry()
primary_registry <- registry[registry$analysis_role == "primary_candidate",
                             , drop = FALSE]
base_terms_for <- function(index) v28_component_model_terms(index, registry)

v28_score_is_valid <- function(score, index) {
  if (!is.numeric(score)) score <- as.numeric(score)
  score[!is.finite(score)] <- NA_real_
  if (index %in% c("CAR", "MHR", "TyG_hsCRP", "hsCRP_BMI", "hsCRP_WC",
                   "SII_BMI", "NLR", "PLR", "SII", "SIRI"))
    score[score <= 0] <- NA_real_
  score
}

# The audit class of an index against an ARBITRARY base, computed by the real
# audit functions (never read from the reference labels).  Used to label each
# arm row with the class of the (index, base) pair that was actually fitted.
audit_class_for_base <- function(index, base_terms) {
  exps <- product_expansions()
  if (index %in% names(exps))
    return(audit_product_in_model(index, base_terms, exps[[index]])$classification)
  t <- derive_additive_terms(index)
  if (is.null(t)) return("outside_supported_class")
  emb <- t$index_terms[intersect(names(t$index_terms), base_terms)]
  audit_index_in_model(t$index_terms, emb, emb,
                       intercept = t$intercept)$classification
}

# ----------------------------- data loading -------------------------------
read_unique <- function(file, cols) {
  p <- file.path(DATA_ROOT, file)
  if (!file.exists(p)) stop("missing official XPT: ", p)
  x <- haven::read_xpt(p)
  if (!"SEQN" %in% names(x) || anyDuplicated(x$SEQN))
    stop("SEQN is not unique: ", file)
  x[, cols, drop = FALSE]
}

# build_period_v33() is byte-identical to the copy in data_build_v33.R.
build_period_v33 <- function(period = c("P", "L")) {
  period <- match.arg(period)
  p <- if (period == "P") {
    list(
      demo   = list(file = "P_DEMO.xpt",   cols = c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH3","SDMVPSU","SDMVSTRA")),
      glu    = list(file = "P_GLU.xpt",    cols = c("SEQN","WTSAFPRP","LBXGLU")),
      ins    = list(file = "P_INS.xpt",    cols = c("SEQN","LBDINSI")),
      trig   = list(file = "P_TRIGLY.xpt", cols = c("SEQN","LBDTRSI","LBDLDLSI")),
      hdl    = list(file = "P_HDL.XPT",    cols = c("SEQN","LBDHDDSI")),
      bmx    = list(file = "P_BMX.xpt",    cols = c("SEQN","BMXBMI","BMXWAIST","BMXHT")),
      biopro = list(file = "P_BIOPRO.xpt", cols = c("SEQN","LBXSAL")),
      hscrp  = list(file = "P_HSCRP.XPT",  cols = c("SEQN","LBXHSCRP")),
      cbc    = list(file = "P_CBC.XPT",    cols = c("SEQN","LBDMONO","LBDNENO","LBDLYMNO","LBXPLTSI","LBXWBCSI")),
      mcq    = list(file = "P_MCQ.xpt",    cols = c("SEQN","MCQ550")),
      diq    = list(file = "P_DIQ.xpt",    cols = c("SEQN","DIQ010")),
      weight = "WTSAFPRP")
  } else {
    list(
      demo   = list(file = "DEMO_L.xpt",   cols = c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH3","SDMVPSU","SDMVSTRA")),
      glu    = list(file = "GLU_L.xpt",    cols = c("SEQN","WTSAF2YR","LBXGLU")),
      ins    = list(file = "INS_L.xpt",    cols = c("SEQN","LBDINSI")),
      trig   = list(file = "TRIGLY_L.xpt", cols = c("SEQN","LBDTRSI","LBDLDLSI")),
      hdl    = list(file = "HDL_L.XPT",    cols = c("SEQN","LBDHDDSI")),
      bmx    = list(file = "BMX_L.xpt",    cols = c("SEQN","BMXBMI","BMXWAIST","BMXHT")),
      biopro = list(file = "BIOPRO_L.xpt", cols = c("SEQN","LBXSAL")),
      hscrp  = list(file = "HSCRP_L.XPT",  cols = c("SEQN","LBXHSCRP")),
      cbc    = list(file = "CBC_L.xpt",    cols = c("SEQN","LBDMONO","LBDNENO","LBDLYMNO","LBXPLTSI","LBXWBCSI")),
      mcq    = list(file = "MCQ_L.xpt",    cols = c("SEQN","MCQ550")),
      diq    = list(file = "DIQ_L.xpt",    cols = c("SEQN","DIQ010")),
      weight = "WTSAF2YR")
  }
  weight_name <- p$weight
  p$weight <- NULL
  dat <- lapply(p, function(z) read_unique(z$file, z$cols))
  d <- Reduce(function(x, y) merge(x, y, by = "SEQN", all.x = TRUE,
                                   sort = FALSE), dat)
  rn <- c(RIAGENDR = "gender", RIDAGEYR = "age", LBXGLU = "FBG",
          LBDINSI = "Insulin", LBDTRSI = "TG", LBDLDLSI = "LDL",
          LBDHDDSI = "HDL", BMXBMI = "BMI", BMXWAIST = "Waist",
          BMXHT = "Height", LBXSAL = "Albumin", LBXHSCRP = "hsCRP",
          LBDMONO = "Monocytes", LBDNENO = "Neutrophils",
          LBDLYMNO = "Lymphocytes", LBXPLTSI = "Platelets")
  names(d)[names(d) %in% names(rn)] <- rn[names(d)[names(d) %in% names(rn)]]
  d$sex <- factor(ifelse(d$gender == 1, "Male",
                         ifelse(d$gender == 2, "Female", NA_character_)),
                  levels = c("Male", "Female"))
  d$dx_gall <- ifelse(d$MCQ550 %in% c(1, 2), as.integer(d$MCQ550 == 1),
                      NA_integer_)
  d$dx_diab <- ifelse(d$DIQ010 %in% c(1, 2), as.integer(d$DIQ010 == 1),
                      NA_integer_)
  d$weight <- d[[weight_name]]
  d$stra <- as.character(d$SDMVSTRA)
  d$psu <- paste(d$stra, as.character(d$SDMVPSU), sep = "_")
  d <- compute_candidate_indices(d)
  d$period <- period
  ok <- is.finite(d$weight) & d$weight > 0 & !is.na(d$SDMVSTRA) &
    !is.na(d$SDMVPSU) & is.finite(d$age) & d$age >= 20
  d[ok, , drop = FALSE]
}

cat("building P and L datasets...\n")
Pd <- build_period_v33("P")
Ld <- build_period_v33("L")

RESULTS <- list()   # the single nested object everything is exported from

flow <- do.call(rbind, lapply(list(Pd, Ld), function(d) {
  data.frame(period = d$period[1], adult_n = nrow(d),
             gallstone_observed = sum(!is.na(d$dx_gall)),
             gallstone_cases = sum(d$dx_gall == 1, na.rm = TRUE),
             diabetes_observed = sum(!is.na(d$dx_diab)),
             diabetes_cases = sum(d$dx_diab == 1, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
print(flow)
RESULTS$flow <- flow
write.csv(flow, file.path(OUT, "v33_flow.csv"), row.names = FALSE)

make_des <- function(d) {
  dd <- d; dd$design_support <- TRUE
  make_full_survey_design(dd, "design_support")
}
P_full <- make_des(Pd)
L_full <- make_des(Ld)
design_df <- data.frame(period = c("P", "L"),
                        design_df = c(degf(P_full), degf(L_full)),
                        n_strata = c(length(unique(Pd$stra)),
                                     length(unique(Ld$stra))),
                        n_psu = c(length(unique(Pd$psu)),
                                  length(unique(Ld$psu))),
                        stringsAsFactors = FALSE)
print(design_df)
RESULTS$design_df <- design_df
write.csv(design_df, file.path(OUT, "v33_design_df.csv"), row.names = FALSE)

# Weighted period statistics (replaces the carried-over v31 table).
period_one <- function(des, period_label) {
  m_age <- svymean(~age, des, na.rm = TRUE)
  f_sex <- svymean(~sex, des, na.rm = TRUE)
  p_dm  <- svymean(~dx_diab, des, na.rm = TRUE)
  p_gs  <- svymean(~dx_gall, des, na.rm = TRUE)
  data.frame(
    period = period_label,
    n = nrow(des$variables),
    weighted_mean_age = as.numeric(m_age),
    women_weighted_pct = 100 * as.numeric(f_sex[["sexFemale"]]),
    diabetes_weighted_pct = 100 * as.numeric(p_dm),
    gallstone_weighted_pct = 100 * as.numeric(p_gs),
    stringsAsFactors = FALSE)
}
period_stats <- rbind(period_one(P_full, "P"), period_one(L_full, "L"))
print(period_stats)
RESULTS$period_stats <- period_stats
write.csv(period_stats, file.path(OUT, "v33_period_stats.csv"),
          row.names = FALSE)

# --------------- shared replicate designs, keyed by domain ----------------
repdes_cache <- new.env()
repdes_for <- function(des) {
  key <- digest::digest(paste(sort(des$variables$SEQN), collapse = ","))
  if (!exists(key, envir = repdes_cache)) {
    assign(key, make_raowu_replicates(des, replicates = BREP, seed = SEED),
           envir = repdes_cache)
  }
  get(key, envir = repdes_cache)
}

outcome_design <- function(des, dx_col) {
  out <- des
  v <- out$variables
  v$dx <- v[[dx_col]]
  out$variables <- v
  out
}

score_vector <- function(d, index, scale_override = NULL) {
  score <- v28_score_is_valid(d[[index]], index)
  scale <- if (is.null(scale_override))
    registry$model_scale[match(index, registry$index)] else scale_override
  if (identical(scale, "log") && index != "HOMA_IR_log")
    score <- ifelse(is.finite(score) & score > 0, log(score), NA_real_)
  score
}

prepare_cand <- function(d, index, base_terms, scale_override = NULL) {
  d$score <- score_vector(d, index, scale_override)
  req <- unique(c("dx", "weight", "stra", "psu", base_terms, "score"))
  d[complete.cases(d[, req, drop = FALSE]), , drop = FALSE]
}

sub_des <- function(des, index, base_terms, scale_override = NULL) {
  v <- des$variables
  v$score <- score_vector(v, index, scale_override)
  req <- unique(c("dx", "weight", "stra", "psu", base_terms, "score"))
  v$.candidate_domain <- complete.cases(v[, req, drop = FALSE])
  out <- des
  out$variables <- v
  subset_analysis_domain(out, ".candidate_domain")
}

# One contrast, evaluated once.
# v33: a failed augmented fit is a FAILURE, not a structural zero.  When the
# augmented fit converged but dropped the score column as aliased, lp1 is
# built from the estimable columns of that augmented fit itself - a genuine
# point of the augmented model space - and the structural zero is then
# VERIFIED by measurement (max |lp1 - lp0| and every replicate delta), not
# written in by hand.
eval_candidate <- function(p, l, p_des, l_des, base_terms, B) {
  bt <- paste(base_terms, collapse = " + ")
  f0 <- stats::as.formula(paste("dx ~", bt))
  f1 <- stats::as.formula(paste("dx ~", bt, "+ score"))
  fit0 <- tryCatch(fit_survey_logistic(p, f0, design = p_des),
                   error = function(e) e)
  if (inherits(fit0, "error"))
    return(list(status = "fit_failed", failure_reason = conditionMessage(fit0),
                P_auc_base = NA_real_, P_auc_aug = NA_real_,
                L_delta_auc = NA_real_, L_lo = NA_real_, L_hi = NA_real_,
                wald_p = NA_real_, wald_p_t15 = NA_real_, valid = 0L,
                wald_note = "base fit failed", aliased = NA,
                fit_failed = TRUE, max_abs_lp_difference = NA_real_))
  mm0 <- stats::model.matrix(stats::delete.response(stats::terms(f0)),
                             data = l)
  c0 <- stats::coef(fit0)
  lp0 <- as.numeric(mm0 %*% c0[match(colnames(mm0), names(c0))])
  y <- as.integer(l$dx); w <- as.numeric(l$weight)
  a0 <- weighted_auc(lp0, y, w)

  fit1 <- tryCatch(fit_survey_logistic(p, f1, design = p_des),
                   error = function(e) e)
  if (inherits(fit1, "error"))
    return(list(status = "fit_failed", failure_reason = conditionMessage(fit1),
                P_auc_base = a0, P_auc_aug = NA_real_,
                L_delta_auc = NA_real_, L_lo = NA_real_, L_hi = NA_real_,
                wald_p = NA_real_, wald_p_t15 = NA_real_, valid = 0L,
                wald_note = "augmented fit failed", aliased = NA,
                fit_failed = TRUE, max_abs_lp_difference = NA_real_))
  c1 <- stats::coef(fit1)
  aliased <- !("score" %in% names(c1)) || !is.finite(c1[["score"]])
  # Rebuilt by hand because predict.svyglm fails after a silently aliased
  # column is dropped.  Only estimable columns contribute.
  mm1 <- stats::model.matrix(stats::delete.response(stats::terms(f1)),
                             data = l)
  estimable <- intersect(colnames(mm1), names(c1)[is.finite(c1)])
  lp1 <- as.numeric(mm1[, estimable, drop = FALSE] %*% c1[estimable])
  a1 <- weighted_auc(lp1, y, w)
  d_auc <- a1 - a0
  repdes <- repdes_for(l_des)
  rw <- weights(repdes, type = "analysis")
  boot <- vapply(seq_len(ncol(rw)), function(b) {
    wb <- as.numeric(rw[, b])
    weighted_auc(lp1, y, wb) - weighted_auc(lp0, y, wb)
  }, numeric(1))
  okb <- is.finite(boot)
  ci <- as.numeric(stats::quantile(boot[okb], c(.025, .975)))
  max_lp <- max(abs(lp1 - lp0))
  verified_zero <- aliased && is.finite(max_lp) && max_lp <= 1e-6 &&
    all(abs(boot[okb]) <= 1e-10)
  if (verified_zero) {
    lo <- 0; hi <- 0; pv <- NA_real_; pvt <- NA_real_
    note <- sprintf(paste0("structural zero VERIFIED: max |lp1-lp0| = %.2e ",
                           "on %d later-period rows; all %d replicate deltas ",
                           "<= 1e-10"),
                    max_lp, length(lp0), sum(okb))
    status <- "structural_zero_verified"
  } else {
    lo <- ci[1]; hi <- ci[2]
    wt <- rao_wu_wald_test(d_auc, boot, repdes$scale, repdes$rscales,
                           mse = repdes$mse, minimum_valid = 1L)
    pv <- wt$p_value
    pvt <- if (is.finite(pv)) 2 * stats::pt(-abs(d_auc / wt$se),
                                            df = survey::degf(l_des))
           else NA_real_
    note <- paste0("Rao-Wu percentile CI on shared domain replicate weights; ",
                   "Wald p from svrVar (valid=", wt$valid, ")")
    if (aliased)
      note <- paste0(note, "; score column aliased yet max |lp1-lp0| = ",
                     sprintf("%.2e", max_lp), " - reported, not certified")
    status <- "estimated"
  }
  list(P_auc_base = a0, P_auc_aug = a1, L_delta_auc = d_auc,
       L_lo = lo, L_hi = hi, wald_p = pv, wald_p_t15 = pvt,
       valid = sum(okb), wald_note = note, aliased = aliased,
       fit_failed = FALSE, max_abs_lp_difference = max_lp,
       status = status)
}

# domain_override: when non-NULL, a list(P = list(df, des), L = list(df, des))
# supplying ONE shared complete-case domain for every contrast in the arm.
run_arm <- function(dx_col, base_fun, arm_label, only = NULL,
                    scale_overrides = list(), B,
                    domain_override = NULL) {
  Pdes <- outcome_design(P_full, dx_col)
  Ldes <- outcome_design(L_full, dx_col)
  rows <- list()
  for (i in seq_len(nrow(primary_registry))) {
    index <- primary_registry$index[i]
    if (!is.null(only) && !(index %in% only)) next
    sc_ov <- scale_overrides[[index]]
    scale_lab <- if (is.null(sc_ov))
      registry$model_scale[match(index, registry$index)] else sc_ov
    bt <- base_fun(index)
    ref_class <- reference_audit_classes()[[index]]
    fit_class <- audit_class_for_base(index, bt)
    if (!is.null(sc_ov)) fit_class <- paste0(fit_class, "_at_", sc_ov,
                                             "_scale")
    if (is.null(domain_override)) {
      p <- prepare_cand(Pdes$variables, index, bt, scale_override = sc_ov)
      l <- prepare_cand(Ldes$variables, index, bt, scale_override = sc_ov)
    } else {
      p <- domain_override$P$df; l <- domain_override$L$df
      p$score <- domain_override$P$scores[[index]]
      l$score <- domain_override$L$scores[[index]]
    }
    gate <- nrow(p) >= 100 && nrow(l) >= 100 &&
      sum(p$dx == 1) >= 20 && sum(l$dx == 1) >= 20 &&
      length(unique(p$psu)) >= 10 && length(unique(l$psu)) >= 10
    if (!gate) {
      rows[[length(rows) + 1L]] <- data.frame(
        arm = arm_label, index = index, audit_class_registry = ref_class,
        audit_class_fitted_base = fit_class, model_scale = scale_lab,
        P_n = nrow(p), P_cases = sum(p$dx == 1),
        L_n = nrow(l), L_cases = sum(l$dx == 1), P_auc_base = NA_real_,
        P_auc_aug = NA_real_, L_delta_auc = NA_real_, L_lo = NA_real_,
        L_hi = NA_real_, wald_p = NA_real_, wald_p_t15 = NA_real_, valid = 0L,
        wald_note = "", aliased = NA, max_abs_lp_difference = NA_real_,
        status = "availability_gate", stringsAsFactors = FALSE)
      next
    }
    if (is.null(domain_override)) {
      l_sub <- sub_des(Ldes, index, bt, scale_override = sc_ov)
      p_sub <- sub_des(Pdes, index, bt, scale_override = sc_ov)
    } else {
      p_sub <- domain_override$P$des; l_sub <- domain_override$L$des
      p_sub$variables$score <- p$score; l_sub$variables$score <- l$score
    }
    z <- tryCatch(
      eval_candidate(p, l, p_sub, l_sub, bt, BREP),
      error = function(e) list(P_auc_base = NA_real_, P_auc_aug = NA_real_,
                               L_delta_auc = NA_real_, L_lo = NA_real_,
                               L_hi = NA_real_, wald_p = NA_real_,
                               wald_p_t15 = NA_real_, valid = 0L,
                               wald_note = paste("fit error:",
                                                 conditionMessage(e)),
                               aliased = NA, fit_failed = TRUE,
                               max_abs_lp_difference = NA_real_,
                               status = "fit_failed"))
    rows[[length(rows) + 1L]] <- data.frame(
      arm = arm_label, index = index, audit_class_registry = ref_class,
      audit_class_fitted_base = fit_class, model_scale = scale_lab,
      P_n = nrow(p), P_cases = sum(p$dx == 1),
      L_n = nrow(l), L_cases = sum(l$dx == 1), P_auc_base = z$P_auc_base,
      P_auc_aug = z$P_auc_aug, L_delta_auc = z$L_delta_auc, L_lo = z$L_lo,
      L_hi = z$L_hi, wald_p = z$wald_p, wald_p_t15 = z$wald_p_t15,
      valid = z$valid, wald_note = z$wald_note,
      aliased = z$aliased, max_abs_lp_difference = z$max_abs_lp_difference,
      status = z$status, stringsAsFactors = FALSE)
    cat(sprintf("[%s] %-12s n=%d/%d cases=%d/%d dAUC=%+.4f (%+.4f, %+.4f) [%s]\n",
                arm_label, index, nrow(p), nrow(l), sum(p$dx == 1),
                sum(l$dx == 1), z$L_delta_auc, z$L_lo, z$L_hi, z$status))
  }
  do.call(rbind, rows)
}

bh <- function(p) {
  p <- as.numeric(p); m <- length(p)
  o <- order(p); ps <- p[o]
  raw <- ps * m / seq_len(m)
  adj <- rev(cummin(rev(raw)))
  q <- pmin(1, adj)
  out <- numeric(m); out[o] <- q; out
}

bh_table <- function(a, family_label, family_id) {
  testable <- a[a$status == "estimated" & is.finite(a$wald_p), , drop = FALSE]
  if (!nrow(testable)) return(NULL)
  testable$bh_q <- round(bh(testable$wald_p), 6)
  testable$family <- family_label
  testable$family_id <- family_id
  testable$family_m <- nrow(testable)
  testable[order(testable$wald_p),
           c("arm", "index", "audit_class_registry", "audit_class_fitted_base",
             "model_scale", "P_n", "L_n",
             "L_delta_auc", "L_lo", "L_hi", "wald_p", "wald_p_t15", "bh_q",
             "family", "family_id", "family_m")]
}

# ------------------------------- Arm A ------------------------------------
cat("\n=== Arm A: diabetes outcome, registry component models ===\n")
armA <- run_arm("dx_diab", base_terms_for, "A_diabetes_registry", B = BREP)
RESULTS$armA <- armA
write.csv(armA, file.path(OUT, "v33_armA_diabetes_registry.csv"),
          row.names = FALSE)

# ------------------ Arm B on ONE shared domain ----------------------------
# The shared complete-case domain: all 24 primary scores, the minimal base,
# the outcome and the design columns non-missing.  Every Arm B contrast (and
# the B0 benchmarks and the metab control below) runs on exactly these rows.
cat("\n=== Arm B: diabetes outcome, minimal base, SHARED domain ===\n")
minimal_base_terms <- c("sex", "log_age", "log_bmi")
minimal_base <- function(index) minimal_base_terms

shared_domain <- function(dx_col) {
  build_side <- function(des_full) {
    des <- outcome_design(des_full, dx_col)
    v <- des$variables
    S <- as.data.frame(lapply(primary_registry$index, function(idx)
      score_vector(v, idx)))
    names(S) <- primary_registry$index
    req <- c("dx", "weight", "stra", "psu", minimal_base_terms)
    mask <- complete.cases(S) &
      complete.cases(v[, req, drop = FALSE])
    # metab is needed by the GLM7 testable prediction; include it in the mask.
    mask <- mask & is.finite(v$metab) & is.finite(v$GLM7)
    vv <- v[mask, , drop = FALSE]
    dd <- des
    dd$variables$.shared_domain <- mask
    sub <- subset_analysis_domain(dd, ".shared_domain")
    list(df = vv, des = sub, scores = S[mask, , drop = FALSE])
  }
  list(P = build_side(P_full), L = build_side(L_full))
}
domB <- shared_domain("dx_diab")
cat("shared Arm B domain n P/L:", nrow(domB$P$df), "/", nrow(domB$L$df),
    " cases:", sum(domB$P$df$dx == 1), "/", sum(domB$L$df$dx == 1), "\n")
stopifnot(nrow(domB$P$df) >= 100, nrow(domB$L$df) >= 100)

armB <- run_arm("dx_diab", minimal_base, "B_diabetes_minimalbase", B = BREP,
                domain_override = domB)
RESULTS$armB <- armB
write.csv(armB, file.path(OUT, "v33_armB_diabetes_minimalbase.csv"),
          row.names = FALSE)
# The shared-domain promise, checked in the output itself.
stopifnot(length(unique(armB$P_n)) == 1L, length(unique(armB$L_n)) == 1L)

# ---------------- Arm B0: single-component benchmarks ---------------------
# On the SAME shared domain: fasting glucose, triglycerides, insulin and HDL
# alone (log-transformed, the same transformations the composite indices
# use).  This is the benchmark the "tautology" critique asks for.
cat("\n=== Arm B0: single-component benchmarks on the shared domain ===\n")
b0_specs <- list(
  list(col = "log_glucose_mg", label = "fasting_glucose_log"),
  list(col = "log_tg_mg",      label = "triglycerides_log"),
  list(col = "log_insulin",    label = "insulin_log"),
  list(col = "log_hdl",        label = "hdl_log"))
b0_rows <- list()
for (spec in b0_specs) {
  p <- domB$P$df; l <- domB$L$df
  p$score <- p[[spec$col]]; l$score <- l[[spec$col]]
  p_sub <- domB$P$des; l_sub <- domB$L$des
  p_sub$variables$score <- p$score; l_sub$variables$score <- l$score
  z <- eval_candidate(p, l, p_sub, l_sub, minimal_base_terms, BREP)
  b0_rows[[length(b0_rows) + 1L]] <- data.frame(
    arm = "B0_single_components", index = spec$label,
    audit_class_registry = "single_component_benchmark",
    audit_class_fitted_base = "rank_adding_exact", model_scale = "log",
    P_n = nrow(p), P_cases = sum(p$dx == 1), L_n = nrow(l),
    L_cases = sum(l$dx == 1), P_auc_base = z$P_auc_base,
    P_auc_aug = z$P_auc_aug, L_delta_auc = z$L_delta_auc, L_lo = z$L_lo,
    L_hi = z$L_hi, wald_p = z$wald_p, wald_p_t15 = z$wald_p_t15,
    valid = z$valid, wald_note = z$wald_note, aliased = z$aliased,
    max_abs_lp_difference = z$max_abs_lp_difference, status = z$status,
    stringsAsFactors = FALSE)
  cat(sprintf("[B0] %-20s dAUC=%+.4f (%+.4f, %+.4f)\n",
              spec$label, z$L_delta_auc, z$L_lo, z$L_hi))
}
armB0 <- do.call(rbind, b0_rows)
RESULTS$armB0 <- armB0
write.csv(armB0, file.path(OUT, "v33_armB0_single_components.csv"),
          row.names = FALSE)

# ---------------- testable predictions ------------------------------------
# (i)  GLM7 against the minimal base MUST equal adding metab alone:
#      GLM7 = metab + log_age + log_bmi, so [base + GLM7] == [base + metab].
# (ii) On the log scale SII*BMI adds exactly log(SII) beyond the base:
#      log(SII*BMI) = log(SII) + ln(BMI) and ln(BMI) is a positive multiple
#      of log10(BMI), so the two augmented spaces coincide.
# (iii) hsCRP*BMI and hsCRP*WC differ, because the base holds ln(BMI) but
#      not log(waist).
cat("\n=== testable predictions ===\n")
pred_rows <- list()
run_single <- function(score_col, label) {
  p <- domB$P$df; l <- domB$L$df
  p$score <- p[[score_col]]; l$score <- l[[score_col]]
  p_sub <- domB$P$des; l_sub <- domB$L$des
  p_sub$variables$score <- p$score; l_sub$variables$score <- l$score
  z <- eval_candidate(p, l, p_sub, l_sub, minimal_base_terms, BREP)
  list(label = label, z = z)
}
g7 <- armB[armB$index == "GLM7", ]
mt <- run_single("metab", "metab alone")
pred_rows[[1]] <- data.frame(
  prediction = "GLM7 augmentation equals metab augmentation on the minimal base",
  lhs = g7$L_delta_auc, rhs = mt$z$L_delta_auc,
  abs_difference = abs(g7$L_delta_auc - mt$z$L_delta_auc),
  tolerance = 1e-8,
  passed = is.finite(g7$L_delta_auc) &&
    abs(g7$L_delta_auc - mt$z$L_delta_auc) < 1e-8 &&
    abs(g7$L_lo - mt$z$L_lo) < 1e-8 && abs(g7$L_hi - mt$z$L_hi) < 1e-8,
  stringsAsFactors = FALSE)
sii <- armB[armB$index == "SII", ]; siib <- armB[armB$index == "SII_BMI", ]
pred_rows[[2]] <- data.frame(
  prediction = "log(SII*BMI) augmentation equals log(SII) augmentation",
  lhs = sii$L_delta_auc, rhs = siib$L_delta_auc,
  abs_difference = abs(sii$L_delta_auc - siib$L_delta_auc),
  tolerance = 1e-8,
  passed = abs(sii$L_delta_auc - siib$L_delta_auc) < 1e-8 &&
    abs(sii$L_lo - siib$L_lo) < 1e-8 && abs(sii$L_hi - siib$L_hi) < 1e-8,
  stringsAsFactors = FALSE)
hb <- armB[armB$index == "hsCRP_BMI", ]; hw <- armB[armB$index == "hsCRP_WC", ]
pred_rows[[3]] <- data.frame(
  prediction = "hsCRP*BMI and hsCRP*WC differ (base holds ln BMI, not log waist)",
  lhs = hb$L_delta_auc, rhs = hw$L_delta_auc,
  abs_difference = abs(hb$L_delta_auc - hw$L_delta_auc),
  tolerance = 1e-10,
  passed = abs(hb$L_delta_auc - hw$L_delta_auc) > 1e-10,
  stringsAsFactors = FALSE)
testable <- do.call(rbind, pred_rows)
print(testable, row.names = FALSE)
stopifnot(all(testable$passed))
RESULTS$testable_predictions <- testable
write.csv(testable, file.path(OUT, "v33_testable_predictions.csv"),
          row.names = FALSE)

# ------------------------------- Arm G ------------------------------------
cat("\n=== Arm G: gallstone worked example, registry component models ===\n")
armG <- run_arm("dx_gall", base_terms_for, "G_gallstone_registry", B = BREP)
RESULTS$armG <- armG
write.csv(armG, file.path(OUT, "v33_armG_gallstone_registry.csv"),
          row.names = FALSE)

# ------------------------------- Arm C ------------------------------------
cat("\n=== Arm C: NLR and CAR at native scale on gallstone ===\n")
armC <- run_arm("dx_gall", base_terms_for, "C_gallstone_native_scale",
                only = c("NLR", "CAR"),
                scale_overrides = list(NLR = "native", CAR = "native"),
                B = BREP)
RESULTS$armC <- armC
write.csv(armC, file.path(OUT, "v33_armC_scale_pairing.csv"),
          row.names = FALSE)

# ------------------------------- Arm F ------------------------------------
cat("\n=== Arm F: GLM7 against base containing metab ===\n")
metab_base <- function(index) c("sex", "log_age", "log_bmi", "metab")
armF_diab <- run_arm("dx_diab", metab_base, "F_diabetes_metab_base",
                     only = "GLM7", B = BREP)
armF_gall <- run_arm("dx_gall", metab_base, "F_gallstone_metab_base",
                     only = "GLM7", B = BREP)
armF <- rbind(armF_diab, armF_gall)
RESULTS$armF <- armF
write.csv(armF, file.path(OUT, "v33_armF_glm7_metab_base.csv"),
          row.names = FALSE)

# ------------------ GLM7 dual-registration table --------------------------
pick <- function(df, index = "GLM7") df[df$index == index, , drop = FALSE]
dual <- rbind(
  cbind(outcome = "dx_diab", record = "GLM7-a", pick(armA)),
  cbind(outcome = "dx_diab", record = "GLM7-b", pick(armF_diab)),
  cbind(outcome = "dx_gall", record = "GLM7-a", pick(armG)),
  cbind(outcome = "dx_gall", record = "GLM7-b", pick(armF_gall)))
dual$registered_component_model <- c(
  rep("sex + log_age + log_bmi", 2), rep("sex + log_age + log_bmi", 2))
dual$registered_component_model[dual$record == "GLM7-b"] <-
  "sex + log_age + log_bmi + metab"
dual$audit_class_for_record <- ifelse(dual$record == "GLM7-a",
                                      "rank_adding_exact",
                                      "rank_preserving_exact")
RESULTS$glm7_dual_registration <- dual
write.csv(dual, file.path(OUT, "v33_glm7_dual_registration.csv"),
          row.names = FALSE)

# ------------------ audit classes (symbolic + numeric) --------------------
audit_rows <- list()
CONFIGS <- list(
  a = list(terms = c("log_age", "log_bmi")),
  b = list(terms = c("log_age", "log_bmi", "metab")))
aud_dat <- Pd[stats::complete.cases(Pd[, c("log_age", "log_bmi", "metab",
                                           "GLM7")]),
              c("log_age", "log_bmi", "metab", "GLM7")]
for (nm in names(CONFIGS)) {
  emb <- setNames(rep(1, length(CONFIGS[[nm]]$terms)), CONFIGS[[nm]]$terms)
  idx <- c(metab = 1, log_age = 1, log_bmi = 1)
  z <- audit_index_in_model(idx, emb, emb, data = aud_dat)
  audit_rows[[nm]] <- data.frame(
    record = paste0("GLM7-", nm),
    registered_component_model = paste(c("sex", CONFIGS[[nm]]$terms),
                                       collapse = " + "),
    audit_class = z$classification,
    rank_gain = z$numerical_diagnostics$rank_gain,
    residual_expression = if (is.na(z$residual_expression)) "" else
      z$residual_expression, stringsAsFactors = FALSE)
}
AUDIT <- do.call(rbind, audit_rows)
print(AUDIT, row.names = FALSE)
RESULTS$glm7_audit_classes <- AUDIT
write.csv(AUDIT, file.path(OUT, "v33_glm7_audit_classes.csv"),
          row.names = FALSE)

# ------------- registry concordance audit (26 records) --------------------
# Every record re-derived from its registered expression plus its full
# registered base through the real audit functions, compared against the
# hand-written reference labels.  The numeric columns use the P-period rows.
cat("\n=== registry concordance audit ===\n")
concordance <- audit_registry_concordance(data = NULL, tol = 1e-10)
print(concordance[, c("record", "expected_class", "actual_class",
                      "concordant")], row.names = FALSE)
stopifnot(nrow(concordance) == 26L, all(concordance$concordant))
RESULTS$registry_concordance <- concordance
write.csv(concordance, file.path(OUT, "v33_registry_audit_concordance.csv"),
          row.names = FALSE)

# ----------------------- Arm D aliased exhibit ----------------------------
cat("\n=== Arm D: aliased coefficient exhibit (gallstone outcome, P) ===\n")
ex_rows <- list()
for (index in c("NLR", "AIP", "SII")) {
  des <- outcome_design(P_full, "dx_gall")
  bt <- base_terms_for(index)
  d_sub <- sub_des(des, index, bt)
  bt_txt <- paste(bt, collapse = " + ")
  f_aug <- stats::as.formula(paste("dx ~", bt_txt, "+ score"))
  f_base <- stats::as.formula(paste("dx ~", bt_txt))
  fit <- survey::svyglm(f_aug, design = d_sub,
                        family = stats::quasibinomial())
  co <- summary(fit)$coefficients
  fit0 <- survey::svyglm(f_base, design = d_sub,
                         family = stats::quasibinomial())
  lp0 <- as.numeric(stats::fitted(fit0))
  lp1 <- as.numeric(stats::fitted(fit))
  for (r in seq_len(nrow(co))) {
    ex_rows[[length(ex_rows) + 1L]] <- data.frame(
      index = index, term = rownames(co)[r],
      estimate = round(co[r, 1], 6), se = round(co[r, 2], 6),
      aliased = is.na(co[r, 1]), note = "", stringsAsFactors = FALSE)
  }
  ex_rows[[length(ex_rows) + 1L]] <- data.frame(
    index = index, term = "max_abs_lp_difference(base_vs_augmented)",
    estimate = signif(max(abs(lp0 - lp1)), 3), se = NA_real_, aliased = FALSE,
    note = "identical linear predictors: coefficient of score is aliased (not identifiable)",
    stringsAsFactors = FALSE)
}
ex <- do.call(rbind, ex_rows)
print(ex)
RESULTS$armD <- ex
write.csv(ex, file.path(OUT, "v33_armD_aliased_exhibit.csv"),
          row.names = FALSE)

# ------------------------- unit scaling -----------------------------------
cat("\n=== Unit scaling: mapped intercept versus frozen intercept ===\n")
req_us <- c("dx_diab", "weight", "stra", "psu", "sex", "log_age", "log_bmi",
            "GLM7", "glucose_mg_dL", "insulin_pmol_L", "tg_mmol_L", "LDL",
            "hdl_mmol_L")
P_us <- Pd[stats::complete.cases(Pd[, req_us]), , drop = FALSE]
L_us <- Ld[stats::complete.cases(Ld[, req_us]), , drop = FALSE]
cat("unit scaling rows P/L:", nrow(P_us), "/", nrow(L_us), "\n")
us <- run_unit_scaling_sensitivity(P_us, L_us, outcome = "dx_diab",
                                   development_design = outcome_design(make_des(P_us), "dx_diab"),
                                   later_design = outcome_design(make_des(L_us), "dx_diab"))
RESULTS$unit_scaling <- us
us_tab <- data.frame(
  quantity = c("unit_shift", "expected_unit_shift",
               "expected_alternate_intercept", "observed_alternate_intercept",
               "intercept_mapping_difference",
               "development_refit_max_lp_difference",
               "later_refit_max_lp_difference",
               "incorrect_transfer_lp_shift",
               "original_later_auc", "original_later_brier",
               "refit_later_auc", "refit_later_brier",
               "incorrect_transfer_later_auc",
               "incorrect_transfer_later_brier"),
  value = c(us$unit_shift, us$expected_unit_shift,
            us$expected_alternate_intercept, us$observed_alternate_intercept,
            us$intercept_mapping_difference,
            us$development_refit_max_lp_difference,
            us$later_refit_max_lp_difference, us$incorrect_transfer_lp_shift,
            us$original_later_metrics[["auc"]],
            us$original_later_metrics[["brier"]],
            us$correctly_reparameterized_later_metrics[["auc"]],
            us$correctly_reparameterized_later_metrics[["brier"]],
            us$incorrectly_transferred_later_metrics[["auc"]],
            us$incorrectly_transferred_later_metrics[["brier"]]),
  stringsAsFactors = FALSE)
print(us_tab, row.names = FALSE)
RESULTS$unit_scaling_table <- us_tab
write.csv(us_tab, file.path(OUT, "v33_unit_scaling.csv"), row.names = FALSE)

# ------------- decision consequences of the mishandled transfer ------------
# The AUC is identical under the constant shift; the decisions are not.
cat("\n=== decision-threshold shift (MIDM demonstration) ===\n")
ts_tab <- threshold_shift_table(
  L_us, us$lp_l_correct, us$lp_l_incorrect_transfer,
  outcome = "dx_diab", weight = "weight",
  thresholds = c(0.05, 0.10, 0.15, 0.20))
print(ts_tab, row.names = FALSE)
RESULTS$threshold_shift <- ts_tab
write.csv(ts_tab, file.path(OUT, "v33_threshold_shift.csv"),
          row.names = FALSE)
ts_curve <- threshold_shift_table(
  L_us, us$lp_l_correct, us$lp_l_incorrect_transfer,
  outcome = "dx_diab", weight = "weight",
  thresholds = round(seq(0.01, 0.50, length.out = 100), 4))
write.csv(ts_curve, file.path(OUT, "v33_threshold_shift_curve.csv"),
          row.names = FALSE)

# ------------------------- BH corrections ---------------------------------
# Families are declared in v33_multiplicity_families.csv; the table below is
# built from the arm results, and every family statement in the manuscript
# cites one of these family_ids.
bhA  <- bh_table(armA, "Arm A registry-base estimable contrasts (10 + GLM7)",
                 "F_armA_11")
bhB3 <- bh_table(armB[armB$index %in% c("TyG", "TyG_BMI", "GLM7"),
                      , drop = FALSE],
                 "Arm B declared primary family", "F_primary3")
armB_all <- armB[armB$status == "estimated" & is.finite(armB$wald_p),
                 , drop = FALSE]
armB_all$bh_q <- round(bh(armB_all$wald_p), 6)
armB_all$family <- "Arm B all evaluated contrasts (sensitivity family)"
armB_all$family_id <- "F_armB_24"
armB_all$family_m <- nrow(armB_all)
armB_all <- armB_all[order(armB_all$wald_p),
  c("arm", "index", "audit_class_registry", "audit_class_fitted_base",
    "model_scale", "P_n", "L_n",
    "L_delta_auc", "L_lo", "L_hi", "wald_p", "wald_p_t15", "bh_q",
    "family", "family_id", "family_m")]
bhG  <- bh_table(armG, "Gallstone registry-base estimable contrasts",
                 "F_gallstone")
bhC  <- bh_table(armC, "Scale-pairing family", "F_scale2")
BH <- do.call(rbind, Filter(Negate(is.null),
                            list(bhA, bhB3, armB_all, bhG, bhC)))
print(BH, row.names = FALSE, digits = 4)
RESULTS$bh <- BH
write.csv(BH, file.path(OUT, "v33_bh_all.csv"), row.names = FALSE)

# ------------------- multiplicity family declaration -----------------------
families <- data.frame(
  family_id = c("F_primary3", "F_armB_24", "F_armA_11",
                "F_gallstone", "F_scale2"),
  label = c("Arm B primary family (predeclared)",
            "Arm B all evaluated contrasts (sensitivity)",
            "Arm A registry-base estimable contrasts",
            "Arm G gallstone estimable contrasts",
            "Arm C scale-pairing contrasts"),
  members = c(
    "TyG; TyG_BMI; GLM7",
    paste(armB_all$index, collapse = "; "),
    paste(bhA$index, collapse = "; "),
    paste(bhG$index, collapse = "; "),
    paste(bhC$index, collapse = "; ")),
  m = c(3L, nrow(armB_all), nrow(bhA), nrow(bhG), nrow(bhC)),
  rationale = c(
    "declared before results were seen; the three prespecified positive contrasts",
    "removes any concern that the 3-member family was chosen after seeing the results",
    "estimable contrasts against registry component bases; verified structural zeros are identities, not hypotheses, and are excluded by construction",
    "null-outcome worked example: no signal is expected",
    "NLR and CAR at the native scale"),
  stringsAsFactors = FALSE)
print(families[, c("family_id", "m")], row.names = FALSE)
RESULTS$multiplicity_families <- families
write.csv(families, file.path(OUT, "v33_multiplicity_families.csv"),
          row.names = FALSE)

# ------------------------- record flow -------------------------------------
# Every one of the 26 audit records, named, with where it is evaluated.
ref <- reference_audit_classes()
flow_rows <- lapply(names(ref), function(id) {
  role <- switch(id,
    GLM7 = "dual-registered: record GLM7-a evaluated in Arms A/B/F; record GLM7-b in Arm F",
    HOMA_IR_log = "scale-sensitivity registration; not one of the 24 primary contrasts",
    "primary candidate; evaluated in Arms A and B")
  evaluated <- switch(id,
    GLM7 = "A; B; F", HOMA_IR_log = "none (scale sensitivity)", "A; B")
  data.frame(record = id, reference_class = unname(ref[[id]]),
             role = role, evaluated_in = evaluated,
             stringsAsFactors = FALSE)
})
record_flow <- do.call(rbind, flow_rows)
record_flow <- rbind(record_flow, data.frame(
  record = "GLM7-b", reference_class = "rank_preserving_exact",
  role = "second GLM7 registration against the full embedded base",
  evaluated_in = "F", stringsAsFactors = FALSE))
n_rp <- sum(record_flow$reference_class == "rank_preserving_exact")
stopifnot(
  nrow(record_flow) == 26L,
  nrow(primary_registry) == 24L,
  n_rp == 15L,
  # 15 rank-preserving = 13 evaluated in Arm A + GLM7-b (Arm F) +
  # HOMA_IR_log (scale sensitivity, not a primary contrast)
  sum(record_flow$reference_class == "rank_preserving_exact" &
        record_flow$evaluated_in == "A; B") == 13L,
  sum(record_flow$reference_class == "reducible_to_interaction") == 5L,
  sum(record_flow$reference_class == "outside_supported_class") == 5L)
RESULTS$record_flow <- record_flow
write.csv(record_flow, file.path(OUT, "v33_record_flow.csv"),
          row.names = FALSE)

# ---------------------- consistency assertions ----------------------------
# v33: Arm A (per-index domain) and Arm B (shared 24-index domain) no longer
# evaluate GLM7 on the same rows, so the v32 byte-identity between the two
# arms is relaxed by design; what MUST still hold is that every contrast is
# internally consistent, the dual-registration table is assembled from the
# arm rows without refitting, and the structural zeros are verified.
gA <- armA[armA$index == "GLM7", ]; gB <- armB[armB$index == "GLM7", ]
gF <- armF_diab[armF_diab$index == "GLM7", ]
gG <- armG[armG$index == "GLM7", ]
gFg <- armF_gall[armF_gall$index == "GLM7", ]
stopifnot(
  is.finite(gA$L_delta_auc), is.finite(gB$L_delta_auc),
  gA$P_n == 3516, gB$P_n == nrow(domB$P$df),
  identical(gF$status, "structural_zero_verified"),
  abs(gF$L_delta_auc) < 1e-8, abs(gF$L_lo) < 1e-12, abs(gF$L_hi) < 1e-12,
  identical(gFg$status, "structural_zero_verified"),
  abs(gFg$L_delta_auc) < 1e-8,
  identical(dual$L_delta_auc[dual$record == "GLM7-a" &
            dual$outcome == "dx_diab"], gA$L_delta_auc),
  identical(dual$L_lo[dual$record == "GLM7-a" & dual$outcome == "dx_gall"],
            gG$L_lo)
)
cat("\nASSERTIONS PASSED: one contrast, one number, one interval per table;\n")
cat("the structural zeros are verified measurements, not assignments.\n")

saveRDS(RESULTS, file.path(OUT, "v33_results.rds"))
cat("\n=== v33 unified results done ===\n")
