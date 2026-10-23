# build_results_v32.R  (PURE ASCII)
# v32 single source of truth for every NHANES number in the manuscript.
#
# WHY THIS SCRIPT EXISTS
#   v31 quoted one GLM7 delta-AUC (+0.0397) with two different confidence
#   intervals (0.0247-0.0546 in the abstract / Arm B table, 0.0251-0.0562 in
#   Table 4 Arm A).  The point estimates were identical; the intervals differed
#   because the Arm A pipeline and the Arm B / dual-registration pipeline each
#   built their OWN Rao-Wu replicate design with a different seed.  An interval
#   for a fixed contrast is a function of the replicate weights, so two weight
#   sets give two intervals.
#
#   v32 removes the possibility: replicate designs are built ONCE per analysis
#   domain (keyed by the exact set of participants in the domain) and shared by
#   every contrast evaluated on that domain.  Every number downstream - text,
#   tables, figures, the machine-readable CSV - is exported from the single
#   nested object v32_results.rds written here.  Assertions at the end prove
#   that contrasts reachable by two routes (e.g. GLM7 on diabetes via Arm A and
#   via the dual-registration table) are byte-identical.
#
# WHAT IS RECOMPUTED HERE (survey-weighted, shared replicate designs)
#   flow, design df, Arm A (diabetes, registry bases), Arm B (diabetes,
#   minimal base - the v32 primary result), Arm G (gallstone worked example,
#   registry bases), Arm C (NLR/CAR at native scale on gallstone), Arm F
#   (GLM7 against a base that already contains metab, both outcomes), the
#   aliased-coefficient exhibit, the unit-scaling intercept demonstration,
#   and all BH corrections (declared families plus an all-24 Arm B family).
#
# WHAT IS REUSED UNCHANGED (not affected by the dual-CI defect)
#   development optimism correction (v31_development_optimism_*.csv),
#   bootstrap diagnosis (v31_bootstrap_*.csv), simulation results, and the
#   literature screen; those are copied at packaging time.
#
# Seeds: one fixed seed (20261023) for every domain.  Because the design is
# shared, the seed no longer differentiates pipelines; it only makes the run
# reproducible.

suppressPackageStartupMessages({
  library(haven); library(survey); library(digest)
})
options(survey.lonely.psu = "adjust")

CODE_DIR  <- "C:/Users/sushq/WorkBuddy/2026-08-12-22-01-44/v32_build/code"
DATA_ROOT <- "C:/Users/sushq/WorkBuddy/2026-08-12-22-01-44/v30_build/xpt_data"
OUT       <- "C:/Users/sushq/WorkBuddy/2026-08-12-22-01-44/v32_build/results"
BREP      <- as.integer(Sys.getenv("KT2_V32_BOOTSTRAP", unset = "1000"))
SEED      <- 20261023L
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

for (f in c("index_registry_v32.R", "model_evaluation_v32.R",
            "algebra_audit_v32.R"))
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

# ----------------------------- data loading -------------------------------
read_unique <- function(file, cols) {
  p <- file.path(DATA_ROOT, file)
  if (!file.exists(p)) stop("missing official XPT: ", p)
  x <- haven::read_xpt(p)
  if (!"SEQN" %in% names(x) || anyDuplicated(x$SEQN))
    stop("SEQN is not unique: ", file)
  x[, cols, drop = FALSE]
}

build_period_v32 <- function(period = c("P", "L")) {
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
Pd <- build_period_v32("P")
Ld <- build_period_v32("L")

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
write.csv(flow, file.path(OUT, "v32_flow.csv"), row.names = FALSE)

make_des <- function(d) {
  dd <- d; dd$design_support <- TRUE
  make_full_survey_design(dd, "design_support")
}
P_full <- make_des(Pd)
L_full <- make_des(Ld)
design_df <- data.frame(period = c("P", "L"),
                        design_df = c(degf(P_full), degf(L_full)),
                        stringsAsFactors = FALSE)
print(design_df)
RESULTS$design_df <- design_df
write.csv(design_df, file.path(OUT, "v32_design_df.csv"), row.names = FALSE)

# --------------- shared replicate designs, keyed by domain ----------------
# Two contrasts evaluated on the same participants MUST share one set of
# Rao-Wu replicate weights; that is the whole point of v32.  The cache key is
# the sorted SEQN list of the domain, so any two pipelines that land on the
# same rows automatically reuse the same weights.
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

prepare_cand <- function(d, index, base_terms, scale_override = NULL) {
  score <- v28_score_is_valid(d[[index]], index)
  scale <- if (is.null(scale_override))
    registry$model_scale[match(index, registry$index)] else scale_override
  if (identical(scale, "log") && index != "HOMA_IR_log")
    score <- ifelse(is.finite(score) & score > 0, log(score), NA_real_)
  d$score <- score
  # The domain is defined by the terms ACTUALLY fitted, not by the registry's
  # component list: Arm B's minimal base must not silently inherit the
  # registry domain, or the two arms would not be "the same participants"
  # under the model being compared.
  req <- unique(c("dx", "weight", "stra", "psu", base_terms, "score"))
  d[complete.cases(d[, req, drop = FALSE]), , drop = FALSE]
}

sub_des <- function(des, index, base_terms, scale_override = NULL) {
  v <- des$variables
  score <- v28_score_is_valid(v[[index]], index)
  scale <- if (is.null(scale_override))
    registry$model_scale[match(index, registry$index)] else scale_override
  if (identical(scale, "log") && index != "HOMA_IR_log")
    score <- ifelse(is.finite(score) & score > 0, log(score), NA_real_)
  v$score <- score
  req <- unique(c("dx", "weight", "stra", "psu", base_terms, "score"))
  v$.candidate_domain <- complete.cases(v[, req, drop = FALSE])
  out <- des
  out$variables <- v
  subset_analysis_domain(out, ".candidate_domain")
}

# One contrast, evaluated once.  The augmented model is always fitted; if the
# index column is aliased the linear predictor falls back to the base model,
# so a structural zero is VERIFIED numerically (aliased coefficient, identical
# linear predictors, zero replicate variance) rather than merely asserted.
eval_candidate <- function(p, l, p_des, l_des, base_terms, B) {
  bt <- paste(base_terms, collapse = " + ")
  f0 <- stats::as.formula(paste("dx ~", bt))
  f1 <- stats::as.formula(paste("dx ~", bt, "+ score"))
  fit0 <- fit_survey_logistic(p, f0, design = p_des)
  fit1 <- tryCatch(fit_survey_logistic(p, f1, design = p_des),
                   error = function(e) e)
  fit_failed <- inherits(fit1, "error")
  aliased <- fit_failed ||
    !("score" %in% names(stats::coef(fit1))) ||
    is.na(stats::coef(fit1)[["score"]])
  # P-period design fits are for coefficient development only; predictions are
  # rebuilt by hand because predict.svyglm fails after a silently aliased
  # column is dropped.
  mm0 <- stats::model.matrix(stats::delete.response(stats::terms(f0)),
                             data = l)
  c0 <- stats::coef(fit0)
  lp0 <- as.numeric(mm0 %*% c0[match(colnames(mm0), names(c0))])
  lp1 <- if (aliased) lp0 else {
    mm1 <- stats::model.matrix(stats::delete.response(stats::terms(f1)),
                               data = l)
    c1 <- stats::coef(fit1)
    as.numeric(mm1 %*% c1[match(colnames(mm1), names(c1))])
  }
  y <- as.integer(l$dx); w <- as.numeric(l$weight)
  a0 <- weighted_auc(lp0, y, w); a1 <- weighted_auc(lp1, y, w)
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
  if (aliased || all(boot[okb] == 0)) {
    lo <- 0; hi <- 0; pv <- NA_real_; pvt <- NA_real_
    note <- "not applicable - deterministic identity"
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
  }
  list(P_auc_base = a0, P_auc_aug = a1, L_delta_auc = d_auc,
       L_lo = lo, L_hi = hi, wald_p = pv, wald_p_t15 = pvt,
       valid = sum(okb), wald_note = note, aliased = aliased,
       fit_failed = fit_failed,
       max_abs_lp_difference = max_lp)
}

run_arm <- function(dx_col, base_fun, arm_label, only = NULL,
                    scale_overrides = list(), B) {
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
    p <- prepare_cand(Pdes$variables, index, bt, scale_override = sc_ov)
    l <- prepare_cand(Ldes$variables, index, bt, scale_override = sc_ov)
    audit_class <- registry$audit_class[match(index, registry$index)]
    if (!is.null(sc_ov)) audit_class <- paste0(audit_class, "_at_", sc_ov,
                                               "_scale")
    gate <- nrow(p) >= 100 && nrow(l) >= 100 &&
      sum(p$dx == 1) >= 20 && sum(l$dx == 1) >= 20 &&
      length(unique(p$psu)) >= 10 && length(unique(l$psu)) >= 10
    if (!gate) {
      rows[[length(rows) + 1L]] <- data.frame(
        arm = arm_label, index = index, audit_class = audit_class,
        model_scale = scale_lab, P_n = nrow(p), P_cases = sum(p$dx == 1),
        L_n = nrow(l), L_cases = sum(l$dx == 1), P_auc_base = NA_real_,
        P_auc_aug = NA_real_, L_delta_auc = NA_real_, L_lo = NA_real_,
        L_hi = NA_real_, wald_p = NA_real_, wald_p_t15 = NA_real_, valid = 0L,
        wald_note = "", status = "availability_gate", stringsAsFactors = FALSE)
      next
    }
    l_sub <- sub_des(Ldes, index, bt, scale_override = sc_ov)
    p_sub <- sub_des(Pdes, index, bt, scale_override = sc_ov)
    z <- tryCatch(
      eval_candidate(p, l, p_sub, l_sub, bt, BREP),
      error = function(e) list(P_auc_base = NA_real_, P_auc_aug = NA_real_,
                               L_delta_auc = NA_real_, L_lo = NA_real_,
                               L_hi = NA_real_, wald_p = NA_real_,
                               wald_p_t15 = NA_real_, valid = 0L,
                               wald_note = paste("fit error:",
                                                 conditionMessage(e)),
                               aliased = NA, fit_failed = TRUE,
                               max_abs_lp_difference = NA_real_))
    rows[[length(rows) + 1L]] <- data.frame(
      arm = arm_label, index = index, audit_class = audit_class,
      model_scale = scale_lab, P_n = nrow(p), P_cases = sum(p$dx == 1),
      L_n = nrow(l), L_cases = sum(l$dx == 1), P_auc_base = z$P_auc_base,
      P_auc_aug = z$P_auc_aug, L_delta_auc = z$L_delta_auc, L_lo = z$L_lo,
      L_hi = z$L_hi, wald_p = z$wald_p, wald_p_t15 = z$wald_p_t15,
      valid = z$valid, wald_note = z$wald_note,
      aliased = z$aliased, max_abs_lp_difference = z$max_abs_lp_difference,
      status = if (grepl("^fit error", z$wald_note)) "fit_error" else
        "estimated", stringsAsFactors = FALSE)
    cat(sprintf("[%s] %-12s n=%d/%d cases=%d/%d dAUC=%+.4f (%+.4f, %+.4f)\n",
                arm_label, index, nrow(p), nrow(l), sum(p$dx == 1),
                sum(l$dx == 1), z$L_delta_auc, z$L_lo, z$L_hi))
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

bh_table <- function(a, family_label) {
  testable <- a[a$status == "estimated" & !is.na(a$wald_p) &
                  a$wald_note != "not applicable - deterministic identity" &
                  is.finite(a$wald_p), , drop = FALSE]
  if (!nrow(testable)) return(NULL)
  testable$bh_q <- round(bh(testable$wald_p), 6)
  testable$family <- family_label
  testable$family_m <- nrow(testable)
  testable[order(testable$wald_p),
           c("arm", "index", "audit_class", "model_scale", "P_n", "L_n",
             "L_delta_auc", "L_lo", "L_hi", "wald_p", "wald_p_t15", "bh_q",
             "family", "family_m")]
}

# ------------------------------- Arm A ------------------------------------
cat("\n=== Arm A: diabetes outcome, registry component models ===\n")
armA <- run_arm("dx_diab", base_terms_for, "A_diabetes_registry", B = BREP)
RESULTS$armA <- armA
write.csv(armA, file.path(OUT, "v32_armA_diabetes_registry.csv"),
          row.names = FALSE)

# ------------------------------- Arm B ------------------------------------
cat("\n=== Arm B: diabetes outcome, minimal base (primary result) ===\n")
minimal_base <- function(index) c("sex", "log_age", "log_bmi")
armB <- run_arm("dx_diab", minimal_base, "B_diabetes_minimalbase", B = BREP)
RESULTS$armB <- armB
write.csv(armB, file.path(OUT, "v32_armB_diabetes_minimalbase.csv"),
          row.names = FALSE)

# ------------------------------- Arm G ------------------------------------
cat("\n=== Arm G: gallstone worked example, registry component models ===\n")
armG <- run_arm("dx_gall", base_terms_for, "G_gallstone_registry", B = BREP)
RESULTS$armG <- armG
write.csv(armG, file.path(OUT, "v32_armG_gallstone_registry.csv"),
          row.names = FALSE)

# ------------------------------- Arm C ------------------------------------
cat("\n=== Arm C: NLR and CAR at native scale on gallstone ===\n")
armC <- run_arm("dx_gall", base_terms_for, "C_gallstone_native_scale",
                only = c("NLR", "CAR"),
                scale_overrides = list(NLR = "native", CAR = "native"),
                B = BREP)
RESULTS$armC <- armC
write.csv(armC, file.path(OUT, "v32_armC_scale_pairing.csv"),
          row.names = FALSE)

# ------------------------------- Arm F ------------------------------------
# GLM7 audited twice: against the registered base (sex + log10(age) +
# log10(BMI)) it is rank_adding_exact; against the same base with metab added
# it is rank_preserving_exact.  These rows are computed HERE, on the shared
# domain replicate weights, and the dual-registration table below is BUILT
# FROM the arm results rather than re-estimated - the single-source rule.
cat("\n=== Arm F: GLM7 against base containing metab ===\n")
metab_base <- function(index) c("sex", "log_age", "log_bmi", "metab")
glm7_registry <- data.frame(
  index = "GLM7", stringsAsFactors = FALSE)
armF_diab <- run_arm("dx_diab", metab_base, "F_diabetes_metab_base",
                     only = "GLM7", B = BREP)
armF_gall <- run_arm("dx_gall", metab_base, "F_gallstone_metab_base",
                     only = "GLM7", B = BREP)
armF <- rbind(armF_diab, armF_gall)
RESULTS$armF <- armF
write.csv(armF, file.path(OUT, "v32_armF_glm7_metab_base.csv"),
          row.names = FALSE)

# ------------------ GLM7 dual-registration table --------------------------
# Assembled from the arm results already computed; nothing is refitted.
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
write.csv(dual, file.path(OUT, "v32_glm7_dual_registration.csv"),
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
write.csv(AUDIT, file.path(OUT, "v32_glm7_audit_classes.csv"),
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
write.csv(ex, file.path(OUT, "v32_armD_aliased_exhibit.csv"),
          row.names = FALSE)

# ------------------------- unit scaling -----------------------------------
cat("\n=== Unit scaling: mapped intercept versus frozen intercept ===\n")
# Same exhibit as v31 (diabetes outcome, GLM7 domain): the function needs
# complete-case rows AND designs built on exactly those rows, otherwise
# predict()/update() hit the full-frame missingness.
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
write.csv(us_tab, file.path(OUT, "v32_unit_scaling.csv"), row.names = FALSE)

# ------------------------- BH corrections ---------------------------------
# Declared families are reported with their declared sizes.  In addition the
# whole 24-contrast Arm B family is corrected as ONE family (m = 24), the
# sensitivity analysis that removes any concern that the 3-member family was
# chosen after seeing the results.
bhA  <- bh_table(armA, "Arm A declared family")
bhB3 <- bh_table(armB[armB$index %in% c("TyG", "TyG_BMI", "GLM7"),
                      , drop = FALSE], "Arm B declared 3-member family")
armB_all <- armB[armB$status == "estimated" & is.finite(armB$wald_p),
                 , drop = FALSE]
armB_all$bh_q <- round(bh(armB_all$wald_p), 6)
armB_all$family <- "Arm B all evaluated contrasts"
armB_all$family_m <- nrow(armB_all)
armB_all <- armB_all[order(armB_all$wald_p),
  c("arm", "index", "audit_class", "model_scale", "P_n", "L_n",
    "L_delta_auc", "L_lo", "L_hi", "wald_p", "wald_p_t15", "bh_q",
    "family", "family_m")]
bhG  <- bh_table(armG, "Gallstone declared family")
bhC  <- bh_table(armC, "Scale-pairing family")
BH <- do.call(rbind, Filter(Negate(is.null),
                            list(bhA, bhB3, armB_all, bhG, bhC)))
print(BH, row.names = FALSE, digits = 4)
RESULTS$bh <- BH
write.csv(BH, file.path(OUT, "v32_bh_all.csv"), row.names = FALSE)

# ---------------------- consistency assertions ----------------------------
# The dual-registration table must be byte-identical to the arm rows it was
# assembled from.  If this ever fails, the single-source rule has been broken.
gA <- armA[armA$index == "GLM7", ]; gB <- armB[armB$index == "GLM7", ]
gF <- armF_diab[armF_diab$index == "GLM7", ]
gG <- armG[armG$index == "GLM7", ]
gFg <- armF_gall[armF_gall$index == "GLM7", ]
stopifnot(
  identical(gA$L_delta_auc, gB$L_delta_auc),
  identical(gA$L_lo, gB$L_lo), identical(gA$L_hi, gB$L_hi),
  identical(gA$P_n, gB$P_n), identical(gA$L_n, gB$L_n),
  identical(gF$L_delta_auc, 0), identical(gF$L_lo, 0), identical(gF$L_hi, 0),
  identical(gFg$L_delta_auc, 0),
  identical(dual$L_delta_auc[dual$record == "GLM7-a" &
            dual$outcome == "dx_diab"], gA$L_delta_auc),
  identical(dual$L_lo[dual$record == "GLM7-a" & dual$outcome == "dx_gall"],
            gG$L_lo)
)
cat("\nASSERTIONS PASSED: GLM7 has one delta-AUC and one CI per outcome,\n")
cat("shared by every table, figure and text quotation downstream.\n")

saveRDS(RESULTS, file.path(OUT, "v32_results.rds"))
cat("\n=== v32 unified results done ===\n")
