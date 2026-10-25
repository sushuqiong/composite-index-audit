# v28 survey-weighted model-evaluation functions.
# AUC is calculated from linear-predictor ranks. Brier and calibration use
# inverse-logit probabilities. Every paired difference uses the same replicate.
#
# v34 addition: threshold_shift_table() quantifies what a mishandled unit
# conversion does to DECISIONS, not to ranking: at a fixed risk threshold the
# number flagged high-risk and the number of true cases missed move, while
# the AUC does not. This is the decision-analytic consequence Figure 4 shows.

suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")

weighted_auc <- function(score, y, w) {
  ok <- is.finite(score) & !is.na(y) & is.finite(w) & w > 0
  score <- score[ok]; y <- as.integer(y[ok]); w <- as.numeric(w[ok])
  if (!length(score) || length(unique(y)) < 2) return(NA_real_)
  w1 <- sum(w[y == 1]); w0 <- sum(w[y == 0])
  if (w1 <= 0 || w0 <= 0) return(NA_real_)
  o <- order(score, y); score <- score[o]; y <- y[o]; w <- w[o]
  neg <- w * (1 - y); gid <- cumsum(c(TRUE, score[-1] != score[-length(score)]))
  group_neg <- rowsum(neg, gid)[, 1]; before <- cumsum(group_neg) - group_neg
  wins <- before[gid] + 0.5 * group_neg[gid]
  sum(wins[y == 1] * w[y == 1]) / (w1 * w0)
}

weighted_brier <- function(prob, y, w) {
  ok <- is.finite(prob) & !is.na(y) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(w[ok] * (prob[ok] - y[ok])^2) / sum(w[ok])
}

make_survey_design <- function(data, weight = "weight", strata = "stra", psu = "psu",
                               allow_zero_weights = FALSE) {
  required <- c(weight, strata, psu)
  if (!all(required %in% names(data))) stop("survey design columns are missing")
  d <- data
  d$.weight <- as.numeric(d[[weight]]); d$.strata <- factor(d[[strata]]); d$.psu <- factor(d[[psu]])
  if (any(!is.finite(d$.weight) |
          if (allow_zero_weights) d$.weight < 0 else d$.weight <= 0)) {
    stop("invalid survey weights")
  }
  survey::svydesign(ids = ~.psu, strata = ~.strata, weights = ~.weight, data = d, nest = TRUE)
}

make_full_survey_design <- function(data, support_indicator = "design_support",
                                    weight = "weight", strata = "stra", psu = "psu") {
  if (!support_indicator %in% names(data)) stop("design-support indicator is missing")
  support <- data[[support_indicator]] %in% TRUE
  if (!any(support)) stop("survey design has no supported records")
  des <- make_survey_design(data[support, , drop = FALSE], weight, strata, psu)
  attr(des, "support_indicator") <- support_indicator
  des
}

subset_analysis_domain <- function(design, indicator) {
  if (!inherits(design, "survey.design2")) stop("design must be a survey.design2 object")
  if (!indicator %in% names(design$variables)) stop("analysis-domain indicator is missing")
  d <- design
  d$variables$.analysis_domain <- d$variables[[indicator]] %in% TRUE
  out <- subset(d, .analysis_domain)
  if (!nrow(out$variables)) stop("analysis domain is empty")
  attr(out, "analysis_domain_indicator") <- indicator
  out
}

fit_survey_logistic <- function(data = NULL, formula, weight = "weight", strata = "stra", psu = "psu",
                                design = NULL) {
  des <- if (is.null(design)) make_survey_design(data, weight = weight, strata = strata, psu = psu) else design
  fit <- survey::svyglm(formula, design = des, family = stats::quasibinomial())
  if (any(!is.finite(stats::coef(fit)))) stop("non-finite survey coefficient")
  fit
}

calibration_fit <- function(data, prob, outcome = "dx", weight = "weight",
                            strata = "stra", psu = "psu", design = NULL) {
  p <- pmin(pmax(as.numeric(prob), 1e-8), 1 - 1e-8)
  d <- data; d$.cal_lp <- qlogis(p)
  des <- if (is.null(design)) make_survey_design(d, weight, strata, psu,
                                                 allow_zero_weights = TRUE) else update(design, .cal_lp = d$.cal_lp)
  fit <- tryCatch(survey::svyglm(stats::as.formula(paste(outcome, "~ .cal_lp")),
                                 design = des, family = stats::quasibinomial()), error = identity)
  if (inherits(fit, "error") || any(!is.finite(stats::coef(fit))))
    return(list(intercept = NA_real_, slope = NA_real_, intercept_se = NA_real_, slope_se = NA_real_,
                intercept_lower = NA_real_, intercept_upper = NA_real_,
                slope_lower = NA_real_, slope_upper = NA_real_, design_df = survey::degf(des), status = "failed",
                reason = if (inherits(fit, "error")) conditionMessage(fit) else "non-finite coefficient"))
  co <- stats::coef(fit)
  covariance <- stats::vcov(fit)
  df <- max(1, survey::degf(des))
  critical <- stats::qt(.975, df = df)
  se <- sqrt(diag(covariance))
  list(intercept = unname(co[[1]]), slope = unname(co[[2]]),
       intercept_se = unname(se[[1]]), slope_se = unname(se[[2]]),
       intercept_lower = unname(co[[1]] - critical * se[[1]]),
       intercept_upper = unname(co[[1]] + critical * se[[1]]),
       slope_lower = unname(co[[2]] - critical * se[[2]]),
       slope_upper = unname(co[[2]] + critical * se[[2]]),
       design_df = df, status = "valid", reason = "")
}

# A replicate estimate uses the same Rao-Wu replicate weights as the survey
# design point estimate.  The coefficient is the weighted quasi-likelihood
# solution; strata and PSU structure enter through the replicate weights, not
# through a pseudo-design that treats individual rows as PSUs.
calibration_fit_replicate <- function(data, prob, outcome = "dx", weight = "weight") {
  p <- pmin(pmax(as.numeric(prob), 1e-8), 1 - 1e-8)
  d <- data.frame(y = as.integer(data[[outcome]]), lp = qlogis(p),
                  w = as.numeric(data[[weight]]))
  ok <- is.finite(d$y) & is.finite(d$lp) & is.finite(d$w) & d$w > 0
  if (sum(ok) < 20L || length(unique(d$y[ok])) < 2L)
    return(list(intercept = NA_real_, slope = NA_real_, status = "failed",
                reason = "insufficient positive replicate weight"))
  lp0 <- d$lp[ok & d$y == 0L]
  lp1 <- d$lp[ok & d$y == 1L]
  separated <- max(lp0) <= min(lp1) || max(lp1) <= min(lp0)
  if (separated)
    return(list(intercept = NA_real_, slope = NA_real_, status = "failed",
                reason = "complete or quasi-complete separation"))
  fit_warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    stats::glm(y ~ lp, data = d[ok, , drop = FALSE],
               weights = w, family = stats::quasibinomial()),
    warning = function(warning_condition) {
      fit_warnings <<- c(fit_warnings, conditionMessage(warning_condition))
      invokeRestart("muffleWarning")
    }), error = identity)
  separation_warning <- any(grepl(
    "did not converge|fitted probabilities numerically 0 or 1 occurred",
    fit_warnings, ignore.case = TRUE))
  covariance <- if (inherits(fit, "error")) NULL else
    tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (inherits(fit, "error") || !isTRUE(fit$converged) ||
      isTRUE(fit$boundary) || separation_warning ||
      any(!is.finite(stats::coef(fit))) || is.null(covariance) ||
      any(!is.finite(covariance)))
    return(list(intercept = NA_real_, slope = NA_real_, status = "failed",
                reason = if (inherits(fit, "error")) conditionMessage(fit) else
                  "non-finite, non-converged or separation-affected replicate calibration"))
  co <- stats::coef(fit)
  list(intercept = unname(co[[1]]), slope = unname(co[[2]]), status = "valid", reason = "")
}

prediction_metrics <- function(data, lp, outcome = "dx", weight = "weight",
                               strata = "stra", psu = "psu", design = NULL) {
  y <- as.integer(data[[outcome]]); w <- as.numeric(data[[weight]]); pr <- stats::plogis(lp)
  cal <- calibration_fit(data, pr, outcome, weight, strata, psu, design = design)
  c(auc = weighted_auc(lp, y, w), brier = weighted_brier(pr, y, w),
    calibration_intercept = cal$intercept, calibration_slope = cal$slope)
}

calibration_uncertainty <- function(data, lp, outcome = "dx", weight = "weight",
                                    strata = "stra", psu = "psu", design = NULL,
                                    model = "model") {
  probability <- stats::plogis(as.numeric(lp))
  fit <- calibration_fit(data, probability, outcome, weight, strata, psu, design)
  data.frame(
    model = model,
    measure = c("calibration_intercept", "calibration_slope"),
    point = c(fit$intercept, fit$slope),
    se = c(fit$intercept_se, fit$slope_se),
    lower = c(fit$intercept_lower, fit$slope_lower),
    upper = c(fit$intercept_upper, fit$slope_upper),
    variance_method = "Taylor linearization",
    design_df = fit$design_df,
    status = fit$status,
    reason = fit$reason,
    stringsAsFactors = FALSE
  )
}

make_raowu_replicates <- function(design, replicates = 1000, seed = 2603L) {
  if (!inherits(design, "survey.design2")) stop("design must be a survey.design2 object")
  set.seed(seed)
  out <- survey::as.svrepdesign(design, type = "subbootstrap", replicates = replicates,
                                mse = TRUE, compress = TRUE)
  actual_replicates <- ncol(stats::weights(out, type = "analysis"))
  expected_scale <- 1 / (actual_replicates - 1)
  if (actual_replicates != replicates ||
      !isTRUE(all.equal(as.numeric(out$scale), expected_scale, tolerance = 1e-12)) ||
      !all(as.numeric(out$rscales) == 1) || !isTRUE(out$mse)) {
    stop("unexpected survey::subbootstrap scaling or replicate count")
  }
  out
}

raowu_specification <- function(design, repdesign) {
  strata <- as.character(design$strata[, 1])
  psu <- as.character(design$cluster[, 1])
  per_stratum <- tapply(psu, strata, function(x) length(unique(x)))
  list(
    implementation = "survey::as.svrepdesign(type='subbootstrap')",
    resampling = "within each stratum, sample n_h-1 PSUs with replacement",
    multiplier = "PSU multiplicity multiplied by n_h/(n_h-1)",
    replicates = ncol(weights(repdesign, type = "analysis")),
    scale = unname(repdesign$scale),
    expected_scale = 1 / (ncol(weights(repdesign, type = "analysis")) - 1),
    scale_matches_expected = isTRUE(all.equal(
      as.numeric(repdesign$scale),
      1 / (ncol(weights(repdesign, type = "analysis")) - 1),
      tolerance = 1e-12)),
    rscales_unique = unique(as.numeric(repdesign$rscales)),
    mse = isTRUE(repdesign$mse),
    lonely_psu_option = getOption("survey.lonely.psu"),
    lonely_strata = sum(per_stratum < 2L),
    min_psu_per_stratum = min(per_stratum),
    max_psu_per_stratum = max(per_stratum),
    interpretation = paste(
      "Prespecified design-based resampling approximation for paired fixed-prediction contrasts;",
      "not a replacement for NCHS-recommended Taylor linearization."
    )
  )
}

compare_taylor_raowu <- function(design, formula, term, replicates = 1000L, seed = 2807L) {
  fit_taylor <- fit_survey_logistic(formula = formula, design = design)
  repdesign <- make_raowu_replicates(design, replicates, seed)
  fit_rw <- survey::svyglm(formula, design = repdesign, family = stats::quasibinomial())
  if (!term %in% names(stats::coef(fit_taylor)) || !term %in% names(stats::coef(fit_rw)))
    stop("variance-method comparison term is absent")
  t_df <- max(1, survey::degf(design)); t_critical <- stats::qt(.975, t_df)
  rw_df <- max(1, survey::degf(repdesign)); rw_critical <- stats::qt(.975, rw_df)
  estimate_t <- unname(stats::coef(fit_taylor)[[term]])
  estimate_r <- unname(stats::coef(fit_rw)[[term]])
  se_t <- unname(sqrt(diag(stats::vcov(fit_taylor)))[[term]])
  se_r <- unname(sqrt(diag(stats::vcov(fit_rw)))[[term]])
  data.frame(
    method = c("Taylor linearization", "Rao-Wu rescaled subbootstrap"),
    term = term,
    estimate = c(estimate_t, estimate_r),
    se = c(se_t, se_r),
    lower = c(estimate_t - t_critical * se_t, estimate_r - rw_critical * se_r),
    upper = c(estimate_t + t_critical * se_t, estimate_r + rw_critical * se_r),
    design_df = c(t_df, rw_df),
    direction = sign(c(estimate_t, estimate_r)),
    se_ratio_to_taylor = c(1, se_r / se_t),
    stringsAsFactors = FALSE
  )
}

rao_wu_wald_test <- function(point, replicate_estimates, scale, rscales,
                             mse = TRUE, minimum_valid = 900L) {
  values <- as.numeric(replicate_estimates)
  ok <- is.finite(values)
  if (sum(ok) < minimum_valid || length(rscales) != length(values) || !is.finite(point)) {
    return(list(se = NA_real_, p_value = NA_real_, valid = sum(ok)))
  }
  variance <- tryCatch(as.numeric(survey::svrVar(
    values[ok], scale = scale, rscales = rscales[ok], mse = mse, coef = point
  )), error = function(e) NA_real_)
  if (!is.finite(variance) || variance < 0) {
    return(list(se = NA_real_, p_value = NA_real_, valid = sum(ok)))
  }
  se <- sqrt(variance)
  p_value <- if (se <= sqrt(.Machine$double.eps)) {
    if (abs(point) <= sqrt(.Machine$double.eps)) 1 else NA_real_
  } else {
    2 * stats::pnorm(-abs(point / se))
  }
  list(se = se, p_value = p_value, valid = sum(ok))
}

development_optimism_bootstrap <- function(data, base_formula, augmented_formula,
                                           B = 1000L, seed = 2605L,
                                           outcome = "dx", weight = "weight",
                                           design = NULL) {
  y <- as.integer(data[[outcome]]); original_weight <- as.numeric(data[[weight]])
  des <- if (is.null(design)) make_survey_design(data, weight = weight) else design
  if (nrow(des$variables) != nrow(data)) stop("optimism design and data rows differ")
  repdes <- make_raowu_replicates(des, B, seed)
  rw <- weights(repdes, type = "analysis")
  status_levels <- c("invalid_replicate_weights", "base_fit_error", "base_nonconvergence",
                     "base_separation_or_boundary", "base_rank_deficiency",
                     "augmented_fit_error", "augmented_nonconvergence",
                     "augmented_separation_or_boundary", "augmented_rank_deficiency",
                     "non_finite_predictions",
                     "non_finite_auc", "valid")
  status <- rep(NA_character_, ncol(rw)); values <- vector("list", ncol(rw))
  detail <- rep("", ncol(rw))
  positive_weight_n <- event_n <- active_psu_n <- integer(ncol(rw))
  weighted_event_rate <- numeric(ncol(rw))
  psu_draw_hash <- character(ncol(rw))
  psu_key <- interaction(data$stra, data$psu, drop = TRUE)
  converged_base <- converged_augmented <- logical(ncol(rw))
  fit_replicate <- function(formula, dd, label) {
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(
      stats::glm(formula, data = dd, weights = .replicate_weight,
                 family = stats::quasibinomial(), control = stats::glm.control(maxit = 50)),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = identity)
    if (inherits(fit, "error"))
      return(list(status = paste0(label, "_fit_error"), detail = conditionMessage(fit), fit = NULL))
    if (!isTRUE(fit$converged))
      return(list(status = paste0(label, "_nonconvergence"), detail = paste(warnings, collapse = "; "), fit = fit))
    separation <- isTRUE(fit$boundary) || any(grepl(
      "fitted probabilities numerically 0 or 1|separation", warnings, ignore.case = TRUE))
    if (separation)
      return(list(status = paste0(label, "_separation_or_boundary"), detail = paste(warnings, collapse = "; "), fit = fit))
    if (fit$rank < length(stats::coef(fit)) || any(!is.finite(stats::coef(fit))))
      return(list(status = paste0(label, "_rank_deficiency"), detail = "rank deficient or non-finite coefficient", fit = fit))
    list(status = "valid", detail = paste(warnings, collapse = "; "), fit = fit)
  }
  for (b in seq_len(ncol(rw))) {
    wb <- as.numeric(rw[, b])
    positive_weight_n[b] <- sum(wb > 0)
    event_n[b] <- sum(y == 1L & wb > 0)
    active_psu_n[b] <- length(unique(psu_key[wb > 0]))
    weighted_event_rate[b] <- if (sum(wb) > 0) sum(wb * y) / sum(wb) else NA_real_
    psu_totals <- tapply(wb, psu_key, sum)
    psu_draw_hash[b] <- digest::digest(round(psu_totals[psu_totals > 0], 10), algo = "sha256")
    if (length(wb) != nrow(data) || any(!is.finite(wb)) || any(wb < 0) ||
        sum(wb) <= 0 || sum(wb[y == 0]) <= 0 || sum(wb[y == 1]) <= 0) {
      status[b] <- "invalid_replicate_weights"; detail[b] <- "invalid or outcome-empty replicate weights"; next
    }
    dd <- data; dd$.replicate_weight <- wb
    bfit <- fit_replicate(base_formula, dd, "base")
    converged_base[b] <- identical(bfit$status, "valid")
    if (!converged_base[b]) { status[b] <- bfit$status; detail[b] <- bfit$detail; next }
    f0 <- bfit$fit
    afit <- fit_replicate(augmented_formula, dd, "augmented")
    converged_augmented[b] <- identical(afit$status, "valid")
    if (!converged_augmented[b]) { status[b] <- afit$status; detail[b] <- afit$detail; next }
    f1 <- afit$fit
    lp0 <- tryCatch(as.numeric(stats::predict(f0, newdata = data, type = "link")), error = function(e) rep(NA_real_, nrow(data)))
    lp1 <- tryCatch(as.numeric(stats::predict(f1, newdata = data, type = "link")), error = function(e) rep(NA_real_, nrow(data)))
    if (any(!is.finite(lp0)) || any(!is.finite(lp1))) { status[b] <- "non_finite_predictions"; detail[b] <- "prediction was non-finite"; next }
    apparent0 <- weighted_auc(lp0, y, wb); apparent1 <- weighted_auc(lp1, y, wb)
    test0 <- weighted_auc(lp0, y, original_weight); test1 <- weighted_auc(lp1, y, original_weight)
    if (any(!is.finite(c(apparent0, apparent1, test0, test1)))) { status[b] <- "non_finite_auc"; detail[b] <- "AUC was non-finite"; next }
    values[[b]] <- c(base = apparent0 - test0, augmented = apparent1 - test1,
                     delta = (apparent1 - apparent0) - (test1 - test0))
    status[b] <- "valid"
  }
  if (anyNA(status)) stop("optimism bootstrap status was not assigned for every replicate")
  counts <- table(factor(status, levels = status_levels)); valid <- sum(status == "valid")
  enough <- valid >= 0.9 * B
  optimism <- if (valid) colMeans(do.call(rbind, values[status == "valid"])) else c(base = NA, augmented = NA, delta = NA)
  f0 <- fit_survey_logistic(formula = base_formula, design = des)
  f1 <- fit_survey_logistic(formula = augmented_formula, design = des)
  a0 <- prediction_metrics(data, predict(f0, type = "link"), outcome = outcome, weight = weight)[["auc"]]
  a1 <- prediction_metrics(data, predict(f1, type = "link"), outcome = outcome, weight = weight)[["auc"]]
  list(requested = B, valid = valid, converged_base = sum(converged_base),
       converged_augmented = sum(converged_augmented),
       jointly_converged = sum(converged_base & converged_augmented),
       failure_counts = as.list(stats::setNames(as.integer(counts), names(counts))),
       replicate_status = data.frame(replicate = seq_along(status), status = status,
                                     detail = detail, positive_weight_n = positive_weight_n,
                                      event_n = event_n, weighted_event_rate = weighted_event_rate,
                                      active_psu_n = active_psu_n,
                                      psu_draw_hash = psu_draw_hash,
                                      same_base_augmented_draw = TRUE,
                                     stringsAsFactors = FALSE),
       stopped = !enough,
       optimism_base = if (enough) unname(optimism[["base"]]) else NA_real_,
       optimism_augmented = if (enough) unname(optimism[["augmented"]]) else NA_real_,
       optimism_delta = if (enough) unname(optimism[["delta"]]) else NA_real_,
       corrected_auc_base = if (enough) a0 - optimism[["base"]] else NA_real_,
       corrected_auc_augmented = if (enough) a1 - optimism[["augmented"]] else NA_real_,
       corrected_delta_auc = if (enough) (a1 - a0) - optimism[["delta"]] else NA_real_,
       raowu_specification = raowu_specification(des, repdes),
       algorithm = "Rao-Wu replicate weights; the same draw is used for base and augmented fits; bootstrap apparent minus original-weight test AUC; output gate 900/1000; failures classified before summarization")
}

summarize_optimism_stability <- function(trace) {
  required <- c("status", "positive_weight_n", "event_n", "weighted_event_rate",
                "active_psu_n", "psu_draw_hash")
  if (!all(required %in% names(trace))) stop("optimism trace lacks stability diagnostics")
  groups <- split(trace, trace$status)
  do.call(rbind, lapply(groups, function(x) {
    pattern_counts <- table(x$psu_draw_hash)
    data.frame(
      status = x$status[[1]],
      replicates = nrow(x),
      event_n_median = stats::median(x$event_n),
      event_n_q1 = unname(stats::quantile(x$event_n, .25)),
      event_n_q3 = unname(stats::quantile(x$event_n, .75)),
      weighted_event_rate_median = stats::median(x$weighted_event_rate, na.rm = TRUE),
      weighted_event_rate_q1 = unname(stats::quantile(x$weighted_event_rate, .25, na.rm = TRUE)),
      weighted_event_rate_q3 = unname(stats::quantile(x$weighted_event_rate, .75, na.rm = TRUE)),
      active_psu_n_min = min(x$active_psu_n),
      active_psu_n_max = max(x$active_psu_n),
      unique_psu_draw_patterns = length(pattern_counts),
      maximum_repeated_pattern_count = max(pattern_counts),
      stringsAsFactors = FALSE
    )
  }))
}

replicate_fixed_prediction_metrics <- function(validation_data, lp_base, lp_augmented,
                                                B = 1000, seed = 2603L,
                                                outcome = "dx", weight = "weight",
                                                strata = "stra", psu = "psu",
                                                minimum_valid = 900L,
                                                validation_design = NULL) {
  stopifnot(length(lp_base) == nrow(validation_data), length(lp_augmented) == nrow(validation_data))
  d <- validation_data; d$.lp0 <- lp_base; d$.lp1 <- lp_augmented
  des <- if (is.null(validation_design)) make_survey_design(d, weight, strata, psu) else validation_design
  if (nrow(des$variables) != nrow(d)) stop("validation design and data rows differ")
  repdes <- make_raowu_replicates(des, B, seed)
  rw <- weights(repdes, type = "analysis")
  y <- as.integer(d[[outcome]])
  one <- function(w) {
    p0 <- plogis(d$.lp0); p1 <- plogis(d$.lp1)
    # Replicate calibration uses the replicate weights inside the same original
    # stratum/PSU structure; each row is not treated as an independent PSU.
    dd <- d; dd$.rw <- as.numeric(w)
    c0 <- calibration_fit_replicate(dd, p0, outcome, ".rw")
    c1 <- calibration_fit_replicate(dd, p1, outcome, ".rw")
    c(auc_base = weighted_auc(d$.lp0, y, w), auc_augmented = weighted_auc(d$.lp1, y, w),
      delta_auc = weighted_auc(d$.lp1, y, w) - weighted_auc(d$.lp0, y, w),
      brier_base = weighted_brier(p0, y, w), brier_augmented = weighted_brier(p1, y, w),
      delta_brier = weighted_brier(p1, y, w) - weighted_brier(p0, y, w),
      calibration_intercept_base = c0$intercept, calibration_intercept_augmented = c1$intercept,
      delta_calibration_intercept = c1$intercept - c0$intercept,
      calibration_slope_base = c0$slope, calibration_slope_augmented = c1$slope,
      delta_calibration_slope = c1$slope - c0$slope)
  }
  vals <- t(apply(rw, 2, one))
  valid_by_measure <- colSums(is.finite(vals))
  enough <- valid_by_measure >= minimum_valid
  intervals <- matrix(NA_real_, nrow = ncol(vals), ncol = 2L,
                      dimnames = list(colnames(vals), c("lower", "upper")))
  single_calibration <- colnames(vals) %in% c(
    "calibration_intercept_base", "calibration_intercept_augmented",
    "calibration_slope_base", "calibration_slope_augmented")
  for (j in seq_len(ncol(vals))) {
    ok <- is.finite(vals[, j])
    if (enough[j] && !single_calibration[j]) {
      bounds <- stats::quantile(vals[ok, j], probs = c(.025, .975), names = FALSE)
      intervals[j, ] <- bounds
    }
  }
  reportable <- enough & !single_calibration
  reasons <- ifelse(!enough, "not estimable under the 900-of-1000 reporting rule",
                    ifelse(single_calibration,
                      "diagnostic only; Taylor linearization is primary for single-model calibration", ""))
  thresholds <- do.call(rbind, lapply(c(800L, 900L, 950L), function(threshold) {
    do.call(rbind, lapply(seq_len(ncol(vals)), function(j) {
      ok <- is.finite(vals[, j]); report <- sum(ok) >= threshold
      bounds <- if (report) stats::quantile(vals[ok, j], c(.025, .975), names = FALSE) else c(NA_real_, NA_real_)
      data.frame(measure = colnames(vals)[j], threshold = threshold,
                 valid = sum(ok), reportable = report, lower = bounds[[1]], upper = bounds[[2]],
                 stringsAsFactors = FALSE)
    }))
  }))
  list(status = if (all(reportable)) "valid" else if (any(reportable)) "partial" else "stopped",
       requested = B, valid_by_measure = valid_by_measure,
       stopped_by_measure = !reportable, failure_by_measure = reasons,
       intervals = intervals, threshold_sensitivity = thresholds, replicates = vals,
       raowu_specification = raowu_specification(des, repdes),
       reason = if (all(reportable)) "" else "one or more metrics failed the replicate-count reporting rule")
}

psu_grouped_cv <- function(data, base_formula, augmented_formula, index,
                           folds = 5L, repeats = 5L, seed = 2601L,
                           outcome = "dx", weight = "weight", strata = "stra", psu = "psu",
                           structural_zero = FALSE, design = NULL) {
  if (!is.null(design) && nrow(design$variables) != nrow(data))
    stop("cross-validation design and data rows differ")
  psu_key <- interaction(data[[strata]], data[[psu]], drop = TRUE)
  psu_levels <- levels(psu_key); out <- list(); row_id <- 1L
  for (r in seq_len(repeats)) {
    set.seed(seed + r - 1L)
    fold_id <- sample(rep(seq_len(folds), length.out = length(psu_levels)))
    names(fold_id) <- psu_levels
    for (k in seq_len(folds)) {
      valid <- fold_id[as.character(psu_key)] == k; train <- !valid
      overlap <- intersect(unique(psu_key[train]), unique(psu_key[valid]))
      status <- "valid"; reason <- ""; mb <- ma <- NULL
      if (length(overlap)) { status <- "failed"; reason <- "PSU leakage" }
      if (status == "valid") {
        train_design <- if (is.null(design)) NULL else {
          z <- design
          z$variables$.cv_train <- train
          subset(z, .cv_train)
        }
        mb <- tryCatch(fit_survey_logistic(data[train, , drop = FALSE], base_formula,
                                           weight, strata, psu, design = train_design),
                       error = identity)
        ma <- if (structural_zero && !inherits(mb, "error")) mb else
          tryCatch(fit_survey_logistic(data[train, , drop = FALSE], augmented_formula,
                                       weight, strata, psu, design = train_design),
                   error = identity)
        if (inherits(mb, "error") || inherits(ma, "error")) {
          status <- "failed"; reason <- paste(c(if (inherits(mb, "error")) conditionMessage(mb),
                                                if (inherits(ma, "error")) conditionMessage(ma)), collapse = "; ")
        }
      }
      met <- rep(NA_real_, 10)
      if (status == "valid") {
        valid_design <- if (is.null(design)) NULL else {
          z <- design
          z$variables$.cv_valid <- valid
          subset(z, .cv_valid)
        }
        lp0 <- as.numeric(stats::predict(mb, newdata = data[valid, , drop = FALSE], type = "link"))
        lp1 <- as.numeric(stats::predict(ma, newdata = data[valid, , drop = FALSE], type = "link"))
        m0 <- prediction_metrics(data[valid, , drop = FALSE], lp0, outcome, weight,
                                 strata, psu, design = valid_design)
        m1 <- prediction_metrics(data[valid, , drop = FALSE], lp1, outcome, weight,
                                 strata, psu, design = valid_design)
        met <- c(m0[c("auc","brier","calibration_intercept","calibration_slope")],
                 m1[c("auc","brier","calibration_intercept","calibration_slope")],
                 delta_auc = m1[["auc"]] - m0[["auc"]], delta_brier = m1[["brier"]] - m0[["brier"]])
      }
      out[[row_id]] <- data.frame(index = index, repetition = r, fold = k, seed = seed + r - 1L,
        n_train = sum(train), n_valid = sum(valid), cases_train = sum(data[[outcome]][train] == 1, na.rm = TRUE),
        cases_valid = sum(data[[outcome]][valid] == 1, na.rm = TRUE),
        psu_train = length(unique(psu_key[train])), psu_valid = length(unique(psu_key[valid])),
        psu_overlap = length(overlap), status = status, failure_reason = reason,
        auc_base = met[1], brier_base = met[2], calibration_intercept_base = met[3],
        calibration_slope_base = met[4], auc_augmented = met[5], brier_augmented = met[6],
        calibration_intercept_augmented = met[7], calibration_slope_augmented = met[8],
        delta_auc = met[9], delta_brier = met[10], stringsAsFactors = FALSE)
      row_id <- row_id + 1L
    }
  }
  do.call(rbind, out)
}

evaluate_incremental_model <- function(development_data, validation_data,
                                       base_formula, augmented_formula, index,
                                       B = 1000L, seed_cv = 2601L, seed_bootstrap = 2603L,
                                       weight = "weight", strata = "stra", psu = "psu", outcome = "dx",
                                       structural_zero = FALSE,
                                       development_design = NULL, validation_design = NULL) {
  base_fit <- fit_survey_logistic(development_data, base_formula, weight, strata, psu,
                                  design = development_design)
  aug_fit <- if (structural_zero) base_fit else
    fit_survey_logistic(development_data, augmented_formula, weight, strata, psu,
                        design = development_design)
  cv <- psu_grouped_cv(development_data, base_formula, augmented_formula, index,
                       folds = 5L, repeats = 5L, seed = seed_cv, outcome = outcome,
                       weight = weight, strata = strata, psu = psu,
                       structural_zero = structural_zero, design = development_design)
  lp0 <- as.numeric(stats::predict(base_fit, newdata = validation_data, type = "link"))
  lp1 <- as.numeric(stats::predict(aug_fit, newdata = validation_data, type = "link"))
  point0 <- prediction_metrics(validation_data, lp0, outcome, weight, strata, psu,
                               design = validation_design)
  point1 <- prediction_metrics(validation_data, lp1, outcome, weight, strata, psu,
                               design = validation_design)
  intervals <- replicate_fixed_prediction_metrics(validation_data, lp0, lp1, B, seed_bootstrap,
                                                    outcome, weight, strata, psu,
                                                    validation_design = validation_design)
  list(index = index, fits = list(base = base_fit, augmented = aug_fit), candidate_cv = cv,
       validation = list(base = point0, augmented = point1,
                         delta_auc = point1[["auc"]] - point0[["auc"]],
                         delta_brier = point1[["brier"]] - point0[["brier"]],
                         delta_calibration_intercept = point1[["calibration_intercept"]] - point0[["calibration_intercept"]],
                         delta_calibration_slope = point1[["calibration_slope"]] - point0[["calibration_slope"]]),
       intervals = intervals)
}

run_glm7_example <- function(development_data, validation_data, ...) {
  evaluate_incremental_model(development_data, validation_data,
    dx ~ sex + log_age + log_bmi,
    dx ~ sex + log_age + log_bmi + GLM7,
    index = "GLM7", ...)
}

# ---------------------------------------------------------------------------
# v34: decision consequences of a mishandled unit conversion.
# lp_correct / lp_bad are two linear predictors for the SAME rows (here the
# later-period predictions before and after an un-mapped unit shift).  At
# each risk threshold t we report, on the survey-weighted population and on
# the raw rows:
#   flagged        predicted risk >= t (referred / treated as high risk)
#   cases_missed   true cases with predicted risk < t
# The AUC is identical for the two predictors; the decision counts are not.
# ---------------------------------------------------------------------------
threshold_shift_table <- function(data, lp_correct, lp_bad, outcome = "dx",
                                  weight = "weight",
                                  thresholds = c(0.05, 0.10, 0.15, 0.20)) {
  stopifnot(length(lp_correct) == nrow(data), length(lp_bad) == nrow(data))
  y <- as.integer(data[[outcome]]); w <- as.numeric(data[[weight]])
  pc <- stats::plogis(as.numeric(lp_correct))
  pb <- stats::plogis(as.numeric(lp_bad))
  ok <- is.finite(pc) & is.finite(pb) & !is.na(y) & is.finite(w) & w > 0
  y <- y[ok]; w <- w[ok]; pc <- pc[ok]; pb <- pb[ok]
  rows <- lapply(thresholds, function(t) {
    fc <- pc >= t; fb <- pb >= t
    data.frame(
      threshold = t,
      n_flagged_correct = sum(fc), n_flagged_bad = sum(fb),
      w_flagged_correct = sum(w[fc]), w_flagged_bad = sum(w[fb]),
      delta_w_flagged = sum(w[fb]) - sum(w[fc]),
      cases_missed_correct = sum(y == 1 & !fc),
      cases_missed_bad = sum(y == 1 & !fb),
      w_cases_missed_correct = sum(w[y == 1 & !fc]),
      w_cases_missed_bad = sum(w[y == 1 & !fb]),
      ppv_correct = if (any(fc)) sum(w[fc & y == 1]) / sum(w[fc]) else NA_real_,
      ppv_bad = if (any(fb)) sum(w[fb & y == 1]) / sum(w[fb]) else NA_real_,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$delta_cases_missed <- out$cases_missed_bad - out$cases_missed_correct
  out$delta_w_cases_missed <- out$w_cases_missed_bad - out$w_cases_missed_correct
  out
}
