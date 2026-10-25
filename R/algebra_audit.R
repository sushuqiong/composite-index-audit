# v32 algebraic audit. The symbolic identity is primary; numerical diagnostics
# are recorded separately and never promoted to proof of an identity.
#
# v32 vocabulary (used by the manuscript, the registry and this file):
#   rank_preserving_exact      C lies in span[Z,E]; the augmented parameter
#                              vector is not uniquely identified and the
#                              increment is zero for any data.
#   rank_adding_exact          C = E + R with R outside span[Z,E]; the map
#                              between the two representations is invertible
#                              and the increment is data-dependent and must be
#                              estimated.
#   reducible_to_interaction   C is not additively decomposable on its native
#                              scale, but expands (distributive law) into a
#                              linear combination of columns already in the
#                              model plus products of pairs of model columns.
#                              Augmenting the model with C adds exactly one
#                              dimension: a fixed-ratio bundle of those
#                              interaction columns. The unrestricted
#                              interaction model y ~ W + I nests the augmented
#                              model y ~ W + C, and the coefficient map is
#                              closed-form.
#   outside_supported_class    C is outside both the additive log-linear class
#                              and the product-expandable class. The audit
#                              reports the boundary honestly and makes no
#                              structural claim.
# The legacy labels exact_reparameterization / not_supported and the v31 label
# not_decomposable are kept as secondary fields so earlier callers keep working.


#' Symbolic decomposition of a composite index against a base model
#'
#' Writes \code{C = E + R}, where \code{E} is the part of the index that the
#' base model already spans and \code{R} is the residual, and checks whether
#' \code{R} adds a dimension.  The numerical rank check verifies the symbolic
#' claim; it never proves it.
#'
#' @param index_terms Named numeric vector: the index as a linear form.
#' @param embedded_terms Named numeric vector: the part present in the base.
#' @param model_terms Named numeric vector: the base model columns.
#' @param data Optional data frame for the numerical layer.
#' @param tol Numerical tolerance.
#' @param approximate_correlation Correlation above which a non-exact match is
#'   flagged as approximate collinearity rather than as a class.
#' @return A list with \code{status}, \code{exact_reparameterization},
#'   \code{residual_terms}, diagnostics and \code{warnings}.
#' @export
#' @examples
#' audit_composite_index(c(a = 1, b = 1), c(a = 1, b = 1))$status
audit_composite_index <- function(index_terms, embedded_terms,
                                   model_terms = embedded_terms, data = NULL,
                                   tol = 1e-10, approximate_correlation = 0.90) {
  stopifnot(is.numeric(index_terms), !is.null(names(index_terms)),
            is.numeric(embedded_terms), !is.null(names(embedded_terms)))
  if (anyDuplicated(names(index_terms)) || anyDuplicated(names(embedded_terms)))
    stop("term names must be unique")
  all_terms <- union(names(index_terms), names(embedded_terms))
  b <- setNames(numeric(length(all_terms)), all_terms); b[names(index_terms)] <- index_terms
  e <- setNames(numeric(length(all_terms)), all_terms); e[names(embedded_terms)] <- embedded_terms
  residual <- b - e
  embedded_match <- all(names(embedded_terms) %in% names(index_terms)) &&
    all(abs(index_terms[names(embedded_terms)] - embedded_terms) <= tol)
  model_match <- setequal(names(model_terms), names(embedded_terms)) &&
    all(abs(model_terms[names(embedded_terms)] - embedded_terms) <= tol)
  exact <- embedded_match && model_match && all(is.finite(residual))
  corr <- NA_real_; rank <- NA_integer_; condition_number <- NA_real_; zero <- character()
  rank_base <- NA_integer_; rank_augmented <- NA_integer_; rank_gain <- NA_integer_
  if (!is.null(data)) {
    missing_terms <- setdiff(all_terms, names(data))
    if (length(missing_terms)) stop("data is missing transformed terms: ", paste(missing_terms, collapse = ", "))
    X <- as.matrix(data[, all_terms, drop = FALSE])
    rank <- qr(X, tol = tol)$rank
    # v31: compare the base space [Z,E] with the augmented space [Z,E,C] on the
    # same rows. A gain of 0 means C lies in the span of the base (structural
    # zero); a gain of 1 means the residual adds one dimension.
    iv_all <- as.numeric(X[, names(index_terms), drop = FALSE] %*% index_terms)
    model_cols <- intersect(names(model_terms), colnames(X))
    if (length(model_cols)) {
      X_base <- cbind(1, as.matrix(data[, model_cols, drop = FALSE]))
      keep <- stats::complete.cases(X_base) & is.finite(iv_all)
      if (sum(keep) > ncol(X_base) + 1L) {
        rank_base <- qr(X_base[keep, , drop = FALSE], tol = tol)$rank
        rank_augmented <- qr(cbind(X_base, iv_all)[keep, , drop = FALSE], tol = tol)$rank
        rank_gain <- rank_augmented - rank_base
      }
    }
    sds <- apply(X, 2, stats::sd, na.rm = TRUE)
    usable <- is.finite(sds) & sds > tol
    zero <- colnames(X)[!usable]
    if (sum(usable) >= 2) condition_number <- tryCatch(
      kappa(scale(X[, usable, drop = FALSE], center = TRUE, scale = sds[usable]), exact = TRUE),
      error = function(e) NA_real_)
    iv <- as.numeric(X[, names(index_terms), drop = FALSE] %*% index_terms)
    ev <- as.numeric(X[, names(embedded_terms), drop = FALSE] %*% embedded_terms)
    corr <- suppressWarnings(cor(iv, ev, use = "complete.obs"))
  }
  approximate <- !exact && is.finite(corr) && abs(corr) >= approximate_correlation
  # v31 classification. The symbolic decomposition says whether C = E + R can be
  # written at all; the numerical rank gain says whether R adds a dimension. The
  # numeric check is a verification of the symbolic claim, never its proof.
  if (!exact) {
    status <- if (approximate) "approximate_collinearity" else "not_decomposable"
  } else if (is.na(rank_gain)) {
    # v34 fix: residual is a fixed-length named vector over all_terms, so its
    # length is never zero. The symbolic question is whether every entry is
    # zero within tolerance (C = E). This makes the no-data branch
    # independent of the sample rank: audit_index_in_model(c(x=1,y=1),
    # c(x=1,y=1), data = NULL) now returns rank_preserving_exact.
    status <- if (all(abs(residual) <= tol)) "rank_preserving_exact" else "rank_adding_exact"
  } else {
    status <- if (rank_gain == 0L) "rank_preserving_exact" else "rank_adding_exact"
  }
  warnings <- character()
  if (!embedded_match) warnings <- c(warnings, "embedded terms or coefficients do not match the index representation")
  if (!model_match) warnings <- c(warnings, "model terms do not match the embedded representation")
  if (approximate) warnings <- c(warnings, "empirical correlation is not an algebraic identity")
  if (length(zero)) warnings <- c(warnings, paste0("zero-variance columns excluded from condition number: ", paste(zero, collapse = ", ")))
  list(status = status, exact_reparameterization = exact,
       decomposable = exact, approximate_collinearity = approximate,
       index_terms = b, embedded_terms = e, residual_terms = residual[abs(residual) > tol],
       model_terms = model_terms, tolerance = tol,
       empirical_index_embedded_correlation = corr, approximate_correlation_threshold = approximate_correlation,
       diagnostics = list(rank = rank, rank_base = rank_base, rank_augmented = rank_augmented,
                          rank_gain = rank_gain, condition_number = condition_number,
                          zero_variance_columns = zero,
                          max_linear_predictor_difference = NA_real_, max_fitted_probability_difference = NA_real_,
                          max_coefficient_difference = NA_real_), warnings = warnings)
}

# ---------------------------------------------------------------------------
# v31 public entry point. Four-field output:
#   classification        rank_preserving_exact | rank_adding_exact |
#                         not_decomposable (plus approximate_collinearity)
#   coefficient_map       the matrix K with gamma = K beta, or NULL together
#                         with blocking_reason when no map exists
#   numerical_diagnostics ranks, rank gain, condition number, correlation
#   failure_states        every state raised by the audit; none downgraded
#
# With base columns [Z,E] and composite C = sum(b_k t_k) where the embedded
# part is E = sum(e_k t_k) and R = C - E, the two augmented representations
# [Z,E,C] and [Z,E,R] are linked by
#   gamma_Z = beta_Z,  gamma_tk = beta_tk + e_k beta_C,  gamma_R = beta_C.
# ---------------------------------------------------------------------------
# v34: intercept is the additive constant of the registered index (for
# example HOMA-IR carries -log(22.5) once on the log scale). It maps to the
# model intercept: gamma_0 = beta_0 + intercept * beta_C.

#' Audit an additive (log-)linear index in a stated base model
#'
#' @param index_terms Named numeric vector giving the index as a linear form.
#' @param embedded_terms Named numeric vector, the embedded part.
#' @param model_terms Named numeric vector, the base model columns.
#' @param data Optional data frame for the numerical verification layer.
#' @param tol Numerical tolerance.
#' @param approximate_correlation Correlation threshold for the approximate
#'   collinearity flag.
#' @param intercept Additive constant of the index; it maps to the model
#'   intercept as \code{gamma_0 = beta_0 + intercept * beta_C}.
#' @return A list with \code{classification}, \code{coefficient_map},
#'   \code{blocking_reason}, \code{residual_expression},
#'   \code{numerical_diagnostics} and \code{failure_states}.
#' @export
#' @examples
#' audit_index_in_model(c(a = 1, b = 1), c(a = 1, b = 1))$classification
audit_index_in_model <- function(index_terms, embedded_terms, model_terms = embedded_terms,
                                 data = NULL, tol = 1e-10,
                                 approximate_correlation = 0.90,
                                 intercept = 0) {
  base <- audit_composite_index(index_terms, embedded_terms, model_terms, data = data,
                                tol = tol, approximate_correlation = approximate_correlation)
  failure_states <- character()
  if (length(base$warnings)) failure_states <- c(failure_states, base$warnings)
  if (!base$decomposable) failure_states <- c(failure_states, "not_decomposable")
  classification <- base$status
  coefficient_map <- NULL; blocking_reason <- NULL
  if (base$decomposable) {
    terms <- union(names(index_terms), names(embedded_terms))
    e_idx <- setNames(numeric(length(terms)), terms)
    e_idx[names(embedded_terms)] <- embedded_terms
    other <- setdiff(names(model_terms), names(embedded_terms))
    # columns of beta: intercept, other model terms, embedded terms, C
    b_names <- c("(Intercept)", other, names(embedded_terms), "C")
    # columns of gamma: intercept, other model terms, embedded terms, R
    g_names <- c("(Intercept)", other, names(embedded_terms), "R")
    K <- diag(length(b_names))
    rownames(K) <- g_names; colnames(K) <- b_names
    for (tk in names(embedded_terms)) K[tk, "C"] <- embedded_terms[[tk]]
    if (intercept != 0) K["(Intercept)", "C"] <- intercept
    coefficient_map <- K
    blocking_reason <- NA_character_
  } else {
    blocking_reason <- if (length(base$residual_terms) == 0L && !base$decomposable)
      "no embedded term of the index appears in the model with the same transformation, unit and scale"
    else "residual is outside the registered additive log-linear expression class"
  }
  list(classification = classification,
       coefficient_map = coefficient_map,
       blocking_reason = blocking_reason,
       residual_expression = if (base$decomposable)
         paste(paste(format(base$residual_terms, digits = 6), names(base$residual_terms), sep = " x "), collapse = " + ")
       else NA_character_,
       numerical_diagnostics = list(
         rank_base = base$diagnostics$rank_base,
         rank_augmented = base$diagnostics$rank_augmented,
         rank_gain = base$diagnostics$rank_gain,
         condition_number = base$diagnostics$condition_number,
         empirical_correlation = base$empirical_index_embedded_correlation,
         tolerance = tol),
       failure_states = failure_states,
       symbolic_layer = "the decomposition C = E + R and the map K follow from the registered definitions",
       numerical_layer = "ranks are computed on the supplied rows and verify, but do not prove, the symbolic claim")
}

design_space_audit <- function(data, tol = 1e-10) {
  required <- c("sex", "log_age", "log_bmi", "GLM7", "metab", "dx")
  if (!all(required %in% names(data))) stop("missing GLM7 design-space columns")
  xb <- stats::model.matrix(~ sex + log_age + log_bmi, data = data)
  xa <- stats::model.matrix(~ sex + log_age + log_bmi + GLM7, data = data)
  xr <- stats::model.matrix(~ sex + log_age + log_bmi + metab, data = data)
  rank_of <- function(x) qr(x, tol = tol)$rank
  fit_check <- tryCatch({
    m1 <- stats::glm(dx ~ sex + log_age + log_bmi + GLM7, data = data, family = stats::binomial())
    m2 <- stats::glm(dx ~ sex + log_age + log_bmi + metab, data = data, family = stats::binomial())
    max(abs(stats::predict(m1, type = "link") - stats::predict(m2, type = "link")), na.rm = TRUE)
  }, error = function(e) NA_real_)
  standardized_kappa <- function(x) {
    z <- x[, colnames(x) != "(Intercept)", drop = FALSE]
    keep <- apply(z, 2, function(v) is.finite(stats::sd(v)) && stats::sd(v) > tol)
    if (!any(keep)) return(NA_real_)
    tryCatch(kappa(scale(z[, keep, drop = FALSE]), exact = TRUE), error = function(e) NA_real_)
  }
  list(rank_base = rank_of(xb), rank_augmented = rank_of(xa), rank_residual = rank_of(xr),
       ncol_base = ncol(xb), ncol_augmented = ncol(xa), ncol_residual = ncol(xr),
       augmented_minus_base = rank_of(xa) - rank_of(xb), residual_minus_base = rank_of(xr) - rank_of(xb),
       condition_number_base = standardized_kappa(xb),
       condition_number_augmented = standardized_kappa(xa),
       condition_number_residual = standardized_kappa(xr),
       condition_number_definition = "SVD 2-norm condition number after removing the intercept and centering/scaling every non-zero-variance column; diagnostic only",
       tolerance = tol, max_linear_predictor_difference = fit_check,
       columns_base = colnames(xb), columns_augmented = colnames(xa), columns_residual = colnames(xr),
       interpretation = "[Z,E,C] and [Z,E,R] are reparameterizations; [Z,E,C] equals [Z,E] alone only when R lies in span([Z,E])")
}

run_diagnostic_benchmark <- function(data, later_delta_auc, audit_status,
                                     tol = 1e-10) {
  required <- c("sex", "log_age", "log_bmi", "GLM7")
  if (!all(required %in% names(data))) stop("benchmark columns are missing")
  embedded_sum <- data$log_age + data$log_bmi
  pearson <- stats::cor(data$GLM7, embedded_sum, method = "pearson", use = "complete.obs")
  spearman <- stats::cor(data$GLM7, embedded_sum, method = "spearman", use = "complete.obs")
  x <- stats::model.matrix(~ sex + log_age + log_bmi + GLM7, data = data)
  z <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  vif <- vapply(seq_len(ncol(z)), function(j) {
    other <- z[, -j, drop = FALSE]
    if (!ncol(other)) return(1)
    fit <- stats::lm.fit(cbind(1, other), z[, j])
    rss <- sum(fit$residuals^2); tss <- sum((z[, j] - mean(z[, j]))^2)
    r2 <- if (tss > tol) 1 - rss/tss else NA_real_
    if (is.finite(r2) && r2 < 1) 1/(1-r2) else Inf
  }, numeric(1))
  condition_number <- design_space_audit(data, tol)$condition_number_augmented
  data.frame(
    method = c("Pearson/Spearman correlation", "VIF/condition number",
               "Later-period AUC-only comparison", "Algebraic-audit workflow"),
    observed_result = c(
      sprintf("Pearson r=%.3f; Spearman rho=%.3f", pearson, spearman),
      sprintf("max VIF=%.3f; standardized SVD kappa=%.3f", max(vif), condition_number),
      sprintf("paired delta AUC=%.4f", later_delta_auc),
      audit_status),
    identifies_exact_reparameterization = c(FALSE, FALSE, FALSE, TRUE),
    interpretation = c(
      "Association only; cannot establish an identity or coefficient map.",
      "Numerical collinearity diagnostic only; not a symbolic test.",
      "Performance contrast only; does not establish model-space equality.",
      "Checks expression, rank, coefficient/covariance map and fitted predictions."),
    stringsAsFactors = FALSE
  )
}

audit_glm7 <- function(data, tol = 1e-10) {
  required <- c("log_age", "log_bmi", "metab", "GLM7")
  if (!all(required %in% names(data))) stop("missing GLM7 audit columns")
  identity_error <- max(abs(data$GLM7 - (data$metab + data$log_age + data$log_bmi)), na.rm = TRUE)
  if (!is.finite(identity_error) || identity_error > tol)
    stop("GLM7 decomposition exceeds tolerance")
  # GLM7 = metab + log10(age) + log10(BMI). Against the base
  # sex + log10(age) + log10(BMI) the residual is metab, which is outside the
  # span of that base, so this record is rank_adding_exact.
  list(
    status = "rank_adding_exact",
    exact_reparameterization = TRUE,
    formula = "GLM7 = metab + log10(age) + log10(BMI)",
    embedded_terms = c(log_age = 1, log_bmi = 1),
    residual_terms = c(metab = 1),
    tolerance = tol,
    identity_error = identity_error,
    diagnostics = list(identity_error = identity_error)
  )
}

check_reparameterization <- function(data, outcome, candidate, residual,
                                     embedded_terms, tol = 1e-10) {
  req <- c(outcome, candidate, residual, embedded_terms)
  if (!all(req %in% names(data))) stop("missing reparameterization columns")
  d <- data[stats::complete.cases(data[, req, drop = FALSE]), , drop = FALSE]
  if (!nrow(d)) stop("no complete rows for reparameterization check")
  base_rhs <- paste(embedded_terms, collapse = " + ")
  f_candidate <- stats::as.formula(paste(outcome, "~", base_rhs, "+", candidate))
  f_residual <- stats::as.formula(paste(outcome, "~", base_rhs, "+", residual))
  m_candidate <- stats::glm(f_candidate, data = d, family = stats::binomial())
  m_residual <- stats::glm(f_residual, data = d, family = stats::binomial())
  lp_candidate <- as.numeric(stats::predict(m_candidate, type = "link"))
  lp_residual <- as.numeric(stats::predict(m_residual, type = "link"))
  lp_diff <- max(abs(lp_candidate - lp_residual), na.rm = TRUE)
  p_diff <- max(abs(stats::plogis(lp_candidate) - stats::plogis(lp_residual)), na.rm = TRUE)
  list(status = if (is.finite(lp_diff) && lp_diff <= tol) "same_model_space" else "diagnostic_mismatch",
       tolerance = tol, n = nrow(d), max_linear_predictor_difference = lp_diff,
       max_fitted_probability_difference = p_diff,
       coefficient_candidate = stats::coef(m_candidate),
       coefficient_residual = stats::coef(m_residual),
       log_likelihood_difference = as.numeric(stats::logLik(m_candidate) - stats::logLik(m_residual)))
}

audit_tolerance_grid <- function(index_terms, embedded_terms, model_terms = embedded_terms,
                                  tolerances = c(1e-8, 1e-10, 1e-12)) {
  do.call(rbind, lapply(tolerances, function(tol) {
    z <- audit_composite_index(index_terms, embedded_terms, model_terms, tol = tol)
    data.frame(tolerance = tol, status = z$status, exact = z$exact_reparameterization,
               residual_terms = paste(names(z$residual_terms), collapse = ";"), stringsAsFactors = FALSE)
  }))
}

# ---------------------------------------------------------------------------
# v32 product-scale audit. Some registered indices are products on their
# native scale (TyG_BMI = TyG * BMI; LAP = (waist - 65 or 58) * TG; HOMA_IR =
# glucose * insulin / 22.5; CVAI is a sex-specific linear form). They are not
# additively decomposable, but the distributive law expands them into a
# linear combination of columns already in the model plus products of pairs
# of model columns. Augmenting y ~ W with such a C therefore spans exactly
# y ~ W + I, where I collects the interaction columns.
#
# expansion is a list with fields:
#   intercept  numeric constant, absorbed by the model intercept (default 0)
#   linear     named numeric: coefficient a_j multiplying model column j
#   products   named numeric: coefficient c_k multiplying the product of two
#              model columns; names use model.matrix syntax "colA:colB". The
#              dummy column "sexFemale" counts as available whenever the
#              factor "sex" is among model_terms.
#
# Classification:
#   reducible_to_interaction  every linear shift targets a model column and
#                             every product involves two available columns
#   outside_supported_class   otherwise; blocking_reason names the offending
#                             expansion terms. No structural claim is made.
#
# Column-space geometry: with base columns W and
#   C = const + sum_j a_j w_j + sum_k c_k I_k,
# the single column C adds one dimension to span[W] (unless C already lies in
# span[W]), namely the fixed-ratio bundle sum_k c_k I_k. The unrestricted
# interaction space span[W, I_1, ..., I_k] has k additional dimensions and
# nests span[W, C]. Hence augmenting with the product index is a 1-df test of
# the bundled interactions, not k separate tests.
#
# Coefficient map (closed form): any fitted coefficient vector beta of
# y ~ W + C maps to the corresponding point of the interaction
# representation:
#   gamma_0  = beta_0 + const * beta_C
#   gamma_j  = beta_j + a_j * beta_C     (columns receiving a linear shift)
#   gamma_Ik = c_k * beta_C              (the bundled interaction columns)
# and X_[W,I] gamma = X_[W,C] beta holds exactly for any beta. The map does
# not make the rank-deficient beta_C itself uniquely estimable.
# ---------------------------------------------------------------------------

#' Audit a product-form index in a stated base model
#'
#' Product indices are handled by expanding them with the distributive law and
#' then asking whether every expansion term is expressible from the base
#' columns.  When it is, the index is a rank-one bundle of interactions and is
#' a 1-df test, not \code{k} separate tests.
#'
#' @param index Character scalar, the index name.
#' @param model_terms Character vector of base model columns.
#' @param expansion A list with \code{intercept}, \code{linear} and
#'   \code{products}, as returned by \code{\link{product_expansions}}.
#' @param data Optional data frame for the numerical verification layer.
#' @param tol Numerical tolerance.
#' @return A list with \code{classification} (reducible_to_interaction or
#'   outside_supported_class), \code{coefficient_map} and
#'   \code{blocking_reason}.
#' @export
#' @examples
#' audit_product_in_model("TyG_BMI", c("sex", "TyG", "BMI"),
#'                        product_expansions()[["TyG_BMI"]])$classification
audit_product_in_model <- function(index, model_terms, expansion,
                                   data = NULL, tol = 1e-10) {
  stopifnot(is.character(index), length(index) == 1L,
            is.character(model_terms), length(model_terms) >= 1L,
            is.list(expansion))
  const <- if (is.null(expansion$intercept)) 0 else as.numeric(expansion$intercept)
  linear <- expansion$linear
  products <- expansion$products
  if (is.null(linear)) linear <- setNames(numeric(0), character(0))
  if (is.null(products)) products <- setNames(numeric(0), character(0))
  stopifnot(is.numeric(linear), is.numeric(products))
  if (anyDuplicated(names(linear)) || anyDuplicated(names(products)))
    stop("expansion term names must be unique")
  # model.matrix-style column availability: factor sex contributes sexFemale.
  model_cols <- ifelse(model_terms == "sex", "sexFemale", model_terms)
  available <- union(model_terms, model_cols)
  split_product <- function(nm) strsplit(nm, ":", fixed = TRUE)[[1]]
  unsupported_products <- character()
  for (nm in names(products)) {
    parts <- split_product(nm)
    if (length(parts) != 2L || !all(parts %in% available))
      unsupported_products <- c(unsupported_products, nm)
  }
  unsupported_linear <- setdiff(names(linear), model_cols)
  reducible <- length(unsupported_products) == 0L && length(unsupported_linear) == 0L
  classification <- if (reducible) "reducible_to_interaction" else "outside_supported_class"
  blocking_reason <- if (reducible) NA_character_ else paste0(
    "expansion term(s) not expressible from the model columns: ",
    paste(c(unsupported_linear, unsupported_products), collapse = ", "))

  coefficient_map <- NULL
  if (reducible) {
    b_names <- c("(Intercept)", model_cols, "C")
    g_names <- c("(Intercept)", model_cols, names(products))
    K <- matrix(0, nrow = length(g_names), ncol = length(b_names),
                dimnames = list(g_names, b_names))
    K["(Intercept)", "(Intercept)"] <- 1
    if (const != 0) K["(Intercept)", "C"] <- const
    for (j in model_cols) K[j, j] <- 1
    for (j in names(linear)) K[j, "C"] <- linear[[j]]
    for (ik in names(products)) K[ik, "C"] <- products[[ik]]
    coefficient_map <- K
  }

  # Numerical verification layer (never the proof): rebuild C from the
  # expansion on the supplied rows and compare the column spaces.
  identity_error <- NA_real_; rank_base <- NA_integer_
  rank_augmented <- NA_integer_; rank_interactions <- NA_integer_
  rank_gain <- NA_integer_; interaction_gain <- NA_integer_
  c_span_residual <- NA_real_
  if (!is.null(data)) {
    needed <- unique(c(model_terms, index))
    missing_cols <- setdiff(needed, names(data))
    if (length(missing_cols))
      stop("data is missing columns: ", paste(missing_cols, collapse = ", "))
    column_of <- function(nm) {
      if (nm %in% names(data)) return(data[[nm]])
      if (nm == "sexFemale") return(as.numeric(data$sex == "Female"))
      stop("cannot materialise column: ", nm)
    }
    W <- stats::model.matrix(
      stats::as.formula(paste("~", paste(model_terms, collapse = " + "))),
      data = data)
    W <- W[, colnames(W) != "(Intercept)", drop = FALSE]
    # The expansion can only be materialised when every term is expressible
    # from the model columns; otherwise the audit stops at the boundary.
    if (reducible) {
      rebuilt <- rep(const, nrow(data))
      for (j in names(linear)) rebuilt <- rebuilt + linear[[j]] * column_of(j)
      for (ik in names(products)) {
        parts <- split_product(ik)
        rebuilt <- rebuilt + products[[ik]] * column_of(parts[1]) * column_of(parts[2])
      }
      identity_error <- suppressWarnings(max(abs(data[[index]] - rebuilt), na.rm = TRUE))
    }
    keep <- stats::complete.cases(data.frame(W, C = data[[index]]))
    rank_fun <- function(x) qr(x, tol = tol)$rank
    if (sum(keep) > ncol(W) + 2L) {
      xb <- cbind(1, W[keep, , drop = FALSE])
      rank_base <- rank_fun(xb)
      rank_augmented <- rank_fun(cbind(xb, data[[index]][keep]))
      rank_gain <- rank_augmented - rank_base
      if (reducible && length(products)) {
        Imat <- sapply(names(products), function(ik) {
          parts <- split_product(ik)
          column_of(parts[1])[keep] * column_of(parts[2])[keep]
        })
        xi <- cbind(xb, Imat)
        rank_interactions <- rank_fun(xi)
        interaction_gain <- rank_interactions - rank_base
        # C must lie (to tolerance) inside span[W, I]: the bundled column is a
        # fixed combination of the interaction columns plus model columns.
        c_span_residual <- max(abs(qr.resid(qr(xi, tol = tol), data[[index]][keep])))
      } else if (reducible) {
        rank_interactions <- rank_base
        interaction_gain <- 0L
        c_span_residual <- max(abs(qr.resid(qr(xb, tol = tol), data[[index]][keep])))
      }
    }
  }

  expansion_expression <- paste(c(
    if (const != 0) format(const, digits = 6),
    paste(format(linear, digits = 6), names(linear), sep = " x "),
    paste(format(products, digits = 6), names(products), sep = " x ")), collapse = " + ")
  failure_states <- character()
  if (!reducible) failure_states <- c(failure_states, "outside_supported_class")
  if (is.finite(identity_error) && identity_error > tol)
    failure_states <- c(failure_states,
      sprintf("expansion does not reproduce the stored column (max |difference| = %.3e)", identity_error))

  list(
    index = index, classification = classification,
    coefficient_map = coefficient_map, blocking_reason = blocking_reason,
    expansion = list(intercept = const, linear = linear, products = products),
    expansion_expression = if (nzchar(expansion_expression)) expansion_expression else "0",
    numerical_diagnostics = list(
      identity_error = identity_error, rank_base = rank_base,
      rank_augmented = rank_augmented, rank_interactions = rank_interactions,
      rank_gain = rank_gain, interaction_gain = interaction_gain,
      c_span_residual = c_span_residual, tolerance = tol),
    failure_states = failure_states,
    symbolic_layer = "the distributive expansion and the map K follow from the registered definitions",
    numerical_layer = "ranks and the identity error are computed on the supplied rows and verify, but do not prove, the symbolic claim")
}
