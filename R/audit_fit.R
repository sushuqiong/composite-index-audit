# audit_fit.R (PURE ASCII)
#
# v34: the entry point a real user reaches for.
#
# Until v33 the only way in was audit_registered_index(name, model_terms),
# which asks the user to transcribe both the coefficient vector of the index
# and the column names of the base model by hand.  A user who already has a
# fitted model object should not have to do that.  audit_fit() reads the
# design matrix, the terms, the aliasing pattern, the family, the offset and
# the weights straight off the fitted object and then delegates to the same
# symbolic layer the rest of the package uses.
#
# What audit_fit() does NOT do: it does not refit anything, it does not look
# at the outcome, and it does not change the classification rule.  The class
# it returns is the class the index has against the base model of this fit,
# which is a property of the pair (index, base).
#
# Term identity is never inferred from a variable name.  A column called
# log_bmi is treated as a bare name with undeclared transformation and units
# unless the caller supplies a NAMED model_terms vector mapping design-matrix
# column -> five-field term key (see term_identity.R).  In bare-name mode
# every comparison is reported as unkeyed and a warning is raised.

AUDIT_FIT_CLASSES <- c("glm", "lm", "svyglm", "coxph", "survreg")

# The four classes the audit can return, plus the states that are explicitly
# NOT classes.  v34 keeps them apart so the text can say "one of four
# classes, or an explicit non-classification state" without overcounting.
AUDIT_FOUR_CLASSES <- c("rank_preserving_exact", "rank_adding_exact",
                        "reducible_to_interaction", "outside_supported_class")
AUDIT_NON_CLASS_STATES <- c("approximate_collinearity", "fit_failed",
                            "not_decomposable")


#' Is a returned label one of the four classes?
#'
#' The audit returns one of four classes, or an explicit non-classification
#' state.  Keeping the two apart is what lets the text say so without
#' overcounting the classes.
#'
#' @param classification Character scalar, the label returned by an audit.
#' @return \code{"class"}, \code{"non_classification_state"} or
#'   \code{"unknown"}.
#' @export
#' @examples
#' audit_state_type("rank_preserving_exact")
#' audit_state_type("approximate_collinearity")
audit_state_type <- function(classification) {
  if (is.null(classification) || length(classification) != 1L ||
      is.na(classification)) return("unknown")
  if (classification %in% AUDIT_FOUR_CLASSES) return("class")
  if (classification %in% AUDIT_NON_CLASS_STATES) return("non_classification_state")
  "unknown"
}

# ---------------------------------------------------------------------------
# Design extraction.  Every supported model class is reduced to the same
# objects: model frame, terms, model matrix, coefficient vector (whose NA
# entries are the aliased columns).
# ---------------------------------------------------------------------------

#' Extract the design of a fitted model for auditing
#'
#' @param fit A fitted model object.
#' @return A list with \code{model_frame}, \code{terms}, \code{model_matrix},
#'   \code{coefficients} and \code{aliased_columns}.
#' @export
#' @examples
#' fit <- glm(dx ~ log_age, data = data.frame(log_age = c(1, 2, 3), dx = c(0, 1, 0)),
#'            family = binomial())
#' names(audit_fit_design(fit))
audit_fit_design <- function(fit) {
  mf <- tryCatch(stats::model.frame(fit), error = function(e) NULL)
  tt <- NULL
  if (!is.null(mf)) tt <- attr(mf, "terms")
  if (is.null(tt)) tt <- tryCatch(stats::terms(fit), error = function(e) NULL)
  X <- NULL
  if (!is.null(tt) && !is.null(mf))
    X <- tryCatch(stats::model.matrix(tt, data = mf), error = function(e) NULL)
  if (is.null(X)) X <- tryCatch(stats::model.matrix(fit), error = function(e) NULL)
  if (is.null(X))
    stop("audit_fit: no model matrix could be built from this fit; ",
         "supported classes are ", paste(AUDIT_FIT_CLASSES, collapse = ", "))
  cf <- tryCatch(stats::coef(fit), error = function(e) NULL)
  aliased <- if (!is.null(cf)) names(cf)[is.na(cf)] else character()
  list(model_frame = mf, terms = tt, model_matrix = X,
       coefficients = cf, aliased_columns = aliased)
}

audit_fit_family <- function(fit) {
  fam <- tryCatch(stats::family(fit), error = function(e) NULL)
  if (is.null(fam)) return(list(family = "unknown", link = "unknown"))
  list(family = fam$family %||% "unknown", link = fam$link %||% "unknown")
}

# Offsets are a genuine confound: if a component of the index enters the
# model as an offset, the index no longer enters a free direction even when
# the column is absent from the design matrix.
audit_fit_offset <- function(fit, design) {
  tt <- design$terms
  off <- if (!is.null(tt)) attr(tt, "offset") else NULL
  vars <- if (!is.null(tt)) attr(tt, "variables") else NULL
  names_off <- character()
  if (!is.null(off) && length(off) && !is.null(vars))
    names_off <- vapply(off, function(i)
      paste(deparse(vars[[i + 1L]]), collapse = " "), character(1))
  mo <- tryCatch(stats::model.offset(design$model_frame), error = function(e) NULL)
  present <- !is.null(mo) && length(mo) > 0L
  if (!length(names_off) && !present)
    return(list(has_offset = FALSE, offset_terms = character(),
                offset_confound = FALSE, note = "no offset in this fit"))
  list(has_offset = TRUE, offset_terms = names_off, offset_confound = TRUE,
       note = paste0("offset present",
         if (length(names_off)) paste0(" (", paste(names_off, collapse = ", "), ")") else "",
         ": the audit reads the design matrix only; an offset is held fixed and",
         " can absorb a component direction"))
}

# Weights do not enter the algebra, but the caller should know which weights
# the fit used, because they determine what the coefficients estimate.
audit_fit_weights <- function(fit) {
  w <- tryCatch(stats::weights(fit), error = function(e) NULL)
  src <- if (is.null(w) || !length(w)) "none"
    else if ("svyglm" %in% class(fit)) "survey design weights (svyglm)"
    else if ("coxph" %in% class(fit)) "survival model case weights (coxph)"
    else if ("survreg" %in% class(fit)) "survival model case weights (survreg)"
    else if ("glm" %in% class(fit)) "prior weights (glm)"
    else "prior weights (lm)"
  list(weights_source = src, n_weights = if (is.null(w)) NA_integer_ else length(w),
       note = "weights do not enter the algebraic classification; they determine what the fitted coefficients estimate")
}

# ---------------------------------------------------------------------------
# Symbolic dispatch against an arbitrary base.  Same three routes the
# registry audit uses, but with the base supplied by the fitted object:
#   product_expandable -> audit_product_in_model
#   additive           -> derive_additive_terms, embedded part = the part of
#                         the index that actually appears in this base
#   neither            -> outside_supported_class, with the parser asked to
#                         refuse the registered formula as a second opinion
# ---------------------------------------------------------------------------

#' Audit an index against an arbitrary stated base
#'
#' The three routes the registry audit uses, with the base supplied by the
#' caller: product expansions, additive decompositions, and the outside class.
#'
#' @param index Character scalar, the index name.
#' @param base_terms Character vector of base model terms.
#' @param data Optional data frame for the numerical layer.
#' @param tol Numerical tolerance.
#' @return The audit list.
#' @export
#' @examples
#' audit_symbolic_for_base("TyG", c("sex", "log_tg_mg", "log_glucose_mg"))$classification
audit_symbolic_for_base <- function(index, base_terms, data = NULL, tol = 1e-10) {
  exps <- product_expansions()
  if (index %in% names(exps))
    return(audit_product_in_model(index, base_terms, exps[[index]],
                                  data = data, tol = tol))
  t <- derive_additive_terms(index)
  if (!is.null(t)) {
    emb <- t$index_terms[intersect(names(t$index_terms), base_terms)]
    return(audit_index_in_model(t$index_terms, emb, emb, data = data,
                                tol = tol, intercept = t$intercept))
  }
  # Outside the supported class.  Derive the reason from the structure, and
  # ask the parser to refuse the registered formula as an independent check.
  reg <- v28_registry()
  f <- if (index %in% reg$index) reg$formula[reg$index == index] else NA_character_
  parser_refused <- if (!is.na(f))
    identical(parse_index_expression(f)$status, "unsupported") else NA
  list(classification = "outside_supported_class",
       coefficient_map = NULL,
       blocking_reason = outside_class_reason(index) %||%
         "the registered expression is outside the additive log-linear and product classes",
       residual_expression = NA_character_,
       numerical_diagnostics = list(rank_base = NA_integer_,
                                    rank_augmented = NA_integer_,
                                    rank_gain = NA_integer_,
                                    condition_number = NA_real_,
                                    empirical_correlation = NA_real_,
                                    tolerance = tol),
       failure_states = "outside_supported_class",
       symbolic_layer = paste0("structure refuses the expression",
         if (isTRUE(parser_refused)) "; the parser refuses the registered formula too"
         else if (isFALSE(parser_refused)) "; the parser accepted the registered formula, so the refusal rests on the algebraic structure"
         else ""),
       numerical_layer = "no decomposition is claimed, so no numerical check applies")
}

# ---------------------------------------------------------------------------
# audit_fit() -- the public entry point.
#
#   fit          a fitted glm / lm / svyglm / coxph / survreg object
#   index        a registered index name, or any name when expression is given
#   expression   optional; when supplied it is parsed and used instead of the
#                registry, so ad-hoc indices work without registration
#   model_terms  optional override: a character vector of design-matrix
#                columns, or a NAMED vector mapping column -> five-field key.
#                Defaults to every non-aliased design-matrix column except
#                the intercept and the index itself.
#   data         optional data frame for the numerical verification layer.
#                When NULL the audit runs in the symbolic branch, which is
#                the honest default: the model frame carries only the model
#                columns, so it can never verify an index whose components
#                lie outside the base.
# ---------------------------------------------------------------------------

#' Audit a composite index against an already-fitted model
#'
#' The entry point a user reaches for.  \code{audit_fit()} reads the design
#' matrix, the terms, the aliased columns, the family, the offset and the
#' weights off the fitted object and then delegates to the same symbolic layer
#' the rest of the package uses.  It does not refit anything and does not look
#' at the outcome: the class it returns is a property of the pair
#' (index, base model).
#'
#' Supported model classes are \code{glm}, \code{lm}, \code{svyglm},
#' \code{coxph} and \code{survreg}.  An offset is reported as
#' \code{"offset_present"} and aliased columns as
#' \code{"aliased_columns_present"}; neither is downgraded.
#'
#' Term identity is never inferred from a variable name: a column called
#' \code{log_bmi} is treated as a bare name with undeclared transformation and
#' units unless the caller passes a NAMED \code{model_terms} vector.
#'
#' @param fit A fitted model object of a supported class.
#' @param index Character scalar, a registered index name (or any name when
#'   \code{expression} is supplied).
#' @param expression Optional index expression, parsed and used instead of the
#'   registry.
#' @param data Optional data frame for the numerical verification layer.  When
#'   \code{NULL} the audit runs in the symbolic branch, which is the honest
#'   default: the model frame carries only the model columns, so it can never
#'   verify an index whose components lie outside the base.
#' @param model_terms Optional override, a character vector of design-matrix
#'   columns or a NAMED vector mapping column to five-field term key.
#' @param tol Numerical tolerance.
#' @param drop_index Drop the index itself from the base model terms.
#' @param keep_aliased Keep aliased columns among the base model terms.
#' @return The audit list, with an added \code{model_object} holding what was
#'   read off the fit, plus \code{state_type} and \code{is_classification}.
#' @export
#' @examples
#' d <- data.frame(sex = factor(c("M", "F", "M", "F")),
#'                 log_tg_mg = c(4.6, 4.8, 4.7, 5.0),
#'                 log_glucose_mg = c(4.3, 4.5, 4.4, 4.6),
#'                 log_age = c(1.6, 1.7, 1.65, 1.66),
#'                 log_bmi = c(1.4, 1.5, 1.45, 1.47),
#'                 dx = c(0, 1, 0, 1))
#' fit <- glm(dx ~ sex + log_age + log_bmi, data = d, family = binomial())
#' audit_fit(fit, "TyG")$classification
audit_fit <- function(fit, index, expression = NULL, data = NULL,
                      model_terms = NULL, tol = 1e-10,
                      drop_index = TRUE, keep_aliased = FALSE) {
  stopifnot(is.character(index), length(index) == 1L, nzchar(index))
  cls <- class(fit)
  if (!any(AUDIT_FIT_CLASSES %in% cls))
    stop("audit_fit: unsupported model class ", paste(cls, collapse = "/"),
         "; supported: ", paste(AUDIT_FIT_CLASSES, collapse = ", "))
  design <- audit_fit_design(fit)
  fam <- audit_fit_family(fit)
  off <- audit_fit_offset(fit, design)
  wts <- audit_fit_weights(fit)

  # Base model TERMS, not design-matrix column names.  A factor term such as
  # sex expands to one or more dummy columns (sexFemale, sexMale, ...), and the
  # symbolic layer is written in terms: product expansions refer to "sexFemale"
  # as the dummy of the term "sex".  attr(model_matrix, "assign") maps every
  # column back to its term, so the audit sees the term the user wrote.
  cols <- colnames(design$model_matrix)
  tl <- attr(design$terms, "term.labels")
  asg <- attr(design$model_matrix, "assign")
  lab <- rep(NA_character_, length(cols))
  if (!is.null(tl) && !is.null(asg) && length(asg) == length(cols))
    lab[asg > 0L] <- tl[asg[asg > 0L]]
  keep <- !is.na(lab) | cols != "(Intercept)"
  keep <- keep & cols != "(Intercept)"
  if (!isTRUE(keep_aliased)) keep <- keep & !(cols %in% design$aliased_columns)
  if (isTRUE(drop_index)) {
    keep <- keep & !(cols %in% index)
    if (any(!is.na(lab))) keep <- keep & !(lab %in% index)
  }
  cn <- unique(lab[keep])
  cn <- cn[!is.na(cn)]
  if (!length(cn)) cn <- unique(cols[keep])   # terms unavailable: fall back

  if (is.null(model_terms)) {
    model_terms <- cn
  } else if (is.null(names(model_terms))) {
    names(model_terms) <- rep("", length(model_terms))
  }

  keyed <- !is.null(names(model_terms)) && all(nzchar(names(model_terms)))
  base_for_symbolic <- if (keyed) names(model_terms) else unname(model_terms)

  res <- if (!is.null(expression)) {
    parsed <- parse_index_expression(expression)
    if (parsed$status == "unsupported")
      stop("expression is outside the supported grammar: ", parsed$expression)
    idx_terms <- c(parsed$linear, parsed$logs)
    emb <- if (keyed) {
      keys <- unname(model_terms)
      idx_keys <- vapply(names(idx_terms), function(nm)
        make_term_key(nm, "identity", "-", "?", "-"), character(1))
      idx_terms[term_key_in(idx_keys, keys) != "mismatch"]
    } else idx_terms[intersect(names(idx_terms), base_for_symbolic)]
    audit_index_in_model(idx_terms, emb, emb, data = data, tol = tol,
                         intercept = parsed$intercept)
  } else {
    audit_symbolic_for_base(index, base_for_symbolic, data = data, tol = tol)
  }

  in_fit <- index %in% colnames(design$model_matrix)
  res$model_object <- list(
    model_class = paste(cls, collapse = "/"),
    family = fam$family, link = fam$link,
    n_obs = nrow(design$model_matrix),
    index_is_in_fit = in_fit,
    index_coefficient = if (in_fit && !is.null(design$coefficients))
      unname(design$coefficients[[index]]) else NA_real_,
    index_coefficient_aliased = if (in_fit)
      isTRUE(is.na(design$coefficients[[index]])) else NA,
    design_matrix_columns = colnames(design$model_matrix),
    aliased_columns = design$aliased_columns,
    aliased_columns_dropped = if (isTRUE(keep_aliased)) character()
      else intersect(design$aliased_columns, colnames(design$model_matrix)),
    base_model_terms = cn,
    design_matrix_terms = unique(lab[!is.na(lab)]),
    term_matching = if (keyed) "keyed: five-field identity supplied by the caller"
      else paste("legacy_bare_names: design-matrix column names only;",
                 "transformation, log base, units and encoding are undeclared",
                 "and are NOT inferred from the column name"),
    weights = wts,
    offset = off)
  res$state_type <- audit_state_type(res$classification)
  res$is_classification <- identical(res$state_type, "class")
  res$failure_states <- unique(c(res$failure_states %||% character(),
    if (isTRUE(off$offset_confound)) "offset_present" else character(),
    if (length(design$aliased_columns)) "aliased_columns_present" else character()))
  res$entry_point <- paste(
    "audit_fit(): base model terms, aliasing, family, offset and weights were",
    "read from the fitted object; the classification rule is unchanged and",
    "the class remains a property of the pair (index, base model)")
  res
}

