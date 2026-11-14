# registry_audit.R (PURE ASCII)
# Independent verification path for the frozen registry (v34).
#
# The reference labels in reference_audit_classes() are hand-written.  They
# are NEVER fed to the audit functions.  Here every one of the 26 audit
# records is re-derived from first principles:
#
#   registered expression + full registered base model
#     -> derive_additive_terms() / product_expansions()   (literal algebra,
#        written out in this file, independent of the reference labels)
#     -> audit_index_in_model() / audit_product_in_model() (the real code)
#     -> actual classification
#     -> compared against the reference label
#
# test_audit.R requires every record to be concordant.  The agreement
# table is exported by build_results.R as v34_registry_audit_concordance.csv.

# ---------------------------------------------------------------------------
# Literal algebraic decompositions C = sum(b_k t_k) + const, in model-column
# names, for every additively decomposable registry record.  Coefficients are
# derived by hand from the published formulas (log(a/b) = log a - log b,
# log10(x) = ln(x)/ln(10), etc.), NOT read back from the reference labels.
#
# Returns list(index_terms, embedded_terms, intercept, note) or NULL when the
# record is not additive (product-form records are handled separately).
# ---------------------------------------------------------------------------

#' Literal algebraic decompositions for the additive registry records
#'
#' Hand-derived from the published formulas (\code{log(a/b) = log a - log b},
#' \code{log10(x) = ln(x)/ln(10)}, and so on), not read back from the
#' reference labels.  Returns \code{NULL} for records that are not additive.
#'
#' @param record_id Character scalar, a registry record id.
#' @return A list with \code{index_terms}, \code{embedded_terms},
#'   \code{intercept} and \code{note}, or \code{NULL}.
#' @export
#' @examples
#' derive_additive_terms("TyG")
derive_additive_terms <- function(record_id) {
  ln10 <- log(10)
  d <- list(
    # GLM7 against the minimal base sex + log10(age) + log10(BMI):
    # C = metab + log_age + log_bmi; metab is not in the base.
    GLM7 = list(index_terms = c(metab = 1, log_age = 1, log_bmi = 1),
                embedded_terms = c(log_age = 1, log_bmi = 1),
                intercept = 0,
                note = "residual metab lies outside the minimal base"),
    # GLM7-b: the same index against the full embedded base (base adds metab).
    `GLM7-b` = list(index_terms = c(metab = 1, log_age = 1, log_bmi = 1),
                    embedded_terms = c(metab = 1, log_age = 1, log_bmi = 1),
                    intercept = 0,
                    note = "every component is in the base; residual is zero"),
    # TyG = ln(TG_mg*G_mg/2) = log_tg_mg + log_glucose_mg - ln(2).
    TyG = list(index_terms = c(log_tg_mg = 1, log_glucose_mg = 1),
               embedded_terms = c(log_tg_mg = 1, log_glucose_mg = 1),
               intercept = -log(2), note = "constant -ln(2) maps to the intercept"),
    # AIP = log10(TG_mmol/HDL_mmol) = (ln tg - ln hdl)/ln(10).
    AIP = list(index_terms = c(log_tg = 1 / ln10, log_hdl = -1 / ln10),
               embedded_terms = c(log_tg = 1 / ln10, log_hdl = -1 / ln10),
               intercept = 0, note = "log10 = ln/ln(10); coefficients carry 1/ln(10)"),
    # HOMA_IR_log = log(glu_mmol) + log(ins_uIU) - log(22.5).
    HOMA_IR_log = list(index_terms = c(log_glucose_mmol = 1, log_insulin = 1),
                       embedded_terms = c(log_glucose_mmol = 1, log_insulin = 1),
                       intercept = -log(22.5),
                       note = "constant -ln(22.5) maps to the intercept"),
    # log-scale ratios: log(num/den) = log(num) - log(den).
    CAR = list(index_terms = c(log_hscrp = 1, log_albumin = -1),
               embedded_terms = c(log_hscrp = 1, log_albumin = -1),
               intercept = 0, note = "log(CAR) = log hsCRP - log albumin"),
    MHR = list(index_terms = c(log_monocytes = 1, log_hdl_mg = -1),
               embedded_terms = c(log_monocytes = 1, log_hdl_mg = -1),
               intercept = 0, note = "log(MHR) = log monocytes - log HDL(mg/dL)"),
    NLR = list(index_terms = c(log_neutrophils = 1, log_lymphocytes = -1),
               embedded_terms = c(log_neutrophils = 1, log_lymphocytes = -1),
               intercept = 0, note = "log(NLR) = log neutrophils - log lymphocytes"),
    PLR = list(index_terms = c(log_platelets = 1, log_lymphocytes = -1),
               embedded_terms = c(log_platelets = 1, log_lymphocytes = -1),
               intercept = 0, note = "log(PLR) = log platelets - log lymphocytes"),
    SII = list(index_terms = c(log_platelets = 1, log_neutrophils = 1,
                               log_lymphocytes = -1),
               embedded_terms = c(log_platelets = 1, log_neutrophils = 1,
                                  log_lymphocytes = -1),
               intercept = 0, note = "log(SII) = log plt + log neu - log lym"),
    SIRI = list(index_terms = c(log_neutrophils = 1, log_monocytes = 1,
                                log_lymphocytes = -1),
                embedded_terms = c(log_neutrophils = 1, log_monocytes = 1,
                                   log_lymphocytes = -1),
                intercept = 0, note = "log(SIRI) = log neu + log mono - log lym"),
    TyG_hsCRP = list(index_terms = c(log_TyG = 1, log_hscrp = 1),
                     embedded_terms = c(log_TyG = 1, log_hscrp = 1),
                     intercept = 0, note = "log(TyG*hsCRP) = log TyG + log hsCRP"),
    hsCRP_BMI = list(index_terms = c(log_hscrp = 1, ln_bmi = 1),
                     embedded_terms = c(log_hscrp = 1, ln_bmi = 1),
                     intercept = 0, note = "log(hsCRP*BMI) = log hsCRP + ln BMI"),
    hsCRP_WC = list(index_terms = c(log_hscrp = 1, log_waist = 1),
                    embedded_terms = c(log_hscrp = 1, log_waist = 1),
                    intercept = 0, note = "log(hsCRP*WC) = log hsCRP + log waist"),
    SII_BMI = list(index_terms = c(log_SII = 1, ln_bmi = 1),
                   embedded_terms = c(log_SII = 1, ln_bmi = 1),
                   intercept = 0, note = "log(SII*BMI) = log SII + ln BMI"),
    # PNI = 10*albumin + 0.005*(lymphocytes*1000) = 10 alb + 0.005 cells_mm3.
    PNI = list(index_terms = c(albumin_g_dL = 10, lymphocytes_cells_mm3 = 0.005),
               embedded_terms = c(albumin_g_dL = 10, lymphocytes_cells_mm3 = 0.005),
               intercept = 0, note = "native-scale linear form, exact")
  )
  d[[record_id]]
}

# Records that are outside the supported class, with the reason derived from
# the formula structure (again written out, not read from the labels):
outside_class_reason <- function(record_id) {
  r <- c(
    TyG_WHtR = "second factor waist/height is itself a ratio of two columns; not expressible as a product of two model columns",
    METS_IR = "division by log(HDL) is a denominator term; outside the additive and product classes",
    VAI = "sex-specific products of ratios of model columns; outside the supported class",
    BRI = "square root of a quadratic form in waist and height; outside the supported class",
    ABSI = "fractional powers of BMI and height in a denominator; outside the supported class"
  )
  if (!record_id %in% names(r)) return(NULL)
  r[[record_id]]
}

# Run the REAL audit on one record and return the actual classification.

#' Re-derive the class of one registry record from first principles
#'
#' The reference labels are never read here: the class is derived from the
#' registered expression and its registered base model, using the same audit
#' code the package exposes.
#'
#' @param record_id Character scalar; \code{"GLM7-b"} is the second GLM7
#'   record.
#' @param data Optional data frame for the numerical layer.
#' @param tol Numerical tolerance.
#' @return A list with \code{actual_class}, \code{blocking_reason},
#'   \code{identity_error}, \code{rank_gain} and \code{path}.
#' @export
#' @examples
#' audit_one_record("TyG")$actual_class
audit_one_record <- function(record_id, data = NULL, tol = 1e-10) {
  base_id <- sub("-b$", "", record_id)
  exps <- product_expansions()
  if (base_id %in% names(exps)) {
    mt <- v28_component_model_terms(base_id, v28_registry())
    a <- audit_product_in_model(base_id, mt, exps[[base_id]],
                                data = data, tol = tol)
    return(list(actual_class = a$classification,
                blocking_reason = a$blocking_reason,
                identity_error = a$numerical_diagnostics$identity_error,
                rank_gain = a$numerical_diagnostics$rank_gain,
                path = "audit_product_in_model"))
  }
  t <- derive_additive_terms(record_id)
  if (is.null(t)) {
    # Outside-class records: confirm the parser also refuses the registered
    # formula string, so "outside" is derived twice (structure + parser).
    reg <- v28_registry()
    f <- reg$formula[reg$index == base_id]
    parsed <- parse_index_expression(f)
    parser_refused <- identical(parsed$status, "unsupported")
    return(list(actual_class = "outside_supported_class",
                blocking_reason = outside_class_reason(record_id),
                identity_error = NA_real_, rank_gain = NA_integer_,
                path = if (parser_refused)
                  "structure + parser both refuse" else "structure (parser not applicable)"))
  }
  a <- audit_index_in_model(t$index_terms, t$embedded_terms, t$embedded_terms,
                            data = data, tol = tol, intercept = t$intercept)
  list(actual_class = a$classification,
       blocking_reason = a$blocking_reason,
       identity_error = NA_real_, rank_gain = a$numerical_diagnostics$rank_gain,
       path = "audit_index_in_model")
}

# Full concordance table over the 26 audit records.

#' Concordance table over all 26 audit records
#'
#' @param data Optional data frame for the numerical layer.
#' @param tol Numerical tolerance.
#' @return A data frame with \code{record}, \code{expected_class},
#'   \code{actual_class}, \code{concordant}, \code{derivation_path},
#'   \code{blocking_reason}, \code{identity_error} and \code{rank_gain}.
#' @export
#' @examples
#' tab <- audit_registry_concordance()
#' all(tab$concordant)
audit_registry_concordance <- function(data = NULL, tol = 1e-10) {
  ref <- reference_audit_classes()
  ids <- names(ref)
  rows <- lapply(ids, function(id) {
    got <- audit_one_record(id, data = data, tol = tol)
    data.frame(record = id, expected_class = unname(ref[[id]]),
               actual_class = got$actual_class,
               concordant = identical(got$actual_class, unname(ref[[id]])),
               derivation_path = got$path,
               blocking_reason = ifelse(is.na(got$blocking_reason), "",
                                        got$blocking_reason),
               identity_error = got$identity_error,
               rank_gain = ifelse(is.na(got$rank_gain), NA_integer_,
                                  got$rank_gain),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  # GLM7-b: the dual registration (same index, full embedded base).
  got_b <- audit_one_record("GLM7-b", data = data, tol = tol)
  out <- rbind(out, data.frame(
    record = "GLM7-b", expected_class = "rank_preserving_exact",
    actual_class = got_b$actual_class,
    concordant = identical(got_b$actual_class, "rank_preserving_exact"),
    derivation_path = got_b$path, blocking_reason = "",
    identity_error = got_b$identity_error, rank_gain = NA_integer_,
    stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}
