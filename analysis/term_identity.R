# term_identity.R (PURE ASCII)
# Five-field identity for a model term, introduced in v34 so that the public
# API no longer matches terms by bare variable name.
#
# A term key is the string
#   variable | transform | base | units | encoding
# Examples:
#   "TG|log10|10|mg/dL|-"         log10-transformed triglycerides in mg/dL
#   "glucose|identity|-|mg/dL|-"  raw fasting glucose in mg/dL
#   "sex|identity|-|-|factor(Male,Female)"
#
# Matching has three outcomes:
#   "match"     all five fields agree (both sides declared units)
#   "unkeyed"   variable, transform, log base and encoding agree, but at
#               least one side left units undeclared ("?"): the identity is
#               plausible, not certified, and the caller must warn
#   "mismatch"  any of variable / transform / base / encoding differs, or
#               both units are declared and differ
# A declared unit never matches a different declared unit; a declared
# transform never matches identity.  Unit metadata therefore fails visibly
# (a warning) rather than silently passing.


#' Build a five-field term identity key
#'
#' @param variable Column name.
#' @param transform Transformation applied, for example \code{"identity"},
#'   \code{"log"} or \code{"log10"}.
#' @param base Log base, or \code{"-"}.
#' @param units Units of the variable, or \code{"-"}.
#' @param encoding Factor coding, or \code{"-"}.
#' @return A character scalar of the form
#'   \code{variable|transform|base|units|encoding}.
#' @export
#' @examples
#' make_term_key("tg_mmol_L", "log10", "10", "mmol/L", "-")
make_term_key <- function(variable, transform = "identity", base = "-",
                          units = "?", encoding = "-") {
  stopifnot(length(variable) == 1L, nzchar(variable))
  transform <- match.arg(transform, c("identity", "log", "log10"))
  base <- if (transform %in% c("log", "log10")) as.character(base) else "-"
  sprintf("%s|%s|%s|%s|%s", variable, transform, base, units, encoding)
}

parse_term_key <- function(key) {
  f <- strsplit(key, "|", fixed = TRUE)[[1]]
  if (length(f) != 5L) stop("malformed term key: ", key)
  list(variable = f[1], transform = f[2], base = f[3], units = f[4],
       encoding = f[5])
}

term_key_compare <- function(a, b) {
  fa <- parse_term_key(a); fb <- parse_term_key(b)
  if (fa$variable != fb$variable) return("mismatch")
  if (fa$transform != fb$transform) return("mismatch")
  if (fa$transform %in% c("log", "log10") && fa$base != fb$base) return("mismatch")
  if (fa$encoding != fb$encoding) return("mismatch")
  if (fa$units == "?" || fb$units == "?") return("unkeyed")
  if (fa$units != fb$units) return("mismatch")
  "match"
}

# Vectorised membership: for each key in `keys`, "match" / "unkeyed" /
# "mismatch" against the best available counterpart in `against`.

#' Compare index term keys against model term keys
#'
#' @param keys Character vector of index term keys.
#' @param against Character vector of model term keys.
#' @return A character vector: \code{"match"}, \code{"mismatch"} or
#'   \code{"unkeyed"}.
#' @export
#' @examples
#' term_key_in(make_term_key("tg", "log", "e", "mmol/L", "-"),
#'             make_term_key("tg", "log", "e", "mmol/L", "-"))
term_key_in <- function(keys, against) {
  vapply(keys, function(k) {
    cmp <- vapply(against, term_key_compare, character(1), b = k)
    if (any(cmp == "match")) "match"
    else if (any(cmp == "unkeyed")) "unkeyed"
    else "mismatch"
  }, character(1))
}

# Convenience: key for a bare column name whose transformation and units the
# caller did not declare.
bare_term_key <- function(variable) {
  make_term_key(variable, transform = "identity", base = "-", units = "?",
                encoding = "-")
}
