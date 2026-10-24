# public_api_v32.R (PURE ASCII)
# Public software layer of the v32 audit tool.  Four entry points:
#
#   parse_index_expression(expr)        mini symbolic parser; classifies an
#                                       R expression as additive_log_linear,
#                                       product_expandable or unsupported
#   register_index(name, expression)    register a user index at runtime;
#                                       returns a registry row usable by
#                                       audit_registered_index()
#   verify_index_column(data, name)     recompute the registered expression
#                                       from the raw columns and compare with
#                                       the stored column (catches unit /
#                                       scale / transcription mistakes)
#   normalize_log_base(terms, ...)      rescale a named coefficient vector
#                                       between natural-log and log10 bases
#
# The parser deliberately supports a small, honest grammar: column names,
# numeric constants, + - * / (division by constants only), log() and log10()
# around single columns at the additive top level.  Anything else (powers,
# square roots, division by columns, logs inside products) is reported as
# unsupported instead of being approximated.

# ---------------------------------------------------------------------------
# Monomial representation: list of list(coef = numeric, cols = character(0..2),
# log_cols = named integer).  A monomial is either a plain product of up to
# two raw columns or a single log-transformed column added at top level.
# ---------------------------------------------------------------------------
.v32_mono_const <- function(k) list(list(coef = k, cols = character(0), log_cols = setNames(integer(0), character(0))))
.v32_mono_col   <- function(nm) list(list(coef = 1, cols = nm, log_cols = setNames(integer(0), character(0))))
.v32_mono_log   <- function(nm, base) list(list(coef = 1, cols = character(0), log_cols = setNames(as.integer(base), nm)))

.v32_mono_add <- function(a, b) c(a, b)
.v32_mono_negate <- function(a) lapply(a, function(m) { m$coef <- -m$coef; m })
.v32_mono_mul <- function(a, b) {
  out <- list()
  for (ma in a) for (mb in b) {
    cols <- c(ma$cols, mb$cols)
    logs <- c(ma$log_cols, mb$log_cols)
    if (length(logs) && (length(cols) || length(ma$log_cols) + length(mb$log_cols) > 1L))
      return("UNSUPPORTED")   # log terms may not enter products or stack
    if (length(cols) > 2L) return("UNSUPPORTED")
    out[[length(out) + 1L]] <- list(coef = ma$coef * mb$coef, cols = cols, log_cols = logs)
  }
  out
}

.v32_walk <- function(e) {
  if (is.numeric(e)) return(.v32_mono_const(as.numeric(e)))
  if (is.name(e)) return(.v32_mono_col(as.character(e)))
  if (!is.call(e)) return("UNSUPPORTED")
  op <- as.character(e[[1]])
  args <- as.list(e[-1])
  if (op == "(") return(.v32_walk(args[[1]]))
  if (op == "+") {
    if (length(args) == 1L) return(.v32_walk(args[[1]]))
    a <- .v32_walk(args[[1]]); b <- .v32_walk(args[[2]])
    if (identical(a, "UNSUPPORTED") || identical(b, "UNSUPPORTED")) return("UNSUPPORTED")
    return(.v32_mono_add(a, b))
  }
  if (op == "-") {
    a <- .v32_walk(args[[1]])
    if (identical(a, "UNSUPPORTED")) return("UNSUPPORTED")
    if (length(args) == 1L) return(.v32_mono_negate(a))
    b <- .v32_walk(args[[2]])
    if (identical(b, "UNSUPPORTED")) return("UNSUPPORTED")
    return(.v32_mono_add(a, .v32_mono_negate(b)))
  }
  if (op == "*") {
    a <- .v32_walk(args[[1]]); b <- .v32_walk(args[[2]])
    if (identical(a, "UNSUPPORTED") || identical(b, "UNSUPPORTED")) return("UNSUPPORTED")
    return(.v32_mono_mul(a, b))
  }
  if (op == "/") {
    a <- .v32_walk(args[[1]])
    if (identical(a, "UNSUPPORTED") || !is.numeric(args[[2]])) return("UNSUPPORTED")
    denom <- as.numeric(args[[2]])
    if (denom == 0) return("UNSUPPORTED")
    return(lapply(a, function(m) { m$coef <- m$coef / denom; m }))
  }
  if (op %in% c("log", "log10")) {
    inner <- args[[1]]
    if (!is.name(inner)) return("UNSUPPORTED")
    return(.v32_mono_log(as.character(inner), if (op == "log") 0L else 10L))
  }
  "UNSUPPORTED"
}

parse_index_expression <- function(expr) {
  e <- if (is.character(expr)) parse(text = expr)[[1]] else substitute(expr)
  monos <- .v32_walk(e)
  if (identical(monos, "UNSUPPORTED"))
    return(list(status = "unsupported", expression = deparse(e),
                linear = numeric(0), products = numeric(0), logs = numeric(0),
                intercept = 0))
  # merge like monomials
  key_of <- function(m) paste0(paste(sort(m$cols), collapse = "*"),
                               if (length(m$log_cols)) paste0("|log", m$log_cols, ":", names(m$log_cols)) else "")
  keys <- vapply(monos, key_of, character(1))
  merged <- list()
  for (i in seq_along(monos)) {
    k <- keys[i]
    if (k %in% names(merged)) merged[[k]]$coef <- merged[[k]]$coef + monos[[i]]$coef
    else merged[[k]] <- monos[[i]]
  }
  monos <- unname(merged)
  linear <- numeric(0); products <- numeric(0); logs <- numeric(0); intercept <- 0
  log_bases <- integer(0)
  for (m in monos) {
    if (length(m$log_cols)) {
      nm <- names(m$log_cols); logs[nm] <- m$coef; log_bases[nm] <- m$log_cols
    } else if (length(m$cols) == 0L) intercept <- intercept + m$coef
    else if (length(m$cols) == 1L) linear[m$cols] <- m$coef
    else products[paste(m$cols, collapse = ":")] <- m$coef
  }
  status <- if (length(products)) "product_expandable"
            else if (length(logs) || length(linear)) "additive_log_linear"
            else "unsupported"
  list(status = status, expression = deparse(e), intercept = intercept,
       linear = linear, products = products, logs = logs, log_bases = log_bases)
}

# ---------------------------------------------------------------------------
.v32_user_registry <- new.env()

register_index <- function(name, expression, units = character(0),
                           overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  if (exists(name, envir = .v32_user_registry) && !overwrite)
    stop("index already registered: ", name, " (use overwrite = TRUE)")
  parsed <- parse_index_expression(expression)
  if (parsed$status == "unsupported")
    stop("expression is outside the supported grammar: ", parsed$expression)
  if (length(units)) {
    stopifnot(is.character(units), !is.null(names(units)))
    parsed$units <- units
  }
  assign(name, parsed, envir = .v32_user_registry)
  invisible(parsed)
}

# Dispatch one registered (or built-in) index to the right audit layer.
#
# v33: term identity is carried through the dispatch.  model_terms may be
#   - a plain character vector of column names (legacy): each name becomes a
#     bare key with undeclared transformation and units, and any reliance on
#     it is reported as "unkeyed" in the output and by warning;
#   - a NAMED character vector mapping data column name -> five-field term
#     key (see term_identity_v33.R), e.g.
#     c(tg_mmol_L = "TG|log10|10|mmol/L|-", glucose_mg_dL = "glucose|log10|10|mg/dL|-").
# A term matches the model only when variable, transformation, log base and
# encoding agree; declared units must agree as well.  log(x) therefore never
# collapses onto x or 1/x, and log10 never silently equals ln.
audit_registered_index <- function(name, model_terms, data = NULL, tol = 1e-10) {
  parsed <- if (exists(name, envir = .v32_user_registry))
    get(name, envir = .v32_user_registry) else {
      exps <- product_expansions()
      if (name %in% names(exps)) {
        e <- exps[[name]]
        list(status = "product_expandable", intercept = e$intercept %||% 0,
             linear = e$linear %||% numeric(0), products = e$products %||% numeric(0),
             logs = numeric(0), log_bases = integer(0))
      } else stop("index is not registered: ", name)
    }
  if (parsed$status == "product_expandable")
    return(audit_product_in_model(
      name, unname(model_terms),
      list(intercept = parsed$intercept, linear = parsed$linear,
           products = parsed$products),
      data = data, tol = tol))
  # additive (log-)linear: C = const + sum_j a_j w_j; the constant is absorbed
  # by the model intercept and enters the coefficient map.
  units <- parsed$units %||% character(0)
  unit_of <- function(nm) if (nm %in% names(units)) units[[nm]] else "?"
  idx_terms <- c(parsed$linear, parsed$logs)
  idx_keys <- c(
    vapply(names(parsed$linear),
           function(nm) make_term_key(nm, "identity", "-", unit_of(nm), "-"),
           character(1)),
    vapply(names(parsed$logs), function(nm) {
      b <- parsed$log_bases[[nm]]  # parser stores 0L for ln, 10L for log10
      make_term_key(nm, if (b == 10L) "log10" else "log",
                    if (b == 10L) "10" else "e", unit_of(nm), "-")
    }, character(1)))
  keyed <- !is.null(names(model_terms)) && all(nzchar(names(model_terms)))
  if (keyed) {
    model_keys <- unname(model_terms)
    names(model_keys) <- names(model_terms)
  } else {
    model_keys <- setNames(vapply(model_terms, bare_term_key, character(1)),
                           model_terms)
  }
  cmp <- term_key_in(idx_keys, unname(model_keys))
  emb <- idx_terms[cmp != "mismatch"]
  unkeyed_terms <- names(idx_terms)[cmp == "unkeyed"]
  out <- audit_index_in_model(idx_terms, emb, emb, data = data, tol = tol,
                              intercept = parsed$intercept)
  out$term_matching <- list(
    mode = if (keyed) "keyed" else "legacy_bare_names",
    unkeyed_terms = unkeyed_terms,
    note = paste("terms matched on the five-field identity",
                 "variable|transform|base|units|encoding;",
                 "'unkeyed' means units were undeclared on at least one side"))
  if (length(unkeyed_terms) || !keyed)
    warning("audit_registered_index: unit metadata undeclared for ",
            paste(c(names(model_terms), unkeyed_terms), collapse = ", "),
            " - declare units via register_index(units = ...) and named model_terms keys")
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
verify_index_column <- function(data, name, expression = NULL, tol = 1e-8) {
  if (!name %in% names(data)) stop("column not in data: ", name)
  parsed <- if (!is.null(expression)) parse_index_expression(expression)
    else if (exists(name, envir = .v32_user_registry)) get(name, envir = .v32_user_registry)
    else stop("no expression supplied and index not registered: ", name)
  if (parsed$status == "unsupported")
    stop("expression is outside the supported grammar: ", parsed$expression)
  rebuilt <- rep(parsed$intercept, nrow(data))
  need <- c(names(parsed$linear), names(parsed$logs))
  for (nm in need) if (!nm %in% names(data)) stop("data is missing column: ", nm)
  for (nm in names(parsed$linear)) rebuilt <- rebuilt + parsed$linear[[nm]] * data[[nm]]
  for (nm in names(parsed$logs)) {
    fun <- if (parsed$log_bases[[nm]] == 10L) log10 else log
    rebuilt <- rebuilt + parsed$logs[[nm]] * fun(data[[nm]])
  }
  for (ik in names(parsed$products)) {
    parts <- strsplit(ik, ":", fixed = TRUE)[[1]]
    for (p in parts) if (!p %in% names(data)) stop("data is missing column: ", p)
    rebuilt <- rebuilt + parsed$products[[ik]] * data[[parts[1]]] * data[[parts[2]]]
  }
  err <- suppressWarnings(max(abs(data[[name]] - rebuilt), na.rm = TRUE))
  list(index = name, max_abs_error = err, ok = is.finite(err) && err <= tol,
       tolerance = tol, expression = parsed$expression)
}

# ---------------------------------------------------------------------------
normalize_log_base <- function(terms, from = c("ln", "log10"), to = c("ln", "log10")) {
  from <- match.arg(from); to <- match.arg(to)
  stopifnot(is.numeric(terms), !is.null(names(terms)))
  if (from == to) return(terms)
  # log10(x) = ln(x) / ln(10): a coefficient fitted against ln(x) is rescaled
  # by ln(10) to be used against log10(x), and vice versa.
  factor <- if (from == "ln") log(10) else 1 / log(10)
  out <- terms * factor
  attr(out, "note") <- sprintf("coefficients rescaled from %s to %s (x %.6f)", from, to, factor)
  out
}
