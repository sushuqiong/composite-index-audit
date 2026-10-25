# reclassify_v34.R (PURE ASCII)
#
# v34: split outside_supported_class into the two things it was conflating.
#
# A review pointed out that "outside the supported class" currently covers two
# very different situations, and that treating them as one boundary is an
# overstatement:
#
#   outside_grammar              the registered expression cannot be written
#                                in the supported grammar at all (a ratio of
#                                two columns, a log of a sum, fractional
#                                powers).  This is a limitation of the parser
#                                and is fixable in future work; it holds for
#                                every base model.
#   outside_no_registered_expr   the registry stores a descriptive label
#                                rather than a closed-form expression, so
#                                there is nothing for the algebra to act on.
#   component_not_in_model       the expression DOES have a known
#                                decomposition, but this particular base
#                                model lacks at least one required column.
#                                This is an algebraic statement about the
#                                pair, not a boundary of the method.
#
# The coarse four-class label is preserved unchanged, and the script asserts
# that collapsing the refined label recovers it, so no existing result moves.

suppressPackageStartupMessages(library(haven))
CD <- Sys.getenv("CIA_CODE_DIR", ".")
OUT <- Sys.getenv("CIA_OUT", ".")
if (!nzchar(CD) || !nzchar(OUT))
  stop("set CIA_CODE_DIR and CIA_OUT before running reclassify_v34.R")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
for (f in c("index_registry.R", "term_identity.R", "algebra_audit.R",
            "public_api.R", "registry_audit.R", "audit_fit.R"))
  source(file.path(CD, f))

# Same rule build_results.R uses: the class of an (index, base) pair, computed
# by the audit functions and never read from the reference labels.
audit_class_for_base <- function(index, base_terms) {
  exps <- product_expansions()
  if (index %in% names(exps))
    return(audit_product_in_model(index, base_terms, exps[[index]])$classification)
  t <- derive_additive_terms(index)
  if (is.null(t)) return("outside_supported_class")
  emb <- t$index_terms[intersect(names(t$index_terms), base_terms)]
  audit_index_in_model(t$index_terms, emb, emb, intercept = t$intercept)$classification
}

# Why this index is outside, for ANY base.
grammar_status <- function(index) {
  exps <- product_expansions()
  if (index %in% names(exps)) return("expansion_supplied")
  if (!is.null(derive_additive_terms(index))) return("additive_decomposition_supplied")
  reg <- v28_registry()
  if (!index %in% reg$index) return("not_in_registry")
  f <- reg$formula[reg$index == index]
  parsed <- tryCatch(parse_index_expression(f)$status, error = function(e) "parse_error")
  if (identical(parsed, "unsupported")) return("parser_refuses_expression")
  # The parser accepted it as a single symbol: the registry stores a label,
  # not a formula.
  "no_closed_form_expression"
}

refine <- function(index, coarse) {
  if (!identical(coarse, "outside_supported_class")) return(coarse)
  gs <- grammar_status(index)
  if (gs %in% c("expansion_supplied", "additive_decomposition_supplied"))
    return("component_not_in_model")
  if (gs == "parser_refuses_expression") return("outside_grammar")
  if (gs == "no_closed_form_expression") return("outside_no_registered_expression")
  "outside_grammar"
}

collapse <- function(refined) {
  ifelse(refined %in% c("outside_grammar", "outside_no_registered_expression",
                        "component_not_in_model"),
         "outside_supported_class", refined)
}

reg <- v28_registry()
MINIMAL_BASE <- c("sex", "log_age", "log_bmi")
rows <- list()
for (i in seq_len(nrow(reg))) {
  idx <- reg$index[i]
  for (b in c("registered_component_base", "minimal_base_W_B")) {
    base <- if (b == "registered_component_base") {
      tryCatch(v28_component_model_terms(idx, reg), error = function(e) character(0))
    } else MINIMAL_BASE
    if (!length(base)) next
    coarse <- audit_class_for_base(idx, base)
    refined <- refine(idx, coarse)
    stopifnot(identical(collapse(refined), coarse))
    rows[[length(rows) + 1L]] <- data.frame(
      index = idx, base = b,
      base_terms = paste(base, collapse = "; "),
      class_coarse = coarse, class_refined = refined,
      grammar_status = grammar_status(idx),
      stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, rows)

cat("\n--- refined classes, registered component base ---\n")
print(table(tab$class_refined[tab$base == "registered_component_base"]))
cat("\n--- refined classes, minimal base W_B ---\n")
print(table(tab$class_refined[tab$base == "minimal_base_W_B"]))

# The split must not move any four-class result.
stopifnot(all(collapse(tab$class_refined) == tab$class_coarse))
for (b in unique(tab$base)) {
  sub <- tab[tab$base == b, ]
  n_out <- sum(sub$class_coarse == "outside_supported_class")
  n_ref <- sum(sub$class_refined %in% c("outside_grammar",
                                        "outside_no_registered_expression",
                                        "component_not_in_model"))
  stopifnot(n_out == n_ref)
}
cat("\nreclassify: coarse labels preserved; outside split into",
    paste(sort(unique(tab$class_refined[tab$class_coarse == "outside_supported_class"])),
          collapse = ", "), "\n")
write.csv(tab, file.path(OUT, "v34_reclassified_classes.csv"), row.names = FALSE)
cat("wrote v34_reclassified_classes.csv\n")
