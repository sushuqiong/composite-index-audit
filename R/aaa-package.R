#' CompositeIndexAudit: pre-fitting algebraic classification of composite
#' indices
#'
#' Composite indices are routinely entered into a regression alongside the very
#' measurements they are computed from.  Whether the index coefficient then
#' estimates anything is a property of the pair (index, base model) and can be
#' decided before fitting, from the registered expression alone.
#'
#' @section Entry points:
#' \itemize{
#'   \item \code{\link{audit_fit}} -- audit an index against an already-fitted
#'         model object.
#'   \item \code{\link{audit_registered_index}} -- the same, with the base
#'         model columns supplied by hand.
#'   \item \code{\link{parse_index_expression}} and
#'         \code{\link{register_index}} -- the symbolic front end.
#'   \item \code{\link{verify_index_column}} -- recompute a stored index column
#'         and catch unit, scale and transcription mistakes.
#' }
#'
#' @section The four classes:
#' \code{rank_preserving_exact}, \code{rank_adding_exact},
#' \code{reducible_to_interaction} and \code{outside_supported_class}.  The
#' audit can also return an explicit non-classification state; see
#' \code{\link{audit_state_type}}.
#'
#' @docType package
#' @name CompositeIndexAudit-package
#' @keywords internal
#' @importFrom stats setNames cor complete.cases
"_PACKAGE"
