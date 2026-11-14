% Generated from the v35 roxygen comments; roxygen2 was not run in this environment.
\name{audit_fit}
\alias{audit_fit}
\title{Audit a composite index against an already-fitted model}
\usage{
audit_fit(fit, index, expression = NULL, data = NULL, model_terms = NULL,
  tol = 1e-10, units = character(0), strict = TRUE, drop_index = TRUE,
  keep_aliased = FALSE)
}
\arguments{
\item{fit}{A fitted model object of a supported class.}
\item{index}{Character scalar, a registered index name, or any name when \code{expression} is supplied.}
\item{expression}{Optional index expression, parsed and used instead of the registry.}
\item{data}{Optional data frame for numerical verification.}
\item{model_terms}{Optional character vector of design terms or named vector mapping columns to five-field term keys.}
\item{tol}{Numerical tolerance.}
\item{units}{Optional named character vector of component units for expression-based certification.}
\item{strict}{Logical; if \code{TRUE} (default), incomplete metadata return \code{not_certified}. Built-in registry names require expression and units for a certified fitted-object audit.}
\item{drop_index}{Drop the index itself from the base model terms.}
\item{keep_aliased}{Keep aliased columns among the base model terms.}
}
\value{The audit list plus model-object metadata, state type and classification indicator.}
\description{Fitted-object convenience interface. It reads design terms, aliasing, family, offset and weights without refitting or using the outcome for classification. For an audit before fitting, use \code{\link{audit_registered_index}} or \code{\link{audit_index_in_model}}.}
\details{Supported classes are \code{glm}, \code{lm}, \code{svyglm}, \code{coxph} and \code{survreg}. Offsets and aliased columns are reported explicitly. Term identity is never inferred from a column name.}
\examples{
d <- data.frame(sex = factor(c("M", "F", "M", "F")),
                log_tg_mg = c(4.6, 4.8, 4.7, 5.0),
                log_glucose_mg = c(4.3, 4.5, 4.4, 4.6),
                log_age = c(1.6, 1.7, 1.65, 1.66),
                log_bmi = c(1.4, 1.5, 1.45, 1.47),
                dx = c(0, 1, 0, 1))
fit <- glm(dx ~ sex + log_age + log_bmi, data = d, family = binomial())
audit_fit(fit, "TyG", strict = FALSE)$classification
}
