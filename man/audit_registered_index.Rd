% Generated from the v35 roxygen comments; roxygen2 was not run in this environment.
\name{audit_registered_index}
\alias{audit_registered_index}
\title{Audit a registered index against a stated base model}
\usage{
audit_registered_index(name, model_terms, data = NULL, tol = 1e-10, strict = TRUE)
}
\arguments{
\item{name}{Character scalar, a registered index name.}
\item{model_terms}{Character vector of base model columns, or a named character vector mapping column name to a five-field term key \code{variable|transform|base|units|encoding}. In strict mode, unnamed or incomplete metadata return \code{not_certified}.}
\item{data}{Optional data frame for the numerical verification layer.}
\item{tol}{Numerical tolerance for algebraic comparisons.}
\item{strict}{Logical; if \code{TRUE} (default), incomplete metadata withhold exact certification.}
}
\value{A list containing classification, coefficient map, blocking reason, residual expression, numerical diagnostics and failure states.}
\description{Pre-fitting audit with the base model terms supplied by the caller. \code{\link{audit_fit}} is the fitted-object convenience interface.}
\examples{
register_index("cia_demo_tyg", "log_tg_mg + log_glucose_mg",
               units = c(log_tg_mg = "ln(mg/dL)", log_glucose_mg = "ln(mg/dL)"))
audit_registered_index("cia_demo_tyg",
  c(log_tg_mg = "log_tg_mg|identity|-|ln(mg/dL)|-",
    log_glucose_mg = "log_glucose_mg|identity|-|ln(mg/dL)|-"))$classification
}
