# Minimal base-R runner for this package's test_that()/expect_*() test files.
# It executes every test block and expectation without loading optional packages.
.cia_harness <- new.env(parent = emptyenv())
.cia_harness$pass <- 0L
.cia_harness$fail <- 0L
.cia_harness$blocks <- 0L

.cia_pass <- function() {
  .cia_harness$pass <- .cia_harness$pass + 1L
  invisible(TRUE)
}

.cia_fail <- function(message) {
  stop(message, call. = FALSE)
}

test_that <- function(desc, code) {
  .cia_harness$blocks <- .cia_harness$blocks + 1L
  expr <- substitute(code)
  tryCatch({
    eval(expr, envir = parent.frame())
    cat("PASS:", desc, "\n")
  }, error = function(e) {
    .cia_harness$fail <- .cia_harness$fail + 1L
    cat("FAIL:", desc, "::", conditionMessage(e), "\n")
  })
  invisible(NULL)
}

expect_true <- function(object, ...) {
  if (!isTRUE(object)) .cia_fail("expected TRUE")
  .cia_pass()
}

expect_false <- function(object, ...) {
  if (!identical(object, FALSE)) .cia_fail("expected FALSE")
  .cia_pass()
}

expect_null <- function(object, ...) {
  if (!is.null(object)) .cia_fail("expected NULL")
  .cia_pass()
}

expect_length <- function(object, n, ...) {
  if (length(object) != n) .cia_fail(paste0("expected length ", n, ", got ", length(object)))
  .cia_pass()
}

expect_equal <- function(object, expected, ..., tolerance = 1e-8) {
  comparison <- all.equal(object, expected, tolerance = tolerance, check.attributes = TRUE)
  if (!isTRUE(comparison)) .cia_fail(paste(comparison, collapse = "; "))
  .cia_pass()
}

expect_setequal <- function(object, expected, ...) {
  a <- unique(object)
  b <- unique(expected)
  same <- length(a) == length(b) && all(a %in% b) && all(b %in% a)
  if (!same) .cia_fail("sets differ")
  .cia_pass()
}

expect_error <- function(object, regexp = NULL, ...) {
  expr <- substitute(object)
  error <- tryCatch({
    eval(expr, envir = parent.frame())
    NULL
  }, error = function(e) e)
  if (is.null(error)) .cia_fail("expected an error, but none occurred")
  if (!is.null(regexp) && !grepl(regexp, conditionMessage(error), ...)) {
    .cia_fail(paste0("error message did not match '", regexp, "': ", conditionMessage(error)))
  }
  .cia_pass()
}
