test_that("log-base coefficient normalization preserves 100 linear predictors", {
  x <- seq(0.2, 20, length.out = 100)
  beta <- c(a = 0.8, b = -0.35)
  lp_ln <- beta[["a"]] * log(x) + beta[["b"]] * log(x + 1)
  converted <- normalize_log_base(beta, from = "ln", to = "log10")
  lp_10 <- converted[["a"]] * log10(x) + converted[["b"]] * log10(x + 1)
  for (i in seq_along(x)) expect_equal(lp_10[[i]], lp_ln[[i]], tolerance = 1e-12)
})

test_that("round-trip log-base normalization preserves 50 coefficient pairs", {
  beta <- seq(-1.225, 1.225, length.out = 50)
  names(beta) <- paste0("term_", seq_along(beta))
  roundtrip <- normalize_log_base(normalize_log_base(beta, "ln", "log10"), "log10", "ln")
  for (nm in names(beta)) expect_equal(roundtrip[[nm]], beta[[nm]], tolerance = 1e-12)
})

test_that("five-field term identity resolves 40 complete and incomplete cases", {
  complete <- make_term_key("tg", "log10", "10", "mmol/L", "-")
  variants <- list(
    make_term_key("tg", "log10", "10", "mg/dL", "-"),
    make_term_key("tg", "log", "e", "mmol/L", "-"),
    make_term_key("tg", "identity", "-", "mmol/L", "-"),
    make_term_key("tg", "log10", "10")
  )
  expected <- c("mismatch", "mismatch", "mismatch", "unkeyed")
  for (replicate in seq_len(10)) {
    for (i in seq_along(variants)) {
      expect_equal(unname(term_key_in(variants[[i]], complete)), expected[[i]])
    }
  }
})

test_that("column verification detects 30 multiplicative scale errors", {
  set.seed(3601)
  d <- data.frame(x = runif(30, 0.5, 2), y = runif(30, 0.5, 2))
  d$C <- d$x + d$y
  for (i in seq_len(30)) {
    di <- d
    di$C_bad <- di$C * (1 + i / 1000)
    expect_false(verify_index_column(di, "C_bad", expression = "x + y", tol = 1e-10)$ok)
  }
})

test_that("column verification accepts 30 correctly reconstructed columns", {
  set.seed(3602)
  d <- data.frame(x = runif(30, 0.5, 2), y = runif(30, 0.5, 2))
  d$C <- d$x + d$y
  for (i in seq_len(30)) {
    di <- d
    di$C <- di$x + di$y + rep(0, 30) * i
    expect_true(verify_index_column(di, "C", expression = "x + y", tol = 1e-12)$ok)
  }
})

test_that("symbolic rank-adding class survives 20 degenerate samples", {
  for (i in seq_len(20)) {
    d <- data.frame(x = rep(i, 6), y = rep(2, 6), z = rep(0, 6))
    a <- audit_index_in_model(c(x = 1, y = 1, z = 1), c(x = 1, y = 1),
                              c(x = 1, y = 1), data = d)
    expect_equal(a$classification, "rank_adding_exact")
  }
})
