test_that("sample degeneration cannot change the symbolic class", {
  d <- data.frame(x = rep(1, 6), y = rep(2, 6), z = rep(0, 6))
  a <- audit_index_in_model(c(x = 1, y = 1, z = 1), c(x = 1, y = 1),
                            c(x = 1, y = 1), data = d)
  expect_equal(a$classification, "rank_adding_exact")
  expect_equal(a$observed_rank_gain, 0L)
  expect_equal(a$verification_status, "sample_degenerate_or_rank_mismatch")
  expect_true("symbolic_numeric_rank_mismatch" %in% a$failure_states)
})

test_that("strict metadata matching withholds exact certification", {
  register_index("cia_v35_sum", "x + y", units = c(x = "u", y = "u"), overwrite = TRUE)
  a <- suppressWarnings(audit_registered_index("cia_v35_sum", c("x", "y")))
  expect_equal(a$classification, "not_certified")
  expect_null(a$coefficient_map)
  keys <- c(x = "x|identity|-|u|-", y = "y|identity|-|u|-")
  b <- audit_registered_index("cia_v35_sum", keys)
  expect_equal(b$classification, "rank_preserving_exact")
})

test_that("an interaction already in the base is not counted as new", {
  a <- audit_product_in_model("xy", c("x", "y", "x:y"),
                              list(products = c("x:y" = 1)))
  expect_equal(a$classification, "rank_preserving_exact")
  expect_equal(a$expected_rank_gain, 0L)
  expect_length(a$effective_expansion$products, 0L)
  expect_equal(unname(a$coefficient_map["x:y", "C"]), 1)
})

test_that("a product expression supplied to audit_fit uses product dispatch", {
  set.seed(3501)
  d <- data.frame(x = runif(60), z = runif(60))
  d$dx <- rbinom(60, 1, plogis(-0.4 + 0.5 * d$x - 0.3 * d$z))
  f <- glm(dx ~ x + z, data = d, family = binomial())
  keys <- c(x = "x|identity|-|u|-", z = "z|identity|-|u|-")
  a <- audit_fit(f, "xz", expression = "x*z", units = c(x = "u", z = "u"),
                 model_terms = keys)
  expect_equal(a$classification, "reducible_to_interaction")
  expect_false(is.null(a$coefficient_map))
})
