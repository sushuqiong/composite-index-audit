test_that("an index whose components are all in the base is rank preserving", {
  a <- audit_index_in_model(c(a = 1, b = 1), c(a = 1, b = 1), c(a = 1, b = 1))
  expect_equal(a$classification, "rank_preserving_exact")
  expect_equal(audit_state_type(a$classification), "class")
  expect_false(is.null(a$coefficient_map))
})

test_that("an index whose components are missing from the base adds a rank", {
  a <- audit_index_in_model(c(a = 1, b = 1), c(a = 1), c(a = 1))
  expect_equal(a$classification, "rank_adding_exact")
})

test_that("the classification is one of four classes or an explicit state", {
  a <- audit_index_in_model(c(a = 1, b = 1), c(a = 1, b = 1), c(a = 1, b = 1))
  expect_equal(audit_state_type(a$classification), "class")
  expect_equal(audit_state_type("approximate_collinearity"),
               "non_classification_state")
  expect_equal(audit_state_type("fit_failed"), "non_classification_state")
  expect_equal(audit_state_type(NA_character_), "unknown")
})

test_that("a registered ad-hoc expression is parsed and checked", {
  register_index("cia_test_sum", "tg_mmol_L + waist_cm", overwrite = TRUE)
  a <- audit_registered_index("cia_test_sum", c("tg_mmol_L", "waist_cm"))
  expect_equal(a$classification, "rank_preserving_exact")
  a2 <- audit_registered_index("cia_test_sum", c("tg_mmol_L"))
  expect_equal(a2$classification, "rank_adding_exact")
})

test_that("verify_index_column reproduces a column from its expression", {
  d <- data.frame(x = c(1, 2, 3), y = c(4, 5, 7))
  d$s <- d$x + d$y
  v <- verify_index_column(d, "s", expression = "x + y")
  expect_true(isTRUE(v$ok))
})
