test_that("the registry carries a reference class for every row", {
  reg <- index_registry()
  expect_true(nrow(reg) > 0L)
  expect_true(all(c("index", "reference_class", "provenance") %in% names(reg)))
  expect_true(all(reg$index %in% names(reference_audit_classes())))
  expect_false(any(is.na(reg$reference_class)))
})

test_that("there are exactly four classes and they are the documented ones", {
  expect_setequal(
    unique(unname(reference_audit_classes())),
    c("rank_preserving_exact", "rank_adding_exact",
      "reducible_to_interaction", "outside_supported_class"))
})

test_that("GLM7 expands to seven routine measures", {
  e <- glm7_full_expression()
  expect_equal(lengths(regmatches(e, gregexpr("log10\\(", e))), 7L)
  for (m in c("glucose", "insulin", "tg", "LDL", "hdl", "age", "BMI"))
    expect_true(grepl(m, e, fixed = TRUE))
})

test_that("product expansions are all named and carry linear parts", {
  exps <- product_expansions()
  expect_true(length(exps) > 0L)
  expect_true(all(nzchar(names(exps))))
  expect_true(is.list(exps))
  expect_true(all(vapply(exps, function(e) is.list(e), logical(1))))
})
