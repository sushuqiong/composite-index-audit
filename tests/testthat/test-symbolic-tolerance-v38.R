test_that("numeric tolerance never erases a nonzero symbolic coefficient", {
  a <- audit_index_in_model(c(x = 1, z = 1e-11), c(x = 1), data = NULL)
  expect_identical(a$classification, "rank_adding_exact")
  expect_identical(a$expected_rank_gain, 1L)
  expect_identical(a$residual_expression, "1e-11 x z")
  b <- audit_composite_index(c(x = 1, z = 1e-11), c(x = 1), data = NULL)
  expect_identical(b$status, "rank_adding_exact")
  expect_identical(b$residual_terms, c(z = 1e-11))
})

test_that("declared coefficient mismatch is not hidden by tolerance", {
  a <- audit_composite_index(c(x = 1), c(x = 1 + 1e-12), data = NULL)
  expect_false(a$exact_reparameterization)
  expect_identical(a$status, "not_decomposable")
  b <- audit_composite_index(c(x = 1, y = 2), c(x = 1, y = 2),
                             c(y = 2, x = 1), data = NULL)
  expect_identical(b$status, "rank_preserving_exact")
})

test_that("small nonzero residual preserves symbolic class despite rank diagnostic", {
  d <- data.frame(x = seq_len(50), z = seq_len(50))
  a <- audit_index_in_model(c(x = 1, z = 1e-11), c(x = 1), data = d)
  expect_identical(a$classification, "rank_adding_exact")
  expect_identical(a$expected_rank_gain, 1L)
  expect_identical(a$observed_rank_gain, 0L)
  expect_identical(a$verification_status, "sample_degenerate_or_rank_mismatch")
  expect_true("symbolic_numeric_rank_mismatch" %in% a$failure_states)
})

test_that("partial interactions and product offsets retain the intended maps", {
  d <- data.frame(sex = factor(rep(c("Male", "Female"), 50), levels = c("Male", "Female")),
                  x = seq_len(100), z = seq_len(100) + 0.5)
  d$LAP <- ifelse(d$sex == "Male", d$x - 65, d$x - 58) * d$z
  d$sexFemale <- as.numeric(d$sex == "Female")
  d$waist_cm <- d$x; d$tg_mmol_L <- d$z
  d$LAP <- ifelse(d$sex == "Male", d$waist_cm - 65, d$waist_cm - 58) * d$tg_mmol_L
  d[["waist_cm:tg_mmol_L"]] <- d$waist_cm * d$tg_mmol_L
  d[["tg_mmol_L:sexFemale"]] <- d$tg_mmol_L * d$sexFemale
  ex <- product_expansions()$LAP
  a <- audit_product_in_model("LAP", c("sex", "waist_cm", "tg_mmol_L", "tg_mmol_L:sexFemale"), ex, data = d)
  expect_identical(a$classification, "reducible_to_interaction")
  expect_identical(a$numerical_diagnostics$rank_gain, 1L)
  expect_identical(a$coefficient_map["tg_mmol_L", "C"], -65)
  b <- audit_product_in_model("LAP", c("sex", "waist_cm", "tg_mmol_L", "waist_cm:tg_mmol_L", "tg_mmol_L:sexFemale"), ex, data = d)
  expect_identical(b$classification, "rank_preserving_exact")
  expect_identical(b$numerical_diagnostics$rank_gain, 0L)
})
