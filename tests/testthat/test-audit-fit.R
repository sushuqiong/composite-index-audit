set.seed(20261025)
n <- 400
d <- data.frame(
  sex = factor(sample(c("Male", "Female"), n, TRUE)),
  log_age = rnorm(n), log_bmi = rnorm(n),
  log_tg_mg = rnorm(n) + 4.7, log_glucose_mg = rnorm(n) + 4.4,
  stringsAsFactors = FALSE)
d$TyG <- d$log_tg_mg + d$log_glucose_mg - log(2)
lp <- -6 + 0.9 * d$log_age + 0.7 * d$log_bmi + 0.4 * d$log_tg_mg
d$dx <- rbinom(n, 1, plogis(lp))

test_that("audit_fit reads a glm and reports base model TERMS", {
  f <- glm(dx ~ sex + log_age + log_bmi, data = d, family = binomial())
  a <- suppressWarnings(audit_fit(f, "TyG", strict = FALSE))
  expect_equal(a$classification, "rank_adding_exact")
  expect_equal(a$model_object$base_model_terms, c("sex", "log_age", "log_bmi"))
  expect_equal(a$model_object$family, "binomial")
  expect_equal(a$model_object$link, "logit")
  expect_true(isTRUE(a$is_classification))
})

test_that("the same index is rank preserving against its component base", {
  f <- glm(dx ~ sex + log_tg_mg + log_glucose_mg, data = d, family = binomial())
  a <- suppressWarnings(audit_fit(f, "TyG", strict = FALSE))
  expect_equal(a$classification, "rank_preserving_exact")
  expect_false(is.null(a$coefficient_map))
})

test_that("audit_fit accepts an lm fit", {
  f <- lm(log_bmi ~ sex + log_tg_mg + log_glucose_mg, data = d)
  a <- suppressWarnings(audit_fit(f, "TyG", strict = FALSE))
  expect_equal(a$classification, "rank_preserving_exact")
})

test_that("an unsupported model class is refused, not silently audited", {
  f <- glm(dx ~ sex + log_age, data = d, family = binomial())
  class(f) <- "somethingelse"
  expect_error(audit_fit(f, "TyG", strict = FALSE), "unsupported model class")
})

test_that("an offset is reported without presuming component redundancy", {
  f <- glm(dx ~ sex + log_age + log_bmi, data = d, family = binomial(),
           offset = rep(0.1, n))
  a <- suppressWarnings(audit_fit(f, "TyG", strict = FALSE))
  expect_false(isTRUE(a$model_object$offset$offset_confound))
  expect_true(isTRUE(a$model_object$offset$offset_review_required))
  expect_identical(a$model_object$offset$fixed_coefficient, 1)
  expect_match(a$model_object$offset$note, "does not prove redundancy")
  expect_true("offset_present" %in% a$failure_states)
  expect_identical(a$classification, "rank_adding_exact")
})

test_that("a fitted-model audit agrees with the registry audit", {
  ref <- reference_audit_classes()
  f <- glm(dx ~ sex + log_tg_mg + log_glucose_mg, data = d, family = binomial())
  a <- suppressWarnings(audit_fit(f, "TyG", strict = FALSE))
  expect_equal(a$classification, unname(ref[["TyG"]]))
})

# ---- v0.40.0 interface contract -------------------------------------------

test_that("v0.40.0: a strict call without a declared base is refused", {
  f <- glm(dx ~ sex + log_age + log_bmi, data = d, family = binomial())
  a <- audit_fit(f, "TyG")                       # strict defaults to TRUE
  expect_equal(a$classification, "not_certified")
  expect_equal(a$state_type, "non_classification_state")
  expect_equal(a$model_object$base_source, "missing")
  expect_true("missing_base_declaration" %in% a$failure_states)
  expect_false(isTRUE(a$is_classification))
  expect_match(a$blocking_reason, "missing_base_declaration")
})

test_that("v0.40.0: an exploratory call tags an inferred base", {
  f <- glm(dx ~ sex + log_age + log_bmi, data = d, family = binomial())
  a <- suppressWarnings(audit_fit(f, "TyG", strict = FALSE))
  expect_equal(a$model_object$base_source, "inferred")
  expect_false(is.na(a$classification))
})

test_that("v0.40.0: a declared common base is invariant to formula order", {
  d$S <- d$log_tg_mg + d$log_glucose_mg
  keys <- c(log_tg_mg = "log_tg_mg|identity|-|mg|-",
            log_glucose_mg = "log_glucose_mg|identity|-|mg|-",
            sex = "sex|identity|-|-|factor(Male,Female)")
  u <- c(log_tg_mg = "mg", log_glucose_mg = "mg")
  f1 <- glm(dx ~ sex + log_tg_mg + log_glucose_mg, data = d, family = binomial())
  f2 <- glm(dx ~ log_glucose_mg + sex + log_tg_mg, data = d, family = binomial())
  a1 <- audit_fit(f1, "S", expression = "log_tg_mg + log_glucose_mg",
                  units = u, model_terms = keys)
  a2 <- audit_fit(f2, "S", expression = "log_tg_mg + log_glucose_mg",
                  units = u, model_terms = keys)
  expect_equal(a1$classification, "rank_preserving_exact")
  expect_equal(a1$classification, a2$classification)
  expect_equal(a1$model_object$base_source, "declared")
})

test_that("v0.40.0: a tampered supplied column keeps the class and withholds certification", {
  dd <- data.frame(x = 1:6, y = c(2, 4, 6, 8, 10, 12))
  dd$C <- dd$x + dd$y
  ok <- verify_index_column(dd, "C", expression = "x + y")
  dd$C[[3]] <- dd$C[[3]] + 0.5
  bad <- verify_index_column(dd, "C", expression = "x + y")
  expect_equal(ok$input_certification, "certified")
  expect_equal(bad$input_certification, "withheld")
  expect_equal(bad$verification_status, "input_column_mismatch")
})
