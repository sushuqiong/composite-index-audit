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
  a <- suppressWarnings(audit_fit(f, "TyG"))
  expect_equal(a$classification, "rank_adding_exact")
  expect_equal(a$model_object$base_model_terms, c("sex", "log_age", "log_bmi"))
  expect_equal(a$model_object$family, "binomial")
  expect_equal(a$model_object$link, "logit")
  expect_true(isTRUE(a$is_classification))
})

test_that("the same index is rank preserving against its component base", {
  f <- glm(dx ~ sex + log_tg_mg + log_glucose_mg, data = d, family = binomial())
  a <- suppressWarnings(audit_fit(f, "TyG"))
  expect_equal(a$classification, "rank_preserving_exact")
  expect_false(is.null(a$coefficient_map))
})

test_that("audit_fit accepts an lm fit", {
  f <- lm(log_bmi ~ sex + log_tg_mg + log_glucose_mg, data = d)
  a <- suppressWarnings(audit_fit(f, "TyG"))
  expect_equal(a$classification, "rank_preserving_exact")
})

test_that("an unsupported model class is refused, not silently audited", {
  f <- glm(dx ~ sex + log_age, data = d, family = binomial())
  class(f) <- "somethingelse"
  expect_error(audit_fit(f, "TyG"), "unsupported model class")
})

test_that("an offset is reported as a confound and never downgraded", {
  f <- glm(dx ~ sex + log_age + log_bmi, data = d, family = binomial(),
           offset = rep(0.1, n))
  a <- suppressWarnings(audit_fit(f, "TyG"))
  expect_true(isTRUE(a$model_object$offset$offset_confound))
  expect_true("offset_present" %in% a$failure_states)
})

test_that("a fitted-model audit agrees with the registry audit", {
  ref <- reference_audit_classes()
  f <- glm(dx ~ sex + log_tg_mg + log_glucose_mg, data = d, family = binomial())
  a <- suppressWarnings(audit_fit(f, "TyG"))
  expect_equal(a$classification, unname(ref[["TyG"]]))
})
