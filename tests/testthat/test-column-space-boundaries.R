# Independent design-matrix oracles for supported and unsupported column spaces.
# Expected ranks are calculated directly with model.matrix()/qr(), never with
# the package's decomposition or coefficient-map routines.
set.seed(39112026)
n <- 240L
x <- rnorm(n); y <- rnorm(n)
d <- data.frame(x=x, y=y, sum_xy=x+y, diff_xy=x-y, C=x, product=x*y,
                sex=rep(c("Female","Male"), each=n/2))

rank_gain <- function(W, C) qr(cbind(W, C=C))$rank - qr(W)$rank

test_that("X is not certified merely because X appears inside X+Y", {
  W <- model.matrix(~ sum_xy, data=d)
  expect_equal(rank_gain(W, d$C), 1L)
  # The registered-expression interface does not solve arbitrary linear
  # combinations of base columns. It therefore declines this certification.
  a <- audit_index_in_model(c(x=1), setNames(numeric(0), character(0)),
                            setNames(numeric(0), character(0)), data=d)
  expect_false(identical(a$classification, "rank_preserving_exact"))
})

test_that("X lies in span(1,X+Y,X-Y) by an independent QR oracle", {
  W <- model.matrix(~ sum_xy + diff_xy, data=d)
  expect_equal(rank_gain(W, d$C), 0L)
  expect_equal(max(abs(d$C - .5*d$sum_xy - .5*d$diff_xy)), 0, tolerance=1e-12)
  # This composite-base expression is outside the present symbolic interface;
  # numerical rank zero is not promoted to a symbolic certificate.
  a <- audit_index_in_model(c(x=1), c(sum_xy=.5,diff_xy=.5),
                            c(sum_xy=.5,diff_xy=.5), data=d)
  expect_true(a$classification %in% c("not_decomposable", "approximate_collinearity"))
})

test_that("no-intercept and categorical encodings use their actual columns", {
  W0 <- model.matrix(~ 0 + x + y, data=d)
  expect_equal(rank_gain(W0, 2*x-y), 0L)
  Wsex <- model.matrix(~ sex, data=d)
  female <- as.numeric(d$sex == "Female")
  expect_equal(rank_gain(Wsex, female), 0L)
  expect_equal(max(abs(female - (1-as.numeric(d$sex == "Male")))), 0, tolerance=1e-12)
})

test_that("an already embedded interaction adds zero but a new product adds one", {
  W_with <- model.matrix(~ x*y, data=d)
  W_without <- model.matrix(~ x+y, data=d)
  expect_equal(rank_gain(W_with, d$product), 0L)
  expect_equal(rank_gain(W_without, d$product), 1L)
})

test_that("centering and scaling are column-space preserving but coefficient mismatch is not certified", {
  shifted <- 3 + 2*(x-mean(x))
  W <- model.matrix(~ x, data=d)
  expect_equal(rank_gain(W, shifted), 0L)
  a <- audit_index_in_model(c(x=1), c(x=2), c(x=2), data=d)
  expect_true(a$classification %in% c("not_decomposable", "approximate_collinearity"))
})

test_that("input mismatch is distinct from symbolic class and model fitting", {
  dd <- data.frame(x=x, y=y, C=x+y+.01*rnorm(n))
  symbolic <- audit_index_in_model(c(x=1,y=1), c(x=1,y=1), data=dd[,c("x","y")])
  checked <- verify_index_column(dd, "C", expression="x+y", tol=1e-12)
  fit <- glm(rbinom(n,1,.2) ~ x+y+C, data=dd, family=binomial())
  expect_identical(symbolic$classification, "rank_preserving_exact")
  expect_identical(checked$input_certification, "withheld")
  expect_identical(checked$verification_status, "input_column_mismatch")
  expect_identical(checked$failure_states, "input_column_mismatch")
  expect_true(inherits(fit, "glm"))
})
