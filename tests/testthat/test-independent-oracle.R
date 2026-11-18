# Independent matrix-algebra oracles: expected values are constructed directly,
# never generated with the package's coefficient-map or decomposition routines.
set.seed(22222026)
n <- 300L
x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
eta <- -0.4 + 0.8*x1 - 0.55*x2 + 0.35*x3
y <- rbinom(n, 1L, plogis(eta))
d <- data.frame(x1=x1, x2=x2, x3=x3, y=y)

test_that("additive E and residual R have an independently verified coefficient map", {
  d$C <- 2*x1 - x2 + x3
  d$R <- x3
  fit_C <- glm(y ~ x1 + x2 + C, data=d, family=binomial())
  fit_R <- glm(y ~ x1 + x2 + R, data=d, family=binomial())
  b <- coef(fit_C)
  expected <- c(b['(Intercept)'], b['x1']+2*b['C'], b['x2']-b['C'], b['C'])
  names(expected) <- c('(Intercept)', 'x1', 'x2', 'R')
  expect_equal(coef(fit_R)[names(expected)], expected, tolerance=1e-8)
  expect_equal(unname(fitted(fit_C)), unname(fitted(fit_R)), tolerance=1e-8)
  expect_equal(logLik(fit_R), logLik(fit_C), tolerance=1e-8)
  X_E <- model.matrix(~ x1 + x2, data=d)
  expect_equal(qr(cbind(X_E, C=d$C))$rank, qr(cbind(X_E, R=d$R))$rank)
  expect_equal(max(abs(d$C - (2*x1-x2+d$R))), 0, tolerance=1e-12)
})

test_that("all embedded components create exact aliasing and zero rank gain", {
  d$C <- 2*x1 - x2
  X0 <- model.matrix(~ x1 + x2, data=d)
  X1 <- cbind(X0, C=d$C)
  independent_residual <- qr.resid(qr(X0), d$C)
  expect_equal(max(abs(independent_residual)), 0, tolerance=1e-10)
  expect_equal(qr(X1)$rank-qr(X0)$rank, 0L)
  a <- audit_index_in_model(c(x1=2, x2=-1), c(x1=2, x2=-1), data=d)
  expect_identical(a$classification, 'rank_preserving_exact')
  expect_identical(a$observed_rank_gain, 0L)
  f <- glm(y ~ x1 + x2 + C, data=d, family=binomial())
  expect_true(is.na(coef(f)['C']) || isTRUE(all.equal(unname(fitted(f)), unname(glm.fit(X0,y,family=binomial())$fitted.values), tolerance=1e-8)))
})

test_that("a constant index shift maps to the independently calculated intercept", {
  constant <- -1.75
  d$C <- 2*x1 - x2 + x3 + constant
  d$R <- x3
  a <- audit_index_in_model(c(x1=2, x2=-1, x3=1), c(x1=2, x2=-1),
                            data=d, intercept=constant)
  beta <- c('(Intercept)'=.7, x1=.2, x2=-.4, C=.8)
  gamma <- as.vector(a$coefficient_map %*% beta)
  names(gamma) <- rownames(a$coefficient_map)
  independent_lp_C <- beta['(Intercept)'] + beta['x1']*x1 + beta['x2']*x2 + beta['C']*d$C
  independent_lp_R <- (beta['(Intercept)']+constant*beta['C']) +
    (beta['x1']+2*beta['C'])*x1 + (beta['x2']-beta['C'])*x2 + beta['C']*x3
  expect_equal(unname(gamma['(Intercept)']), unname(beta['(Intercept)']+constant*beta['C']), tolerance=1e-12)
  expect_equal(max(abs(independent_lp_C-independent_lp_R)), 0, tolerance=1e-12)
})

test_that("native product and explicit interaction have identical independent fits", {
  d$product <- (2*x1-1.5)*(x2+0.25)
  expansion <- list(intercept=-.375, linear=c(x1=.5, x2=-1.5), products=c('x1:x2'=2))
  fit_C <- glm(y ~ x1 + x2 + product, data=d, family=binomial())
  fit_I <- glm(y ~ x1*x2, data=d, family=binomial())
  b <- coef(fit_C)
  expected <- c(b['(Intercept)']-.375*b['product'], b['x1']+.5*b['product'],
                b['x2']-1.5*b['product'], 2*b['product'])
  names(expected) <- c('(Intercept)', 'x1', 'x2', 'x1:x2')
  expect_equal(coef(fit_I)[names(expected)], expected, tolerance=1e-8)
  expect_equal(unname(fitted(fit_C)), unname(fitted(fit_I)), tolerance=1e-8)
  X_I <- model.matrix(~ x1*x2, data=d)
  expect_equal(max(abs(qr.resid(qr(X_I), d$product))), 0, tolerance=1e-8)
  a <- audit_product_in_model('product', c('x1','x2'), expansion, data=d)
  expect_identical(a$classification, 'reducible_to_interaction')
  expect_identical(a$numerical_diagnostics$rank_gain, 1L)
})

test_that("random positive definite covariance gives an exact orthogonal decomposition", {
  for (i in seq_len(12)) {
    A <- matrix(rnorm(36), nrow=6L)
    V <- crossprod(A)+diag(.5, 6L)
    fixed <- c(1, rep(0,5))
    b <- rnorm(6)
    ee <- eigen(V, symmetric=TRUE)
    Wi <- ee$vectors %*% (t(ee$vectors)/sqrt(ee$values))
    Vi <- solve(V)
    precision_c <- as.numeric(crossprod(fixed, Vi %*% fixed))
    whitened_c <- as.vector(Wi %*% fixed)
    Q <- qr.Q(qr(t(matrix(whitened_c, nrow=1))), complete=TRUE)
    D <- t(Q[,2:6,drop=FALSE]) %*% Wi
    direction <- as.vector(Vi %*% fixed/sqrt(precision_c))
    q_constraint <- sum((D %*% b)^2)
    q_direction <- as.numeric(crossprod(direction,b))^2
    q_total <- as.numeric(crossprod(b,Vi %*% b))
    expect_equal(max(abs(D %*% fixed)), 0, tolerance=1e-10)
    expect_equal(max(abs(D %*% V %*% t(D)-diag(5))), 0, tolerance=1e-10)
    expect_equal(q_constraint+q_direction, q_total, tolerance=1e-10)
    expect_equal(q_direction, as.numeric(crossprod(fixed,Vi %*% b))^2/precision_c, tolerance=1e-10)
  }
})
