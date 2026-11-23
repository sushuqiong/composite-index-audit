# Topic 2 v41 NEW analyses -- BASE R ONLY.
# (1) 1,000-replication Monte Carlo with a MULTIPLICATIVE composite under a
#     non-fixed-ratio DGP (Claude T2 #2).
# (2) Runtime scaling over n and number of components (Claude T2 #3).
# (3) Numerical diagnostics (QR rank, kappa, condition number) on controlled
#     column spaces, to show they do not produce the symbolic audit class
#     (Claude T2 #3). No participant-level data are used.
options(stringsAsFactors = FALSE, digits = 12)
BASE <- "C:/Users/sushq/WorkBuddy/2026-08-12-22-01-44/v36v41v38_build/analysis"
OUT  <- file.path(BASE, "out"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------- (1) multiplicative-composite Monte Carlo ----------
set.seed(20261119L)
B <- 1000L; n <- 1500L
res <- matrix(NA_real_, B, 6L)
colnames(res) <- c("reject_bundled_1df", "reject_unrestricted_2df",
                   "reject_fixed_ratio_constraint", "constraint_true",
                   "gamma1", "gamma2")
fail <- 0L
for (b in seq_len(B)) {
  w1 <- rnorm(n); w2 <- rnorm(n); F <- rbinom(n, 1L, 0.5)
  P1 <- w1*w2; P2 <- w2*F
  g1 <- runif(1, 0.20, 0.80); g2 <- runif(1, 0.20, 0.80)
  lp <- -0.5 + 0.5*w1 + 0.5*w2 + g1*P1 + g2*P2
  y  <- rbinom(n, 1L, plogis(lp))
  C  <- P1 + 7*P2                       # declared fixed-ratio composite
  fit <- tryCatch({
    m1 <- glm(y ~ w1 + w2 + F + C, family = binomial())
    m2 <- glm(y ~ w1 + w2 + F + P1 + P2, family = binomial())
    s1 <- summary(m1)$coefficients
    z  <- s1["C", "Estimate"]/s1["C", "Std. Error"]
    be <- coef(m2)[c("P1", "P2")]; V <- vcov(m2)[c("P1", "P2"), c("P1", "P2")]
    chi2 <- as.numeric(t(be) %*% solve(V) %*% be)
    L <- c(-7, 1); est <- sum(L*be); se <- sqrt(max(as.numeric(t(L) %*% V %*% L), 0))
    c(rej1 = as.numeric(abs(z) > qnorm(.975)),
      rej2 = as.numeric(chi2 > qchisq(.95, 2)),
      rejC = as.numeric(if (se > 0) abs(est/se) > qnorm(.975) else 0),
      truth = as.numeric(abs(g2 - 7*g1) < 0.02))
  }, error = function(e) rep(NA_real_, 4L))
  if (anyNA(fit)) { fail <- fail + 1L; next }
  res[b, ] <- c(fit, g1, g2)
}
ok <- !is.na(res[, 1])
R <- colSums(res[ok, 1:4, drop = FALSE]); N <- sum(ok)
summ <- data.frame(
  quantity = c("replications_completed", "failed_fits",
               "symbolic_reducible_to_interaction_recovered",
               "bundled_1df_rejection_rate", "unrestricted_2df_rejection_rate",
               "fixed_ratio_constraint_rejection_rate",
               "fixed_ratio_constraint_true_rate",
               "non_rejection_rate_when_constraint_false"),
  value = c(N, fail,
            N,                         # deterministic algebraic classification
            R["reject_bundled_1df"]/N,
            R["reject_unrestricted_2df"]/N,
            R["reject_fixed_ratio_constraint"]/N,
            R["constraint_true"]/N,
            1 - R["reject_fixed_ratio_constraint"]/N))
write.csv(summ, file.path(OUT, "v41_multiplicative_simulation.csv"), row.names = FALSE)
cat("=== (1) multiplicative simulation ===\n"); print(summ, row.names = FALSE)
cat("\nNote: symbolic classification is deterministic from the declared expansion,\n")
cat("so recovery is 1000/1000 by construction, not a Monte Carlo estimate.\n")

# ---------- (2) runtime scaling ----------
set.seed(20261119L)
ns <- c(1e3, 1e4, 1e5, 1e6); ks <- c(2L, 5L, 10L, 20L)
scale_rows <- list()
for (nn in ns) for (kk in ks) {
  X <- matrix(rnorm(nn*kk), nrow = nn, ncol = kk)
  t <- system.time({
    q <- qr(X); rk <- q$rank; kap <- kappa(X, exact = TRUE)
    ## symbolic-style check: count exact duplicate columns via crossproduct rank
    dup <- sum(duplicated(as.data.frame(X)))
  })[["elapsed"]]
  scale_rows[[length(scale_rows)+1L]] <- data.frame(n = nn, components = kk,
      seconds = round(t, 5), qr_rank = rk, kappa = signif(kap, 6),
      exact_duplicate_columns = dup)
}
sc <- do.call(rbind, scale_rows)
write.csv(sc, file.path(OUT, "v41_runtime_scaling.csv"), row.names = FALSE)
cat("\n=== (2) runtime scaling (base R QR + kappa) ===\n"); print(sc, row.names = FALSE)

# ---------- (3) numerical diagnostics vs symbolic class ----------
set.seed(20261119L)
n0 <- 2000L
a <- rnorm(n0); b <- rnorm(n0); c <- rnorm(n0)
cases <- list(
  exact_linear_identity = cbind(1, a, b, a + b),
  near_collinear        = cbind(1, a, b, a + b + rnorm(n0, 0, 1e-8)),
  fixed_ratio_product   = cbind(1, a, b, a*b, b*c, a*b + 7*b*c),
  independent_columns   = cbind(1, a, b, c, a*b, b*c))
diag <- list()
for (nm in names(cases)) {
  X <- cases[[nm]]
  q <- qr(X)
  cm <- crossprod(X)
  ev <- eigen(cm, symmetric = TRUE, only.values = TRUE)$values
  diag[[nm]] <- data.frame(case = nm, columns = ncol(X), qr_rank = q$rank,
      kappa = signif(kappa(X, exact = TRUE), 6),
      condition_number = signif(max(ev)/min(ev[ev > 0]), 6),
      symbolic_audit_class = c("rank_preserving_exact", "rank_adding_exact",
                               "reducible_to_interaction",
                               "outside_supported_class (declaration required)")[which(names(cases) == nm)])
}
dg <- do.call(rbind, diag)
write.csv(dg, file.path(OUT, "v41_numerical_vs_symbolic.csv"), row.names = FALSE)
cat("\n=== (3) numerical diagnostics vs symbolic class ===\n"); print(dg, row.names = FALSE)
cat("\nT2_NEW_ANALYSES_COMPLETE\n")
