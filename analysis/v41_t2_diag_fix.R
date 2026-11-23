# Topic 2 v41 -- recompute the numeric-vs-symbolic diagnostic table CLEANLY.
# The earlier condition_number used max(ev)/min(ev[ev>0]) on crossprod(X), which
# for rank-deficient X discards the true near-zero eigenvalue and returns a
# spuriously small number. Here condition_number is the standard 2-norm ratio
# of singular values (equal to kappa(X, exact=TRUE)); rank deficiency is explicit.
options(stringsAsFactors = FALSE, digits = 12)
OUT <- "C:/Users/sushq/WorkBuddy/2026-08-12-22-01-44/v36v41v38_build/analysis/out"
set.seed(20261119L)
n0 <- 2000L
a <- rnorm(n0); b <- rnorm(n0); c <- rnorm(n0)
cases <- list(
  exact_linear_identity = cbind(1, a, b, a + b),
  near_collinear        = cbind(1, a, b, a + b + rnorm(n0, 0, 1e-8)),
  fixed_ratio_product   = cbind(1, a, b, a*b, b*c, a*b + 7*b*c),
  independent_columns   = cbind(1, a, b, c, a*b, b*c))
classes <- c("rank_preserving_exact", "rank_adding_exact",
             "reducible_to_interaction",
             "outside_supported_class (declaration required)")
dg <- do.call(rbind, lapply(seq_along(cases), function(ii) {
  X <- cases[[ii]]; nm <- names(cases)[ii]
  q <- qr(X); d <- svd(X)$d
  data.frame(case = nm, columns = ncol(X), qr_rank = q$rank,
             rank_deficiency = ncol(X) - q$rank,
             condition_number_2norm = signif(d[1]/d[length(d)], 6),
             symbolic_audit_class = classes[ii])
}))
write.csv(dg, file.path(OUT, "v41_numerical_vs_symbolic.csv"), row.names = FALSE)
print(dg, row.names = FALSE)
cat("T2_DIAG_FIX_COMPLETE\n")
