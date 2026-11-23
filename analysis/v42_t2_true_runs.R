# T2 v42 NEW analyses -- BASE R ONLY, source the package R/ files directly.
# (1) QR rank at three tolerances on the near-collinear matrix (Claude #2)
# (2) Column-order test: y ~ x1+x2+C vs y ~ C+x1+x2 (GPT #1)
# (3) Strict-mode product certification + input-column tampering (Claude #3, GPT #3)
# (4) Multiplicative simulation DGP re-verification (Claude #1 ADEMP)
# (5) Calibration intercept CI covariance decomposition (Claude a)
options(stringsAsFactors=FALSE, digits=12)
PK <- "C:/Users/sushq/WorkBuddy/2026-08-12-22-01-44/v37v42v39_build"
OUT <- "C:/Users/sushq/WorkBuddy/2026-08-12-22-01-44/v37v42v39_build/analysis/out"
for (f in c("index_registry.R","term_identity.R","algebra_audit.R","audit_fit.R","registry_audit.R","public_api.R"))
  source(file.path(PK,"rpkg_ascii",f))

cat("=== (1) QR rank at three tolerances ===\n")
set.seed(20261120L); n <- 1500L
x1 <- rnorm(n); x2 <- rnorm(n)
C_exact   <- x1 + x2
# near-collinear: C = x1+x2 + tiny noise -> condition number ~3e8 (per published)
eps <- 3e-9
C_near <- x1 + x2 + eps*rnorm(n)
M_exact <- cbind(x1,x2,C_exact); M_near <- cbind(x1,x2,C_near)
k2 <- function(X) { s<-svd(X)$d; max(s)/min(s) }
rows <- list()
mats <- list(exact_identity=M_exact, near_collinear=M_near)
for (lab in names(mats)) {
  M <- mats[[lab]]
  for (tol in c(1e-8,1e-10,1e-12)) {
    r <- qr(M, tol=tol)$rank
    rows[[length(rows)+1L]] <- data.frame(matrix=lab, tolerance=tol, qr_rank=r,
               condition_number_2norm=signif(k2(M),6))
    cat(sprintf("  %-16s tol=%.0e -> rank=%d  (2-norm kappa=%.3g)\n", lab, tol, r, k2(M)))
  }
}
qr_tab <- do.call(rbind, rows); write.csv(qr_tab, file.path(OUT,"v42_qr_tolerance_rank.csv"), row.names=FALSE)

cat("\n=== (2) Column-order test ===\n")
set.seed(20261120L); n<-1500L
w1<-rnorm(n); w2<-rnorm(n); C<-w1+w2
y <- rbinom(n,1,plogis(-0.5+0.5*w1+0.5*w2))
f1 <- glm(y ~ w1 + w2 + C, family=binomial())
f2 <- glm(y ~ C + w1 + w2, family=binomial())
a1 <- audit_fit(f1,"C",expression="w1 + w2",data=data.frame(y,w1,w2,C),strict=FALSE)
a2 <- audit_fit(f2,"C",expression="w1 + w2",data=data.frame(y,w1,w2,C),strict=FALSE)
rank1<-qr(model.matrix(f1))$rank; rank2<-qr(model.matrix(f2))$rank
base<-glm(y ~ w1 + w2, family=binomial()); rbase<-qr(model.matrix(base))$rank
cat(sprintf("  y~w1+w2+C : class=%s  rank=%d  aliased=%s  fitted values max diff=%.3g\n",
    a1$classification, rank1, paste(names(coef(f1))[is.na(coef(f1))],collapse=",") ,max(abs(fitted(f1)-fitted(f2)))))
cat(sprintf("  y~C+w1+w2 : class=%s  rank=%d  aliased=%s\n",
    a2$classification, rank2, paste(names(coef(f2))[is.na(coef(f2))],collapse=",")))
cat(sprintf("  base (no C) rank=%d ; increment both=%d\n", rbase, rank1-rbase))
cat(sprintf("  same symbolic class? %s ; same fitted values? %s\n",
    identical(a1$classification,a2$classification), isTRUE(all.equal(fitted(f1),fitted(f2)))))
ord <- data.frame(formula=c("y~w1+w2+C","y~C+w1+w2"),
  classification=c(a1$classification,a2$classification),
  qr_rank=c(rank1,rank2), rank_increment_vs_base=c(rank1-rbase,rank2-rbase),
  aliased_column=c(paste(names(coef(f1))[is.na(coef(f1))],collapse=","),
                   paste(names(coef(f2))[is.na(coef(f2))],collapse=",")),
  max_abs_fitted_diff=max(abs(fitted(f1)-fitted(f2))))
write.csv(ord, file.path(OUT,"v42_column_order_test.csv"), row.names=FALSE)

cat("\n=== (3) Strict certification + input-column tampering ===\n")
set.seed(20261120L); n<-600L
d <- data.frame(w1=rnorm(n), w2=rnorm(n))
d$C <- d$w1 + d$w2
register_index("DemoAdd","w1 + w2",units=c(w1="-",w2="-"),overwrite=TRUE)
v_ok <- verify_index_column(d,"C",expression="w1 + w2")
d2 <- d; d2$C[1] <- d2$C[1] + 0.5
v_bad <- verify_index_column(d2,"C",expression="w1 + w2")
cat(sprintf("  intact  : ok=%s max_abs_error=%.3g input_certification=%s verification_status=%s failure_states=%s\n",
    v_ok$ok, v_ok$max_abs_error, v_ok$input_certification, v_ok$verification_status,
    if(length(v_ok$failure_states)) v_ok$failure_states else "(none)"))
cat(sprintf("  tampered: ok=%s max_abs_error=%.3g input_certification=%s verification_status=%s failure_states=%s\n",
    v_bad$ok, v_bad$max_abs_error, v_bad$input_certification, v_bad$verification_status,
    paste(v_bad$failure_states,collapse=",")))
cert <- data.frame(case=c("intact","one_value_plus_0.5"),
  max_abs_error=c(v_ok$max_abs_error,v_bad$max_abs_error),
  tolerance=c(v_ok$tolerance,v_bad$tolerance),
  input_certification=c(v_ok$input_certification,v_bad$input_certification),
  verification_status=c(v_ok$verification_status,v_bad$verification_status),
  failure_states=c("none",paste(v_bad$failure_states,collapse=",")))
write.csv(cert, file.path(OUT,"v42_strict_certification_tamper.csv"), row.names=FALSE)

cat("\n=== (4) Multiplicative simulation DGP re-verification (B=1000) ===\n")
set.seed(20261119L); B<-1000L; n<-1500L
res<-matrix(NA_real_,B,4L); fail<-0L
for (b in seq_len(B)) {
  w1<-rnorm(n); w2<-rnorm(n); F<-rbinom(n,1L,.5)
  P1<-w1*w2; P2<-w2*F
  g1<-runif(1,.20,.80); g2<-runif(1,.20,.80)
  lp<- -0.5+0.5*w1+0.5*w2+g1*P1+g2*P2
  y<-rbinom(n,1L,plogis(lp)); C<-P1+7*P2
  fit<-tryCatch({
    m1<-glm(y~w1+w2+F+C,family=binomial()); m2<-glm(y~w1+w2+F+P1+P2,family=binomial())
    s1<-summary(m1)$coefficients; z<-s1["C","Estimate"]/s1["C","Std. Error"]
    be<-coef(m2)[c("P1","P2")]; V<-vcov(m2)[c("P1","P2"),c("P1","P2")]
    chi2<-as.numeric(t(be)%*%solve(V)%*%be)
    L<-c(-7,1); est<-sum(L*be); se<-sqrt(max(as.numeric(t(L)%*%V%*%L),0))
    c(as.numeric(abs(z)>qnorm(.975)),as.numeric(chi2>qchisq(.95,2)),
      as.numeric(if(se>0) abs(est/se)>qnorm(.975) else 0),as.numeric(abs(g2-7*g1)<0.02))
  },error=function(e) rep(NA_real_,4L))
  if(anyNA(fit)){fail<-fail+1L;next}; res[b,]<-fit
}
ok<-!is.na(res[,1]); N<-sum(ok); R<-colSums(res[ok,,drop=FALSE])
cat(sprintf("  completed=%d failed=%d ; bundled 1df reject=%.1f%% ; unrestricted 2df reject=%.1f%% ; fixed-ratio reject=%.1f%% ; constraint true=%.1f%%\n",
   N,fail,100*R[1]/N,100*R[2]/N,100*R[3]/N,100*R[4]/N))
sim<-data.frame(quantity=c("replications_completed","failed_fits","bundled_1df_rejection_rate",
  "unrestricted_2df_rejection_rate","fixed_ratio_constraint_rejection_rate",
  "fixed_ratio_constraint_true_rate","non_rejection_when_false"),
  value=c(N,fail,R[1]/N,R[2]/N,R[3]/N,R[4]/N,1-R[3]/N))
write.csv(sim,file.path(OUT,"v42_multiplicative_simulation_verify.csv"),row.names=FALSE)

cat("\n=== (5) Calibration intercept CI covariance (published values) ===\n")
# published: correctly_mapped a=-0.1597 SE=.1176 CI[-.4104,.0910]; slope b=.939
# omitted intercept -1.6996 SE=.0940 CI[-1.8999,-1.4993]; intercept shift = -1.5405
# Under a + shift: Var(a+shift)=Var(a) (shift fixed). Observed SE shrinks -> delta method
# with cov(a,b). Report only the algebraic consistency, not an invented number.
a_se<-.1176; b_se<-(function(){ 0.0940 })()  # placeholder: see note
cat("  NOTE: full joint covariance is not stored in the released CSV; only the\n")
cat("  algebraic identity Var(a + k*b) = Var(a) + k^2 Var(b) + 2k Cov(a,b) applies.\n")
cat("  A footnote can state this WITHOUT inventing a numeric covariance.\n")
cat("\nT2_V42_TRUE_RUNS_COMPLETE\n")
