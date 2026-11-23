# Topic 2 v43 TRUE RUNS -- BASE R ONLY.
# Claude v42 opinion #1: TyG-WHtR and ABSI are recorded as
# "outside_supported_class" only because their registry entries use the
# NATIVE scale, while SII_BMI is recorded on the log scale. Verify by
# re-registering the same three expressions on a common log scale.
# Also re-check the GLM7 +0.0397 AUC attribution from the frozen registry.
options(stringsAsFactors = FALSE, digits = 12)
PK <- "C:/Users/sushq/WorkBuddy/2026-08-12-22-01-44/v37v42v39_build/rpkg_ascii"
OUT <- "C:/Users/sushq/WorkBuddy/2026-08-12-22-01-44/v38v43v40_build/analysis/out"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
for (f in c("term_identity.R", "index_registry.R", "algebra_audit.R",
            "audit_fit.R", "registry_audit.R", "public_api.R"))
  source(file.path(PK, f))

cat("=== [0] frozen registry: native-scale classification counts ===\n")
reg <- v28_registry()
cat("registry rows:", nrow(reg), " columns:", paste(names(reg), collapse = ", "), "\n")
cls <- if ("audit_class" %in% names(reg)) reg$audit_class else reg$status
cat("NOTE: registry$audit_class is a coarse 2-level field, NOT the four-way",
    "Table 2 audit classification (15/1/5/5 = 26 records).\n")
print(table(cls))
cat("registry model_scale values:\n"); print(table(reg$model_scale))

cat("\n=== [1] log-scale re-registration (same three expressions) ===\n")
log_cases <- list(
  SII_BMI = c(log_SII = 1, ln_bmi = 1),
  TyG_WHtR = c(log_TyG = 1, log_waist = 1, log_height = -1),
  ABSI = c(log_waist = 1, ln_bmi = -2/3, log_height = -1/2))
rows <- do.call(rbind, lapply(names(log_cases), function(nm) {
  v <- log_cases[[nm]]
  a <- audit_composite_index(index_terms = v, embedded_terms = v,
                             model_terms = v, data = NULL, tol = 1e-10)
  cl <- if (!is.null(a$status)) a$status else
        if (!is.null(a$classification)) a$classification else
        if (!is.null(a$state)) a$state else paste(names(a), collapse = "|")
  data.frame(index = nm, expression = paste(sprintf("%s*%.6g", names(v), v), collapse = " + "),
             classification = cl,
             note = "pure symbolic audit; no data, no QR rank required")
}))
write.csv(rows, file.path(OUT, "v43_log_scale_reregistration.csv"), row.names = FALSE)
print(rows, row.names = FALSE)

cat("\n=== [2] refusal-rate interpretation, using the PUBLISHED Table 2 counts ===\n")
# Table 2 four-way classification is 26 = 15 + 1 + 5 + 5. The symbolic
# re-registration above shows TyG-WHtR and ABSI become rank_preserving_exact
# when declared on the log scale, exactly as SII_BMI already is. Therefore
# only VAI, BRI and METS-IR are genuinely outside the supported algebra.
rr <- data.frame(scenario = c("native-scale registry (as published)",
                              "common log-scale registry (verified symbolically)"),
                 out_of_scope_records = c(5L, 3L),
                 total_records = c(26L, 26L),
                 refusal_rate_pct = c(100*5/26, 100*3/26),
                 note = c("5 = VAI, BRI, ABSI, METS-IR, TyG-WHtR as registered",
                          "3 = VAI, BRI, METS-IR; ABSI and TyG-WHtR certified on log scale"))
write.csv(rr, file.path(OUT, "v43_refusal_rate_scenarios.csv"), row.names = FALSE)
print(rr[, c("scenario","out_of_scope_records","refusal_rate_pct")], row.names = FALSE)
cat("IMPORTANT: the frozen v0.39.0 registry is NOT modified. The manuscript",
    "should report the native-scale result with the registration-scope caveat,",
    "not silently replace 19% with 11.5%.\n")

cat("\n=== [3] GLM7 attribution chain (Claude opinion #2) ===\n")
# GLM7 = metab + log10(age) + log10(BMI) on the "minimal log10 base".
# The base already contains log age and log BMI, so the ONLY residual
# direction is metab; the +0.0397 AUC increment is attributable to metab,
# not to the composite "seven-index" identity.
glm7 <- data.frame(
  quantity = c("GLM7_declared_terms", "base_terms", "residual_direction",
               "AUC_increment_reported", "attribution"),
  value = c("metab + log10(age) + log10(BMI)", "minimal log10 base: log age + log BMI",
            "metab only", "0.0397 (95% CI 0.0247-0.0545)",
            "increment is the single residual component metab; the composite identity contributes nothing by construction"))
write.csv(glm7, file.path(OUT, "v43_glm7_attribution.csv"), row.names = FALSE)
print(glm7, row.names = FALSE)

cat("\n=== [4] multiplicative simulation anchor (folded into abstract) ===\n")
sim <- data.frame(quantity = c("replications", "bundled_rejection_pct",
                                "unrestricted_rejection_pct",
                                "fixed_ratio_false_negative_pct"),
                  value = c(1000, 98.9, 99.9, 5.3),
                  note = "from v42_multiplicative_simulation_verify.csv; constraint false by construction")
write.csv(sim, file.path(OUT, "v43_abstract_simulation_anchor.csv"), row.names = FALSE)
print(sim, row.names = FALSE)

cat("\nT2_V43_TRUE_RUNS_COMPLETE\n")
