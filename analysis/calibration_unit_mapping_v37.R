# Topic 2 v37: calibration consequences of omitting the unit-map intercept.
# This is a focused rerun; it does not repeat the Arm A/B ranking analyses.
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")
need_env <- function(x) {
  v <- Sys.getenv(x, unset = "")
  if (!nzchar(v)) stop("set ", x)
  v
}
CODE_DIR <- need_env("CIA_CODE_DIR")
OUT <- need_env("CIA_OUT")
FIG_OUT <- need_env("CIA_FIG_OUT")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)
for (f in c("index_registry.R", "model_evaluation.R", "algebra_audit.R", "data_build.R"))
  source(file.path(CODE_DIR, f))
req <- c("dx_diab", "weight", "stra", "psu", "sex", "log_age", "log_bmi",
         "GLM7", "glucose_mg_dL", "insulin_pmol_L", "tg_mmol_L", "LDL", "hdl_mmol_L")
P <- build_period("P"); L <- build_period("L")
P <- P[complete.cases(P[, req]), ]; L <- L[complete.cases(L[, req]), ]
Pdes <- make_survey_design(P); Ldes <- make_survey_design(L)
us <- run_unit_scaling_sensitivity(P, L, outcome = "dx_diab",
                                  development_design = Pdes, later_design = Ldes)
stopifnot(nrow(L) == 2653L, sum(L$dx_diab) == 368L)
models <- list(correctly_mapped = us$lp_l_correct,
               intercept_omitted = us$lp_l_incorrect_transfer)
metrics <- do.call(rbind, lapply(names(models), function(nm) {
  m <- prediction_metrics(L, models[[nm]], "dx_diab", design = Ldes)
  data.frame(model = nm, n = nrow(L), cases = sum(L$dx_diab),
             survey_weighted = TRUE, AUC = unname(m[["auc"]]),
             Brier = unname(m[["brier"]]),
             calibration_intercept = unname(m[["calibration_intercept"]]),
             calibration_slope = unname(m[["calibration_slope"]]))
}))
uncertainty <- do.call(rbind, lapply(names(models), function(nm)
  calibration_uncertainty(L, models[[nm]], "dx_diab", design = Ldes, model = nm)))
write.csv(metrics, file.path(OUT, "nhanes_unit_calibration_metrics_v37.csv"), row.names = FALSE)
write.csv(uncertainty, file.path(OUT, "nhanes_unit_calibration_uncertainty_v37.csv"), row.names = FALSE)

# Survey-weighted LOESS is a descriptive curve; Taylor-linearized intervals
# for calibration intercept/slope are supplied in the separate result table.
smooth_curve <- function(lp) {
  dd <- data.frame(p = plogis(lp), y = L$dx_diab, w = L$weight/mean(L$weight))
  fit <- loess(y ~ p, data = dd, weights = w, span = 0.75, degree = 1,
               control = loess.control(surface = "direct"))
  grid <- seq(min(dd$p), max(dd$p), length.out = 500)
  data.frame(predicted = grid, observed = pmax(0, pmin(1, predict(fit, newdata = data.frame(p = grid)))))
}
curves <- lapply(models, smooth_curve)
draw <- function() {
  par(mar = c(4.6, 4.8, 3.2, 1.0), las = 1)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Predicted probability",
       ylab = "Observed proportion", main = "GLM7 later-period calibration")
  abline(0, 1, lty = 3, col = "grey45")
  lines(curves$correctly_mapped, col = "#1261A0", lwd = 2.6)
  lines(curves$intercept_omitted, col = "#C43C39", lwd = 2.6, lty = 2)
  legend("topleft", bty = "n", cex = 0.85,
         legend = c("Correctly mapped intercept", "Intercept omitted", "Ideal"),
         col = c("#1261A0", "#C43C39", "grey45"), lty = c(1, 2, 3), lwd = c(2.6, 2.6, 1))
  mtext("Survey-weighted LOESS (span = 0.75); axes truncated to [0,1]", side = 1, line = 3.5, cex = 0.72)
}
png(file.path(FIG_OUT, "Figure5.png"), width = 1800, height = 1600, res = 220)
draw(); dev.off()
pdf(file.path(FIG_OUT, "Figure5.pdf"), width = 7.2, height = 6.4)
draw(); dev.off()
curve_out <- do.call(rbind, lapply(names(curves), function(nm)
  data.frame(model = nm, curves[[nm]])))
write.csv(curve_out, file.path(OUT, "nhanes_unit_calibration_curve_v37.csv"), row.names = FALSE)
# Survey-weighted net benefit as a sensitivity analysis on the same predictions.
# Net benefit is a population-standardized weighted proportion; it is not a
# clinical recommendation and does not add development-coefficient uncertainty.
nb_rows <- list(); ii <- 0L
for (threshold in c(0.05, 0.10, 0.15, 0.20)) {
  for (nm in names(models)) {
    p <- plogis(models[[nm]]); y <- L$dx_diab; w <- L$weight
    flag <- p >= threshold; odds <- threshold/(1-threshold)
    nb <- sum(w*(flag*y - flag*(1-y)*odds))/sum(w)
    treat_all <- sum(w*(y - (1-y)*odds))/sum(w)
    ii <- ii + 1L
    nb_rows[[ii]] <- data.frame(threshold=threshold, model=nm, n=nrow(L),
      weighted_net_benefit=nb, weighted_treat_all_net_benefit=treat_all,
      weighted_flagged_percent=100*sum(w*flag)/sum(w),
      weighted_detected_cases_per1000=1000*sum(w*flag*y)/sum(w),
      weighted_missed_cases_per1000=1000*sum(w*(!flag)*y)/sum(w),
      estimand="Survey-weighted net-benefit sensitivity on frozen predictions")
  }
}
nb <- do.call(rbind, nb_rows)
write.csv(nb, file.path(OUT, "nhanes_unit_weighted_net_benefit_v39.csv"), row.names=FALSE)
print(nb, row.names=FALSE)
print(metrics, row.names = FALSE)
print(uncertainty, row.names = FALSE)
cat("UNIT_CALIBRATION_V39_WITH_WEIGHTED_NB_COMPLETE\n")
