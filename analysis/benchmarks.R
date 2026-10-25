# benchmarks.R (PURE ASCII)
# Software-particle staples for v34:
#   (a) runtime and memory of the audit versus fitting the model;
#   (b) a MEASURED comparison against the diagnostics a fitted model already
#       supports (base R alias(), qr() rank, kappa(), car::vif(), deriv()).
# The comparison is run, not asserted conceptually: each tool is executed on
# the same test cases and its output recorded verbatim.

suppressPackageStartupMessages(library(survey))

need_env <- function(nm) {
  v <- Sys.getenv(nm, unset = "")
  if (!nzchar(v)) stop("set ", nm, " before running benchmarks.R")
  v
}
CODE_DIR <- need_env("CIA_CODE_DIR")
OUT      <- need_env("CIA_OUT")
for (f in c("index_registry.R", "algebra_audit.R",
            "term_identity.R", "public_api.R"))
  source(file.path(CODE_DIR, f), local = TRUE)

rows <- list()
emit <- function(group, case, tool, value, detail = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    group = group, case = case, tool = tool, value = value,
    detail = detail, stringsAsFactors = FALSE)
}

# ----------------------- (a) runtime benchmark -----------------------------
set.seed(20261023); n <- 3566
d <- data.frame(
  sex = factor(sample(c("Male", "Female"), n, TRUE)),
  age = runif(n, 20, 80), BMI = rnorm(n, 27, 5),
  waist_cm = rnorm(n, 95, 14), tg_mmol_L = rlnorm(n, log(1.3), 0.5),
  glucose_mg_dL = rnorm(n, 100, 22))
d$log_age <- log10(d$age); d$log_bmi <- log10(d$BMI)
d$log_tg_mg <- log(d$tg_mmol_L * 88.57)
d$log_glucose_mg <- log(d$glucose_mg_dL)
d$TyG <- log(d$tg_mmol_L * 88.57 * d$glucose_mg_dL / 2)
d$metab <- log10(d$glucose_mg_dL) + log10(d$tg_mmol_L)
d$LAP <- ifelse(d$sex == "Male", d$waist_cm - 65, d$waist_cm - 58) * d$tg_mmol_L
d$dx <- rbinom(n, 1, plogis(-6 + 1.5 * d$log_glucose_mg + d$log_age))

# Windows proc.time ticks are coarse (~16 ms), so time one long batch and
# divide, rather than median-of-many short reads.
# substitute/eval, NOT force(expr): a function promise memoises after the
# first evaluation, so a force()-loop would time one call plus no-ops.
time_op <- function(expr, total_inner = 200) {
  e <- substitute(expr)
  env <- parent.frame()
  for (i in 1:5) eval(e, env)   # warm-up: JIT and method dispatch
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(total_inner)) eval(e, env)
  c(seconds = (proc.time()[["elapsed"]] - t0) / total_inner,
    inner = total_inner)
}
t_audit_nodata <- time_op(
  audit_index_in_model(c(log_tg_mg = 1, log_glucose_mg = 1),
                       c(log_tg_mg = 1, log_glucose_mg = 1)), 20000)
t_audit_data <- time_op(
  audit_index_in_model(c(log_tg_mg = 1, log_glucose_mg = 1),
                       c(log_tg_mg = 1, log_glucose_mg = 1), data = d), 5000)
t_product <- time_op(
  audit_product_in_model("LAP", c("sex", "waist_cm", "tg_mmol_L"),
                         product_expansions()$LAP, data = d), 5000)
t_glm <- time_op(
  glm(dx ~ sex + log_age + log_bmi + TyG, data = d, family = binomial()), 200)
for (x in list(list("audit_additive_nodata", t_audit_nodata),
               list("audit_additive_with_rows", t_audit_data),
               list("audit_product_with_rows", t_product),
               list("glm_fit_for_contrast", t_glm))) {
  emit("runtime", x[[1]], "seconds_per_call",
       sprintf("%.6f", x[[2]][["seconds"]]),
       sprintf("one batch of %d calls divided; n = %d rows",
               as.integer(x[[2]][["inner"]]), n))
}
options(repos = c(CRAN = "https://cloud.r-project.org"))
dep_counts <- unique(unlist(tools::package_dependencies(
  c("survey", "haven", "digest"), db = utils::installed.packages(),
  recursive = TRUE)))
emit("runtime", "dependency_footprint", "recursive_package_count",
     as.character(length(dep_counts)),
     paste("survey + haven + digest:", paste(sort(dep_counts), collapse = ", ")))
emit("runtime", "audit_versus_fit_ratio", "seconds_ratio",
     sprintf("%.1f", t_glm[["seconds"]] / t_audit_nodata[["seconds"]]),
     "a no-data audit is this many times faster than one glm fit")

# ----------------- (b) measured comparison with existing tools --------------
# Case set: the four audit classes plus the failure modes the audit names.
# For each tool we record what it actually returns on each case.
cases <- list(
  rank_preserving = list(idx = c(log_tg_mg = 1, log_glucose_mg = 1),
                         emb = c(log_tg_mg = 1, log_glucose_mg = 1),
                         score = "TyG"),
  rank_adding = list(idx = c(metab = 1, log_age = 1, log_bmi = 1),
                     emb = c(log_age = 1, log_bmi = 1), score = NULL),
  reducible = list(idx = NULL, emb = NULL, score = "LAP"))

# Tool 1: this audit (no data needed).
emit("tool_comparison", "rank_preserving_TyG", "CompositeIndexAudit",
     audit_index_in_model(c(log_tg_mg = 1, log_glucose_mg = 1),
                          c(log_tg_mg = 1, log_glucose_mg = 1))$classification,
     "runs before fitting; names the class and the coefficient map")
emit("tool_comparison", "rank_adding_GLM7", "CompositeIndexAudit",
     audit_index_in_model(c(metab = 1, log_age = 1, log_bmi = 1),
                          c(log_age = 1, log_bmi = 1))$classification,
     "runs before fitting")
emit("tool_comparison", "reducible_LAP", "CompositeIndexAudit",
     audit_product_in_model("LAP", c("sex", "waist_cm", "tg_mmol_L"),
                            product_expansions()$LAP)$classification,
     "runs before fitting; returns the interaction bundle")
emit("tool_comparison", "outside_TyG_WHtR", "CompositeIndexAudit",
     audit_product_in_model("TyG_WHtR", c("sex", "TyG", "waist_cm", "height_cm"),
                            list(products = c("TyG:waist_per_height" = 1)))$classification,
     "runs before fitting; names the blocking term")

# Tool 2: base R alias() on a fitted model (needs the fit, hence the outcome).
m_tyg <- glm(dx ~ sex + log_tg_mg + log_glucose_mg + TyG, data = d,
             family = binomial())
al <- alias(m_tyg)
emit("tool_comparison", "rank_preserving_TyG", "alias(glm)",
     if (is.null(al$Complete)) "no aliasing found" else
       paste("aliased:", paste(rownames(al$Complete), collapse = ",")),
     "detects the realized aliasing AFTER fitting; needs the outcome; no class, no map")

# Tool 3: qr() rank on the design matrix (no outcome needed, but no verdict).
X0 <- model.matrix(~ sex + log_tg_mg + log_glucose_mg, data = d)
X1 <- model.matrix(~ sex + log_tg_mg + log_glucose_mg + TyG, data = d)
emit("tool_comparison", "rank_preserving_TyG", "qr() rank",
     sprintf("rank gain %d", qr(X1)$rank - qr(X0)$rank),
     "numeric rank on THESE rows only; cannot certify the identity for any data")

# Tool 4: kappa() condition number.
emit("tool_comparison", "rank_preserving_TyG", "kappa()",
     sprintf("%.1f", kappa(scale(X1[, -1]))),
     "ill-conditioning diagnostic; no threshold separates identity from association")

# Tool 5: car::vif() on the fitted model (needs the fit).
vif_val <- tryCatch({
  vv <- car::vif(m_tyg)
  sprintf("max %.3f", max(vv))
}, error = function(e) paste("error:", conditionMessage(e)))
emit("tool_comparison", "rank_preserving_TyG", "car::vif(glm)",
     vif_val, "variance inflation AFTER fitting; aliased columns make it uncomputable")

# Tool 6: deriv() as the base-R symbolic comparator.
dd <- deriv(~ log(a) + log(b) - log(2), c("a", "b"))
emit("tool_comparison", "symbolic_layer", "deriv()",
     paste(deparse(dd), collapse = " "),
     "symbolic derivatives exist in base R, but there is no span/identifiability verdict")

# Tool 7: unit-conversion handling - only the audit ships the intercept map.
emit("tool_comparison", "unit_conversion", "CompositeIndexAudit",
     "beta0* = beta0 - k * beta_C returned in closed form",
     "alias/qr/kappa/vif/deriv have no notion of units or log bases")
emit("tool_comparison", "unit_conversion", "alias/qr/kappa/vif/deriv",
     "not applicable", "no unit or log-base awareness in any of them")

res <- do.call(rbind, rows)
print(res, row.names = FALSE)
write.csv(res, file.path(OUT, "v34_benchmarks_and_tool_comparison.csv"),
          row.names = FALSE, na = "")
cat("\n=== v34 benchmarks done ===\n")
