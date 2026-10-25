# verify_unit_factor.R (PURE ASCII)
#
# v34: the manuscript states that converting released insulin from pmol/L to
# uIU/mL and lipid measures from mmol/L to mg/dL "adds a constant k = 1.1691
# to GLM7".  A reviewer pointed out that the registered GLM7 expression does
# not show insulin at all, and asked for the composition of k: which
# conversion factors, which log base, and why the lipid terms do not all
# contribute.  This script derives k term by term from the seven measures,
# checks it against the shift measured on the data, and reports the SI
# alternative that the code does NOT use.
#
# Run: CIA_OUT=... Rscript verify_unit_factor.R   (no NHANES data needed;
# the empirical check is read back from v34_unit_scaling.csv when present.)

OUT <- Sys.getenv("CIA_OUT", unset = "")
if (!nzchar(OUT)) stop("set CIA_OUT before running verify_unit_factor.R")

# ---------------------------------------------------------------------------
# The seven routine measures and the unit each is released in.
#   measure        released unit        coefficient in GLM7   conversion
#   glucose        mg/dL                +1                    none
#   insulin        pmol/L               +1                    / 6.00
#   triglycerides  mmol/L               +1                    x 88.57
#   LDL            mmol/L               +1                    x 38.67
#   HDL            mmol/L               -1                    x 38.67
#   age            years (log10)        +1                    none
#   BMI            kg/m2 (log10)        +1                    none
# ---------------------------------------------------------------------------
GLM7_TERMS <- data.frame(
  measure = c("glucose", "insulin", "triglycerides", "LDL", "HDL", "age", "BMI"),
  released_unit = c("mg/dL", "pmol/L", "mmol/L", "mmol/L", "mmol/L", "years", "kg/m2"),
  coefficient = c(1, 1, 1, 1, -1, 1, 1),
  conversion_to_mg_dL_or_uIU = c(1, 1 / 6, 88.57, 38.67, 38.67, 1, 1),
  stringsAsFactors = FALSE)
GLM7_TERMS$log10_shift <- GLM7_TERMS$coefficient *
  log10(GLM7_TERMS$conversion_to_mg_dL_or_uIU)

k_analytic <- sum(GLM7_TERMS$log10_shift)
k_closed   <- log10(88.57 / 6)
# The cholesterol factor 38.67 enters once with +1 (LDL) and once with -1
# (HDL), so it cancels exactly and does not appear in k.
cholesterol_cancels <- isTRUE(all.equal(
  sum(GLM7_TERMS$log10_shift[GLM7_TERMS$measure %in% c("LDL", "HDL")]), 0))
stopifnot(cholesterol_cancels)
stopifnot(isTRUE(all.equal(k_analytic, k_closed, tolerance = 1e-12)))

# The SI conversion for insulin is 1 uIU/mL = 6.945 pmol/L, not 6.00 pmol/L.
# 6.00 is the factor implied by the released NHANES insulin variable
# (LBXIN: uIU/mL; LBDINSI: pmol/L, where 1 uIU/mL = 6.00 pmol/L).  Reporting
# both makes the choice auditable instead of implicit.
k_si <- log10(88.57 / 6.945)

print(GLM7_TERMS, row.names = FALSE)
cat("\nk = sum of the seven per-term log10 shifts =", format(k_analytic, digits = 15), "\n")
cat("k in closed form                          = log10(88.57 / 6) =",
    format(k_closed, digits = 15), "\n")
cat("k if insulin were converted with 6.945     = log10(88.57 / 6.945) =",
    format(k_si, digits = 15), "\n")

# ---- empirical cross-check against the measured shift ----------------------
emp <- NA_real_
p <- file.path(OUT, "v34_unit_scaling.csv")
if (file.exists(p)) {
  u <- read.csv(p, stringsAsFactors = FALSE)
  row <- u[u$quantity == "unit_shift", ]
  if (nrow(row)) emp <- as.numeric(row$value[1])
}
if (is.finite(emp)) {
  cat("measured shift on the data                =", format(emp, digits = 15), "\n")
  stopifnot(abs(emp - k_analytic) < 1e-9)
  cat("analytic k reproduces the measured shift\n")
} else {
  cat("v34_unit_scaling.csv not found; analytic derivation only\n")
}

out <- data.frame(
  quantity = c("n_measures_in_GLM7", "k_analytic", "k_closed_form",
               "k_closed_form_expression", "k_with_SI_insulin_factor_6.945",
               "cholesterol_factor_cancels", "insulin_factor_used",
               "triglyceride_factor_used", "cholesterol_factor_used",
               "log_base", "glucose_shift", "insulin_shift",
               "triglyceride_shift", "ldl_shift", "hdl_shift",
               "age_shift", "bmi_shift", "measured_shift_on_data",
               "analytic_minus_measured"),
  value = c("7", format(k_analytic, digits = 15), format(k_closed, digits = 15),
            "log10(88.57 / 6)", format(k_si, digits = 15),
            "TRUE (LDL +1 and HDL -1 both use 38.67, so the term is exactly 0)",
            "6.00 pmol/L per uIU/mL (NHANES LBDINSI convention)",
            "88.57 mg/dL per mmol/L", "38.67 mg/dL per mmol/L",
            "10", format(GLM7_TERMS$log10_shift[1], digits = 15),
            format(GLM7_TERMS$log10_shift[2], digits = 15),
            format(GLM7_TERMS$log10_shift[3], digits = 15),
            format(GLM7_TERMS$log10_shift[4], digits = 15),
            format(GLM7_TERMS$log10_shift[5], digits = 15),
            format(GLM7_TERMS$log10_shift[6], digits = 15),
            format(GLM7_TERMS$log10_shift[7], digits = 15),
            if (is.finite(emp)) format(emp, digits = 15) else "not measured",
            if (is.finite(emp)) format(k_analytic - emp, digits = 3) else "NA"),
  stringsAsFactors = FALSE)
write.csv(out, file.path(OUT, "v34_unit_factor_derivation.csv"), row.names = FALSE)
cat("wrote v34_unit_factor_derivation.csv\n")
