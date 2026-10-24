# vignette_v33.R  (PURE ASCII)  -  runs in a clean session without NHANES data
#
# Box 1 of the v33 manuscript, verbatim and self-contained: the example builds
# its own synthetic data, so it runs from a fresh directory with
#   Rscript vignette_v33.R
# Ten minutes, four audit classes, the coefficient map, the aliased
# coefficient a rank-preserving index produces when fitted anyway, the
# no-data audit, the keyed term identity, and the public API
# (parse / register / verify / normalize).
#
# Expected output (abridged):
#   [1] "rank_preserving_exact" "rank_adding_exact" "reducible_to_interaction"
#       "outside_supported_class"
#   [1] "rank_preserving_exact"        # the no-data audit of TyG
#   GLM7 coefficient NA (aliased)
#   LAP: X_[W,I] K beta reproduces the linear predictor exactly
#   verify_index_column flags the deliberately miss-scaled column
#   normalize_log_base rescales by ln(10) = 2.302585

source("index_registry_v33.R")
source("algebra_audit_v33.R")
source("term_identity_v33.R")
source("public_api_v33.R")

set.seed(20261023); n <- 1500
d <- data.frame(
  sex = factor(sample(c("Male", "Female"), n, TRUE), levels = c("Male", "Female")),
  age = runif(n, 20, 80),
  BMI = rnorm(n, 27, 5),
  waist_cm = rnorm(n, 95, 14),
  tg_mmol_L = rlnorm(n, log(1.3), 0.5),
  glucose_mg_dL = rnorm(n, 100, 22))
d$sexFemale <- as.numeric(d$sex == "Female")
d$log_tg_mg <- log(d$tg_mmol_L * 88.57)
d$log_glucose_mg <- log(d$glucose_mg_dL)
d$TyG <- log(d$tg_mmol_L * 88.57 * d$glucose_mg_dL / 2)          # ln(TG*G/2)
# metab is the worked-example residual of GLM7 = metab + log10(age) +
# log10(BMI).  In this synthetic example it is a stand-in built from two
# measured columns; the registered GLM7 used on NHANES takes the five
# published components (see index_registry_v33.R).
d$metab <- log10(d$glucose_mg_dL) + log10(d$tg_mmol_L)
d$log_age <- log10(d$age); d$log_bmi <- log10(d$BMI)
d$GLM7 <- d$metab + d$log_age + d$log_bmi
d$LAP <- ifelse(d$sex == "Male", d$waist_cm - 65, d$waist_cm - 58) * d$tg_mmol_L

# --- the four classes --------------------------------------------------------
reg <- v33_registry(); exps <- product_expansions()

a_tyg <- audit_index_in_model(c(log_tg_mg = 1, log_glucose_mg = 1),
                              c(log_tg_mg = 1, log_glucose_mg = 1), data = d)
a_glm7 <- audit_index_in_model(c(metab = 1, log_age = 1, log_bmi = 1),
                               c(log_age = 1, log_bmi = 1), data = d)
a_lap <- audit_product_in_model("LAP", v28_component_model_terms("LAP", reg),
                                exps$LAP, data = d)
# VAI multiplies THREE quantities; a three-factor product is outside the
# supported class.  data = NULL: the symbolic verdict needs no rows.
a_vai <- audit_product_in_model("VAI", v28_component_model_terms("VAI", reg),
                                list(products = c("waist_cm:tg_mmol_L:hdl_mmol_L" = 1)))
classes <- c(a_tyg$classification, a_glm7$classification,
             a_lap$classification, a_vai$classification)
print(classes)
stopifnot(identical(unname(classes),
  c("rank_preserving_exact", "rank_adding_exact",
    "reducible_to_interaction", "outside_supported_class")))

# --- the no-data audit (v33): the verdict is symbolic, not a sample property -
# TyG against its own components, with no rows at all.
a_nodata <- audit_index_in_model(c(log_tg_mg = 1, log_glucose_mg = 1),
                                 c(log_tg_mg = 1, log_glucose_mg = 1),
                                 data = NULL)
print(a_nodata$classification)
stopifnot(a_nodata$classification == "rank_preserving_exact")

# --- coefficient map and the aliased coefficient -----------------------------
print(a_glm7$coefficient_map)          # gamma_log_age = beta_log_age + beta_C
d$dx <- rbinom(n, 1, plogis(-3 + 0.3 * d$metab + 0.2 * d$log_age))
fit <- glm(dx ~ log_age + log_bmi + metab + GLM7, data = d, family = binomial())
print(coef(fit)["GLM7"])               # NA: rank-preserving, aliased
stopifnot(is.na(coef(fit)[["GLM7"]]))

# --- reducible_to_interaction: the bundled coefficient is exact --------------
d$dx2 <- rbinom(n, 1, plogis(-3 + 0.02 * d$LAP + 0.01 * d$waist_cm))
m_aug <- glm(dx2 ~ sex + waist_cm + tg_mmol_L + LAP, data = d, family = binomial())
d$w_tg <- d$waist_cm * d$tg_mmol_L
d$tg_f <- d$tg_mmol_L * d$sexFemale
K <- a_lap$coefficient_map
b <- coef(m_aug)
mm_c <- model.matrix(~ sex + waist_cm + tg_mmol_L + LAP, d)
lp_aug <- as.numeric(mm_c %*% b[colnames(mm_c)])
names(b)[names(b) == "LAP"] <- "C"
lp_map <- as.numeric(model.matrix(~ sex + waist_cm + tg_mmol_L + w_tg + tg_f, d) %*%
                       (K %*% b[colnames(K)]))
cat("LAP max |lp difference| under gamma = K beta:",
    max(abs(lp_aug - lp_map)), "\n")
stopifnot(max(abs(lp_aug - lp_map)) < 1e-8)

# --- public API ---------------------------------------------------------------
p <- parse_index_expression("TyG * BMI")
print(p$status)                        # "product_expandable"
stopifnot(p$status == "product_expandable", p$products[["TyG:BMI"]] == 1)

register_index("DemoTyGBMI", "TyG * BMI")
d$DemoTyGBMI <- d$TyG * d$BMI
a_demo <- audit_registered_index("DemoTyGBMI", c("sex", "TyG", "BMI"), data = d)
print(a_demo$classification)           # "reducible_to_interaction"
stopifnot(a_demo$classification == "reducible_to_interaction")

# Keyed model terms: the same call, but with the five-field term identity
# declared explicitly, so units and log bases are part of the match.
register_index("DemoLogTG", "log(tg_mmol_L)", units = c(tg_mmol_L = "mmol/L"),
               overwrite = TRUE)
d$DemoLogTG <- log(d$tg_mmol_L)
a_keyed <- audit_registered_index(
  "DemoLogTG",
  c(log_tg = "tg_mmol_L|log|e|mmol/L|-"),
  data = d)
print(a_keyed$classification)          # "rank_preserving_exact"
stopifnot(a_keyed$classification == "rank_preserving_exact",
          a_keyed$term_matching$mode == "keyed")

v_ok <- verify_index_column(d, "DemoTyGBMI")
print(v_ok$ok)                         # TRUE
d$DemoTyGBMI_bad <- d$TyG * d$BMI * 1.01   # a 1% scale slip, as in unit mistakes
v_bad <- verify_index_column(d, "DemoTyGBMI_bad", expression = "TyG * BMI")
print(c(v_bad$ok, signif(v_bad$max_abs_error, 3)))
stopifnot(v_ok$ok, !v_bad$ok)

cf_ln <- c(log_tg_mg = 0.8, log_glucose_mg = 1.1)   # fitted against ln(x)
cf_10 <- normalize_log_base(cf_ln, from = "ln", to = "log10")
print(round(cf_10, 4))                 # 1.8421, 2.5328 (x ln(10) = 2.302585)
stopifnot(abs(cf_10[["log_tg_mg"]] - 0.8 * log(10)) < 1e-12)

cat("=== vignette_v33 done: classes, no-data audit, keyed identity, map, aliasing and API verified ===\n")
