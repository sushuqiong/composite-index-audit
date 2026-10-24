# v28 pre-specified candidate registry and formula engine.
# Scope: documented composite indices with an explicit component-conditional model.
# Raw NHANES files are supplied by the user through GLM7_T2_V28_DATA and are not
# copied into the submission package.

v28_safe_log <- function(x) ifelse(is.finite(x) & x > 0, log(x), NA_real_)
v28_safe_log10 <- function(x) ifelse(is.finite(x) & x > 0, log10(x), NA_real_)

v28_registry <- function() {
  # Raw HOMA_IR is one of the 24 primary candidates but cannot support an exact
  # additive identity on its product scale. HOMA_IR_log is a separate scale
  # sensitivity and is not used to select a preferred index.
  data.frame(
    index = c("GLM7", "TyG", "TyG_BMI", "TyG_WC", "TyG_WHtR", "METS_IR",
      "LAP", "VAI", "CVAI", "BRI", "ABSI", "AIP", "HOMA_IR",
      "CAR", "MHR", "TyG_hsCRP", "hsCRP_BMI", "hsCRP_WC", "SII_BMI",
      "NLR", "PLR", "SII", "SIRI", "PNI", "HOMA_IR_log"),
    family = c("GLM7_worked_example", "TyG", "TyG_body", "TyG_body", "TyG_body",
      "metabolic_body", "metabolic_body", "metabolic_body", "metabolic_body",
      "metabolic_body", "metabolic_body", "metabolic", "metabolic",
      "inflammation_metabolic", "inflammation_metabolic", "inflammation_metabolic",
      "inflammation_body", "inflammation_body", "inflammation_body",
      "inflammation_hematologic", "inflammation_hematologic", "inflammation_hematologic",
      "inflammation_hematologic", "inflammation_nutrition", "metabolic_sensitivity"),
    analysis_role = c(rep("primary_candidate", 24), "scale_sensitivity"),
    status = c(rep("published_formula", 3), rep("exploratory_constructed", 2),
      rep("published_formula", 8), rep("exploratory_constructed", 8),
      rep("published_formula", 4)),
    formula = c(
      "metab + log10(age) + log10(BMI)",
      "ln(TG_mg_dL*glucose_mg_dL/2)", "TyG*BMI", "TyG*waist_cm", "TyG*(waist_cm/height_cm)",
      "ln(2*glucose_mg_dL+TG_mg_dL)*BMI/ln(HDL_mg_dL)",
      "(waist_cm-65_or_58)*TG_mmol_L", "sex_specific_VAI(WC,BMI,TG,HDL)",
      "sex_specific_CVAI(age,BMI,waist_cm,log10(TG_mmol_L),HDL_mmol_L)", "BRI_from_waist_and_height",
      "(waist_m)/(BMI^(2/3)*height_m^(1/2))", "log10(TG_mmol_L/HDL_mmol_L)",
      "glucose_mmol_L*insulin_uIU_mL/22.5", "hsCRP_mg_L/albumin_g_dL",
      "monocytes/HDL_mg_dL", "TyG*hsCRP_mg_L", "hsCRP_mg_L*BMI", "hsCRP_mg_L*waist_cm",
      "SII*BMI", "neutrophils/lymphocytes", "platelets/lymphocytes",
      "platelets*neutrophils/lymphocytes", "neutrophils*monocytes/lymphocytes",
      "10*albumin_g_dL+0.005*(lymphocytes_1000_cells_uL*1000)", "log(glucose_mmol_L)+log(insulin_uIU_mL)-log(22.5)"),
    model_scale = c("native", "native", "native", "native", "native", "native", "native",
      "native", "native", "native", "native", "native", "native", "log", "log", "log", "log",
      "log", "log", "log", "log", "log", "log", "native", "log"),
    component_terms = c(
      "sex;log10(age);log10(BMI)", "sex;log(TG_mg_dL);log(glucose_mg_dL)",
      "sex;TyG;BMI", "sex;TyG;waist_cm", "sex;TyG;waist_cm;height_cm",
      "sex;log(2*glucose_mg_dL+TG_mg_dL);log(HDL_mg_dL);BMI", "sex;waist_cm;TG_mmol_L",
      "sex;waist_cm;BMI;TG_mmol_L;HDL_mmol_L", "sex;age;BMI;waist_cm;log10(TG_mmol_L);HDL_mmol_L",
      "sex;waist_cm;height_cm", "sex;waist_cm;BMI;height_cm", "sex;log(TG_mmol_L);log(HDL_mmol_L)",
      "sex;glucose_mmol_L;insulin_uIU_mL", "sex;log(hsCRP_mg_L);log(albumin_g_dL)", "sex;log(monocytes);log(HDL_mg_dL)",
      "sex;log(TyG);log(hsCRP_mg_L)", "sex;log(hsCRP_mg_L);log(BMI)", "sex;log(hsCRP_mg_L);log(waist_cm)", "sex;log(SII);log(BMI)",
      "sex;log(neutrophils);log(lymphocytes)", "sex;log(platelets);log(lymphocytes)", "sex;log(platelets);log(neutrophils);log(lymphocytes)",
      "sex;log(neutrophils);log(monocytes);log(lymphocytes)", "sex;albumin_g_dL;lymphocytes_cells_mm3",
      "sex;log(glucose_mmol_L);log(insulin_uIU_mL)"),
    audit_class = c(
      "exact_reparameterization", "exact_reparameterization",
      "not_supported", "not_supported", "not_supported", "not_supported",
      "not_supported", "not_supported", "not_supported", "not_supported",
      "not_supported", "exact_reparameterization", "not_supported",
      "exact_reparameterization", "exact_reparameterization", "exact_reparameterization",
      "exact_reparameterization", "exact_reparameterization", "exact_reparameterization",
      "exact_reparameterization", "exact_reparameterization", "exact_reparameterization",
      "exact_reparameterization", "exact_reparameterization", "exact_reparameterization"),
    component_increment_role = c(
      "residual_adds_dimension", "structural_zero",
      rep("estimated_nonadditive", 9), "structural_zero", "estimated_nonadditive",
      rep("structural_zero", 12)),
    citation_key = c("ref8_Wang2025GLM7", "ref4_GuerreroRomero2010TyG", "ref9_Er2016TyGBMI",
      "exploratory_TyG_times_waist", "exploratory_TyG_times_WHtR", "ref7_BelloChavolla2018METSIR",
      "ref6_Kahn2005LAP", "ref10_Amato2010VAI", "ref11_Xia2016CVAI", "ref12_Thomas2013BRI",
      "ref13_Krakauer2012ABSI", "ref14_Dobiasova2001AIP", "ref5_Matthews1985HOMA",
      "exploratory_CRP_albumin_ratio", "exploratory_monocyte_HDL_ratio", "exploratory_TyG_times_hsCRP",
      "exploratory_hsCRP_times_BMI", "exploratory_hsCRP_times_waist", "exploratory_SII_times_BMI",
      "standard_NLR_definition", "standard_PLR_definition", "ref15_Hu2014SII", "ref16_Qi2016SIRI",
      "ref17_Onodera1984PNI", "ref5_Matthews1985HOMA"),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# v32 symbolic layer: native-scale product expansions (distributive law).
# Each expansion is expressed in model-column names (the term_map targets),
# so "log10_tg" means log10(TG_mmol_L) and "sexFemale" is the dummy column of
# the factor term sex. Coefficients follow the published formulas exactly.
#
#   TyG_BMI = TyG * BMI
#   TyG_WC  = TyG * waist_cm
#   LAP     = (waist - 65) * TG [male] / (waist - 58) * TG [female]
#           = waist_cm:tg_mmol_L - 65 * tg_mmol_L + 7 * tg_mmol_L:sexFemale
#   CVAI    = f_M + (f_F - f_M) * sexFemale, with
#             f_F - f_M = 80.61 + 1.03 age + 4.20 BMI - 2.88 waist
#                         + 17.76 log10(TG) + 4.66 HDL
#   HOMA_IR = (1/22.5) * glucose_mmol_L : insulin_uIU_mL
# ---------------------------------------------------------------------------
# v33: this is part of the public API (see public_api_v33.R); the v32 name
# v32_product_expansions is kept below as a deprecated alias.
product_expansions <- function() {
  list(
    TyG_BMI = list(products = c("TyG:BMI" = 1)),
    TyG_WC = list(products = c("TyG:waist_cm" = 1)),
    LAP = list(linear = c(tg_mmol_L = -65),
               products = c("waist_cm:tg_mmol_L" = 1, "tg_mmol_L:sexFemale" = 7)),
    CVAI = list(intercept = -267.93,
                linear = c(age = 0.68, BMI = 0.03, waist_cm = 4.00,
                           log10_tg = 22.00, hdl_mmol_L = -16.32,
                           sexFemale = 80.61),
                products = c("age:sexFemale" = 1.03, "BMI:sexFemale" = 4.20,
                             "waist_cm:sexFemale" = -2.88,
                             "log10_tg:sexFemale" = 17.76,
                             "hdl_mmol_L:sexFemale" = 4.66)),
    HOMA_IR = list(products = c("glucose_mmol_L:insulin_uIU_mL" = 1 / 22.5))
  )
}

# v33: the four-class labels below are a HAND-WRITTEN reference answer, one
# per registry row. They are NOT the algorithm's output. registry_audit_v33.R
# derives the classification independently from each registered expression
# plus its full registered base model and compares the two; test_audit_v33.R
# requires the agreement to hold record by record. GLM7 is additionally
# dual-registered (its second record, against the full embedded base, is
# rank_preserving_exact), so the 25 rows yield 26 audit records:
# 15 rank_preserving_exact + 1 rank_adding_exact + 5 reducible_to_interaction
# + 5 outside_supported_class.
reference_audit_classes <- function() {
  c(GLM7 = "rank_adding_exact", TyG = "rank_preserving_exact",
    TyG_BMI = "reducible_to_interaction", TyG_WC = "reducible_to_interaction",
    TyG_WHtR = "outside_supported_class", METS_IR = "outside_supported_class",
    LAP = "reducible_to_interaction", VAI = "outside_supported_class",
    CVAI = "reducible_to_interaction", BRI = "outside_supported_class",
    ABSI = "outside_supported_class", AIP = "rank_preserving_exact",
    HOMA_IR = "reducible_to_interaction", CAR = "rank_preserving_exact",
    MHR = "rank_preserving_exact", TyG_hsCRP = "rank_preserving_exact",
    hsCRP_BMI = "rank_preserving_exact", hsCRP_WC = "rank_preserving_exact",
    SII_BMI = "rank_preserving_exact", NLR = "rank_preserving_exact",
    PLR = "rank_preserving_exact", SII = "rank_preserving_exact",
    SIRI = "rank_preserving_exact", PNI = "rank_preserving_exact",
    HOMA_IR_log = "rank_preserving_exact")
}

v33_registry <- function() {
  reg <- v28_registry()
  cls <- reference_audit_classes()
  stopifnot(all(reg$index %in% names(cls)))
  reg$reference_class <- unname(cls[reg$index])
  exps <- product_expansions()
  reg$added_columns <- vapply(reg$index, function(nm) {
    e <- exps[[nm]]
    if (is.null(e) || is.null(e$products)) return("")
    paste(names(e$products), collapse = ";")
  }, character(1))
  reg
}

# Deprecated v32 aliases, kept so earlier calling code keeps working.
v32_product_expansions <- product_expansions
v32_audit_classes <- reference_audit_classes
v32_registry <- v33_registry

compute_candidate_indices <- function(d) {
  # The caller must supply official-unit columns. Values outside each formula's
  # support domain are set to NA rather than imputed or winsorized.
  d$sex_male <- as.integer(as.character(d$sex) %in% c("Male", "1"))
  d$glucose_mg_dL <- d$FBG
  d$glucose_mmol_L <- d$FBG / 18
  d$tg_mmol_L <- d$TG
  d$tg_mg_dL <- d$TG * 88.57
  d$hdl_mmol_L <- d$HDL
  d$hdl_mg_dL <- d$HDL * 38.67
  d$insulin_pmol_L <- d$Insulin
  # GLM7 uses the released pmol/L field. HOMA-IR requires the source uU/mL
  # scale, recovered with the official factor of 6.0.
  d$insulin_uIU_mL <- d$Insulin / 6.0
  d$height_cm <- d$Height
  d$waist_cm <- d$Waist
  d$albumin_g_dL <- d$Albumin
  d$hsCRP_mg_L <- d$hsCRP
  d$monocytes <- d$Monocytes
  d$neutrophils <- d$Neutrophils
  d$lymphocytes <- d$Lymphocytes
  d$platelets <- d$Platelets
  d$SII <- d$platelets * d$neutrophils / d$lymphocytes
  d$log_age <- v28_safe_log10(d$age)
  d$log_bmi <- v28_safe_log10(d$BMI)
  d$metab <- v28_safe_log10(d$glucose_mg_dL) + v28_safe_log10(d$insulin_pmol_L) +
    v28_safe_log10(d$tg_mmol_L) + v28_safe_log10(d$LDL) - v28_safe_log10(d$hdl_mmol_L)
  d$GLM7 <- d$metab + d$log_age + d$log_bmi
  tyg_product <- d$tg_mg_dL * d$glucose_mg_dL / 2
  d$TyG <- ifelse(is.finite(tyg_product) & tyg_product > 0, log(tyg_product), NA_real_)
  d$log_TyG <- v28_safe_log(d$TyG)
  d$TyG_BMI <- d$TyG * d$BMI
  d$TyG_WC <- d$TyG * d$waist_cm
  d$TyG_WHtR <- d$TyG * (d$waist_cm / d$height_cm)
  d$METS_IR <- log(2 * d$glucose_mg_dL + d$tg_mg_dL) * d$BMI / log(d$hdl_mg_dL)
  d$LAP <- ifelse(d$sex_male == 1, d$waist_cm - 65, d$waist_cm - 58) * d$tg_mmol_L
  # Published formula implementations use the stated sex-specific constants
  # and the original measurement units. Unsupported rows remain NA.
  d$VAI <- ifelse(d$sex_male == 1,
                  (d$waist_cm / (39.68 + 1.88 * d$BMI)) *
                    (d$tg_mmol_L / 1.03) * (1.31 / d$hdl_mmol_L),
                  (d$waist_cm / (36.58 + 1.89 * d$BMI)) *
                    (d$tg_mmol_L / 0.81) * (1.52 / d$hdl_mmol_L))
  d$CVAI <- ifelse(d$sex_male == 1,
                   -267.93 + 0.68 * d$age + 0.03 * d$BMI + 4.00 * d$waist_cm +
                     22.00 * log10(d$tg_mmol_L) - 16.32 * d$hdl_mmol_L,
                   -187.32 + 1.71 * d$age + 4.23 * d$BMI + 1.12 * d$waist_cm +
                     39.76 * log10(d$tg_mmol_L) - 11.66 * d$hdl_mmol_L)
  ratio <- d$waist_cm / (2 * pi)
  bri_radicand <- 1 - (ratio^2 / (0.5 * d$height_cm)^2)
  d$BRI <- ifelse(is.finite(bri_radicand) & bri_radicand >= 0,
                  364.2 - 365.5 * sqrt(bri_radicand), NA_real_)
  d$ABSI <- (d$waist_cm / 100) /
    (d$BMI^(2/3) * (d$height_cm / 100)^(1/2))
  d$AIP <- log10(d$tg_mmol_L / d$hdl_mmol_L)
  d$HOMA_IR <- d$glucose_mmol_L * d$insulin_uIU_mL / 22.5
  d$HOMA_IR_log <- v28_safe_log(d$HOMA_IR)
  d$CAR <- d$hsCRP_mg_L / d$albumin_g_dL
  d$MHR <- d$monocytes / d$hdl_mg_dL
  d$TyG_hsCRP <- d$TyG * d$hsCRP_mg_L
  d$hsCRP_BMI <- d$hsCRP_mg_L * d$BMI
  d$hsCRP_WC <- d$hsCRP_mg_L * d$waist_cm
  d$SII_BMI <- d$SII * d$BMI
  d$NLR <- d$neutrophils / d$lymphocytes
  d$PLR <- d$platelets / d$lymphocytes
  d$SIRI <- d$neutrophils * d$monocytes / d$lymphocytes
  d$lymphocytes_cells_mm3 <- d$lymphocytes * 1000
  d$PNI <- 10 * d$albumin_g_dL + 0.005 * d$lymphocytes_cells_mm3
  d$log_glucose <- v28_safe_log(d$glucose_mg_dL)
  d$log_glucose_mmol <- v28_safe_log(d$glucose_mmol_L)
  d$log_tg <- v28_safe_log(d$tg_mmol_L)
  d$log10_tg <- v28_safe_log10(d$tg_mmol_L)
  d$log_tg_mg <- v28_safe_log(d$tg_mg_dL)
  d$log_glucose_mg <- v28_safe_log(d$glucose_mg_dL)
  d$log_hdl <- v28_safe_log(d$hdl_mmol_L)
  d$log_hdl_mg <- v28_safe_log(d$hdl_mg_dL)
  d$log_insulin <- v28_safe_log(d$insulin_uIU_mL)
  d$log_mets_sum <- v28_safe_log(2 * d$glucose_mg_dL + d$tg_mg_dL)
  d$log_waist <- v28_safe_log(d$waist_cm)
  d$log_height <- v28_safe_log(d$height_cm)
  d$ln_bmi <- v28_safe_log(d$BMI)
  d$log_SII <- v28_safe_log(d$SII)
  d$log_hscrp <- v28_safe_log(d$hsCRP_mg_L)
  d$log_albumin <- v28_safe_log(d$albumin_g_dL)
  d$log_monocytes <- v28_safe_log(d$monocytes)
  d$log_platelets <- v28_safe_log(d$platelets)
  d$log_neutrophils <- v28_safe_log(d$neutrophils)
  d$log_lymphocytes <- v28_safe_log(d$lymphocytes)
  d
}

audit_candidate_registry <- function(registry = v28_registry(), data = NULL,
                                     tolerances = c(1e-8, 1e-10, 1e-12)) {
  stopifnot(is.data.frame(registry), all(c("index", "formula", "audit_class") %in% names(registry)))
  out <- registry
  out$rank_base <- NA_integer_; out$rank_augmented <- NA_integer_; out$rank_residual <- NA_integer_
  out$rank_gain <- NA_integer_; out$condition_number <- NA_real_; out$available_P <- NA; out$available_L <- NA
  out$tol_1e8 <- NA_character_; out$tol_1e10 <- NA_character_; out$tol_1e12 <- NA_character_
  if (!is.null(data)) for (i in seq_len(nrow(out))) {
    nm <- out$index[i]
    out$available_P[i] <- nm %in% names(data)
    out$available_L[i] <- nm %in% names(data)
    if (!isTRUE(out$available_P[i])) next
    term_value <- function(term) {
      if (term == "sex") return(data$sex)
      if (grepl("^log10\\(", term)) {
        inner <- sub("^log10\\((.*)\\)$", "\\1", term)
        z <- switch(inner, age = data$age, BMI = data$BMI,
                     TG_mmol_L = data$tg_mmol_L,
                     stop("unknown audit log10 term: ", term))
        return(v28_safe_log10(z))
      }
      if (grepl("^log\\(", term)) {
        inner <- sub("^log\\((.*)\\)$", "\\1", term)
        z <- switch(inner,
          TG_mmol_L = data$tg_mmol_L, HDL_mmol_L = data$hdl_mmol_L,
          TG_mg_dL = data$tg_mg_dL, glucose_mg_dL = data$glucose_mg_dL,
          glucose_mmol_L = data$glucose_mmol_L, insulin_uIU_mL = data$insulin_uIU_mL,
          hsCRP_mg_L = data$hsCRP_mg_L, albumin_g_dL = data$albumin_g_dL,
          monocytes = data$monocytes, HDL_mg_dL = data$hdl_mg_dL,
          neutrophils = data$neutrophils, lymphocytes = data$lymphocytes,
          platelets = data$platelets, SII = data$SII, BMI = data$BMI,
          TyG = data$TyG, waist_cm = data$waist_cm, height_cm = data$height_cm,
          `2*glucose_mg_dL+TG_mg_dL` = 2 * data$glucose_mg_dL + data$tg_mg_dL,
          stop("unknown audit log term: ", term))
        return(v28_safe_log(z))
      }
      switch(term,
        log10_age = data$log_age, log10_BMI = data$log_bmi,
        TyG = data$TyG, SII = data$SII,
        age = data$age, BMI = data$BMI, waist_cm = data$waist_cm,
        height_cm = data$height_cm, glucose_mmol_L = data$glucose_mmol_L,
        glucose_mg_dL = data$glucose_mg_dL, TG_mmol_L = data$tg_mmol_L,
        TG_mg_dL = data$tg_mg_dL, HDL_mmol_L = data$hdl_mmol_L,
        HDL_mg_dL = data$hdl_mg_dL,
        insulin_uIU_mL = data$insulin_uIU_mL, hsCRP_mg_L = data$hsCRP_mg_L,
        albumin_g_dL = data$albumin_g_dL, lymphocytes_cells_mm3 = data$lymphocytes_cells_mm3,
        stop("unknown audit term: ", term))
    }
    terms <- unique(strsplit(out$component_terms[i], ";", fixed = TRUE)[[1]])
    vals <- lapply(terms, term_value); names(vals) <- paste0("z", seq_along(vals))
    score <- data[[nm]]
    if (identical(out$model_scale[i], "log") && nm != "HOMA_IR_log") score <- v28_safe_log(score)
    vals$composite <- score
    ok <- complete.cases(as.data.frame(vals))
    if (sum(ok) < 20) next
    base_df <- data.frame(vals[setdiff(names(vals), "composite")], check.names = FALSE)
    aug_df <- data.frame(vals, check.names = FALSE)
    xb <- stats::model.matrix(~ ., data = base_df[ok, , drop = FALSE])
    xa <- stats::model.matrix(~ ., data = aug_df[ok, , drop = FALSE])
    rank_fun <- function(x, tol) qr(x, tol = tol)$rank
    out$rank_base[i] <- rank_fun(xb, tolerances[[2]])
    out$rank_augmented[i] <- rank_fun(xa, tolerances[[2]])
    out$rank_gain[i] <- out$rank_augmented[i] - out$rank_base[i]
    scaled <- if (ncol(xb) > 1L) scale(xb[, -1, drop = FALSE], center = TRUE, scale = TRUE) else xb
    out$condition_number[i] <- tryCatch(kappa(scaled), error = function(e) NA_real_)
    out$tol_1e8[i] <- as.character(rank_fun(xa, tolerances[[1]]) - rank_fun(xb, tolerances[[1]]))
    out$tol_1e10[i] <- as.character(rank_fun(xa, tolerances[[2]]) - rank_fun(xb, tolerances[[2]]))
    out$tol_1e12[i] <- as.character(rank_fun(xa, tolerances[[3]]) - rank_fun(xb, tolerances[[3]]))
  }
  out
}

v28_component_model_terms <- function(index, registry = v28_registry()) {
  row <- registry[registry$index == index, , drop = FALSE]
  if (nrow(row) != 1L) stop("index not found or duplicated in registry: ", index)
  symbolic <- unique(strsplit(row$component_terms[[1]], ";", fixed = TRUE)[[1]])
  term_map <- c(
    sex = "sex", `log10(age)` = "log_age", `log10(BMI)` = "log_bmi",
    `log(TG_mg_dL)` = "log_tg_mg", `log(glucose_mg_dL)` = "log_glucose_mg",
    TyG = "TyG", BMI = "BMI", waist_cm = "waist_cm", height_cm = "height_cm",
    `log(2*glucose_mg_dL+TG_mg_dL)` = "log_mets_sum",
    `log(HDL_mg_dL)` = "log_hdl_mg", TG_mmol_L = "tg_mmol_L",
    HDL_mmol_L = "hdl_mmol_L", age = "age", `log10(TG_mmol_L)` = "log10_tg",
    glucose_mmol_L = "glucose_mmol_L", insulin_uIU_mL = "insulin_uIU_mL",
    `log(TG_mmol_L)` = "log_tg", `log(HDL_mmol_L)` = "log_hdl",
    `log(hsCRP_mg_L)` = "log_hscrp", `log(albumin_g_dL)` = "log_albumin",
    `log(monocytes)` = "log_monocytes", `log(TyG)` = "log_TyG",
    `log(BMI)` = "ln_bmi", `log(waist_cm)` = "log_waist", `log(SII)` = "log_SII",
    `log(neutrophils)` = "log_neutrophils", `log(lymphocytes)` = "log_lymphocytes",
    `log(platelets)` = "log_platelets", albumin_g_dL = "albumin_g_dL",
    lymphocytes_cells_mm3 = "lymphocytes_cells_mm3",
    `log(glucose_mmol_L)` = "log_glucose_mmol", `log(insulin_uIU_mL)` = "log_insulin"
  )
  missing <- setdiff(symbolic, names(term_map))
  if (length(missing)) stop("unmapped registry component term(s) for ", index, ": ", paste(missing, collapse = ", "))
  unname(term_map[symbolic])
}

v28_component_increment_is_structural_zero <- function(index, registry = v28_registry()) {
  row <- registry[registry$index == index, , drop = FALSE]
  if (nrow(row) != 1L) stop("index not found or duplicated in registry: ", index)
  identical(row$component_increment_role[[1]], "structural_zero")
}

build_component_conditional_model <- function(index, registry = v28_registry(), outcome = "dx") {
  row <- registry[registry$index == index, , drop = FALSE]
  if (!nrow(row)) stop("index not found in registry: ", index)
  terms <- v28_component_model_terms(index, registry)
  rhs <- paste(c(terms, index), collapse = " + ")
  list(index = index, base_formula = stats::as.formula(paste(outcome, "~", paste(terms, collapse = " + "))),
       augmented_formula = stats::as.formula(paste(outcome, "~", rhs)),
       component_terms = terms, analysis_role = row$analysis_role[[1]],
       component_increment_role = row$component_increment_role[[1]],
       structural_zero = v28_component_increment_is_structural_zero(index, registry))
}

evaluate_candidate_family <- function(index, data, design, registry = v28_registry(), ...) {
  spec <- build_component_conditional_model(index, registry = registry)
  if (!index %in% names(data)) stop("computed index is missing: ", index)
  list(index = index, family = registry$family[match(index, registry$index)],
       base_formula = deparse(spec$base_formula), augmented_formula = deparse(spec$augmented_formula),
       status = "delegated_to_evaluate_incremental_model", analysis_role = spec$analysis_role)
}
