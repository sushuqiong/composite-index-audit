# data_build.R (PURE ASCII)
# Canonical NHANES P/L loader for every v34 script.  The body of
# build_period() is byte-identical to the copy inside
# build_results.R; the package audit checks the two copies match.
#
# v34: no absolute paths are hard-coded.  Point CIA_DATA_ROOT at the
# directory containing the official NHANES XPT files before sourcing; the
# loader fails closed if the variable is unset.

suppressPackageStartupMessages(library(haven))

DATA_ROOT <- Sys.getenv("CIA_DATA_ROOT", unset = "")
if (!nzchar(DATA_ROOT))
  stop("set CIA_DATA_ROOT to the directory containing the NHANES XPT files")
if (!dir.exists(DATA_ROOT))
  stop("CIA_DATA_ROOT does not exist: ", DATA_ROOT)

read_unique <- function(file, cols) {
  p <- file.path(DATA_ROOT, file)
  if (!file.exists(p)) stop("missing official XPT: ", p)
  x <- haven::read_xpt(p)
  if (!"SEQN" %in% names(x) || anyDuplicated(x$SEQN))
    stop("SEQN is not unique: ", file)
  x[, cols, drop = FALSE]
}

build_period <- function(period = c("P", "L")) {
  period <- match.arg(period)
  p <- if (period == "P") {
    list(
      demo   = list(file = "P_DEMO.xpt",   cols = c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH3","SDMVPSU","SDMVSTRA")),
      glu    = list(file = "P_GLU.xpt",    cols = c("SEQN","WTSAFPRP","LBXGLU")),
      ins    = list(file = "P_INS.xpt",    cols = c("SEQN","LBDINSI")),
      trig   = list(file = "P_TRIGLY.xpt", cols = c("SEQN","LBDTRSI","LBDLDLSI")),
      hdl    = list(file = "P_HDL.XPT",    cols = c("SEQN","LBDHDDSI")),
      bmx    = list(file = "P_BMX.xpt",    cols = c("SEQN","BMXBMI","BMXWAIST","BMXHT")),
      biopro = list(file = "P_BIOPRO.xpt", cols = c("SEQN","LBXSAL")),
      hscrp  = list(file = "P_HSCRP.XPT",  cols = c("SEQN","LBXHSCRP")),
      cbc    = list(file = "P_CBC.XPT",    cols = c("SEQN","LBDMONO","LBDNENO","LBDLYMNO","LBXPLTSI","LBXWBCSI")),
      mcq    = list(file = "P_MCQ.xpt",    cols = c("SEQN","MCQ550")),
      diq    = list(file = "P_DIQ.xpt",    cols = c("SEQN","DIQ010")),
      weight = "WTSAFPRP")
  } else {
    list(
      demo   = list(file = "DEMO_L.xpt",   cols = c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH3","SDMVPSU","SDMVSTRA")),
      glu    = list(file = "GLU_L.xpt",    cols = c("SEQN","WTSAF2YR","LBXGLU")),
      ins    = list(file = "INS_L.xpt",    cols = c("SEQN","LBDINSI")),
      trig   = list(file = "TRIGLY_L.xpt", cols = c("SEQN","LBDTRSI","LBDLDLSI")),
      hdl    = list(file = "HDL_L.XPT",    cols = c("SEQN","LBDHDDSI")),
      bmx    = list(file = "BMX_L.xpt",    cols = c("SEQN","BMXBMI","BMXWAIST","BMXHT")),
      biopro = list(file = "BIOPRO_L.xpt", cols = c("SEQN","LBXSAL")),
      hscrp  = list(file = "HSCRP_L.XPT",  cols = c("SEQN","LBXHSCRP")),
      cbc    = list(file = "CBC_L.xpt",    cols = c("SEQN","LBDMONO","LBDNENO","LBDLYMNO","LBXPLTSI","LBXWBCSI")),
      mcq    = list(file = "MCQ_L.xpt",    cols = c("SEQN","MCQ550")),
      diq    = list(file = "DIQ_L.xpt",    cols = c("SEQN","DIQ010")),
      weight = "WTSAF2YR")
  }
  weight_name <- p$weight
  p$weight <- NULL
  dat <- lapply(p, function(z) read_unique(z$file, z$cols))
  d <- Reduce(function(x, y) merge(x, y, by = "SEQN", all.x = TRUE,
                                   sort = FALSE), dat)
  rn <- c(RIAGENDR = "gender", RIDAGEYR = "age", LBXGLU = "FBG",
          LBDINSI = "Insulin", LBDTRSI = "TG", LBDLDLSI = "LDL",
          LBDHDDSI = "HDL", BMXBMI = "BMI", BMXWAIST = "Waist",
          BMXHT = "Height", LBXSAL = "Albumin", LBXHSCRP = "hsCRP",
          LBDMONO = "Monocytes", LBDNENO = "Neutrophils",
          LBDLYMNO = "Lymphocytes", LBXPLTSI = "Platelets")
  names(d)[names(d) %in% names(rn)] <- rn[names(d)[names(d) %in% names(rn)]]
  d$sex <- factor(ifelse(d$gender == 1, "Male",
                         ifelse(d$gender == 2, "Female", NA_character_)),
                  levels = c("Male", "Female"))
  d$dx_gall <- ifelse(d$MCQ550 %in% c(1, 2), as.integer(d$MCQ550 == 1),
                      NA_integer_)
  d$dx_diab <- ifelse(d$DIQ010 %in% c(1, 2), as.integer(d$DIQ010 == 1),
                      NA_integer_)
  d$weight <- d[[weight_name]]
  d$stra <- as.character(d$SDMVSTRA)
  d$psu <- paste(d$stra, as.character(d$SDMVPSU), sep = "_")
  d <- compute_candidate_indices(d)
  d$period <- period
  ok <- is.finite(d$weight) & d$weight > 0 & !is.na(d$SDMVSTRA) &
    !is.na(d$SDMVPSU) & is.finite(d$age) & d$age >= 20
  d[ok, , drop = FALSE]
}
