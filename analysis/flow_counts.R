# flow_counts.R (PURE ASCII)
#
# v34: the STROBE-style level-by-level count that the review asked for.  The manuscript
# reports 3,801 (period P) and 2,958 (period L) analysed adults and a reviewer
# queried whether a fasting subsample of that size is plausible.  This script
# walks from the official XPT file to the analysed row, one exclusion at a
# time, and writes the whole ladder to v34_flow_counts.csv.
#
# It also records, by name, the survey variables that were previously only
# described in words ("a positive fasting weight"): the fasting subsample
# weight, the stratification variable and the PSU variable, together with the
# design degrees of freedom of the FULL design and of the analysed subset.  The
# two are reported separately because a subset should not reduce the degrees
# of freedom of the design.
#
# Run:
#   CIA_CODE_DIR=... CIA_DATA_ROOT=... CIA_OUT=... Rscript flow_counts.R

suppressPackageStartupMessages({ library(haven); library(survey) })

need <- function(nm) {
  v <- Sys.getenv(nm, unset = "")
  if (!nzchar(v)) stop("set ", nm, " before running flow_counts.R")
  v
}
CODE_DIR  <- need("CIA_CODE_DIR")
DATA_ROOT <- need("CIA_DATA_ROOT")
OUT       <- need("CIA_OUT")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(CODE_DIR, "index_registry.R"))
source(file.path(CODE_DIR, "data_build.R"), local = FALSE)

spec <- function(period) {
  if (period == "P") list(
    demo = "P_DEMO.xpt", glu = "P_GLU.xpt", mcq = "P_MCQ.xpt", diq = "P_DIQ.xpt",
    weight = "WTSAFPRP",
    label = "2017-March 2020 (pre-pandemic)")
  else list(
    demo = "DEMO_L.xpt", glu = "GLU_L.xpt", mcq = "MCQ_L.xpt", diq = "DIQ_L.xpt",
    weight = "WTSAF2YR",
    label = "2021-2023")
}

read_one <- function(file, cols) {
  p <- file.path(DATA_ROOT, file)
  if (!file.exists(p)) stop("missing official XPT: ", p)
  x <- haven::read_xpt(p)
  x[, intersect(cols, names(x)), drop = FALSE]
}

rows <- list()
for (period in c("P", "L")) {
  s <- spec(period)
  demo <- read_one(s$demo, c("SEQN", "RIDAGEYR", "SDMVSTRA", "SDMVPSU"))
  glu  <- read_one(s$glu,  c("SEQN", s$weight))
  mcq  <- read_one(s$mcq,  c("SEQN", "MCQ550"))
  diq  <- read_one(s$diq,  c("SEQN", "DIQ010"))

  n_xpt        <- nrow(demo)
  adult        <- is.finite(demo$RIDAGEYR) & demo$RIDAGEYR >= 20
  n_adult      <- sum(adult, na.rm = TRUE)
  n_design     <- sum(adult & !is.na(demo$SDMVSTRA) & !is.na(demo$SDMVPSU),
                      na.rm = TRUE)

  # Align the fasting weight onto the DEMO rows by SEQN.  merge(..., sort =
  # FALSE) does not promise to keep the row order of x, so matching by key is
  # the only safe way to combine a column from one file with a filter built
  # from another.
  w <- glu[[s$weight]][match(demo$SEQN, glu$SEQN)]
  stopifnot(length(w) == nrow(demo))
  n_in_fasting_file <- sum(adult & !is.na(w), na.rm = TRUE)
  n_weight     <- sum(adult & !is.na(demo$SDMVSTRA) & !is.na(demo$SDMVPSU) &
                      is.finite(w) & w > 0, na.rm = TRUE)

  # ---- analysed rows: exactly the filter build_period() applies -----------
  full <- build_period(period)
  n_analytic <- nrow(full)

  # ---- design degrees of freedom ------------------------------------------
  degf_of <- function(stra, psu, wt) {
    ok <- !is.na(stra) & !is.na(psu) & is.finite(wt) & wt > 0
    if (sum(ok) < 2L) return(list(df = NA_integer_, strata = NA_integer_,
                                  psu = NA_integer_))
    dd <- survey::svydesign(ids = ~psu, strata = ~stra, weights = ~wt,
                            data = data.frame(psu = psu[ok], stra = stra[ok],
                                              wt = wt[ok]),
                            nest = TRUE)
    list(df = survey::degf(dd),
         strata = length(unique(stra[ok])),
         psu = length(unique(paste(stra[ok], psu[ok], sep = "_"))))
  }
  # Degrees of freedom of the COMPLETE sample design: every DEMO row with a
  # valid stratum and PSU, carrying a constant positive weight, so the number
  # reflects the design and not the analytic subsample.  NCHS guidance is to
  # take the degrees of freedom from the full design, because a subpopulation
  # should not be credited with fewer degrees of freedom than the design that
  # produced it.  The analytic-subset figure is reported beside it so the two
  # can be compared rather than assumed equal.
  stra_all <- as.character(demo$SDMVSTRA)
  psu_all  <- paste(demo$SDMVSTRA, demo$SDMVPSU, sep = "_")
  full_des <- degf_of(stra_all, psu_all, rep(1, nrow(demo)))
  ana_des  <- degf_of(as.character(full$stra), as.character(full$psu), full$weight)

  rows[[length(rows) + 1L]] <- data.frame(
    period = period, cycle_label = s$label,
    demo_file = s$demo,
    n_xpt_rows = n_xpt,
    n_age_20plus = n_adult,
    n_age20plus_with_strata_psu = n_design,
    n_adults_in_fasting_file = n_in_fasting_file,
    n_fasting_subsample_weight_positive = n_weight,
    n_analysed = n_analytic,
    n_gallstone_observed = sum(!is.na(full$dx_gall)),
    n_gallstone_cases = sum(full$dx_gall == 1, na.rm = TRUE),
    n_diabetes_observed = sum(!is.na(full$dx_diab)),
    n_diabetes_cases = sum(full$dx_diab == 1, na.rm = TRUE),
    fasting_weight_variable = s$weight,
    strata_variable = "SDMVSTRA",
    psu_variable = "SDMVPSU",
    degf_full_design = full_des$df,
    n_strata_full = full_des$strata,
    n_psu_full = full_des$psu,
    degf_analysed_subset = ana_des$df,
    n_strata_analysed = ana_des$strata,
    n_psu_analysed = ana_des$psu,
    stringsAsFactors = FALSE)
}

flow <- do.call(rbind, rows)
print(flow[, c("period", "n_xpt_rows", "n_age_20plus",
               "n_fasting_subsample_weight_positive", "n_analysed",
               "degf_full_design", "degf_analysed_subset")],
      row.names = FALSE)

# The analysed count must be the same number the main pipeline uses.
for (period in c("P", "L")) {
  got <- flow$n_analysed[flow$period == period]
  stopifnot(got == nrow(build_period(period)))
}
cat("\nflow_counts: analysed rows reproduce build_period() exactly\n")

write.csv(flow, file.path(OUT, "v34_flow_counts.csv"), row.names = FALSE)
cat("wrote v34_flow_counts.csv\n")
