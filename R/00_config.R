#' Configuration and shared settings
#'
#' Single place where the NCA business rules are declared. Every engine in this
#' repo (PKNCA, NonCompart, Phoenix WinNonlin) is configured from these values so
#' that any observed difference is attributable to the engine, not the settings.
#'
#' This is the whole point of the exercise: NCA engines agree on arithmetic and
#' disagree on conventions. Pin the conventions and the comparison becomes
#' meaningful.

# ---- Paths ------------------------------------------------------------------
PATHS <- list(
  root       = normalizePath(".", mustWork = FALSE),
  derived    = file.path("data", "derived"),
  outputs    = file.path("outputs"),
  phoenix    = file.path("phoenix"),
  phoenix_in = file.path("phoenix", "input"),
  phoenix_rf = file.path("phoenix", "reference")
)

for (p in PATHS) dir.create(p, showWarnings = FALSE, recursive = TRUE)

# ---- NCA business rules (the settings contract) -----------------------------
#
# These map onto the equivalent controls in each engine:
#
#   Rule                 PKNCA                       NonCompart      Phoenix WinNonlin
#   -------------------  --------------------------  --------------  ---------------------------
#   AUC method           auc.method                  down=           Calculation method
#   lambda-z selection   adj.r.squared.factor        R2ADJ=          Lambda Z / Best Fit
#   min lambda-z points  min.hl.points               (fixed 3)       Min. no. of points
#   Tmax in lambda-z     allow.tmax.in.half.life     (excluded)      Exclude Tmax
#   AUC%extrap flag      max.aucinf.pext             (reported)      Acceptance criteria
#
NCA_RULES <- list(
  # "lin up/log down" is the modern default and matches WinNonlin's
  # "Linear Log Trapezoidal". Use "linear" to match WinNonlin's "Linear
  # Trapezoidal (Linear/Log interpolation)" with linear down.
  auc_method        = "lin up/log down",

  # NonCompart's equivalent switch. "Linear" or "Log".
  noncompart_down   = "Log",

  # Minimum adjusted R-squared improvement to accept an extra terminal point.
  # WinNonlin's Best Fit uses 0.0001; PKNCA's default is the same idea.
  adj_r2_factor     = 0.0001,

  # NonCompart's R2ADJ threshold: lambda-z is reported only if adjusted R^2
  # exceeds this. 0.7 is the package default; 0.8 is stricter.
  r2adj_threshold   = 0.7,

  # Minimum number of points in the terminal regression.
  min_hl_points     = 3L,

  # Whether Cmax may be included in the lambda-z regression. WinNonlin's Best
  # Fit excludes it by default; keep both engines aligned.
  allow_tmax_in_hl  = FALSE,

  # Acceptance threshold for percent of AUCinf that was extrapolated.
  max_aucinf_pext   = 20,

  # BLQ handling, by region of the profile.
  #   "drop" | "keep" | a number
  #
  # CAUTION: PKNCA classifies any concentration of exactly 0 as BLQ. In
  # datasets::Theoph, 9 of 12 subjects have conc = 0 at time zero -- genuine
  # pre-dose measurements, not censored values. Setting blq_first = "drop"
  # deletes those records, which pushes the first measurement to 0.25 h, makes
  # the interval start (0) precede it, and causes PKNCA to return NA for every
  # AUC-derived parameter on those 9 subjects. WinNonlin keeps the zero and
  # integrates the 0 -> 0.25 h trapezoid.
  #
  # "keep" is both PKNCA's default and the WinNonlin-matching choice.
  blq_first         = "keep",
  blq_middle        = 0,
  blq_last          = "drop",

  # Numerical tolerance for declaring two engines equivalent. 0.001 = 0.1%,
  # which is the tolerance used in the published NonCompart/WinNonlin
  # validation.
  tolerance         = 0.001
)

# ---- Theophylline dose convention -------------------------------------------
#
# datasets::Theoph carries Dose in mg/kg and Wt in kg, so the true total dose
# varies slightly by subject (318-320 mg). The published WinNonlin reference
# output that this repo compares against was generated with a flat 320 mg for
# every subject. To keep the comparison valid we adopt the same convention.
#
# This is a comparability choice, not a PK recommendation. For real analyses,
# use the subject-specific dose.
THEOPH_DOSE_MG <- 320

# ---- External reference data ------------------------------------------------
#
# Phoenix WinNonlin 6.3/7.0 output for datasets::Theoph, published as part of
# the NonCompart validation package by Sungpil Han.
#
# Fetched at runtime and gitignored rather than vendored into this repository.
# That is deliberate: as of checking, the source repository carries no LICENSE
# file, so redistribution terms are not declared. Fetching keeps provenance
# with the source and avoids redistributing someone else's files. See
# README "Data provenance and licensing".
#
#   Kim H, Han S, Cho YS, Yoon SK, Bae KS. Development of R packages
#   'NonCompart' and 'ncar' for noncompartmental analysis (NCA).
#   Transl Clin Pharmacol. 2018;26(1):10-15.
#
WINNONLIN_REFERENCE <- list(
  theoph_linear = paste0(
    "https://raw.githubusercontent.com/asancpt/NonCompart-tests/master/",
    "Final_Parameters_Pivoted_Theoph_Linear.csv"
  ),
  theoph_log = paste0(
    "https://raw.githubusercontent.com/asancpt/NonCompart-tests/master/",
    "Final_Parameters_Pivoted_Theoph_Log.csv"
  ),
  indometh_bolus_linear = paste0(
    "https://raw.githubusercontent.com/asancpt/NonCompart-tests/master/",
    "Final_Parameters_Pivoted_Indometh_Linear.csv"
  )
)

# ---- Parameter name crosswalk -----------------------------------------------
#
# The three engines name the same quantities differently. CDISC PP domain
# PPTESTCD terms are used as the canonical key.
#
PARAM_MAP <- data.frame(
  PPTESTCD  = c("CMAX", "TMAX", "CLST", "TLST", "AUCLST", "AUCIFO",
                "AUCPEO", "LAMZ", "LAMZHL", "LAMZNPT", "R2ADJ",
                "CLFO", "VZFO", "AUMCLST", "AUMCIFO", "MRTEVIFO"),
  label     = c("Cmax", "Tmax", "Clast", "Tlast", "AUClast", "AUCinf(obs)",
                "AUC %extrap(obs)", "Lambda_z", "t1/2", "N points lambda_z",
                "R2 adjusted", "CL/F(obs)", "Vz/F(obs)", "AUMClast",
                "AUMCinf(obs)", "MRTinf(obs)"),
  pknca     = c("cmax", "tmax", "clast.obs", "tlast", "auclast", "aucinf.obs",
                "aucpext.obs", "lambda.z", "half.life", "lambda.z.n.points",
                "adj.r.squared", "cl.obs", "vz.obs", "aumclast", "aumcinf.obs",
                "mrt.obs"),
  noncompart = c("CMAX", "TMAX", "CLST", "TLST", "AUCLST", "AUCIFO",
                 "AUCPEO", "LAMZ", "LAMZHL", "LAMZNPT", "R2ADJ",
                 "CLFO", "VZFO", "AUMCLST", "AUMCIFO", "MRTEVIFO"),
  winnonlin = c("Cmax", "Tmax", "Clast", "Tlast", "AUClast", "AUCINF_obs",
                "AUC_.Extrap_obs", "Lambda_z", "HL_Lambda_z",
                "No_points_lambda_z", "Rsq_adjusted", "Cl_F_obs", "Vz_F_obs",
                "AUMClast", "AUMCINF_obs", "MRTINF_obs"),
  stringsAsFactors = FALSE
)

# ---- Session capture --------------------------------------------------------
capture_provenance <- function(file = file.path(PATHS$outputs, "provenance.txt")) {
  pkgs <- c("PKNCA", "NonCompart", "dplyr", "tidyr", "ggplot2")
  vers <- vapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE))
      as.character(utils::packageVersion(p)) else "not installed"
  }, character(1))

  lines <- c(
    "NCA engine comparison - provenance record",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("R: ", R.version.string),
    paste0("Platform: ", R.version$platform),
    "",
    "Package versions:",
    paste0("  ", names(vers), ": ", vers),
    "",
    "NCA rules in force:",
    paste0("  auc_method: ", NCA_RULES$auc_method),
    paste0("  noncompart_down: ", NCA_RULES$noncompart_down),
    paste0("  min_hl_points: ", NCA_RULES$min_hl_points),
    paste0("  allow_tmax_in_hl: ", NCA_RULES$allow_tmax_in_hl),
    paste0("  r2adj_threshold: ", NCA_RULES$r2adj_threshold),
    paste0("  tolerance: ", NCA_RULES$tolerance),
    paste0("  theoph_dose_mg: ", THEOPH_DOSE_MG)
  )
  writeLines(lines, file)
  invisible(lines)
}

message("Config loaded. AUC method: ", NCA_RULES$auc_method,
        " | tolerance: ", NCA_RULES$tolerance)
