#' NonCompart engine
#'
#' NonCompart is written to reproduce Phoenix WinNonlin's NCA conventions --
#' its lambda-z selection criterion, its AUC methods, and its CDISC PP domain
#' parameter names. A published validation package reports agreement with
#' WinNonlin 6.3 and 7.0 to within 0.1% on Theophylline and Indomethacin.
#'
#' That makes it the useful middle term in this comparison: it stands in for
#' WinNonlin when no license is available, and it is the thing to check first
#' when PKNCA and Phoenix disagree.
#'
#' Reference:
#'   Kim H, Han S, Cho YS, Yoon SK, Bae KS. Development of R packages
#'   'NonCompart' and 'ncar' for noncompartmental analysis (NCA).
#'   Transl Clin Pharmacol. 2018;26(1):10-15.

source(file.path("R", "00_config.R"))

run_noncompart <- function(conc,
                           dose_mg,
                           subject_col = "Subject",
                           time_col    = "Time",
                           conc_col    = "conc",
                           adm         = "Extravascular",
                           dur         = 0,
                           rules       = NCA_RULES) {

  if (!requireNamespace("NonCompart", quietly = TRUE))
    stop("NonCompart is required: install.packages('NonCompart')")

  res <- NonCompart::tblNCA(
    concData  = as.data.frame(conc),
    key       = subject_col,
    colTime   = time_col,
    colConc   = conc_col,
    dose      = dose_mg,
    adm       = adm,
    dur       = dur,
    doseUnit  = "mg",
    timeUnit  = "h",
    concUnit  = "mg/L",
    down      = rules$noncompart_down,
    R2ADJ     = rules$r2adj_threshold
  )

  res <- as.data.frame(res, stringsAsFactors = FALSE)

  # tblNCA returns wide, with PPTESTCD terms as column names already.
  id <- res[[subject_col]]
  param_cols <- intersect(names(res), PARAM_MAP$noncompart)

  out <- do.call(rbind, lapply(param_cols, function(p) {
    data.frame(
      Subject  = id,
      PPTESTCD = PARAM_MAP$PPTESTCD[match(p, PARAM_MAP$noncompart)],
      value    = suppressWarnings(as.numeric(res[[p]])),
      engine   = "NonCompart",
      stringsAsFactors = FALSE
    )
  }))

  out$Subject <- suppressWarnings(as.numeric(as.character(out$Subject)))
  rownames(out) <- NULL
  out
}

#' Map the AUC-method vocabulary between engines
#'
#' PKNCA and NonCompart use different words for the same two integration rules.
#' Getting this wrong is the single most common cause of a spurious AUC
#' difference, so it gets its own function rather than living as a magic string.
auc_method_crosswalk <- function(pknca_method) {
  switch(
    pknca_method,
    "linear"           = "Linear",   # linear up, linear down
    "lin up/log down"  = "Log",      # linear up, log down
    "lin-log"          = "Log",
    stop("Unmapped AUC method: ", pknca_method)
  )
}

#' Assert the two engines are configured consistently before comparing them
assert_rules_consistent <- function(rules = NCA_RULES) {
  expected <- auc_method_crosswalk(rules$auc_method)
  if (!identical(expected, rules$noncompart_down)) {
    stop(sprintf(
      paste0("AUC method mismatch: PKNCA '%s' maps to NonCompart down='%s', ",
             "but config says down='%s'. Fix NCA_RULES in 00_config.R."),
      rules$auc_method, expected, rules$noncompart_down))
  }
  invisible(TRUE)
}
