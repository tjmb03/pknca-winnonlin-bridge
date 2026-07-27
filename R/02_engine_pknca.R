#' PKNCA engine
#'
#' Runs NCA through PKNCA with the settings declared in 00_config.R and returns
#' a tidy per-subject frame keyed on CDISC PPTESTCD terms, so it can be lined up
#' against NonCompart and Phoenix WinNonlin output.

source(file.path("R", "00_config.R"))

run_pknca <- function(conc, dose,
                      subject_col = "Subject",
                      time_col    = "Time",
                      conc_col    = "conc",
                      dose_col    = "Dose",
                      rules       = NCA_RULES) {

  if (!requireNamespace("PKNCA", quietly = TRUE))
    stop("PKNCA is required: install.packages('PKNCA')")

  # Options are passed as a plain named list of overrides.
  #
  # Do NOT build this with PKNCA.options(...): calling it with arguments sets
  # the options *globally* as a side effect, and returns the complete 17-entry
  # option set including NULL-valued entries (keep_interval_cols, debug) that
  # fail PKNCA's own re-validation when handed back to PKNCAdata(). Both
  # behaviours are wrong here -- a comparison harness must not leak state
  # between engine runs.
  opts <- list(
    auc.method              = rules$auc_method,
    adj.r.squared.factor    = rules$adj_r2_factor,
    min.hl.points           = rules$min_hl_points,
    allow.tmax.in.half.life = rules$allow_tmax_in_hl,
    max.aucinf.pext         = rules$max_aucinf_pext,
    conc.na                 = "drop",
    conc.blq                = list(first  = rules$blq_first,
                                   middle = rules$blq_middle,
                                   last   = rules$blq_last)
  )

  # Formula interfaces. Documented signature is data-first, formula-second.
  f_conc <- stats::as.formula(
    sprintf("%s ~ %s | %s", conc_col, time_col, subject_col))
  f_dose <- stats::as.formula(
    sprintf("%s ~ %s | %s", dose_col, time_col, subject_col))

  o_conc <- PKNCA::PKNCAconc(as.data.frame(conc), f_conc)
  o_dose <- PKNCA::PKNCAdose(as.data.frame(dose), f_dose)

  # Ask for the full single-dose parameter set. PKNCA derives prerequisites
  # automatically (lambda.z is computed because aucinf.obs needs it), but
  # naming them explicitly makes the request self-documenting.
  intervals <- data.frame(
    start              = 0,
    end                = Inf,
    cmax               = TRUE,
    tmax               = TRUE,
    clast.obs          = TRUE,
    tlast              = TRUE,
    auclast            = TRUE,
    aucinf.obs         = TRUE,
    aucpext.obs        = TRUE,
    lambda.z           = TRUE,
    half.life          = TRUE,
    lambda.z.n.points  = TRUE,
    adj.r.squared      = TRUE,
    cl.obs             = TRUE,
    vz.obs             = TRUE,
    aumclast           = TRUE,
    aumcinf.obs        = TRUE,
    mrt.obs            = TRUE
  )

  # Parameter availability shifts a little across PKNCA versions. Ask the
  # package what it actually supports rather than assuming, so a version bump
  # produces a clear message instead of an opaque error deep in pk.nca().
  known <- tryCatch(rownames(PKNCA::get.interval.cols()),
                    error = function(e) NULL)
  if (!is.null(known)) {
    requested <- setdiff(names(intervals), c("start", "end"))
    missing   <- setdiff(requested, known)
    if (length(missing) > 0) {
      warning("PKNCA ", utils::packageVersion("PKNCA"),
              " does not expose: ", paste(missing, collapse = ", "),
              " -- dropping from the request.")
      intervals <- intervals[, !names(intervals) %in% missing, drop = FALSE]
    }
  }

  o_data <- PKNCA::PKNCAdata(o_conc, o_dose,
                             intervals = intervals,
                             options   = opts)
  o_res  <- PKNCA::pk.nca(o_data)

  long <- as.data.frame(o_res$result)

  # PKNCA returns long format: one row per subject per parameter.
  # Translate its parameter names into CDISC PPTESTCD via the crosswalk.
  keep <- long$PPTESTCD %in% PARAM_MAP$pknca
  long <- long[keep, , drop = FALSE]
  long$PPTESTCD_CDISC <- PARAM_MAP$PPTESTCD[
    match(long$PPTESTCD, PARAM_MAP$pknca)]

  out <- data.frame(
    Subject  = long[[subject_col]],
    PPTESTCD = long$PPTESTCD_CDISC,
    value    = suppressWarnings(as.numeric(long$PPORRES)),
    engine   = "PKNCA",
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$PPTESTCD), ]
  rownames(out) <- NULL

  attr(out, "pknca_object") <- o_res
  out
}

#' Pivot the tidy engine output to one row per subject
widen_engine <- function(tidy) {
  w <- stats::reshape(tidy[, c("Subject", "PPTESTCD", "value")],
                      idvar = "Subject", timevar = "PPTESTCD",
                      direction = "wide")
  names(w) <- sub("^value\\.", "", names(w))
  w[order(w$Subject), , drop = FALSE]
}
