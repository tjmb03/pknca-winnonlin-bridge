#' Prepare the Theophylline validation dataset
#'
#' datasets::Theoph is a single-dose oral theophylline study in 12 subjects,
#' shipped with base R. It is the anchor dataset for this repo because a
#' Phoenix WinNonlin 6.3/7.0 output for it has been published, giving a real
#' commercial-software reference that anyone can check against without owning
#' a Phoenix license.
#'
#' Produces:
#'   data/derived/theoph_conc.csv  - concentration records
#'   data/derived/theoph_dose.csv  - dosing records
#'   data/derived/theoph_phoenix.csv - Phoenix-ready flat file

source(file.path("R", "00_config.R"))

prepare_theoph <- function(dose_mg = THEOPH_DOSE_MG) {

  d <- datasets::Theoph

  # Subject arrives as an ordered factor with levels sorted by peak
  # concentration, not by subject number. Coerce through character to get the
  # printed label, not the level index -- this is a classic silent bug and it
  # would silently misalign every subject against the WinNonlin reference.
  d$Subject <- as.numeric(as.character(d$Subject))

  conc <- data.frame(
    Subject = d$Subject,
    Time    = d$Time,
    conc    = d$conc,
    Wt      = d$Wt,
    stringsAsFactors = FALSE
  )
  conc <- conc[order(conc$Subject, conc$Time), ]

  # Flat dose convention -- see note in 00_config.R
  dose <- data.frame(
    Subject = sort(unique(conc$Subject)),
    Time    = 0,
    Dose    = dose_mg,
    stringsAsFactors = FALSE
  )

  # Phoenix-ready flat file: one row per observation, dose carried on the
  # subject. Phoenix's NCA object takes dose either from a separate column or
  # from a dosing worksheet; a carried column is the simplest mapping.
  phoenix <- merge(conc, dose[, c("Subject", "Dose")], by = "Subject")
  phoenix <- phoenix[order(phoenix$Subject, phoenix$Time),
                     c("Subject", "Time", "conc", "Dose", "Wt")]
  names(phoenix) <- c("Subject", "Time", "Conc", "Dose", "Weight")

  list(conc = conc, dose = dose, phoenix = phoenix)
}

# ---- Sanity checks ----------------------------------------------------------
validate_theoph <- function(d) {
  stopifnot(
    "expected 12 subjects"        = length(unique(d$conc$Subject)) == 12,
    "expected 132 concentrations" = nrow(d$conc) == 132,
    "no missing concentrations"   = !any(is.na(d$conc$conc)),
    "times are non-negative"      = all(d$conc$Time >= 0),
    "one dose row per subject"    = nrow(d$dose) == length(unique(d$conc$Subject))
  )

  # Each subject must have a pre-dose or near-zero first sample and a
  # descending tail; otherwise lambda-z selection is meaningless.
  by_subj <- split(d$conc, d$conc$Subject)
  n_short <- sum(vapply(by_subj, function(x) nrow(x) < 5L, logical(1)))
  if (n_short > 0)
    warning(n_short, " subject(s) have fewer than 5 samples")

  invisible(TRUE)
}

# ---- Run --------------------------------------------------------------------
# sys.nframe() is 0 only when this file is executed directly (Rscript).
# source()ing it from run_all.R adds a frame, so the block below is skipped and
# the caller drives the pipeline.
if (sys.nframe() == 0L) {
  theoph <- prepare_theoph()
  validate_theoph(theoph)

  write.csv(theoph$conc,
            file.path(PATHS$derived, "theoph_conc.csv"), row.names = FALSE)
  write.csv(theoph$dose,
            file.path(PATHS$derived, "theoph_dose.csv"), row.names = FALSE)
  write.csv(theoph$phoenix,
            file.path(PATHS$derived, "theoph_phoenix.csv"), row.names = FALSE)

  message("Theophylline prepared: ",
          length(unique(theoph$conc$Subject)), " subjects, ",
          nrow(theoph$conc), " concentration records")
  message("Flat dose: ", THEOPH_DOSE_MG, " mg (matches published WinNonlin reference)")
}
