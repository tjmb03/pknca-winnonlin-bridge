#' CDISC track: SDTM PC/EX -> NCA -> PP domain
#'
#' The Theophylline track proves the engines agree. This track shows the
#' workflow in the shape a regulatory submission actually takes: start from
#' SDTM domains, run NCA, and emit a CDISC PP-style parameter dataset.
#'
#' Source data is the CDISC pilot study (xanomeline transdermal) distributed in
#' the {pharmaversesdtm} package -- public, CDISC-standard, and the same data
#' the pharmaverse ADPC/ADPPK templates use. Those templates stop at the
#' analysis dataset; this picks up where they leave off and does the NCA.
#'
#' Note this is a multiple-dose study, so the interval specification matters
#' more than it does for single-dose Theophylline. Steady-state parameters
#' (AUCtau, Ctrough, accumulation) are where PKNCA's interval machinery earns
#' its keep and where naive per-subject loops fall apart.

source(file.path("R", "00_config.R"))
source(file.path("R", "02_engine_pknca.R"))

prepare_cdisc_pk <- function() {

  if (!requireNamespace("pharmaversesdtm", quietly = TRUE))
    stop("Install pharmaversesdtm: install.packages('pharmaversesdtm')")
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Install dplyr: install.packages('dplyr')")

  # Package data is lazy-loaded into a data environment, not the namespace, so
  # get(..., envir = asNamespace(...)) fails. Load explicitly into a local env.
  e <- new.env(parent = emptyenv())
  utils::data("pc", package = "pharmaversesdtm", envir = e)
  utils::data("ex", package = "pharmaversesdtm", envir = e)

  pc <- as.data.frame(get("pc", envir = e), stringsAsFactors = FALSE)
  ex <- as.data.frame(get("ex", envir = e), stringsAsFactors = FALSE)

  # ---- Concentrations from PC ----------------------------------------------
  # PCSTRESN is the standardised numeric result; PCSTRESC carries the "<BLQ"
  # text flag. Keep both so BLQ handling is explicit rather than implied by a
  # missing value.
  conc <- pc[, intersect(
    c("STUDYID", "USUBJID", "PCTESTCD", "PCTEST", "PCSTRESN", "PCSTRESC",
      "PCLLOQ", "PCTPT", "PCTPTNUM", "VISITDY", "PCDTC"),
    names(pc)), drop = FALSE]

  # Nominal time from first dose. Pre-dose timepoints carry negative PCTPTNUM;
  # clamp to zero the way the CDISC NCA implementation guide does.
  conc$NFRLT <- ifelse(conc$PCTPTNUM < 0, 0, conc$PCTPTNUM)
  if ("VISITDY" %in% names(conc))
    conc$NFRLT <- conc$NFRLT + 24 * pmax(conc$VISITDY - 1, 0)

  conc$BLQ  <- !is.na(conc$PCSTRESC) & grepl("BLQ|BLOQ|<", conc$PCSTRESC)
  conc$CONC <- conc$PCSTRESN

  # BLQ policy, stated once and applied here rather than scattered:
  #   before first measurable -> 0 (drug genuinely absent)
  #   after                   -> LLOQ/2
  # This mirrors the pharmaverse ADPC template so the two are comparable.
  lloq <- if ("PCLLOQ" %in% names(conc)) conc$PCLLOQ else NA_real_
  conc$CONC <- ifelse(conc$BLQ & conc$NFRLT == 0, 0,
               ifelse(conc$BLQ & conc$NFRLT >  0, 0.5 * lloq, conc$CONC))

  # ---- Doses from EX --------------------------------------------------------
  dose <- ex[ex$EXDOSE > 0, intersect(
    c("STUDYID", "USUBJID", "EXTRT", "EXDOSE", "EXDOSU", "EXDOSFRQ",
      "VISITDY", "EXSTDTC"), names(ex)), drop = FALSE]
  dose$NFRLT <- 24 * pmax(dose$VISITDY - 1, 0)

  # One representative dose per subject for the single-dose interval.
  first_dose <- do.call(rbind, lapply(split(dose, dose$USUBJID), function(d) {
    d[which.min(d$NFRLT), , drop = FALSE]
  }))

  list(conc = conc, dose = dose, first_dose = first_dose)
}

#' Run NCA on the CDISC data with an explicit interval specification
run_cdisc_nca <- function(prepared,
                          analyte = NULL,
                          interval_end = 24) {

  conc <- prepared$conc
  if (!is.null(analyte)) conc <- conc[conc$PCTESTCD == analyte, ]

  conc <- conc[!is.na(conc$CONC) & !is.na(conc$NFRLT), ]
  conc <- conc[order(conc$USUBJID, conc$NFRLT), ]

  dose <- prepared$first_dose[, c("USUBJID", "NFRLT", "EXDOSE")]
  names(dose) <- c("USUBJID", "TIME", "DOSE")

  conc_in <- data.frame(
    USUBJID = conc$USUBJID,
    TIME    = conc$NFRLT,
    CONC    = conc$CONC,
    stringsAsFactors = FALSE
  )

  # Subjects with too few points cannot support a terminal fit; drop them here
  # with a message rather than letting the engine emit 30 warnings.
  n_per <- table(conc_in$USUBJID)
  too_few <- names(n_per)[n_per < NCA_RULES$min_hl_points + 1]
  if (length(too_few) > 0) {
    message("Dropping ", length(too_few), " subject(s) with < ",
            NCA_RULES$min_hl_points + 1, " samples")
    conc_in <- conc_in[!conc_in$USUBJID %in% too_few, ]
    dose    <- dose[!dose$USUBJID %in% too_few, ]
  }

  run_pknca(conc = conc_in, dose = dose,
            subject_col = "USUBJID",
            time_col    = "TIME",
            conc_col    = "CONC",
            dose_col    = "DOSE")
}

#' Emit a CDISC PP-style parameter dataset
#'
#' Not a full PP domain -- that needs PPSEQ, PPGRPID, PPSPEC, and study-level
#' metadata this repo does not have. It is the analysis-ready core that a
#' programmer would join study metadata onto.
to_pp_domain <- function(tidy_nca, studyid = "CDISCPILOT01",
                         domain = "PP", ppcat = "PHARMACOKINETIC") {

  units <- c(
    CMAX = "ug/mL", TMAX = "h", AUCLST = "h*ug/mL", AUCIFO = "h*ug/mL",
    AUCPEO = "%", LAMZ = "1/h", LAMZHL = "h", CLFO = "L/h", VZFO = "L",
    CLST = "ug/mL", TLST = "h", AUMCLST = "h*h*ug/mL",
    AUMCIFO = "h*h*ug/mL", MRTEVIFO = "h", LAMZNPT = "", R2ADJ = ""
  )

  pp <- data.frame(
    STUDYID  = studyid,
    DOMAIN   = domain,
    USUBJID  = tidy_nca$Subject,
    PPTESTCD = tidy_nca$PPTESTCD,
    PPTEST   = PARAM_MAP$label[match(tidy_nca$PPTESTCD, PARAM_MAP$PPTESTCD)],
    PPCAT    = ppcat,
    PPORRES  = tidy_nca$value,
    PPORRESU = unname(units[tidy_nca$PPTESTCD]),
    PPSTRESN = tidy_nca$value,
    PPSTRESU = unname(units[tidy_nca$PPTESTCD]),
    stringsAsFactors = FALSE
  )

  pp <- pp[order(pp$USUBJID, pp$PPTESTCD), ]
  pp$PPSEQ <- stats::ave(seq_len(nrow(pp)), pp$USUBJID, FUN = seq_along)
  rownames(pp) <- NULL

  # Report terminal-fit coverage. Parameters that depend on lambda-z are only
  # available for subjects whose profile supports a terminal regression; in a
  # real study that is routinely well under 100%. Stating it here keeps the
  # same discipline the engine comparison uses -- a gap is a finding, not
  # something to discover later from a column of NAs.
  n_subj <- length(unique(pp$USUBJID))
  n_lamz <- sum(!is.na(pp$PPSTRESN[pp$PPTESTCD == "LAMZHL"]))
  message(sprintf(
    "Terminal-phase fit obtained for %d/%d subjects (%.0f%%). %s",
    n_lamz, n_subj, 100 * n_lamz / n_subj,
    "AUCinf, CL/F, Vz/F and MRT are NA for the remainder."))

  if (n_lamz / n_subj < 0.8)
    warning(sprintf(
      paste0("Only %.0f%% of subjects support a terminal fit. Check sampling ",
             "duration, BLQ handling, and NCA_RULES$min_hl_points (%d) before ",
             "interpreting AUCinf-derived parameters."),
      100 * n_lamz / n_subj, NCA_RULES$min_hl_points))

  pp[, c("STUDYID", "DOMAIN", "USUBJID", "PPSEQ", "PPTESTCD", "PPTEST",
         "PPCAT", "PPORRES", "PPORRESU", "PPSTRESN", "PPSTRESU")]
}
