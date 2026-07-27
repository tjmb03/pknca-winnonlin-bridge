#' Three-way engine comparison
#'
#' Takes tidy output from any combination of PKNCA, NonCompart and Phoenix
#' WinNonlin and produces a per-subject, per-parameter agreement table plus
#' summary statistics and diagnostic plots.
#'
#' The comparison is deliberately pairwise-against-a-reference rather than
#' all-against-all: pick the engine you are treating as the system of record
#' and measure everything else against it.

source(file.path("R", "00_config.R"))

#' Combine engine outputs into a comparison table
#'
#' @param ... tidy frames from run_pknca / run_noncompart / read_winnonlin_csv
#' @param reference engine name to treat as the denominator
compare_engines <- function(..., reference = "WinNonlin") {

  frames <- Filter(Negate(is.null), list(...))
  if (length(frames) < 2)
    stop("Need at least two engines to compare; got ", length(frames))

  all <- do.call(rbind, frames)
  engines <- unique(all$engine)

  if (!reference %in% engines) {
    reference <- engines[1]
    message("Reference engine not available; using '", reference, "'")
  }

  wide <- stats::reshape(all,
                         idvar     = c("Subject", "PPTESTCD"),
                         timevar   = "engine",
                         direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))

  others <- setdiff(engines, reference)

  for (e in others) {
    ref <- wide[[reference]]
    tst <- wide[[e]]
    # Percent difference relative to the reference. Guard against a zero
    # denominator, which legitimately occurs for Tlag and some Tmax values.
    pct <- ifelse(is.na(ref) | is.na(tst), NA_real_,
                  ifelse(abs(ref) < .Machine$double.eps,
                         ifelse(abs(tst) < .Machine$double.eps, 0, Inf),
                         100 * (tst - ref) / ref))
    wide[[paste0("pct_diff_", e)]] <- pct

    # A missing value is NOT a pass. Treating NA as agreement is how a
    # validation harness reports success by not testing -- the single most
    # dangerous failure mode here. Missing comparisons are tracked separately
    # and surfaced in the verdict.
    wide[[paste0("pass_", e)]] <-
      !is.na(pct) & abs(pct) <= 100 * NCA_RULES$tolerance
    wide[[paste0("missing_", e)]] <-
      is.na(tst) | is.na(ref)
  }

  attr(wide, "reference") <- reference
  attr(wide, "engines")   <- engines
  wide[order(wide$PPTESTCD, wide$Subject), , drop = FALSE]
}

#' Summarise agreement by parameter
summarise_agreement <- function(cmp) {
  reference <- attr(cmp, "reference")
  engines   <- setdiff(attr(cmp, "engines"), reference)

  do.call(rbind, lapply(engines, function(e) {
    pct  <- cmp[[paste0("pct_diff_", e)]]
    pass <- cmp[[paste0("pass_", e)]]
    miss <- cmp[[paste0("missing_", e)]]

    params <- sort(unique(cmp$PPTESTCD))
    do.call(rbind, lapply(params, function(p) {
      idx <- cmp$PPTESTCD == p
      pv  <- pct[idx]
      pa  <- pass[idx]
      mi  <- miss[idx]
      data.frame(
        engine      = e,
        reference   = reference,
        PPTESTCD    = p,
        label       = PARAM_MAP$label[match(p, PARAM_MAP$PPTESTCD)],
        n_subjects  = sum(idx),
        n_compared  = sum(!mi),
        n_missing   = sum(mi),
        max_abs_pct = if (all(is.na(pv))) NA_real_
                      else max(abs(pv[is.finite(pv)]), na.rm = TRUE),
        n_pass      = sum(pa, na.rm = TRUE),
        n_fail      = sum(!pa & !mi, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

#' Print a compact verdict
report_verdict <- function(agree, tolerance = NCA_RULES$tolerance) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat("ENGINE AGREEMENT (tolerance = ", 100 * tolerance, "%)\n", sep = "")
  cat(strrep("=", 72), "\n", sep = "")

  clean <- TRUE

  for (e in unique(agree$engine)) {
    sub <- agree[agree$engine == e, ]
    n_expected <- sum(sub$n_subjects)
    n_compared <- sum(sub$n_compared)
    n_missing  <- sum(sub$n_missing)
    n_fail     <- sum(sub$n_fail)

    cat("\n", e, " vs ", sub$reference[1], "\n", sep = "")
    cat(sprintf("  coverage : %d/%d comparisons made", n_compared, n_expected))
    if (n_missing > 0) cat(sprintf("  (%d MISSING)", n_missing))
    cat("\n")
    cat(sprintf("  agreement: %d/%d within tolerance\n",
                n_compared - n_fail, n_compared))

    # Coverage gaps come first: a parameter that was never compared is a
    # bigger problem than one that was compared and differed slightly.
    gaps <- sub[sub$n_missing > 0, ]
    if (nrow(gaps) > 0) {
      clean <- FALSE
      cat("  NOT COMPARED (engine returned no value):\n")
      for (i in seq_len(nrow(gaps)))
        cat(sprintf("    %-10s %-20s %d/%d subjects\n",
                    gaps$PPTESTCD[i], gaps$label[i],
                    gaps$n_missing[i], gaps$n_subjects[i]))
    }

    failing <- sub[sub$n_fail > 0, ]
    if (nrow(failing) > 0) {
      clean <- FALSE
      cat("  EXCEEDS TOLERANCE:\n")
      for (i in seq_len(nrow(failing)))
        cat(sprintf("    %-10s %-20s %d/%d subjects, max |diff| %.4f%%\n",
                    failing$PPTESTCD[i], failing$label[i],
                    failing$n_fail[i], failing$n_compared[i],
                    failing$max_abs_pct[i]))
    }

    if (nrow(gaps) == 0 && nrow(failing) == 0)
      cat("  All parameters compared and in agreement.\n")
  }

  cat("\n", strrep("=", 72), "\n", sep = "")
  if (!clean)
    cat("VERDICT: NOT CLEAN -- investigate gaps and/or differences above.\n")
  else
    cat("VERDICT: CLEAN -- full coverage, all within tolerance.\n")
  cat(strrep("=", 72), "\n", sep = "")

  invisible(list(summary = agree, clean = clean))
}

#' Scatter plots, one facet per parameter
plot_agreement <- function(cmp, engine, reference = attr(cmp, "reference"),
                           path = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed; skipping plots")
    return(invisible(NULL))
  }

  d <- data.frame(
    PPTESTCD = cmp$PPTESTCD,
    ref      = cmp[[reference]],
    test     = cmp[[engine]],
    stringsAsFactors = FALSE
  )
  d <- d[is.finite(d$ref) & is.finite(d$test), ]
  d$label <- PARAM_MAP$label[match(d$PPTESTCD, PARAM_MAP$PPTESTCD)]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = ref, y = test)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey55") +
    ggplot2::geom_point(size = 1.8, alpha = 0.75) +
    ggplot2::facet_wrap(~ label, scales = "free") +
    ggplot2::labs(
      title = paste0(engine, " vs ", reference),
      subtitle = "Dashed line is identity; departures indicate convention differences",
      x = reference, y = engine) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey92"))

  if (!is.null(path)) {
    ggplot2::ggsave(path, p, width = 11, height = 8, dpi = 120)
    message("Plot written: ", path)
  }
  p
}

#' Attribute an observed difference to reference precision vs engine error
#'
#' The headline agreement number is easy to misread. A comparison against a
#' reference file that stores values to a fixed number of decimal places is
#' bounded below by that file's rounding, not by the engines' accuracy. For
#' small-magnitude parameters (lambda-z especially) the rounding term dominates
#' entirely.
#'
#' This computes, per parameter, the largest relative difference attributable
#' purely to the reference's stored precision, so it can be set against the
#' observed difference. If observed <= predicted, the comparison is precision-
#' limited and says nothing about engine error.
#'
#' @param cmp output of compare_engines()
#' @param ref_path path to the reference CSV, read as text to recover the
#'   decimal places actually written (read.csv would discard them)
attribute_precision <- function(cmp, ref_path, engine, subject_col = "Subject") {

  if (!file.exists(ref_path)) return(NULL)

  raw  <- readLines(ref_path)
  hdr  <- strsplit(raw[1], ",")[[1]]
  body <- lapply(raw[-1], function(l) strsplit(l, ",")[[1]])

  decimals <- function(x) {
    x <- trimws(x)
    if (!grepl("\\.", x)) return(0L)
    nchar(sub("^[^.]*\\.", "", x))
  }

  reference <- attr(cmp, "reference")
  params <- sort(unique(cmp$PPTESTCD))

  do.call(rbind, lapply(params, function(p) {
    wcol <- PARAM_MAP$winnonlin[match(p, PARAM_MAP$PPTESTCD)]
    # read.csv mangles '%' to '.', so match on the sanitised name too
    j <- which(hdr == wcol | make.names(hdr) == make.names(wcol))
    if (length(j) == 0) return(NULL)

    idx <- cmp$PPTESTCD == p
    vals <- cmp[[reference]][idx]
    dp <- vapply(body, function(r) decimals(r[j[1]]), integer(1))
    dp <- dp[seq_along(vals)]

    # Largest relative error a half-ulp at that decimal place could produce
    pred <- ifelse(abs(vals) < .Machine$double.eps, NA_real_,
                   100 * (0.5 * 10^(-dp)) / abs(vals))
    obs  <- abs(cmp[[paste0("pct_diff_", engine)]][idx])

    data.frame(
      PPTESTCD          = p,
      label             = PARAM_MAP$label[match(p, PARAM_MAP$PPTESTCD)],
      ref_decimals      = if (length(unique(dp)) == 1) as.character(dp[1])
                          else paste0(min(dp), "-", max(dp)),
      max_observed_pct  = suppressWarnings(max(obs,  na.rm = TRUE)),
      max_rounding_pct  = suppressWarnings(max(pred, na.rm = TRUE)),
      precision_limited = suppressWarnings(
        max(obs, na.rm = TRUE) <= max(pred, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
}

#' Report engine-vs-engine agreement, free of any reference-file rounding
engine_vs_engine <- function(tidy_a, tidy_b) {
  cmp <- compare_engines(tidy_a, tidy_b,
                         reference = unique(tidy_a$engine)[1])
  other <- setdiff(attr(cmp, "engines"), attr(cmp, "reference"))
  d <- abs(cmp[[paste0("pct_diff_", other)]])
  list(
    n              = nrow(cmp),
    n_bit_identical = sum(d == 0, na.rm = TRUE),
    max_pct_diff   = suppressWarnings(max(d, na.rm = TRUE)),
    double_eps_pct = 100 * .Machine$double.eps
  )
}
